/*
 * GlobalPrrController.cc
 *
 *  Created on: Nov 25, 2025
 *      Author: yychen
 */


#include "GlobalPrrController.h"

#include "inet/physicallayer/wireless/common/medium/RadioMedium.h"
#include "inet/physicallayer/wireless/common/signal/WirelessSignal.h"
#include "inet/common/ModuleAccess.h"
#include "inet/mobility/contract/IMobility.h"
#include "inet/common/packet/Packet.h"
#include "inet/networklayer/common/L3AddressResolver.h"
#include "inet/common/Units.h"

#include "inet/common/INETUtils.h"
#include "inet/common/ProtocolTag_m.h"
#include "inet/common/Simsignals.h"
#include "inet/common/packet/chunk/FieldsChunk.h"
#include "inet/linklayer/common/MacAddressTag_m.h"
#include "inet/networklayer/common/NetworkInterface.h"
#include "inet/networklayer/contract/IInterfaceTable.h"
#include "inet/physicallayer/wireless/common/analogmodel/bitlevel/LayeredTransmission.h"
#include "inet/physicallayer/wireless/common/radio/packetlevel/Radio.h"
#include "inet/physicallayer/wireless/common/signal/Interference.h"
#include "inet/networklayer/common/InterfaceTable.h"


void GlobalPrrController::initialize(int stage)
{
    cSimpleModule::initialize(stage);
    if (stage == INITSTAGE_LOCAL) {
        const char *modeStr = par("mode");
        if (!strcmp(modeStr, "fast"))
            mode = FAST;
        else if (!strcmp(modeStr, "physical"))
            mode = PHYSICAL;
        else
            throw cRuntimeError("Unknown mode parameter: %s", modeStr);

        // find radioMedium (common path: system module has radioMedium)
        cModule *system = getSystemModule();
        cModule *rm = system->getSubmodule("radioMedium");
        if (!rm) {
            // fallback: try find by glob
            rm = getModuleByPath("*.radioMedium");
        }
        if (!rm)
            throw cRuntimeError("radioMedium not found in network");
        radioMedium = check_and_cast<IRadioMedium*>(rm);

        // try to get commonly used submodules (may need casts)
        if (rm->getSubmodule("pathLoss"))
            pathLoss = check_and_cast<IPathLoss*>(rm->getSubmodule("pathLoss"));
        if (rm->getSubmodule("obstacleLoss"))
            obstacleLoss = dynamic_cast<IObstacleLoss*>(rm->getSubmodule("obstacleLoss"));
        if (rm->getSubmodule("backgroundNoise"))
            backgroundNoise = dynamic_cast<IBackgroundNoise*>(rm->getSubmodule("backgroundNoise"));
        if (rm->getSubmodule("propagation"))
            propagation = check_and_cast<IPropagation*>(rm->getSubmodule("propagation"));
        if (rm->getSubmodule("mediumLimitCache"))
            mediumLimitCache = check_and_cast<IMediumLimitCache*>(rm->getSubmodule("mediumLimitCache"));

        collectRadios();

        // Register signals (MUST match NED names)
        reachableSignal = registerSignal("reachable");
        receiverCountSignal = registerSignal("receiverCount");
        theoreticalPRRSignal = registerSignal("theoreticalPRR");

        EV_ERROR << "GlobalPrrController initialized. mode=" << (mode==FAST ? "FAST" : "PHYSICAL")
                << ", radios=" << radios.size() << "\n";
    }
    if (stage == INITSTAGE_LAST) {
        // 1) register signals
        //reachableSignal = registerSignal("reachable");
        //receiverCountSignal = registerSignal("receiverCount");
        //theoreticalPRRSignal = registerSignal("theoreticalPRR");

        // 2) parent is network,GlobalPrrController and node[...] at the same level,through getParentModule() get network
        cModule *parent = getParentModule();
        if (!parent) {
            throw cRuntimeError("GlobalPrrController: parent module not found");
        }

        // 3) get the size of node[] (assume named "node")
        //    use getVectorSize("node")

        // try reading parameter named "numNodes" from parent; fallback to probing
        int numNodes = -1;
        if (parent->hasPar("numNodes"))
            numNodes = parent->par("numNodes").intValue();

        nodes.clear();
        if (numNodes >= 0) {
            EV_ERROR<<"GlobalPrrController::initialize numNodes = "<<numNodes<<"\n";
            for (int i = 0; i < numNodes; ++i) {
                cModule *n = parent->getSubmodule("node", i);
                if (!n) {
                    EV_WARN << "GlobalPrrController: expected node[" << i << "] but got nullptr\n";
                    continue;
                }
                nodes.push_back(n);
            }
        } else {
            // fallback probe if parameter not present
            EV_ERROR<<"GlobalPrrController::initialize numNodes fallback\n";
            int i = 0;
            while (true) {
                cModule *n = parent->getSubmodule("node", i);
                if (!n) break;
                nodes.push_back(n);
                ++i;
            }
        }
        nodeVectorSize = (int)nodes.size();


        //nodeVectorSize = parent->getVectorSize("node");
        EV_ERROR << "GlobalPrrController: found node vector size = " << nodeVectorSize << "\n";

        // 4) cache the ptr of node Module
        /*nodes.clear();
        for (int i = 0; i < nodeVectorSize; ++i) {
            cModule *n = parent->getSubmodule("node", i);
            if (n == nullptr) {
                EV_WARN << "GlobalPrrController: parent->getSubmodule(\"node\", " << i << ") returned nullptr\n";
                continue;
            }
            nodes.push_back(n);
            EV_DEBUG << "GlobalPrrController: cached node: " << n->getFullPath() << "\n";
        }*/

        // 5) if radios list is needed,do as follow
        //
        // radios.clear();
        // for (cModule* n : nodes) {
        //     cModule *wlan = n->getSubmodule("wlan", 0);
        //     if (wlan) {
        //         cModule *r = wlan->getSubmodule("radio");
        //         if (r) radios.push_back(check_and_cast<Radio*>(r));
        //     } else {
        //         cModule *r2 = n->getSubmodule("radio");
        //         if (r2) radios.push_back(check_and_cast<Radio*>(r2));
        //     }
        // }
    }
}


// return node module containing the radio (or nullptr)
static cModule* getNodeOfRadio(const IRadio *radio)
{
    if (!radio)
        return nullptr;
    // radio is a module derived from cModule in INET; getContainingNode expects cModule*
    const cModule *radioModule = dynamic_cast<const cModule*>(radio);
    if (!radioModule)
        return nullptr;
    return getContainingNode(const_cast<cModule*>(radioModule));
}

cModule* GlobalPrrController::findNodeByAddress(const L3Address &addr)
{
    if (addr.isUnspecified())
        return nullptr;

    for (cModule *node : nodes) {
        // get IInterfaceTable of node
        IInterfaceTable *ift = L3AddressResolver().findInterfaceTableOf(node);
        if (!ift)
            continue;

        // NetworkInterface*
        NetworkInterface *nic = ift->findInterfaceByAddress(addr);

        if (nic != nullptr)
            return node;
    }
    return nullptr;
}


// FAST: use pathLoss + obstacleLoss + backgroundNoise + receiver sensitivity
bool GlobalPrrController::canRadioReachRadio_Fast(Radio *txRadio, Radio *rxRadio, const PacketMeta &meta)
{
    //EV_ERROR<<"GlobalPrrController::canRadioReachRadio_Fast\n";
    // get positions
    cModule *txNode = getContainingNode(txRadio);
    cModule *rxNode = getContainingNode(rxRadio);
    if (!txNode || !rxNode) return false;

    IMobility *mtx = nullptr;
    IMobility *mrx = nullptr;
    if (txNode->getSubmodule("mobility"))
        mtx = check_and_cast<IMobility*>(txNode->getSubmodule("mobility"));
    if (rxNode->getSubmodule("mobility"))
        mrx = check_and_cast<IMobility*>(rxNode->getSubmodule("mobility"));
    Coord txPos = mtx ? mtx->getCurrentPosition() : Coord::ZERO;
    Coord rxPos = mrx ? mrx->getCurrentPosition() : Coord::ZERO;

    m distance = m(txPos.distance(rxPos));
    Hz freq = meta.centerFrequency;
    if (std::isnan(freq.get()))
        freq = Hz(2.4e9); // fallback

    mps propagationSpeed = propagation ? propagation->getPropagationSpeed() : mps(299792458.0);

    double pathLossFraction = 1.0;
    if (pathLoss)
        pathLossFraction = pathLoss->computePathLoss(propagationSpeed, freq, distance);

    double obstacleLoss_dB = 0.0;
    if (obstacleLoss)
        obstacleLoss_dB = obstacleLoss->computeObstacleLoss(freq, txPos, rxPos);
    double obstacleLinear = pow(10.0, -obstacleLoss_dB/10.0);

    // tx power
    double Pt_dBm = meta.txPower_dBm;
    if (std::isnan(Pt_dBm)) {
        // attempt to read param from transmitter node
        cModule *node = getContainingNode(txRadio);
        if (node && node->hasPar("txPower")) Pt_dBm = node->par("txPower").doubleValue();
        else Pt_dBm = 0.0; // fallback
    }
    double Pt_mW = math::dBmW2mW(Pt_dBm);

    // antenna gains: skipping detailed model  assume 1.0 unless you read antenna objects
    double txGain = 1.0, rxGain = 1.0;

    double Pr_mW = Pt_mW * pathLossFraction * txGain * rxGain * obstacleLinear;
    double Pr_dBm = math::mW2dBmW(Pr_mW);

    // receiver sensitivity heuristic: read param "minReceptionPower" if exists
    double rxSens_dBm = -85.0; // fallback
    if (rxRadio->hasPar("minReceptionPower"))
        rxSens_dBm = rxRadio->par("minReceptionPower").doubleValue();

    return Pr_dBm >= rxSens_dBm;
}

// PHYSICAL: use RadioMedium::testReception (the wrapper you added). Build a transient transmission for txRadio using its transmitter.
// WARNING: ownership of tmp packet/transmission depends on INET; we try to be conservative.
bool GlobalPrrController::canRadioReachRadio_Physical(Radio *txRadio, Radio *rxRadio, const PacketMeta &meta)
{
    RadioMedium *rm = dynamic_cast<RadioMedium*>(radioMedium);
    if (!rm) {
        // fallback
        return canRadioReachRadio_Fast(txRadio, rxRadio, meta);
    }

    Packet *tmpPkt = new Packet("globalprr-phy-snap");
    // optionally add tags if transmitter expects them (e.g., center frequency)
    const ITransmission *transmission = nullptr;
    try {
        if (txRadio->getTransmitter())
            transmission = txRadio->getTransmitter()->createTransmission(txRadio, tmpPkt, simTime());
    } catch(...) {
        // creation failed: fallback
        delete tmpPkt;
        return canRadioReachRadio_Fast(txRadio, rxRadio, meta);
    }
    if (!transmission) {
        delete tmpPkt;
        return canRadioReachRadio_Fast(txRadio, rxRadio, meta);
    }

    bool ok = false;
    try {
        ok = rm->testReception(txRadio, transmission, rxRadio); // your wrapper
    } catch(...) {
        ok = false;
    }

    // cleanup; conservative: delete tmpPkt (but check your createTransmission semantics)
    delete tmpPkt;
    // do NOT delete transmission here unless your createTransmission contract says so.
    return ok;
}

bool GlobalPrrController::canNodeReachNode(cModule *srcNode, cModule *dstNode, const PacketMeta &meta)
{
    //EV_ERROR<<"GlobalPrrController::canNodeReachNode\n";
    if (!srcNode || !dstNode) return false;
    // gather radios of src and dst
    std::vector<Radio*> srcRadios;
    std::vector<Radio*> dstRadios;
    for (Radio *r : radios) {
        cModule *n = getContainingNode(r);
        if (n == srcNode) srcRadios.push_back(r);
        if (n == dstNode) dstRadios.push_back(r);
    }
    if (srcRadios.empty() || dstRadios.empty()) return false;

    for (auto txR : srcRadios) {
        for (auto rxR : dstRadios) {
            bool ok = false;
            if(mode == PHYSICAL){
                ok = canRadioReachRadio_Physical(txR, rxR, meta);
                //EV_ERROR<<"GlobalPrrController::canNodeReachNode,PHYSICAL\n";
            }
            else if(mode == FAST){
                ok = canRadioReachRadio_Fast(txR, rxR, meta);
                //EV_ERROR<<"GlobalPrrController::canNodeReachNode,fast\n";
            }
            //bool ok = (mode == PHYSICAL) ? canRadioReachRadio_Physical(txR, rxR, meta) : canRadioReachRadio_Fast(txR, rxR, meta);
            if (ok) return true;
        }
    }
    return false;
}

/*bool GlobalPrrController::isDestinationReachableFromSource(cModule *srcNode, cModule *dstNode, const PacketMeta &meta)
{
    // build adjacency list per node
    std::vector<cModule*> nodes;
    std::unordered_map<cModule*, int> nodeIndex;
    cModule *system = getSystemModule();
    for (cModule::SubmoduleIterator it(system); !it.end(); ++it) {
        cModule *m = *it;
        // skip radioMedium / controller etc  heuristics: nodes that have "interfaceTable" are hosts
        if (m->getSubmodule("interfaceTable")) {
            nodeIndex[m] = nodes.size();
            nodes.push_back(m);
        }
    }
    int N = nodes.size();
    if (nodeIndex.find(srcNode) == nodeIndex.end() || nodeIndex.find(dstNode) == nodeIndex.end())
        return false;

    std::vector<std::vector<int>> adj(N);
    // fill edges (directed)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            cModule *ni = nodes[i];
            cModule *nj = nodes[j];
            // check if ni -> nj is possible
            bool edge = canNodeReachNode(ni, nj, meta);
            if (edge) adj[i].push_back(j);
        }
    }

    // BFS from src to dst
    int s = nodeIndex[srcNode];
    int t = nodeIndex[dstNode];
    std::vector<char> vis(N, 0);
    std::deque<int> dq;
    vis[s] = 1;
    dq.push_back(s);
    while (!dq.empty()) {
        int u = dq.front(); dq.pop_front();
        if (u == t) return true;
        for (int v : adj[u]) {
            if (!vis[v]) { vis[v]=1; dq.push_back(v); }
        }
    }
    return false;
}*/
bool GlobalPrrController::isDestinationReachableFromSource(cModule *srcNode,cModule *dstNode,const PacketMeta &meta)
{
    if (!srcNode || !dstNode)
    return false;
    //EV_ERROR<<"GlobalPrrController::isDestinationReachableFromSource begin\n";
    int N = nodes.size();

    // build nodeIndex
    std::unordered_map<cModule*, int> nodeIndex;
    nodeIndex.reserve(N);
    for (int i = 0; i < N; ++i)
        nodeIndex[nodes[i]] = i;

    //  source or destination not at nodes
    if (!nodeIndex.count(srcNode) || !nodeIndex.count(dstNode)){
        //EV_ERROR<<"GlobalPrrController::isDestinationReachableFromSource error: source or destination not at nodes\n";
        return false;
    }

    int srcIdx = nodeIndex[srcNode];
    int dstIdx = nodeIndex[dstNode];

    //
    std::vector<std::vector<int>> adj(N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            cModule *ni = nodes[i];
            cModule *nj = nodes[j];

            if (canNodeReachNode(ni, nj, meta)) {
                adj[i].push_back(j);
            }
        }
    }

    // BFS
    std::vector<char> visited(N, 0);
    std::deque<int> dq;

    visited[srcIdx] = 1;
    dq.push_back(srcIdx);

    while (!dq.empty()) {
        int u = dq.front();
        dq.pop_front();

        if (u == dstIdx)
            return true;    //

        for (int v : adj[u]) {
            if (!visited[v]) {
                visited[v] = 1;
                dq.push_back(v);
            }
        }
    }

    return false; //
}
void GlobalPrrController::collectRadios()
{
    radios.clear();
    cModule *system = getSystemModule();
    // common host layout: each node may have wlan[0].radio or radio
    for (cModule::SubmoduleIterator it(system); !it.end(); ++it) {
        cModule *m = *it;
        cModule *wlan = m->getSubmodule("wlan", 0);
        if (wlan) {
            cModule *r = wlan->getSubmodule("radio");
            if (r)
                radios.push_back(check_and_cast<Radio*>(r));
            continue;
        }
        cModule *r2 = m->getSubmodule("radio");
        if (r2)
            radios.push_back(check_and_cast<Radio*>(r2));
    }
}

Radio *GlobalPrrController::findRadioForNode(cModule *node)
{
    if (!node) return nullptr;
    // try wlan[0].radio
    cModule *wlan = node->getSubmodule("wlan", 0);
    if (wlan) {
        cModule *r = wlan->getSubmodule("radio");
        if (r) return check_and_cast<Radio*>(r);
    }
    cModule *r = node->getSubmodule("radio");
    if (r) return check_and_cast<Radio*>(r);
    // fallback: search in cached radios by ancestor
    for (auto rad : radios) {
        cModule *radNode = getContainingNode(rad);
        if (radNode == node) return rad;
    }
    return nullptr;
}

/*void GlobalPrrController::notifySend(cModule *senderNode, L3Address destAddr, const PacketMeta& meta)
{
    Enter_Method("notifySend");
    totalTx++;
    Radio *txRadio = findRadioForNode(senderNode);
    if (!txRadio) {
        EV_WARN << "notifySend: radio for node " << senderNode->getFullPath() << " not found\n";
        return;
    }
    if (mode == FAST)
        handleNotifyFast(txRadio, senderNode, destAddr, meta);
    else
        handleNotifyPhysical(txRadio, senderNode, destAddr, meta);
}*/

void GlobalPrrController::notifySend(cModule *senderNode, L3Address destAddr, const PacketMeta& meta) {
    Enter_Method("notifySend");
    totalTx++;
    // find destination node by address
    cModule *destNode = findNodeByAddress(destAddr);
    if (!destNode) {
        EV_ERROR << "notifySend: destination node for address " << destAddr << " not found. Marking as unreachable.\n";
        // record scalar or emit unreachable signal
        // update stats accordingly
        emit(reachableSignal, false);
        emit(receiverCountSignal, 0);
        emit(theoreticalPRRSignal, 0.0);
        return;
    }
    bool reachable = isDestinationReachableFromSource(senderNode, destNode, meta);
    if (reachable) {
        EV_ERROR<<"packet reachable\n";
        packetsWithAtLeastOneTheoreticalRx++; // or more appropriate name e.g., packetsThatReachDestination
        theoreticalTotalReceives += 1; // or different metric for end-to-end (maybe count per-packet)
    }
    else{
        EV_ERROR<<"packet not reachable\n";
    }
    // emit signals if declared in NED
    emit(reachableSignal, reachable ? true : false);
    emit(receiverCountSignal, reachable ? 1 : 0); // receiverCount for end-to-end maybe 1 if reachable else 0
    emit(theoreticalPRRSignal, reachable ? 1.0 : 0.0);
}

// FAST :use pathLoss/obstacleLoss/backgroundNoise and compute directly
/*void GlobalPrrController::handleNotifyFast(Radio *txRadio, cModule *senderNode, L3Address destAddr, const PacketMeta &meta)
{
    // 1) gather tx params
    double Pt_dBm = meta.txPower_dBm;
    if (std::isnan(Pt_dBm)) {
        // try to read from transmitter if available
        if (txRadio && txRadio->getTransmitter()) {
            // Many transmitter classes expose getPower() or getGain(); try to read conservatively
            try {
                // This is a heuristic: some Transmitter implementations have getPower() returning W
                // Replace with actual API if different in your version
                const ITransmitter *tx = txRadio->getTransmitter();
                // if transmitter exposes getPower(), use it; else keep NAN
                // (No compile-time guarantee here; if your transmitter has different API, adapt.)
                // Example fallback: read parameter from ancestor node
            } catch (...) {
                // ignore
            }
        }
        // fallback: try to read node parameter
        if (senderNode && senderNode->hasPar("txPower"))
            Pt_dBm = senderNode->par("txPower").doubleValue();
    }

    if (std::isnan(Pt_dBm)) {
        EV_WARN << "handleNotifyFast: txPower unknown, assuming 0 dBm for estimate\n";
        Pt_dBm = 0.0;
    }

    // center frequency
    Hz freq = meta.centerFrequency;
    if (std::isnan(freq.get())) {
        // try read from txRadio parameters or ancestor
        if (txRadio && txRadio->hasPar("centerFrequency"))
            freq = Hz(txRadio->par("centerFrequency"));
        else {
            EV_WARN << "handleNotifyFast: centerFrequency unknown, using 2.4GHz default\n";
            freq = Hz(2.4e9);
        }
    }

    // propagation speed
    mps propagationSpeed = propagation ? propagation->getPropagationSpeed() : mps(299792458.0);

    // prepare candidate list: ideally use mediumLimitCache to prune
    std::vector<Radio*> candidates;
    if (mediumLimitCache) {
        // If mediumLimitCache provides API for getting candidate radios by tx position, use it.
        // Fallback: use full radios list
        for (auto r : radios)
            if (r != txRadio)
                candidates.push_back(r);
    }
    else {
        for (auto r : radios)
            if (r != txRadio)
                candidates.push_back(r);
    }

    // tx position
    Coord txPos;
    {
        cModule *txNode = getContainingNode(txRadio);
        IMobility *mtx = nullptr;
        if (txNode && txNode->getSubmodule("mobility"))
            mtx = check_and_cast<IMobility*>(txNode->getSubmodule("mobility"));
        if (mtx)
            txPos = mtx->getCurrentPosition();
        else
            txPos = Coord::ZERO;
    }

    double Pt_mW = math::dBmW2mW(Pt_dBm);

    bool anyReceiver = false;
    int perPacketReceives = 0;

    for (auto rxRadio : candidates) {
        // rx pos
        cModule *rxNode = getContainingNode(rxRadio);
        if (!rxNode) continue;
        IMobility *mrx = nullptr;
        if (rxNode->getSubmodule("mobility"))
            mrx = check_and_cast<IMobility*>(rxNode->getSubmodule("mobility"));
        Coord rxPos = mrx ? mrx->getCurrentPosition() : Coord::ZERO;

        // distance
        m dist = m(txPos.distance(rxPos));

        // path loss (FreeSpacePathLoss::computePathLoss returns linear fraction in INET 4.x)
        double pathLossFraction = 1.0;
        if (pathLoss) {
            // computePathLoss(propagationSpeed, freq, distance)
            pathLossFraction = pathLoss->computePathLoss(propagationSpeed, freq, dist);
            if (pathLossFraction <= 0.0)
                pathLossFraction = 1e-300; // avoid zero
        }

        // obstacle loss (dB)
        double obstacleLoss_dB = 0.0;
        if (obstacleLoss) {
            obstacleLoss_dB = obstacleLoss->computeObstacleLoss(freq, txPos, rxPos);
        }
        double obstacleLinear = math::dB2fraction(-obstacleLoss_dB); // careful: computeObstacleLoss returns dB loss; subtracting it already applied later
        // Above: if computeObstacleLoss returns positive dB loss (e.g., 10 dB), we need multiply Pr by 10^{-loss/10}
        obstacleLinear = pow(10.0, -obstacleLoss_dB/10.0);

        // simple antenna gain handling: try to read simple antenna gains if available (fallback 1.0)
        double txAntennaGain = 1.0;
        double rxAntennaGain = 1.0;
        // many radio/antenna modules have parameters "antennaGain" or antenna objects; skip detailed lookup for brevity

        // Received power in mW (linear)
        double Pr_mW = Pt_mW * pathLossFraction * txAntennaGain * rxAntennaGain * obstacleLinear;
        double Pr_dBm = math::mW2dBmW(Pr_mW);

        // get receiver sensitivity: try receiver object (IReceiver), else fallback to parameter
        double rxSensitivity_dBm = -85.0;
        const IReceiver *receiver = rxRadio->getReceiver();
        if (receiver) {
            // Many IReceiver implementations give minReceptionPower in W, but no guaranteed method name.
            // Try typical API: getMinReceptionPower() -> returns W
            try {
                // This is not a compile-time guaranteed API; adapt if necessary.
                // We'll try to dynamic_cast to a type that has getMinReceptionPower (if known).
                // Fallback: read module parameter
            } catch (...) {}
        }
        // fallback read param
        if (rxRadio->hasPar("minReceptionPower"))
            rxSensitivity_dBm = rxRadio->par("minReceptionPower").doubleValue();

        // final check
        if (Pr_dBm >= rxSensitivity_dBm) {
            anyReceiver = true;
            perPacketReceives++;
        }
    } // end candidates

    if (anyReceiver) packetsWithAtLeastOneTheoreticalRx++;
    theoreticalTotalReceives += perPacketReceives;
}
// PHYSICAL :use new RadioMedium::testReception(...)wrapper
void GlobalPrrController::handleNotifyPhysical(Radio *txRadio, cModule *senderNode, L3Address destAddr, const PacketMeta &meta)
{
    // We rely on RadioMedium::testReception(transmitter, transmission, receiver) which you added.
    RadioMedium *rm = dynamic_cast<RadioMedium*>(radioMedium);
    if (!rm) {
        EV_WARN << "handleNotifyPhysical: radioMedium is not RadioMedium*; falling back to FAST\n";
        handleNotifyFast(txRadio, senderNode, destAddr, meta);
        return;
    }

    // create a dummy packet to form a transmission. Ownership semantics vary by INET version.
    Packet *tmpPkt = new Packet("globalprr-theory-snapshot");
    // Optionally add tags to tmpPkt if transmitter requires them (e.g., center frequency)
    // e.g., tmpPkt->addTag<...>()

    // create transmission via transmitter; typical signature: createTransmission(radio, packet, simTime());
    const ITransmission *transmission = nullptr;
    try {
        if (txRadio->getTransmitter())
            transmission = txRadio->getTransmitter()->createTransmission(txRadio, tmpPkt, simTime());
    } catch (std::exception &e) {
        EV_WARN << "handleNotifyPhysical: createTransmission threw: " << e.what() << ", fallback to FAST\n";
        delete tmpPkt;
        handleNotifyFast(txRadio, senderNode, destAddr, meta);
        return;
    }
    if (!transmission) {
        EV_WARN << "handleNotifyPhysical: transmission creation failed, fallback to FAST\n";
        delete tmpPkt;
        handleNotifyFast(txRadio, senderNode, destAddr, meta);
        return;
    }

    bool anyReceiver = false;
    int perPacketReceives = 0;

    for (auto r : radios) {
        if (r == txRadio) continue;
        // Use the wrapper you added to RadioMedium
        bool ok = false;
        try {
            ok = rm->testReception(txRadio, transmission, r);
        } catch (...) {
            EV_WARN << "handleNotifyPhysical: testReception threw exception; treating as not received\n";
            ok = false;
        }
        if (ok) {
            anyReceiver = true;
            perPacketReceives++;
        }
    }
    // CLEANUP: we created tmpPkt and transmission; memory ownership semantics differ across INET versions:
    // - In many INET versions createTransmission does NOT delete the packet: caller should delete the packet.
    // - The transmission object might be owned by medium after addTransmission; in our testReception wrapper we transiently add and remove,
    //   so medium doesn't take permanent ownership. The safe option is to delete tmpPkt and let transmission object be destroyed if needed.
    // But deleting transmission object here may be invalid if it is managed internally. To be conservative:
    // - delete tmpPkt (we created it)
    // - do NOT delete transmission (let implementation manage it)  if you observe leaks, you must inspect createTransmission ownership.
    delete tmpPkt;

    if (anyReceiver) packetsWithAtLeastOneTheoreticalRx++;
    theoreticalTotalReceives += perPacketReceives;
}*/


void GlobalPrrController::finish()
{
    cSimpleModule::finish();
    recordScalar("total_tx", totalTx);
    recordScalar("theoretical_packets_with_at_least_one_rx", packetsWithAtLeastOneTheoreticalRx);
    recordScalar("theoretical_total_receiver_events", theoreticalTotalReceives);
    if (totalTx > 0) {
        double ratio = (double)packetsWithAtLeastOneTheoreticalRx / (double)totalTx;
        recordScalar("theoretical_at_least_one_rx_ratio", ratio);
    }
    EV_ERROR << "GlobalPrrController finished: totalTx=" << totalTx << " packetsWithAtLeastOneTheoreticalRx=" << packetsWithAtLeastOneTheoreticalRx << "\n";
}


