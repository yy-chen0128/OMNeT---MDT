#include "../../mysrc/routing/MDTRouting.h"

#include "inet/common/Ptr.h"

#include "inet/common/IProtocolRegistrationListener.h"
#include "inet/common/ModuleAccess.h"
#include "inet/common/ProtocolTag_m.h"
#include "inet/common/packet/Packet.h"
#include "inet/common/stlutils.h"
#include "inet/linklayer/common/InterfaceTag_m.h"
#include "inet/networklayer/common/HopLimitTag_m.h"
#include "inet/networklayer/common/L3AddressResolver.h"
#include "inet/networklayer/common/L3AddressTag_m.h"
#include "inet/networklayer/common/L3Tools.h"
#include "inet/networklayer/ipv4/IcmpHeader.h"
#include "inet/networklayer/ipv4/Ipv4Header_m.h"
#include "inet/networklayer/ipv4/Ipv4Route.h"
#include "inet/transportlayer/common/L4PortTag_m.h"
#include "inet/transportlayer/contract/udp/UdpControlInfo.h"

#include "inet/common/ModuleAccess.h"
#include "inet/mobility/contract/IMobility.h"
#include "../../mysrc/routing/MdtRelayTag_m.h"
#include "../../mysrc/routing/MdtRelayChunk_m.h"

#include "../../mysrc/routing/Neighbor.h"
#include "ForwardProtocol.h"
#include "inet/linklayer/common/InterfaceTag_m.h"
#include "JoinProtocol.h"
#include "MaintenanceProtocol.h"
#include "inet/networklayer/common/NextHopAddressTag_m.h"

#include "Neighbor.h"

#include <iomanip>
#include "inet/common/packet/chunk/Chunk.h"
namespace mysrc {
namespace routing {
Define_Module(MDTRouting);

const int KIND_DELAYEDSEND = 100;//delay
const int KIND_DELAYEDSEND_NS = 101;//delay

//KA
static double dBmToW(double dBm) {
    return pow(10.0, (dBm - 30.0) / 10.0);
}
static double wToDbm(double W) {
    return 10.0 * log10(std::max(W, 1e-300)) + 30.0;
}
static double theoreticalReceivedPower_dBm(double Ptx_dBm, double freq_Hz, double d_m, double alpha) {
    const double c = 3e8;
    double lambda = c / freq_Hz;
    double PL1m_dB = 20.0 * log10(4.0 * M_PI / lambda); // FSPL at 1m (dB)
    double PLd_dB = PL1m_dB + 10.0 * alpha * log10(std::max(d_m, 1e-6));
    double Pr_dBm = Ptx_dBm - PLd_dB;
    return Pr_dBm;
}
static double coordDistance(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.empty() || b.empty()) return 1e300;
    double s = 0;
    size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) {
        double d = a[i] - b[i];
        s += d * d;
    }
    return sqrt(s);
}

void MDTRouting::initialize(int stage){
    if (stage == INITSTAGE_ROUTING_PROTOCOLS){
        addressType = getSelfIPAddress().getAddressType(); // needed for handleStartOperation()
    }
    RoutingProtocolBase::initialize(stage);
    if (stage == INITSTAGE_LOCAL){
        //first token
        if(getContainingNode(this)->getId() == 7){//TODO
            attached = true;
            EV_INFO<<"node[7] hold the first token";
        }

        host = getContainingNode(this);
        routingTable.reference(this, "routingTableModule", true);
        interfaceTable.reference(this, "interfaceTableModule", true);
        networkProtocol.reference(this, "networkProtocolModule", true);

        activeRouteTimeout = par("activeRouteTimeout");

        mdtUDPPort = par("udpPort").intValue();

        //cModule *host = getContainingNode(this);
        mobilityModule = host->getSubmodule("mobility");
        mobility = check_and_cast<inet::IMobility *>(mobilityModule);

        wlan = host->getSubmodule("wlan",0);
        if(!wlan)
            EV_ERROR<<"no wlan!\n";
        else
            radioModule = wlan->getSubmodule("radio");

        if (radioModule) {
            // transmitter.power is typically a unit; careful: par() might return string or unit type
            if (radioModule->hasPar("transmitter.power")) {
                // par could be unit; doubleValue expects SI units (W)
                double ptx_W = radioModule->par("transmitter.power").doubleValue(); // W
                radioTxPower_dBm = 10.0 * log10(ptx_W) + 30.0;
            }
            else{
                EV_ERROR<<"radioModule->hasPar(transmitter.power) fail\n";
            }
            if (radioModule->hasPar("transmitter.centerFrequency"))
                radioFreq_Hz = radioModule->par("transmitter.centerFrequency").doubleValue();
            else{
                EV_ERROR<<"radioModule->hasPar(transmitter.centerFrequency) fail\n";
            }
        }
        // try to read pathloss alpha from medium params
        cModule *medium = getSimulation()->getSystemModule()->getSubmodule("radioMedium"); // may need path adjust
        if (medium && medium->hasPar("pathLoss.alpha")) {
            pathlossAlpha = medium->par("pathLoss.alpha").doubleValue();
        }
        else{
            EV_ERROR<<"medium && medium->hasPar(pathLoss.alpha) fail\n";
        }
        // compute strict threshold using 225m
        if (!std::isnan(radioTxPower_dBm)){
            kaThreshold_dBm_strict = theoreticalReceivedPower_dBm(radioTxPower_dBm, radioFreq_Hz, 225.0, pathlossAlpha);
            EV_ERROR<<"compute kaThreshold_dBm_strict local\n";
        }
        else{
            kaThreshold_dBm_strict = kaThreshold_dBm_default;
            EV_ERROR<<"get kaThreshold_dBm_strict from default\n";
        }

        // allow overriding via module parameter in NED/.ini
        if (hasPar("kaThreshold_dBm")) {
            kaThreshold_dBm_strict = par("kaThreshold_dBm").doubleValue(); // user can set exact dBm
        }
        else{
            EV_ERROR<<"no hasPar(kaThreshold_dBm)\n";
        }
        EV_ERROR << "KA strict threshold (dBm) = " << kaThreshold_dBm_strict
                << ", default fallback = " << kaThreshold_dBm_default
                << ", radioTx_dBm=" << radioTxPower_dBm << "\n";

        joinProtocol = check_and_cast<JoinProtocol*>(getModuleByPath("^.joinProtocol"));
        forwardProtocol = check_and_cast<ForwardProtocol*>(getModuleByPath("^.forwardProtocol"));
        maintenanceProtocol = check_and_cast<MaintenanceProtocol*>(getModuleByPath("^.maintenanceProtocol"));
        neighbor = check_and_cast<Neighbor*>(getModuleByPath("^.neighbor"));

        routeTimer = new cMessage("core-route-timer");
        scheduleAt(simTime() + uniform(1, activeRouteTimeout), routeTimer);

        //rpc test
        RPC_EV("INFO", "node connect, node=" << getContainingNode(this)->getId() );
    }
    else if (stage == INITSTAGE_ROUTING_PROTOCOLS) {
        networkProtocol->registerHook(0, this);
        host->subscribe(linkBrokenSignal, this);
        host->subscribe(receptionStartedSignal, this);
        host->subscribe(receptionEndedSignal, this);
    }
}
L3Address MDTRouting::getSelfIPAddress() const
{
    return routingTable->getRouterIdAsGeneric();
}
bool MDTRouting::getSelfCoordinate(std::vector<double> &coordVec) const
{
    coordVec.clear();
    if (!mobility) {
        EV_WARN << "getSelfCoordinate: mobility pointer is null\n";
        return false;
    }

    inet::Coord pos = mobility->getCurrentPosition();
    coordVec.resize(3);
    coordVec[0] = pos.x;
    coordVec[1] = pos.y;
    coordVec[2] = pos.z;
    return true;
}
void MDTRouting::socketDataArrived(UdpSocket *socket, Packet *packet)
{
    //Enter_Method("socketDataArrived");
    EV_INFO<< "socketDataArrived: pkt=" << (packet ? packet->getName() : std::string("nullptr")) << "\n";

    if (!packet) {
        EV_ERROR << "socketDataArrived: null packet\n";
        return;
    }
    // process incoming packet
    processPacket(packet);
}

void MDTRouting::socketErrorArrived(UdpSocket *socket, Indication *indication)
{
    EV_WARN << "Ignoring UDP error report " << indication->getName() << endl;
    delete indication;
}

void MDTRouting::socketClosed(UdpSocket *socket)
{
    if (operationalState == State::STOPPING_OPERATION)
        startActiveOperationExtraTimeOrFinish(par("stopOperationExtraTime"));
}
void MDTRouting::handleMessageWhenUp(cMessage *msg){
  //handle selfmessage
    EV_INFO<<"handle message\n";
    if (msg->isSelfMessage()){
        EV_INFO<<"self message\n";
        if(msg == routeTimer){
            refreshRoutetuple();
            scheduleAt(simTime() + (activeRouteTimeout / 2), routeTimer);
        }
        else if (msg->getKind() == KIND_DELAYEDSEND_NS) {
            auto timer = check_and_cast<PacketHolderForNSMessage *>(msg);
            L3Address dest = timer->getNsrqDestAddr();
            L3Address src = timer->getNsrqSrcAddr();
            Packet *p = timer->removeOwnedPacket();

            //Packet*p = maintenanceProtocol->createNeighborSetRequestTo(entry.addr,getSelfIPAddress(),L3Address());
            if(src == L3Address()){
                EV_ERROR<<"nsrq send to pneighbor: "<<dest<<"\n";
                sendPacketTo(p,dest,0);
            }
            else{
                EV_ERROR<<"new neighbor and nsrq send to the sender: "<<src<<"\n";
                sendPacketTo(p,src,0);
            }
            maintenanceProtocol->AddwaitingReplyNodes(dest,simTime());
            delete timer;
        }
        else if (msg->getKind() == KIND_DELAYEDSEND){
            auto timer = check_and_cast<PacketHolderMessage *>(msg);
            L3Address dest = timer->getDestAddr();
            Packet *p = timer->removeOwnedPacket();
            sendPacketTo(p,dest,0);
            delete timer;
        }
        else{
            throw cRuntimeError("Unknown self message");
        }
    }
    else{
        EV_INFO<<"receive message from outside\n";
        socket.processMessage(msg);
    }
}
void MDTRouting::checkIpVersionAndPacketTypeCompatibility(MDTControlPacketType packetType){
  //use Ipv6 or not,now assume using Ipv4 only,leave this method empty
}

void MDTRouting::processPacket(Packet *packet)
{
    //distribute different pachet to their specific handle method,but our protocol consist of
    //four parts:forward,join,maintenance and initialization.MDTRouting.cc responses for
    //forward part.what type of msg it should handle?
    Enter_Method("processPacket");
    if (!packet) {
        EV_ERROR << "processPacket: null packet\n";
        return;
    }
    // extract source address if available
    L3Address sourceAddr;//last jump node,not origin node
    if (auto l3ind = packet->findTag<inet::L3AddressInd>())
        sourceAddr = l3ind->getSrcAddress();
    // Pop the control chunk from the front (chunk API)
    const auto& ctrl = packet->popAtFront<MDTControlPacket>();
    if (!ctrl) {
        EV_WARN << "processPacket: no MDTControlPacket chunk at front, dropping packet\n";
        delete packet;
        return;
    }
    EV_ERROR<<"packet id:"<<packet->getId()<<"\n";

    bool weakLink = false;
    simtime_t arrival = simTime(); // time when this handler runs; you observed it's nearly equal to reception end
    std::optional<size_t> bestIdx;
    double bestScore = 1e9;
    // iterate recent deque to find match
    for (size_t i = 0; i < recentRxDeque.size(); ++i) {
        const RecentRxEntry &e = recentRxDeque[i];
        // time difference
        double dt = fabs((arrival - e.endTime).dbl());
        // primary: match txid if available -> prefer
        double score = dt;
        if (dt < matchTimeEps) {
            if (score < bestScore) {
                bestScore = score;
                bestIdx = i;
            }
        }
    }
    if (bestIdx.has_value()) {
        RecentRxEntry matched = recentRxDeque[*bestIdx];
        // remove matched entry
        recentRxDeque.erase(recentRxDeque.begin() + *bestIdx);

        EV_ERROR << "Matched packet " << packet->getName() << " to rx power " << matched.power_W
               << ")\n";
        // use matched.power_W to drive KA decision / neighbor update etc.
        double rx_dBm = 10.0 * log10(matched.power_W) + 30.0;
        double rxMinPow_dBm = rx_dBm;

        EV_ERROR << "Packet " << packet->getName()
                << " power=" << matched.power_W << " W (" << rxMinPow_dBm << " dBm)\n";

        //size_t neighCount = getPNeighborNodes().size();

        double selectedThreshold_dBm = kaThreshold_dBm_strict; //strict
        /*if (neighCount <= 3) {
            EV_ERROR << "Neighbor count <= 3 (" << neighCount << "), accepting KA regardless of power\n";
            selectedThreshold_dBm = -1e9; // effectively accept any power
        }*/

        if (rxMinPow_dBm < selectedThreshold_dBm) {
            EV_ERROR << "packet " << " considered on BOUNDARY/WEAK ("
                    << rxMinPow_dBm << " dBm < threshold " << selectedThreshold_dBm << " dBm)\n";
            // do not add to neighbor
            // pass a copy of the chunk to Neighbor (avoid ownership transfer of Packet)
            weakLink = true;
        } else {
            EV_ERROR << "packet " << " considered STRONG ("
                    << rxMinPow_dBm << " dBm >= threshold " << selectedThreshold_dBm << " dBm)\n";
            // accept
            // pass a copy of the chunk to Neighbor (avoid ownership transfer of Packet)
            weakLink = false;
        }
    }
    else {
        EV_ERROR << "No recent RX entry matched for packet " << packet->getName() << " arrival=" << arrival << "\n";
        // fallback behavior: handle packet without power info
    }


    auto packetType = (MDTControlPacketType)ctrl->getPacketType();
    switch (packetType) {
        case MDTControlPacketType::KA: { // KeepAlive
            // dynamic cast to concrete chunk type
            const auto ka = CHK(dynamicPtrCast<KeepAlive>(ctrl->dupShared()));
            if (!ka) {
                EV_ERROR << "processPacket: KA chunk dynamic cast failed\n";
                delete packet;
                return;
            }
            EV_ERROR << "processPacket: received KeepAlive src=" << sourceAddr
                    << " ts=" << simTime() << " attached=" << ka->getAttached() << "\n";

            if (weakLink) {
                EV_ERROR << "KA packet " << " considered on BOUNDARY/WEAK (";
                // do not add to neighbor
                // pass a copy of the chunk to Neighbor (avoid ownership transfer of Packet)
                if (neighbor) {
                    neighbor->processKeepAliveWeakLinkChunk(ka, sourceAddr);
                }
                else {
                    EV_WARN << "processPacket: no neighbor module to handle KeepAlive\n";
                }
            } else {
                EV_ERROR << "KA packet " << " considered STRONG (";
                // accept
                // pass a copy of the chunk to Neighbor (avoid ownership transfer of Packet)
                if (neighbor) {
                    neighbor->processKeepAliveChunk(ka, sourceAddr);
                }
                else {
                    EV_WARN << "processPacket: no neighbor module to handle KeepAlive\n";
                }
            }
            delete packet; // we popped the chunk already, safe to delete packet
            return;
        }
        case MDTControlPacketType::KAN: {
            const auto ka = CHK(dynamicPtrCast<KeepAliveNotification>(ctrl->dupShared()));
            if (!ka) {
                EV_ERROR << "processPacket: KAN chunk dynamic cast failed\n";
                delete packet;
                return;
            }
            EV_ERROR << "processPacket: received KAN src=" << sourceAddr
                    << " ts=" << simTime() << " attached=" << ka->getAttached() << "\n";

            if (weakLink) {
                EV_ERROR << "KAN packet " << " considered on BOUNDARY/WEAK ( and do not handle it";
                // do not add to neighbor
                // pass a copy of the chunk to Neighbor (avoid ownership transfer of Packet)
                if (neighbor) {
                    neighbor->processKeepAliveNotificationChunk(ka, sourceAddr,true);
                }
                else {
                    EV_WARN << "processPacket: no neighbor module to handle KAN\n";
                }
            } else {
                EV_ERROR << "KAN packet " << " considered STRONG (";
                // accept
                // pass a copy of the chunk to Neighbor (avoid ownership transfer of Packet)
                if (neighbor) {
                    neighbor->processKeepAliveNotificationChunk(ka, sourceAddr);
                }
                else {
                    EV_WARN << "processPacket: no neighbor module to handle KAN\n";
                }
            }


            delete packet; // we popped the chunk already, safe to delete packet
            return;
        }
        case MDTControlPacketType::JRQ: {
            const auto jrq = CHK(dynamicPtrCast<JoinRequest>(ctrl->dupShared()));
            if(!jrq){
                EV_ERROR << "processPacket: JRQ chunk dynamic cast failed\n";
                delete packet;
                return;
            }


            EV_ERROR << "processPacket: received JoinRequest src=" << jrq->getSourceAddr()
                    << " dest=" << jrq->getDestAddr() << " pred=" << jrq->getPredAddr() <<"\n";

            if(!attached){//attached nodes only
                EV_ERROR << "this node is not attached,can not proccess JRQ\n";
                delete packet;
                return;
            }
            if(joinProtocol){
                joinProtocol->processJoinRequestChunk(jrq);
            }
            else{
                EV_WARN << "processPacket: no join module to handle JoinRequest\n";
            }
            delete packet;
            return;
        }

        case MDTControlPacketType::JRP: {
            const auto jrp = CHK(dynamicPtrCast<JoinReply>(ctrl->dupShared()));
            if(!jrp){
                EV_ERROR << "processPacket: JRP chunk dynamic cast failed\n";
                delete packet;
                return;
            }
            EV_ERROR << "processPacket: received JoinReply src=" << jrp->getSourceAddr()
                    << " dest=" << jrp->getDestAddr() << " pred=" << jrp->getPredAddr() <<"\n";

            if(!attached && !initializing){//attached or target node only
                EV_ERROR << "this node is not attached,can not proccess JRP\n";
                delete packet;
                return;
            }
            if(joinProtocol){
                joinProtocol->processJoinReplyChunk(jrp);
            }
            else{
                EV_WARN << "processPacket: no join module to handle JoinReply\n";
            }
            delete packet;
            return;
        }

        case MDTControlPacketType::NSRQ:{
            const auto nsrq = CHK(dynamicPtrCast<NeighborSetRequest>(ctrl->dupShared()));
            if(!nsrq){
                EV_ERROR << "processPacket: NSRQ chunk dynamic cast failed\n";
                delete packet;
                return;
            }
            EV_ERROR << "processPacket: received NeighborSetRequest src=" << nsrq->getSourceAddr()
                    << " dest=" << nsrq->getDestAddr() << " pred=" << nsrq->getPredAddr() <<"\n";
            if(!attached){//attached nodes only
                EV_ERROR << "this node is not attached,can not proccess NSRQ\n";
                delete packet;
                return;
            }
            if(maintenanceProtocol){
                maintenanceProtocol->processNeighborSetRequestChunk(nsrq);
            }
            else{
                EV_WARN << "processPacket: no join module to handle NeighborSetRequest\n";
            }
            delete packet;
            return;
        }
        case MDTControlPacketType::NSRP:{
            const auto nsrp = CHK(dynamicPtrCast<NeighborSetReply>(ctrl->dupShared()));
            if(!nsrp){
                EV_ERROR << "processPacket: NSRP chunk dynamic cast failed\n";
                delete packet;
                return;
            }
            EV_ERROR << "processPacket: received NeighborSetReply src=" << nsrp->getSourceAddr()
                    << " dest=" << nsrp->getDestAddr() << " pred=" << nsrp->getPredAddr() <<"\n";
            if(!attached && !initializing){//attached nodes or target node only
                EV_ERROR << "this node is not attached,can not proccess NSRP\n";
                delete packet;
                return;
            }
            if(maintenanceProtocol){
                maintenanceProtocol->processNeighborSetReplyChunk(nsrp);
            }
            else{
                EV_WARN << "processPacket: no join module to handle NeighborSetReply\n";
            }
            delete packet;
            return;
        }
        case MDTControlPacketType::NSN: {
            const auto nsn = CHK(dynamicPtrCast<NeighborSetNotification>(ctrl->dupShared()));
            if(!nsn){
                EV_ERROR << "processPacket: NSN chunk dynamic cast failed\n";
                delete packet;
                return;
            }
            EV_ERROR << "processPacket: received NeighborSetNotification src=" << nsn->getSourceAddr()
                        << " dest=" << nsn->getDestAddr() << " pred=" << nsn->getPredAddr()<<"\n";
            if(maintenanceProtocol){
                maintenanceProtocol->processNeighborSetNotificationChunk(nsn);
            }
            else{
                EV_WARN << "processPacket: no join module to handle NeighborSetNotification\n";
            }
            delete packet;
            return;
        }
        case MDTControlPacketType::CDQ:{
            const auto cd = CHK(dynamicPtrCast<CoordsDiscoverRequest>(ctrl->dupShared()));
            if(!cd){
                EV_ERROR << "processPacket: CDQ chunk dynamic cast failed\n";
                delete packet;
                return;
            }
            EV_ERROR << "processPacket: received CoordsDiscoverRequest src=" << cd->getSourceAddr()
                         << " tar=" << cd->getTargetAddr()<<"\n";
            processCoordsDiscoverRequestChunk(cd);
            delete packet;
            return;
        }
        case MDTControlPacketType::CDP:{
            const auto cd = CHK(dynamicPtrCast<CoordsDiscoverReply>(ctrl->dupShared()));
            if(!cd){
                EV_ERROR << "processPacket: CDP chunk dynamic cast failed\n";
                delete packet;
                return;
            }
            EV_ERROR << "processPacket: received CoordsDiscoverReply src=" << cd->getSourceAddr()
                        << " dest=" << cd->getDestAddr() << " tar=" << cd->getTargetAddr()<<"\n";
            processCoordsDiscoverReplyChunk(cd);
            delete packet;
            return;
        }
        // add other cases (IT, LN, FN, etc.) similar way:
        default:
            EV_WARN << "processPacket: unknown MDT control packet type=" << (int)packetType << ", dropping\n";
            delete packet;
            return;
    }
}
void MDTRouting::sendPacketTo(Packet *packet, const L3Address& destAddr, int destPort,double delay)
{
    Enter_Method("sendPacketTo");

    if (!packet) {
        EV_ERROR << "sendPacketTo: called with null packet\n";
        return;
    }

    take(packet);

    int port = (destPort == 0 ? (int)mdtUDPPort : destPort);
    EV_INFO << "sendPacketTo: sending to " << destAddr << " port=" << port
            << " originalPacket='" << packet->getName() << "'\n";

    auto setStandardTags = [&](Packet *p) {
        int interfaceId = CHK(interfaceTable->findInterfaceByName(par("interface")))->getInterfaceId();
        p->addTagIfAbsent<InterfaceReq>()->setInterfaceId(interfaceId);
        //???
        p->addTagIfAbsent<HopLimitReq>()->setHopLimit(64);
        //
        p->addTagIfAbsent<L3AddressReq>()->setDestAddress(destAddr);
        p->addTagIfAbsent<L4PortReq>()->setDestPort(port);
    };

    EV_INFO << "sendPacketTo: FALLBACK sending original packet with tags\n";
    setStandardTags(packet);
    EV_ERROR<<"packet id:"<<packet->getId()<<"\n";
    //if (auto l4 = packet->findTag<L4PortReq>()) EV_INFO << "  L4 destPort tag = " << l4->getDestPort() << "\n";
    //if (auto l3 = packet->findTag<L3AddressReq>()) EV_INFO << "  L3 dest = " << l3->getDestAddress() << "\n";

    //socket.sendTo(packet, destAddr, port);
    if (delay == 0) {
            socket.sendTo(packet, destAddr, port);
            EV_INFO << "sendPacketTo: sent immediately\n";
    } else {
        auto *timer = new PacketHolderMessage("mdt-delay-send", KIND_DELAYEDSEND);
        timer->setOwnedPacket(packet);
        timer->setDestAddr(destAddr);
        scheduleAfter(delay, timer);
        EV_INFO << "sendPacketTo: scheduled for delayed send after " << delay << "\n";
    }
}

void MDTRouting::sendJoinRequest(const L3Address& destAddr){
    std::vector<double> selfcoords;
    getSelfCoordinate(selfcoords);
    Packet*p = joinProtocol->createJoinRequestTo(destAddr,getSelfIPAddress(),L3Address(),selfcoords);
    //L3Address dest = L3Address(inet::Ipv4Address::ALLONES_ADDRESS);
    EV_INFO<<"send first jrq\n";
    sendPacketTo(p,destAddr,0);
}
void MDTRouting::sendCoordsDiscoverRequest(const L3Address& targetAddr){
    L3Address broadcastAddr = L3Address(inet::Ipv4Address::ALLONES_ADDRESS);
    auto cd = makeShared<CoordsDiscoverRequest>();

    cd->setPacketType(MDTControlPacketType::CDQ);
    cd->setSourceAddr(getSelfAddress());

    cd->setTargetAddr(targetAddr);
    cd->setHopcount(0);
    int bytes = 0;
    bytes += 1;                // packet type
    bytes += 4;                // src IPv4
    bytes += 4;                 //pred
    bytes += 4;                 //hopcount
    cd->setChunkLength(inet::B(bytes)); // set chunk length explicitly

    Packet* p = new Packet("mdt-coordsdiscoverrequest",cd);
    sendPacketTo(p,broadcastAddr,0);
}
void MDTRouting::processCoordsDiscoverRequestChunk(const Ptr<CoordsDiscoverRequest>& cd){
    L3Address target = cd->getTargetAddr();
    std::map<L3Address, NeighborEntry> known = getKnownNodes();
    for(auto& [key,value]:known){
        if(key == target){
            sendCoordsDiscoverReply(cd->getSourceAddr(),target,value);
            EV_INFO<<"cdrq find the target node: "<<target<<"\n";
        }
    }
}
void MDTRouting::sendCoordsDiscoverReply(const L3Address& destAddr,const L3Address& targetAddr,const NeighborEntry& entry){
    auto cd = makeShared<CoordsDiscoverReply>();

    cd->setPacketType(MDTControlPacketType::CDP);
    cd->setSourceAddr(getSelfAddress());
    cd->setDestAddr(destAddr);
    cd->setTargetAddr(targetAddr);
    cd->setNeighbor(entry);

    int bytes = 0;
    bytes += 1;                // packet type
    bytes += 4;                // src IPv4
    bytes += 4;                 //dest
    bytes += 4;                 //tar

    bytes += 4;            // addr IPv4
    bytes += 1;            // attached flag (1 byte)
    int ndim = entry.dim;
    if (ndim < 0) ndim = 0;
    bytes += 4;            // dim (uint32)
    bytes += 8 * ndim;     // coords: ndim * 8 bytes (double as uint64)
    bytes += 8;            // timeout (uint64 ms)

    cd->setChunkLength(inet::B(bytes)); // set chunk length explicitly

    Packet* p = new Packet("mdt-coordsdiscoverreply",cd);
    sendPacketTo(p,destAddr,0);
}
void MDTRouting::processCoordsDiscoverReplyChunk(const Ptr<CoordsDiscoverReply>& cd){
    if(cd->getDestAddr() == getSelfAddress()){
        L3Address target = cd->getTargetAddr();
        EV_INFO<<"receive CoordsDiscoverReply target = "<<target<<"\n";
        NeighborEntry entry = cd->getNeighbor();
        std::vector<NeighborEntry> entrys;
        entrys.push_back(entry);
        addOrUpdateKnownNode(entrys);
        completeRouteDiscovery(target);
    }
}
void MDTRouting::sendNSRQ(const std::vector<NeighborEntry>& entrys,const L3Address& srcAddr){//srcAddr is who sent the entrys
    Enter_Method_Silent("sendNSRQ");
    if(maintenanceProtocol->getStage() == MaintenanceState::INACTIVE){
        //change stage
        EV_INFO<<"maintenanceProtocol change stage to JOINSTAGE\n";
        maintenanceProtocol->setStage(MaintenanceState::JOINSTAGE);
        //srcAddr has been inquired
        maintenanceProtocol->AddinquiredNodes(srcAddr);
    }

    for(auto entry:entrys){
        //don't repeat sending
        if(maintenanceProtocol->findNodeinInquiredNodes(entry.addr)){
            continue;
        }

        Packet*p = maintenanceProtocol->createNeighborSetRequestTo(entry.addr,getSelfIPAddress(),getSelfIPAddress());//TODO attention: choose which addr as pred?
        //auto *timer = new PacketHolderForNSMessage("mdt-sendnsrq-jitter", KIND_DELAYEDSEND_NS);
        //double delay = 0;//uniform(0,0.1);       //TODO
        //timer->setOwnedPacket(p);
        //timer->setNsrqDestAddr(entry.addr);
        EV_ERROR<<"send to new neighbor: "<<entry.addr<<"\n";
        maintenanceProtocol->AddwaitingReplyNodes(entry.addr,simTime());
        if(neighbor->isPhysicalNeighbor(entry.addr)){
            //timer->setNsrqSrcAddr(L3Address());
            //scheduleAfter(delay, timer);
            sendPacketTo(p,entry.addr,0);
        }
        else{//MDT forward
            //timer->setNsrqSrcAddr(srcAddr);
            //scheduleAfter(delay, timer);
            if(maintenanceProtocol->getStage() == MaintenanceState::JOINSTAGE){
                if(srcAddr == L3Address()){
                    EV_ERROR<<"src is empty send NSRQ error\n";
                    sendPacketTo(p,entry.addr,0);
                }
                else{
                    sendPacketTo(p,srcAddr,0);
                }
            }
            else{
                sendPacketTo(p,entry.addr,0);
            }
        }
        maintenanceProtocol->AddinquiredNodes(entry.addr);
    }


    //else if(maintenanceProtocol->getStage() == MaintenanceState::JOINSTAGE){

   // }
}
// core methods:
//when this node is source node

INetfilter::IHook::Result MDTRouting::ensureRouteForDatagram(Packet *datagram,int state)
{
    Enter_Method("ensureRouteForDatagram");

    const auto& networkHeader = getNetworkProtocolHeader(datagram);
    const L3Address& destAddr = networkHeader->getDestinationAddress();
    const L3Address& srcAddr = networkHeader->getSourceAddress();//source addr not the last previous one

    EV_ERROR<<" ensureRouteForDatagram called, packet=" << datagram <<" self ip addr = "<<getSelfAddress()<<"\n";

    if (!datagram) {
        EV_ERROR << "ensureRouteForDatagram: packet is nullptr\n";
        return INetfilter::IHook::ACCEPT;
    }

    if (destAddr.isUnspecified()) {
        EV_WARN << "ensureRouteForDatagram: cannot determine dest address -> ACCEPT\n";
        // debug helper to print tags/chunks may be helpful here
        // printPacketDebugInfo(datagram);
        return INetfilter::IHook::ACCEPT;
    }

    EV_ERROR << "ensureRouteForDatagram: src=" << srcAddr
             << ", dest=" << destAddr << "\n";

    // basic checks replicating your original logic
    if (destAddr.isBroadcast() || destAddr.isMulticast()) {
        EV_ERROR << "Destination is broadcast/multicast -> ACCEPT\n";
        return INetfilter::IHook::ACCEPT;
    }
    if (routingTable->isLocalAddress(destAddr)) {
        EV_ERROR << "ensureRouteForDatagram: dest is local -> ACCEPT\n";
        return INetfilter::IHook::ACCEPT;
    }

    IRoute *route = nullptr;
    try {
        route = routingTable->findBestMatchingRoute(destAddr);
        if(route){
            EV_ERROR <<"findBestMatchingRoute :dest :"<< route->getDestinationAsGeneric()
                   << " NextHop: " << route->getNextHopAsGeneric() << endl;
        }
    } catch (...) {
        route = nullptr;
    }


    std::string pkname = datagram->getName();
    std::string lowerName = pkname;
    std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    if (lowerName.find("icmp") != std::string::npos) {//mdt packet or icmp error
        return INetfilter::IHook::ACCEPT;
    }
    bool is_mdt = false;
    if(lowerName.find("mdt") != std::string::npos){
        is_mdt = true;
    }
    /*auto ctrl = datagram->peekAtFront<MDTControlPacket>();
    if (ctrl) {
        return INetfilter::IHook::ACCEPT;
    }*/

    //add relay tag
    Ptr<MdtRelayChunk> relayChunk;
    L3Address predaddr = L3Address();

    if (state == 1) {//forward
        //Ptr<const Chunk> raw = datagram->popAtBack(B(9));
        //Ptr<const MdtRelayChunk> c = datagram->popAtBack<MdtRelayChunk>(B(9));
        //printLast9Bytes(datagram);

        EV_ERROR <<"cheak mdtrelay chunk\n";
        EV_ERROR << "RECV: packet ptr=" << datagram
                << " dataLen=" << datagram->getDataLength() << "\n";
        //EV_ERROR << "RECV: back class=" << (datagram->peekAtBack() ? datagram->peekAtBack()->getClassName() : std::string("NULL")) << "\n";

        auto c = datagram->peekAtBack<MdtRelayChunk>(B(9));
        EV_ERROR << "c name:"<< c->getClassName() << "\n";
        if (c) {
            predaddr = c->getPredAddr();

            relayChunk = makeShared<MdtRelayChunk>();
            relayChunk->setHasRelay(c->getHasRelay());
            relayChunk->setRelayAddr(c->getRelayAddr());
            relayChunk->setPredAddr(getSelfIPAddress());
            relayChunk->setChunkLength(B(9));

            datagram->popAtBack<MdtRelayChunk>(B(9));

            //debug
            std::string pkname = datagram->getName();
            std::string lowerName = pkname;
            std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            EV_ERROR <<"find predaddr\n";
            //RPC_EV("INFO", "note predaddr = "<< predaddr << "at simtime: "<<simTime()<< "at node: "<< getSelfIPAddress());
            // not content "mdt"
            if (lowerName.find("mdt") == std::string::npos) {//app packet
                //RPC_EV("INFO", "predaddr = "<< predaddr << "at node: "<< getSelfIPAddress());
            }
            if(predaddr.isUnspecified())
                RPC_EV("INFO", "predaddr isUnspecified!"<<"dest = "<< destAddr << "at node: "<< getSelfIPAddress());
        } else {
            //datagram->insertAtBack(raw);
            EV_ERROR << "fail to get MdtRelayChunk\n";
        }
    }
    else {//local out
        relayChunk = makeShared<MdtRelayChunk>();
        relayChunk->setHasRelay(false);
        relayChunk->setRelayAddr(L3Address());
        relayChunk->setPredAddr(getSelfIPAddress());
        relayChunk->setChunkLength(B(9));

        EV_ERROR <<"packet name: "<< datagram->getName();
        EV_ERROR << " ensureRouteForDatagram: added MdtRelayTag to packet\n";
    }

    simtime_t time = simTime();
    if (routingTable) {
        auto P = getPNeighborNodes();
        if(P.empty()){//no pneighbor
            return INetfilter::IHook::DROP;
        }
        if(P.find(destAddr) != P.end()){//rule one pneighbor
            //note that predaddr
            L3Address nxt = L3Address();
            for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
                IRoute* rt = routingTable->getRoute(i);
                if (!rt) continue;
                if(rt->getDestinationAsGeneric() == destAddr){//rule one pneighbor
                    nxt = rt->getNextHopAsGeneric();
                }
            }
            std::map<L3Address, NeighborEntry> weaklinkneighbors = getweaklinkNeighbors();
            if(weaklinkneighbors.count(nxt) != 0){//greedy for weaklink
                auto K = getKnownNodes();
                bool flag = true;
                for(auto &[v,value] : K){
                    if(v == nxt){
                        if(value.timeout <= simTime()){
                            flag = false;
                        }
                    }
                }
                std::vector<double> destCoord;
                std::vector<double> nxtCoord;
                if(getKnownNodeCoordinate(destAddr,destCoord) && getKnownNodeCoordinate(nxt,nxtCoord)){
                    double bestDist = coordDistance(nxtCoord, destCoord);
                    L3Address best = L3Address();
                    for (auto &[v,value] : P){
                        std::vector<double> vcoord;
                        if(value.timeout <= simTime() || v == getSelfIPAddress() || v == predaddr || !value.attachedToDT)continue;
                        if(weaklinkneighbors.count(v))continue;
                        if(!getDTNeighborNodeCoordinate(v, vcoord))continue;
                        double d = coordDistance(vcoord, destCoord);
                        if(d < bestDist){
                            best = v;
                            bestDist = d;
                        }
                    }
                    if (!best.isUnspecified() && best != getSelfIPAddress() && flag){
                        nxt = best;
                        Routetuple tuple("physics",getSelfAddress(),L3Address(),nxt,destAddr,-1,-1,simTime() + activeRouteTimeout);
                        ensureRoute(destAddr,nxt,tuple);
                    }
                }
            }
            if(!predaddr.isUnspecified() && nxt == predaddr){//loop
                //not likely
            }

            if (relayChunk){
                //datagram->insertAtBack(relayChunk);
                //auto ch = datagram->peekAtBack<MdtRelayChunk>(B(9));

                //printLast9Bytes(datagram);
                EV_ERROR << "back name:"<< datagram << "\n";
                //addnewchunk(datagram,relayChunk);
                auto ipv4Hdr = datagram->removeAtFront<Ipv4Header>();
                datagram->insertAtBack(relayChunk);
                ipv4Hdr->setTotalLengthField(ipv4Hdr->getChunkLength() +
                                             datagram->getDataLength());   // relayChunk
                ipv4Hdr->setCrc(0);
                ipv4Hdr->updateCrc();
                datagram->insertAtFront(ipv4Hdr);

                EV_ERROR<<"addnewchunk,datagram = "<<datagram<<"\n";
            }
            return INetfilter::IHook::ACCEPT;
        }
        L3Address dest = destAddr;
        /*for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
            IRoute* rt = routingTable->getRoute(i);
            if (!rt) continue;
            if(rt->getNextHopAsGeneric() == dest){//rule one pneighbor
                return INetfilter::IHook::ACCEPT;
            }
        }*/

        std::vector<double> destCoord;
        if(!getKnownNodeCoordinate(dest,destCoord)){//unknown node
            delayDatagram(datagram);
            if(!hasOngoingRouteDiscovery(dest)){
                sendCoordsDiscoverRequest(dest);
            }
            if (relayChunk){
                //datagram->insertAtBack(relayChunk);
                //auto ch = datagram->peekAtBack<MdtRelayChunk>(B(9));
                //printLast9Bytes(datagram);
                //EV_ERROR << "back name:"<< ch->getClassName() << "\n";
                EV_ERROR << "back name:"<< datagram << "\n";
                //addnewchunk(datagram,relayChunk);
                auto ipv4Hdr = datagram->removeAtFront<Ipv4Header>();
                datagram->insertAtBack(relayChunk);
                ipv4Hdr->setTotalLengthField(ipv4Hdr->getChunkLength() +
                                             datagram->getDataLength());   // relayChunk
                ipv4Hdr->setCrc(0);
                ipv4Hdr->updateCrc();
                datagram->insertAtFront(ipv4Hdr);
                EV_ERROR<<"addnewchunk,datagram = "<<datagram;
            }
            return INetfilter::IHook::QUEUE;
        }

        L3Address relaydest = L3Address();
        //auto tag = datagram->addTagIfAbsent<MdtRelayTag>();
        if (relayChunk->getHasRelay()) {
            relaydest = relayChunk->getRelayAddr();
        }
        if (relaydest == getSelfIPAddress()) {
            relaydest = L3Address();
            relayChunk->setHasRelay(false);
        }

        //auto P = getPNeighborNodes();
        L3Address self = getSelfIPAddress();
        L3Address nxthop = L3Address();
        if(!relaydest.isUnspecified()){//MDT
            //short cut first
            if(forwardProtocol->shortcutforMDTrelay(dest,relaydest,predaddr)){
                RPC_EV("INFO", "find new relay as short cut: "<<relaydest<<"at: "<<time);
                EV_ERROR<<"find new relay as short cut: "<<relaydest<<"\n";
                relayChunk->setRelayAddr(relaydest);
            }
            if(forwardProtocol->findnexthop(relaydest, nxthop)){//greedy

            }
            else{//MDT
                IRoute* r = routingTable->findBestMatchingRoute(nxthop);
                if(r != nullptr){
                    nxthop = r->getNextHopAsGeneric();
                }
                else{
                    L3Address relaythop = L3Address();
                    forwardProtocol->findnexthop(nxthop, relaythop);//TODO now we assume the result is greedy
                    nxthop = relaythop;
                }
                /*if(!nxthop.isUnspecified() && nxthop != self  && P.find(nxthop) != P.end()){
                    Routetuple tuple("physics",getSelfAddress(),L3Address(),nxthop,dest,-1,-1,simTime() + activeRouteTimeout);
                    ensureRoute(dest,nxthop,tuple);
                }
                return INetfilter::IHook::ACCEPT;*/
            }
            if(!nxthop.isUnspecified() && nxthop != self  && P.find(nxthop) != P.end()){
                if(nxthop == predaddr){//loop!
                    //RPC_EV("INFO", "MDT find loop, nextHop=" << predaddr<<"dest = "<< dest << "at node: "<< self);
                    for (int i = 0; i < routingTable->getNumRoutes(); ++i) {//delete first
                        IRoute *rt = routingTable->getRoute(i);
                        if(rt->getNextHopAsGeneric() == predaddr && rt->getDestinationAsGeneric() == relaydest){
                            routingTable->deleteRoute(rt);
                            break;
                        }
                    }

                    L3Address ban = nxthop;
                    nxthop = L3Address();
                    if(forwardProtocol->findnexthop(relaydest, nxthop,ban)){

                    }
                    else{
                        IRoute* r = routingTable->findBestMatchingRoute(nxthop);
                        if(r != nullptr){
                            nxthop = r->getNextHopAsGeneric();
                        }
                        else{
                            L3Address relaythop = L3Address();
                            forwardProtocol->findnexthop(nxthop, relaythop,ban);//TODO now we assume the result is greedy
                            nxthop = relaythop;
                        }
                    }
                    if(!nxthop.isUnspecified()){
                        //RPC_EV("INFO", "find new nextHop=" << nxthop);
                    }
                    else{
                        RPC_EV("INFO", "mdt can not find new nextHop at: "<<time<<" relaydest = "<<relaydest<<" predaddr: "<<predaddr<<" chunk pred: "<<relayChunk->getPredAddr());
                    }
                }
                if(!nxthop.isUnspecified()){
                    Routetuple tuple("physics",getSelfAddress(),L3Address(),nxthop,dest,-1,-1,simTime() + activeRouteTimeout);
                    ensureRoute(dest,nxthop,tuple);
                }
            }
            if (relayChunk){
                //datagram->insertAtBack(relayChunk);
                //auto ch = datagram->peekAtBack<MdtRelayChunk>(B(9));
                //printLast9Bytes(datagram);
                //EV_ERROR << "back name:"<< ch->getClassName() << "\n";
                EV_ERROR << "back name:"<< datagram << "\n";
                //addnewchunk(datagram,relayChunk);
                auto ipv4Hdr = datagram->removeAtFront<Ipv4Header>();
                datagram->insertAtBack(relayChunk);
                ipv4Hdr->setTotalLengthField(ipv4Hdr->getChunkLength() +
                                             datagram->getDataLength());   // relayChunk
                ipv4Hdr->setCrc(0);
                ipv4Hdr->updateCrc();
                datagram->insertAtFront(ipv4Hdr);
                EV_ERROR<<"addnewchunk,datagram = "<<datagram;
            }
            return INetfilter::IHook::ACCEPT;
        }
        else{//greedy
            if(forwardProtocol->findnexthop(dest, nxthop)){//greedy
                /*if(!nxthop.isUnspecified() && nxthop != self  && P.find(nxthop) != P.end()){
                    Routetuple tuple("physics",getSelfAddress(),L3Address(),nxthop,dest,-1,-1,simTime() + activeRouteTimeout);
                    ensureRoute(dest,nxthop,tuple);
                }
                return INetfilter::IHook::ACCEPT;*/
            }
            else{//MDT
                //change tag
                L3Address relay = nxthop;
                relayChunk->setHasRelay(true);
                relayChunk->setRelayAddr(relay);
                //find nxt
                IRoute* r = routingTable->findBestMatchingRoute(relay);
                if(r != nullptr){
                    nxthop = r->getNextHopAsGeneric();
                }
                else{
                    L3Address relaythop = L3Address();
                    forwardProtocol->findnexthop(relay, relaythop);//TODO now we assume the result is greedy
                    nxthop = relaythop;
                }
                /*if(!nxthop.isUnspecified() && nxthop != self  && P.find(nxthop) != P.end()){
                    Routetuple tuple("physics",getSelfAddress(),L3Address(),nxthop,dest,-1,-1,simTime() + activeRouteTimeout);
                    ensureRoute(dest,nxthop,tuple);
                }
                return INetfilter::IHook::ACCEPT;*/
            }
            if(!nxthop.isUnspecified() && nxthop != self  && P.find(nxthop) != P.end()){
                if(nxthop == predaddr){//loop!
                    if(!is_mdt)RPC_EV("INFO", "greedy find loop, nextHop=" << predaddr<<"dest = "<< dest << "at time: "<< time);
                    for (int i = 0; i < routingTable->getNumRoutes(); ++i) {//delete first
                        IRoute *rt = routingTable->getRoute(i);
                        if(rt->getNextHopAsGeneric() == predaddr && rt->getDestinationAsGeneric() == destAddr){
                            routingTable->deleteRoute(rt);
                            if(!is_mdt)RPC_EV("INFO", "delete route for loop");
                            break;
                        }
                    }

                    L3Address ban = nxthop;
                    nxthop = L3Address();
                    if(forwardProtocol->findnexthop(dest, nxthop,ban)){
                        if(!is_mdt)RPC_EV("INFO","greedy can find nxthop: "<< nxthop<<" at time: "<< time);
                    }
                    else{
                        //change tag
                        L3Address relay = nxthop;
                        relayChunk->setHasRelay(true);
                        relayChunk->setRelayAddr(relay);
                        if(!is_mdt)RPC_EV("INFO","mdt can find nxthop: "<< nxthop<<" at time: "<< time);
                        //find nxt
                        IRoute* r = routingTable->findBestMatchingRoute(relay);
                        if(r != nullptr){
                            nxthop = r->getNextHopAsGeneric();
                        }
                        else{
                            L3Address relaythop = L3Address();
                            forwardProtocol->findnexthop(relay, relaythop,ban);//TODO now we assume the result is greedy
                            nxthop = relaythop;
                        }
                    }
                    if(!nxthop.isUnspecified()){
                        //RPC_EV("INFO", "find new nextHop=" << nxthop);
                    }
                    else{
                        if(!is_mdt)RPC_EV("INFO", "greedy can not find new nextHop at: "<<time<<" dest = "<<dest <<" predaddr: "<<predaddr<<" chunk pred: "<<relayChunk->getPredAddr());
                    }
                }
                if(!nxthop.isUnspecified()){
                    Routetuple tuple("physics",getSelfAddress(),L3Address(),nxthop,dest,-1,-1,simTime() + activeRouteTimeout);
                    ensureRoute(dest,nxthop,tuple);
                }
            }
            if (relayChunk){
                //datagram->insertAtBack(relayChunk);
                //auto ch = datagram->peekAtBack<MdtRelayChunk>(B(9));
                //printLast9Bytes(datagram);
                //EV_ERROR << "back name:"<< ch->getClassName() << "\n";
                EV_ERROR << "back name:"<< datagram << "\n";
                //addnewchunk(datagram,relayChunk);
                auto ipv4Hdr = datagram->removeAtFront<Ipv4Header>();
                datagram->insertAtBack(relayChunk);
                ipv4Hdr->setTotalLengthField(ipv4Hdr->getChunkLength() +
                                             datagram->getDataLength());   // relayChunk
                ipv4Hdr->setCrc(0);
                ipv4Hdr->updateCrc();
                datagram->insertAtFront(ipv4Hdr);
                EV_ERROR<<"addnewchunk,datagram = "<<datagram;
            }
            return INetfilter::IHook::ACCEPT;
        }
    }

    else{
        EV_ERROR<<"deleteexpiredneighbors: no routetable\n";
    }
    if (relayChunk){
        //datagram->insertAtBack(relayChunk);
        //auto ch = datagram->peekAtBack<MdtRelayChunk>(B(9));
        //printLast9Bytes(datagram);
        //EV_ERROR << "back name:"<< ch->getClassName() << "\n";
        EV_ERROR << "back name:"<< datagram << "\n";
        //addnewchunk(datagram,relayChunk);
        auto ipv4Hdr = datagram->removeAtFront<Ipv4Header>();
        datagram->insertAtBack(relayChunk);
        ipv4Hdr->setTotalLengthField(ipv4Hdr->getChunkLength() +
                                     datagram->getDataLength());   // relayChunk
        ipv4Hdr->setCrc(0);
        ipv4Hdr->updateCrc();
        datagram->insertAtFront(ipv4Hdr);
        EV_ERROR<<"addnewchunk,datagram = "<<datagram;
    }
    return INetfilter::IHook::ACCEPT;


}

//INetfilter::IHook::Result MDTRouting::datagramForwardHook(Packet *datagram){
//}
void MDTRouting::delayDatagram(Packet *datagram)
{
    const auto& networkHeader = getNetworkProtocolHeader(datagram);
    EV_INFO << "Queuing datagram, source " << networkHeader->getSourceAddress() << ", destination " << networkHeader->getDestinationAddress() << endl;
    const L3Address& target = networkHeader->getDestinationAddress();
    targetAddressToDelayedPackets.insert(std::pair<L3Address, Packet *>(target, datagram));
}
bool MDTRouting::hasOngoingRouteDiscovery(const L3Address& target) {
    auto lt = targetAddressToDelayedPackets.lower_bound(target);
    auto ut = targetAddressToDelayedPackets.upper_bound(target);
    return lt != ut;
}
void MDTRouting::completeRouteDiscovery(const L3Address& target){
    EV_INFO << "Completing route discovery, originator " << getSelfIPAddress() << ", target " << target << endl;
    //ASSERT(hasOngoingRouteDiscovery(target));

    auto lt = targetAddressToDelayedPackets.lower_bound(target);
    auto ut = targetAddressToDelayedPackets.upper_bound(target);

    // reinject the delayed datagrams
    for (auto it = lt; it != ut; it++) {
        Packet *datagram = it->second;
        const auto& networkHeader = getNetworkProtocolHeader(datagram);
        EV_ERROR << "Sending queued datagram: source " << networkHeader->getSourceAddress() << ", destination " << networkHeader->getDestinationAddress() << endl;
        networkProtocol->reinjectQueuedDatagram(datagram);
    }

    // clear the multimap
    targetAddressToDelayedPackets.erase(lt, ut);

    // we have a route for the destination, thus we must cancel the WaitForRREPTimer events
    //auto waitRREPIter = waitForRREPTimers.find(target);
    //ASSERT(waitRREPIter != waitForRREPTimers.end());
    //cancelAndDelete(waitRREPIter->second);
    //waitForRREPTimers.erase(waitRREPIter);
}

bool MDTRouting::updateValidRouteLifeTime(const L3Address& destAddr, simtime_t lifetime){//not use
    //TODO,example:
    //IRoute *route = routingTable->findBestMatchingRoute(destAddr);
    //if (route && route->getSource() == this) {
    //    AodvRouteData *routeData = check_and_cast<AodvRouteData *>(route->getProtocolData());
    //    if (routeData->isActive()) {
    //        simtime_t newLifeTime = std::max(routeData->getLifeTime(), lifetime);
    //        EV_DETAIL << "Updating " << route << " lifetime to " << newLifeTime << endl;
    //        routeData->setLifeTime(newLifeTime);
    //        return true;
    //    }
    //}
    //return false;
}
MDTRouting::MDTRouting()
{
}
void MDTRouting::sendmessage(){//add Parameters and for different messages there are different
                          //send methods for them
//call method:sendMDTPacket to send message
}
//accordingly there is a method to create the message:
void MDTRouting::createmessage(){
  //add Parameters and maybe several similar methods
}
void MDTRouting::handlemessage(){

  //add Parameters and maybe several similar methods
  //may updateRoutingTable
}
void MDTRouting::updateRoutingTable(IRoute *route,const L3Address& nextHop,const Routetuple& routetuple){
    EV_DETAIL << "Updating existing route: " << route << endl;

    route->setNextHop(nextHop);
    //route->setMetric();

    MDTData *routingData = check_and_cast<MDTData *>(route->getProtocolData());
    ASSERT(routingData != nullptr);

    routingData->setroutetuple(routetuple);

    EV_DETAIL << "Route updated: " << route << endl;

    //scheduleExpungeRoutes();//refresh timer state
}
void MDTRouting::refreshRoutetuple(){
    std::vector<IRoute *> toDelete;
    for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
        IRoute *rt = routingTable->getRoute(i);
        if (!rt) continue;
        auto *mdt = dynamic_cast<MDTData *>(rt->getProtocolData());
        if (mdt && mdt->getroutetuple().timeout < simTime())
            toDelete.push_back(rt);
    }
    for (IRoute *rt : toDelete){
        EV_INFO<<"delete expired route to:"<<rt->getDestinationAsGeneric()<<"\n";
        routingTable->deleteRoute(rt);
    }
    logRoutingTable(this,0,routingTable);

}
IRoute *MDTRouting::createRoute(const L3Address& destAddr,const L3Address& nextHop,const Routetuple& routetuple){
    // create a new route
    IRoute *newRoute = routingTable->createRoute();

    // adding generic fields
    newRoute->setDestination(destAddr);
    newRoute->setNextHop(nextHop);
    newRoute->setPrefixLength(addressType->getMaxPrefixLength()); // TODO
    //newRoute->setMetric();
    NetworkInterface *ifEntry = interfaceTable->findInterfaceByName(par("interface")); // TODO IMPLEMENT: multiple interfaces
    if (ifEntry)
        newRoute->setInterface(ifEntry);
    newRoute->setSourceType(IRoute::MDTRouting);
    newRoute->setSource(this);

    // A route towards a destination that has a routing table entry
    // that is marked as valid.  Only active routes can be used to
    // forward data packets.

    // adding protocol-specific fields
    MDTData *newProtocolData = new MDTData();
    newProtocolData->setroutetuple(routetuple);

    newRoute->setProtocolData(newProtocolData);

    EV_DETAIL << "Adding new route " << newRoute << endl;

    routingTable->addRoute(newRoute);

    return newRoute;
}
IRoute *MDTRouting::ensureRoute(const L3Address& dest,const L3Address& nextHop,const Routetuple& rt){
    if(!neighbor->isPhysicalNeighbor(nextHop)){
        return nullptr;
    }
    if(dest == getSelfIPAddress()){
        Routetuple tuple("physics",dest,L3Address(),L3Address(),dest,-1,-1,simTime() + activeRouteTimeout);
        IRoute * r = createRoute(dest, dest, tuple);
        return r;
    }
    for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
        IRoute *route = routingTable->getRoute(i);
        if (route->getSourceType() == IRoute::MDTRouting &&
                route->getDestinationAsGeneric() == dest) {
                updateRoutingTable(route, nextHop, rt);
                logRoutingTable(this,0,routingTable);
                return route;
            }
    }
    IRoute * r = createRoute(dest, nextHop, rt);
    logRoutingTable(this,0,routingTable);
    return r;
}
void MDTRouting::findroutefornewdt(const std::vector<L3Address> newdtset){
    std::set<L3Address> routeDestinations;
    if (routingTable) {
        for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
            IRoute* rt = routingTable->getRoute(i);
            if (!rt) continue;

            MDTData* routeDestData = dynamic_cast<MDTData*>(rt->getProtocolData());
            if (routeDestData) {
                Routetuple tuple = routeDestData->getroutetuple();
                routeDestinations.insert(tuple.dest);
            }
        }
        std::vector<NeighborEntry> result;
        for (const auto& key : newdtset) {
            if (routeDestinations.find(key) == routeDestinations.end()) {
                std::map<L3Address, NeighborEntry> dtNeighbors = neighbor->getDTNeighbors();
                auto it = dtNeighbors.find(key);
                if (it != dtNeighbors.end()) {
                    result.push_back(it->second);
                }
            }
        }
        if(!result.empty()){
            EV_ERROR<<"new dt neighbor and use maintenanceProtocol to discover route\n";
            maintenanceProtocol->findroutefornewdt(result);
        }
    }
    else{
        EV_ERROR<<"no routetable and can not find new route\n";
        return;
    }
}

void MDTRouting::deleteexpiredneighbors(std::vector<L3Address> nodes){
    std::vector<double> selfCoord;
    getSelfCoordinate(selfCoord);
    L3Address self = getSelfIPAddress();

    auto P = getPNeighborNodes();
    if (routingTable) {// add nodes to expired set if it is NextHop in routetable
        for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
            IRoute* rt = routingTable->getRoute(i);
            if (!rt) continue;
            L3Address nextHop = rt->getNextHopAsGeneric();
            if (P.find(nextHop) == P.end()) {
                if (std::find(nodes.begin(), nodes.end(), nextHop) == nodes.end()) {
                    nodes.push_back(nextHop);
                }
            }
        }
    }
    if(nodes.empty()){
        return;
    }
    std::map<L3Address, std::vector<L3Address>> nodeToDeleteSet;
    for (auto node : nodes) {
        std::vector<L3Address> deleteSet;
        std::vector<IRoute*> routesToRemove;

        if (routingTable) {
            for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
                IRoute* rt = routingTable->getRoute(i);
                if (!rt) continue;
                if (rt->getNextHopAsGeneric() == node) {
                    deleteSet.push_back(rt->getDestinationAsGeneric());
                    routesToRemove.push_back(rt);
                }
            }
        }
        else {
            EV_ERROR << "deleteExpiredNeighbors: no routing table\n";
        }
        if (!deleteSet.empty()) {
            nodeToDeleteSet[node] = deleteSet;
        }
        for (auto rt : routesToRemove) {
            EV_ERROR << "delete route has expired node info, dest = "
                     << rt->getDestinationAsGeneric()
                     << " nxt = " << rt->getNextHopAsGeneric() << "\n";
            routingTable->deleteRoute(rt);
        }
    }

    for(auto node :nodes){
        L3Address nxthop = L3Address();
        L3Address dest = node;
        std::vector<double> destCoord;
        if(getKnownNodeCoordinate(dest,destCoord)){
            L3Address nxthop = L3Address();
            if(forwardProtocol->findnexthop(dest, nxthop)){//greedy
                if(!nxthop.isUnspecified() && nxthop != self && nxthop != node && P.find(nxthop) != P.end()){
                    for(auto destaddr : nodeToDeleteSet[node]){
                        Routetuple tuple("physics",getSelfAddress(),L3Address(),nxthop,destaddr,-1,-1,simTime() + activeRouteTimeout);
                        ensureRoute(destaddr,nxthop,tuple);
                    }
                }
            }
            else{//MDT
                //change tag
                L3Address relay = nxthop;
                //find nxt
                IRoute* r = routingTable->findBestMatchingRoute(relay);
                if(r != nullptr){
                    nxthop = r->getNextHopAsGeneric();
                }
                else{
                    L3Address relaythop = L3Address();
                    if(forwardProtocol->findnexthop(relay, relaythop))//TODO can not handle
                        nxthop = relaythop;
                }
                if(!nxthop.isUnspecified() && nxthop != self && nxthop != node && P.find(nxthop) != P.end()){
                    for(auto destaddr : nodeToDeleteSet[node]){
                        Routetuple tuple("physics",getSelfAddress(),L3Address(),nxthop,destaddr,-1,-1,simTime() + activeRouteTimeout);
                        ensureRoute(destaddr,nxthop,tuple);
                    }
                }
                else{
                    EV_ERROR<<"no nxthop\n";
                }
            }
        }
        else{//unlikely
            EV_ERROR<<"deleteexpiredneighbors: can not find the coords of expired node"<<dest<<"\n";
        }

    }

}
void MDTRouting::sendMDTPacket(const Ptr<MDTControlPacket>& MDTPacket){//not use
  //send message and it is an essential interface
  Packet *packet = new Packet();
  //TODO:add necessary labels to this packet
  socket.send(packet);

}

void MDTRouting::receiveSignal(cComponent *source, simsignal_t signalID, cObject *obj, cObject *details){
  //use host->subscribe(linkBrokenSignal, this); to subscribe linkBrokenSignal
  //this method is temporarily useless
    Enter_Method("%s", cComponent::getSignalName(signalID));

    if (signalID == linkBrokenSignal) {//not use
        EV_DETAIL << "Received link break signal" << endl;
        Packet *datagram = check_and_cast<Packet *>(obj);
        const auto& networkHeader = findNetworkProtocolHeader(datagram);
        if (networkHeader != nullptr) {
            L3Address unreachableAddr = networkHeader->getDestinationAddress();
            if (unreachableAddr.getAddressType() == addressType) {
                // A node initiates processing for a RERR message in three situations:
                //
                //   (i)     if it detects a link break for the next hop of an active
                //           route in its routing table while transmitting data (and
                //           route repair, if attempted, was unsuccessful), or

                // TODO Implement: local repair

                //IRoute *route = routingTable->findBestMatchingRoute(unreachableAddr);

                //if (route && route->getSource() == this)
                    //handleLinkBreakSendRERR(route->getNextHopAsGeneric());
            }
        }
    }
    else if(signalID == receptionEndedSignal){
        EV_INFO << "receiveEndSignal from " << (source ? source->getFullPath() : "null")
                << " objClass=" << (obj ? obj->getClassName() : "null") << "\n";

        if (!obj) return;
        if (auto *sr = dynamic_cast<const inet::physicallayer::ScalarReception *>(obj)) {
            RecentRxEntry e;
            e.power_W = sr->getPower().get();
            e.endTime = simTime(); // reception end handler invoked at end time
            // try to get tx id
            /*if (const inet::physicallayer::ITransmission *tx = sr->getTransmission()) {
                e.txid = tx->getId();
            }*/
            // try transmitter id via sr->getReceiver()/getTransmitter() chain or details
            /*if (const inet::physicallayer::WirelessSignal *ws = dynamic_cast<const inet::physicallayer::WirelessSignal *>(obj)) {
                // unlikely, since obj is ScalarReception, but safe check
            }*/
            // If sr has methods to get startPosition or transmitter, collect them:
            // e.startPosition = sr->getStartPosition(); // if exists
            // if available: e.transmitterPath = ... (via sr->getTransmission()->getTransmitter()->getModule()->getFullPath())

            // Best-effort: try to get transmitter module and its id
            /*if (const inet::physicallayer::ITransmission *txptr = sr->getTransmission()) {
                if (const inet::physicallayer::IRadio *txRadio = txptr->getTransmitter()) {
                    try {
                        const cModule *txMod = dynamic_cast<const cModule *>(txRadio->getModule());
                        if (txMod) {
                            e.transmitterPath = txMod->getFullPath();
                            e.transmitterId = txRadio->getId();
                        }
                    } catch (...) { }
                }
            }*/

            // push into deque (maintain capacity)
            recentRxDeque.push_front(e);
            while (recentRxDeque.size() > recentRxCapacity) recentRxDeque.pop_back();

            EV_ERROR << "Stored recent RX entry: "<< " power=" << e.power_W
                     << " W endTime=" << e.endTime << "\n";

            // clean expired entries (optional)
            while (!recentRxDeque.empty() && simTime() - recentRxDeque.back().endTime > recentEntryLifetime)
                recentRxDeque.pop_back();

            return;
        }
    }
    else if(signalID == receptionStartedSignal){
        EV_INFO << "receiveSignal from " << (source ? source->getFullPath() : "null")
                << " objClass=" << (obj ? obj->getClassName() : "null") << "\n";

        if (!obj) return;
        if (auto *sr = dynamic_cast<const ScalarReception *>(obj)) {
            // 1
            double rxPower_W = sr->getPower().get(); // W

            // 2 get transmission id
            int txid = -1;
            if (const ITransmission *tx = sr->getTransmission()) {
                txid = tx->getId();
            }
            else {
                EV_ERROR<<"ITransmission *tx = sr->getTransmission() fail\n";
            }

            EV_INFO << "Direct ScalarReception power W: " << rxPower_W << ", txid=" << txid << "\n";

            // 3 to map
            if (txid != -1) {
                txPowerMap[txid] = { rxPower_W, simTime() };
            } else {
                //
                EV_WARN << "ScalarReception without transmission id  cannot map power to txid\n";
            }

            // 4 delete expired entry
            for (auto it = txPowerMap.begin(); it != txPowerMap.end(); ) {
                if (simTime() - it->second.ts > txPowerEntryLifetime)
                    it = txPowerMap.erase(it);
                else
                    ++it;
            }

            return;
        }
    }
}
//these follow three methods are indispensable

void MDTRouting::handleStartOperation(LifecycleOperation *operation)
{
    rebootTime = simTime();

    socket.setOutputGate(gate("socketOut"));
    socket.setCallback(this);
    EV_INFO << "binding to port " << mdtUDPPort << endl;
    socket.bind(L3Address(), mdtUDPPort);
    socket.setBroadcast(true);

}

void MDTRouting::handleStopOperation(LifecycleOperation *operation)
{
    socket.close();
    clearState();
}

void MDTRouting::handleCrashOperation(LifecycleOperation *operation)
{
    socket.destroy();
    clearState();
}

void MDTRouting::clearState()
{
    if(routeTimer)cancelAndDelete(routeTimer);
    targetAddressToDelayedPackets.clear();
}
IRoute *MDTRouting::findBestMatchingRoute(const L3Address& dest){
    return routingTable->findBestMatchingRoute(dest);
}
void MDTRouting::updateRoutetuple(const L3Address& sourceAddr,const L3Address& predAddr,const L3Address& succAddr,const L3Address& destAddr,
                            double cost,double error){
    for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
        IRoute *rt = routingTable->getRoute(i);
        if(rt->getDestinationAsGeneric() == destAddr){
            MDTData *routeData = rt ? dynamic_cast<MDTData *>(rt->getProtocolData()) : nullptr;
            if(routeData){
                Routetuple tuple = routeData->getroutetuple();
                if(tuple.source == sourceAddr && tuple.dest == destAddr){
                    EV_INFO<<"update route dest = "<<destAddr<<"\n";
                    routeData->setroutetuplepath(sourceAddr,predAddr,succAddr,destAddr);
                    routeData->settimeout(simTime() + activeRouteTimeout);
                    rt->setNextHop(succAddr);
                    return;
                }
            }
        }
    }
    EV_INFO<<"add new route: src = "<<sourceAddr<<"pred = "<<predAddr<< "succ = "<<succAddr<<"dest = "<<destAddr<<"\n";
    Routetuple tuple("physics",sourceAddr,predAddr,succAddr,destAddr,-1,-1,simTime() + activeRouteTimeout);
    ensureRoute(destAddr,succAddr,tuple);
}
bool MDTRouting::isDTNeighborOfTarget(const std::vector<double>& targetCoord, const L3Address& addr){
    return neighbor->isDTNeighborOfTarget(targetCoord,addr);
}
std::vector<NeighborEntry> MDTRouting::computeDTNeighborsForCoord(const std::vector<double>& targetCoord, const L3Address& addr){
    return neighbor->computeDTNeighborsForCoord(targetCoord,addr);
}
void MDTRouting::updateDTNeighbors(){
    neighbor->updateDTNeighbors(true);
}
bool MDTRouting::addOrUpdateKnownNode(const std::vector<NeighborEntry>& entrys){
    if(neighbor){
        return neighbor->addOrUpdateKnownNode(entrys);
    }
    else{
        EV_ERROR<<"no neighbor for addOrUpdateKnownNode";
        return false;
    }
}
bool MDTRouting::addOrUpdateDTNeighbors(const std::vector<NeighborEntry>& entrys){
    if(neighbor){
        return neighbor->addOrUpdateDTNeighbors(entrys);
    }
    else{
        EV_ERROR<<"no neighbor for addOrUpdateDTNeighbors";
        return false;
    }
}
std::map<L3Address, NeighborEntry> MDTRouting::getDTNeighbors() const{
    if (neighbor)
        return neighbor->getDTNeighbors();
    else
        return std::map<L3Address, NeighborEntry>();
}
std::map<L3Address, NeighborEntry> MDTRouting::getMinSimplexNodes() const{
    if (neighbor)
        return neighbor->getMinSimplexNodes();
    else
        return std::map<L3Address, NeighborEntry>();
}
std::map<L3Address, NeighborEntry> MDTRouting::getPNeighborNodes() const{
    if (neighbor)
        return neighbor->getPNeighborNodes();
    else
        return std::map<L3Address, NeighborEntry>();
}
std::map<L3Address, NeighborEntry> MDTRouting::getweaklinkNeighbors() const{
    if (neighbor)
        return neighbor->getweaklinkNeighbors();
    else
        return std::map<L3Address, NeighborEntry>();
}
L3Address MDTRouting::getSelfAddress() const
{
    return getSelfIPAddress();
}



bool MDTRouting::findTupleByDestandUpdate(const L3Address& dest, Routetuple &outTuple) const
{
    for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
        IRoute *rt = routingTable->getRoute(i);
        MDTData *routeData = rt ? dynamic_cast<MDTData *>(rt->getProtocolData()) : nullptr;
        if(routeData){
            routeData->settimeout(simTime() + activeRouteTimeout);
            Routetuple routetuple = routeData->getroutetuple();
            if(routetuple.dest == dest){
                outTuple = routetuple;
                return true;
            }
        }
    }
    return false;
}

bool MDTRouting::getNeighborNodeCoordinate(const L3Address& node, std::vector<double>& outCoords) const
{
    if (neighbor)
        return neighbor->getNeighborCoords(node, outCoords);
    return false;
}
bool MDTRouting::getDTNeighborNodeCoordinate(const L3Address& node, std::vector<double>& outCoords) const
{
    if (neighbor)
        return neighbor->getDTNeighborCoords(node, outCoords);
    return false;
}
bool MDTRouting::getKnownNodeCoordinate(const L3Address& node, std::vector<double>& outCoords) const{
    if (neighbor)
        return neighbor->getKnownNodeCoords(node, outCoords);
    return false;
}
std::map<L3Address, NeighborEntry> MDTRouting::getKnownNodes() const {
    return neighbor->getKnownNodes();
}

//rpc
static int addressToNodeId(const std::string &addrStr) {
    auto pos = addrStr.find_last_of('.');
    if (pos == std::string::npos) return -1;
    try {
        int last = std::stoi(addrStr.substr(pos+1));
        return last - 1; // mapping rule
    } catch (...) {
        return -1;
    }
}
std::vector<int> MDTRouting::getMsgEventInfo(Packet *datagram,std::string type){
    L3Address origSrc = L3Address(); // unspecified by default
    L3Address origDst = L3Address();
    bool haveNetworkHeader = false;

    if (auto ipv4hdr = datagram->peekAtFront<Ipv4Header>()) {
        origSrc = ipv4hdr->getSrcAddress();
        origDst = ipv4hdr->getDestAddress();
        haveNetworkHeader = true;
    }

    L3Address prevHop = L3Address();
    if (auto l3ind = datagram->findTag<L3AddressInd>()) {
        prevHop = l3ind->getSrcAddress();
    }

    if (!haveNetworkHeader) {
        if (auto l3ind2 = datagram->findTag<L3AddressInd>()) {
            origSrc = l3ind2->getSrcAddress();
            origDst = l3ind2->getDestAddress();
        } else if (auto l3req = datagram->findTag<L3AddressReq>()) {
            origDst = l3req->getDestAddress();
            origSrc = l3req->getSrcAddress();
        }
    }

    L3Address nextHop = L3Address(); // unspecified
    if (!origDst.isUnspecified()) {
        IRoute *rt = routingTable->findBestMatchingRoute(origDst);
        if (rt) {
            try {
                nextHop = rt->getNextHopAsGeneric();
            } catch (...) {
                // TODO fallback
            }
        }
        if (nextHop.isUnspecified()) {
            if (auto nhReq = datagram->findTag<NextHopAddressReq>()) {
                nextHop = nhReq->getNextHopAddress();
            }
        }
    }

    std::string srcStr = origSrc.str();
    std::string prevStr = prevHop.str();
    std::string nextStr = nextHop.str();
    std::string destStr = origDst.str();
    // map
    int srcNodeId = addressToNodeId(srcStr);
    int prevNodeId = addressToNodeId(prevStr);
    int nextNodeId = addressToNodeId(nextStr);
    int destNodeId = addressToNodeId(destStr);

    std::vector<int> res = {-1,-1,-1,-1};
    if(type == "IN"){
        res[0] = srcNodeId;
        res[1] = prevNodeId;
        res[3] = destNodeId;
        return res;
    }
    else if(type == "FORWARD"){
        res[0] = srcNodeId;
        res[1] = prevNodeId;
        res[2] = nextNodeId;
        res[3] = destNodeId;
        return res;
    }
    else if(type == "OUT"){
        res[0] = srcNodeId;
        res[2] = nextNodeId;
        res[3] = destNodeId;
        return res;
    }
    else{
        return res;
    }
}
void MDTRouting::routetablelog(){
    logRoutingTable(this, 0, routingTable);
}
void MDTRouting::finish() {
    EV_ERROR<<"route table:"<<"\n";
    for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
        IRoute *rt = routingTable->getRoute(i);
        EV_ERROR <<"dest :"<< rt->getDestinationAsGeneric() << "PrefixLength: " << rt->getPrefixLength()
           << " NextHop: " << rt->getNextHopAsGeneric() << endl;
        MDTData *routeDestData = rt ? dynamic_cast<MDTData *>(rt->getProtocolData()) : nullptr;
        if(routeDestData){
            Routetuple tuple = routeDestData->getroutetuple();
            EV_ERROR <<"src = "<<tuple.source<<"dest = "<<tuple.dest<<"pred = "<<tuple.pred<<"succ = "<<tuple.succ<<"\n";
        }
    }
}


MDTRouting::~MDTRouting()
{
    clearState();
    //delete source;
}
}
}







