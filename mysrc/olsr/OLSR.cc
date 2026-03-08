/*
 * OLSR.cc
 *
 *  Created on: Feb 9, 2026
 *      Author: yychen
 */
#include "../olsr/OLSR.h"

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

#include "inet/linklayer/common/InterfaceTag_m.h"
#include "inet/networklayer/common/NextHopAddressTag_m.h"

#include <iomanip>
#include <cstdint>
#include "inet/common/packet/chunk/Chunk.h"

#include "inet/transportlayer/udp/UdpHeader_m.h"

#include "inet/networklayer/common/L3AddressResolver.h"

namespace mysrc {
namespace olsr {
using namespace inet;

Define_Module(OLSR);

const int KIND_DELAYEDSEND = 100;//delay

void OLSR::initialize(int stage){
    if (stage == INITSTAGE_ROUTING_PROTOCOLS){
        addressType = getSelfIPAddress().getAddressType(); // needed for handleStartOperation()
    }
    RoutingProtocolBase::initialize(stage);
    if (stage == INITSTAGE_LOCAL){
        host = getContainingNode(this);
        routingTable.reference(this, "routingTableModule", true);
        interfaceTable.reference(this, "interfaceTableModule", true);
        networkProtocol.reference(this, "networkProtocolModule", true);

        helloMsgTimer = new cMessage("core-route-timer");
        helloInterval = par("helloInterval").doubleValue();
        scheduleAt(simTime() + uniform(0, helloInterval), helloMsgTimer);

        tcCacheTimer = new cMessage("core-route-timer");
        tcCacheInterval = par("tcCacheInterval").doubleValue();
        scheduleAt(simTime() + uniform(0, tcCacheInterval), tcCacheTimer);

        tcGroupTimer = new cMessage("core-route-timer");
        tcInterval = par("tcInterval").doubleValue();
        scheduleAt(simTime() + uniform(0, tcInterval), tcGroupTimer);
    }
    else if (stage == INITSTAGE_ROUTING_PROTOCOLS) {
        networkProtocol->registerHook(0, this);
        //host->subscribe(linkBrokenSignal, this);
    }
}
OLSR::OLSR(){

}
L3Address OLSR::getSelfIPAddress() const
{
    return routingTable->getRouterIdAsGeneric();
}
void OLSR::handleMessageWhenUp(cMessage *msg){
    if (msg->isSelfMessage()){
        EV_INFO<<"self message\n";
        if(msg == helloMsgTimer){
            sendHelloPacket();
            scheduleAt(simTime() + helloInterval, helloMsgTimer);
        }
        else if(msg == tcGroupTimer){
            constructLocalTc();
            sendTcPacket();
            scheduleAt(simTime() + tcInterval, tcGroupTimer);
        }
        else if(msg == tcCacheTimer){
            sendTcPacket();
            scheduleAt(simTime() + tcCacheInterval, tcCacheTimer);
        }
        else if (msg->getKind() == KIND_DELAYEDSEND){
            auto timer = check_and_cast<OlsrPacketHolderMessage *>(msg);
            Packet *p = timer->removeOwnedPacket();
            sendPacket(p,olsrUDPPort,0);
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
void OLSR::socketDataArrived(UdpSocket *socket, Packet *packet)
{
    // process incoming packet
    processPacket(packet);
}

void OLSR::socketErrorArrived(UdpSocket *socket, Indication *indication)
{
    EV_WARN << "Ignoring UDP error report " << indication->getName() << endl;
    delete indication;
}

void OLSR::socketClosed(UdpSocket *socket)
{
    if (operationalState == State::STOPPING_OPERATION)
        startActiveOperationExtraTimeOrFinish(par("stopOperationExtraTime"));
}

void OLSR::handleStartOperation(LifecycleOperation *operation)
{
    //rebootTime = simTime();

    socket.setOutputGate(gate("socketOut"));
    socket.setCallback(this);
    EV_INFO << "binding to port " << olsrUDPPort << endl;
    socket.bind(L3Address(), olsrUDPPort);
    socket.setBroadcast(true);

    const char *mcast = nullptr;//par("multicastGroup").stringValue(); // add this param in NED / omnetpp.ini
    if (!mcast || strlen(mcast) == 0)
        mcast = "224.0.0.109"; // default (LL-MANET-Routers), or use 224.0.0.103 for OLSRv1 if you prefer
    multicastGroup = L3AddressResolver().resolve(mcast);
    // find interface id
    auto *it = CHK(interfaceTable->findInterfaceByName(par("interface")));
    ifId = it->getInterfaceId();
    socket.joinMulticastGroup(multicastGroup, ifId);
}

void OLSR::handleStopOperation(LifecycleOperation *operation)
{
    socket.close();
    //clearState();
}

void OLSR::handleCrashOperation(LifecycleOperation *operation)
{
    socket.destroy();
    //clearState();
}

void OLSR::processPacket(Packet *packet){
    L3Address sourceAddr = packet->getTag<L3AddressInd>()->getSrcAddress();
    unsigned int arrivalPacketTTL = packet->getTag<HopLimitInd>()->getHopLimit() - 1;
    const auto& olsrPacket = packet->popAtFront<OlsrControlPacket>();

    auto packetType = olsrPacket->getPacketType();
    switch (packetType) {
        case OlsrPktType::Hello:
            processHello(CHK(dynamicPtrCast<OlsrHello>(olsrPacket->dupShared())));
            RoutingTableComputation();
            //delete expireTuples
            state.expireTuples(simTime());
            delete packet;
            return;
        case OlsrPktType::Tc:
            processTc(CHK(dynamicPtrCast<OlsrTcGroup>(olsrPacket->dupShared())), sourceAddr);
            RoutingTableComputation();
            if(getContainingNode(this)->getId() == 10 || getContainingNode(this)->getId() == 11){//TODO
                printRouteTable();
                state.printOlsrState();
            }
            delete packet;
            return;
        default:
            throw cRuntimeError("OLSR Control Packet arrived with undefined packet type: %d", packetType);
    }
}
void OLSR::processHello(const Ptr<OlsrHello>& hello){
    L3Address sender = hello->getOriginator();

    LinkSensing(sender, hello);

    PopulateNeighborSet(sender);

    PopulateTwoHopNeighborSet(sender, hello);

    MprComputation();

    PopulateMprSelectorSet(sender, hello);
}
void OLSR::LinkSensing(const L3Address& sender, const Ptr<OlsrHello>& hello)
{
    SimTimeType now = simTime();
    L3Address local = getSelfIPAddress();

    // neighborHoldTime
    SimTimeType holdTime = neighborHoldTime;

    // get LinkTuple
    OlsrLinkTuple *lt = state.getOrCreateLink(local, sender);

    // update asymTime
    lt->asymTime = now + holdTime;

    bool symmetric = false;

    int blocks = hello->getLinkCodesArraySize();
    int offset = 0;

    for (int i = 0; i < blocks; ++i)
    {
        uint8_t linkCode = hello->getLinkCodes(i);
        uint8_t linkType = linkCode & 0x03;
        uint8_t neighborType = (linkCode >> 2) & 0x03;

        int count = hello->getLinkCounts(i);

        for (int j = 0; j < count; ++j)
        {
            L3Address addr = hello->getNbIfaceAddrs(offset + j);

            if (addr != local)
                continue;

            if (linkType == 3) // LOST_LINK
                lt->symTime = now;
            else
                symmetric = true;
        }

        offset += count;
    }

    //  symTime
    if (symmetric)
        lt->symTime = now + holdTime;

    //
    lt->time = std::max(lt->asymTime, lt->symTime);

    //
    EV_INFO << "OLSR: LinkSensing link " << lt->localIfaceAddr.str() << " <-> "
              << lt->neighborIfaceAddr.str()
              << " symTime=" << lt->symTime
              << " asymTime=" << lt->asymTime
              << " time=" << lt->time << "\n";
}
void OLSR::PopulateNeighborSet(const L3Address& sender)
{
    SimTime now = simTime();

    OlsrLinkTuple* lt = state.findLink(getSelfIPAddress(),sender);

    if (!lt)
        return;

    uint8_t status = (lt->symTime >= now) ? 2 : 1;

    state.addOrUpdateNeighbor(
        sender,
        willingness,
        lt->time
    )->status = status;
}
void OLSR::PopulateTwoHopNeighborSet(
        const L3Address& sender,
        const Ptr<OlsrHello>& hello)
{
    EV_ERROR<<"PopulateTwoHopNeighborSet "<<"sender: "<<sender<<"\n";
    SimTimeType now = simTime();

    int offset = 0;

    int blocks = hello->getLinkCodesArraySize();
    EV_ERROR<<"blocks = "<<blocks<<"\n";
    for (int i = 0; i < blocks; ++i)
    {
        uint8_t linkCode = hello->getLinkCodes(i);

        uint8_t linkType = linkCode & 0x03;
        uint8_t neighborType = (linkCode >> 2) & 0x03;
        EV_ERROR<<"linkType: "<<(int)linkType<<"neighborType: "<<(int)neighborType<<"\n";

        int count = hello->getLinkCounts(i);

        // only symmetric neighbors
        if (neighborType != 1 && neighborType != 2) // SYM_NEIGH or MPR_NEIGH
        {
            offset += count;
            continue;
        }

        if (linkType != 2) // SYM_LINK
        {
            offset += count;
            continue;
        }

        for (int j = 0; j < count; ++j)
        {
            L3Address addr = hello->getNbIfaceAddrs(offset + j);

            if (addr == getSelfIPAddress())
                continue;

            if (state.findNeighbor(addr))
                continue;

            state.addOrUpdateTwoHop(
                    getSelfIPAddress(),
                    sender,
                    addr,
                    now + neighborHoldTime);
            EV_ERROR<<"new 2-hop neighbor: "<<addr<<"\n";
        }

        offset += count;
    }
}
void OLSR::MprComputation()
{
    std::set<L3Address> mpr;
    std::set<L3Address> uncovered;
    auto nb2sets = state.getnb2HopSet();
    auto nset = state.getNbSet();

    // fill uncovered
    for (auto t : nb2sets)
        uncovered.insert(t->nb2hopAddr);

    // find one cover neighbor first
    for (auto t2 : nb2sets)
    {
        // find 2-hop neighbor cover set
        std::vector<L3Address> candidates;
        for (auto nb : nset)
        {
            if (nb->status != 2 || nb->willingness == 0)
                continue;

            for (auto t : nb2sets)
                if (t->nbMainAddr == nb->mainAddr && t->nb2hopAddr == t2->nb2hopAddr)
                    candidates.push_back(nb->mainAddr);
        }

        if (candidates.size() == 1)
        {
            L3Address uniqueNB = candidates[0];
            mpr.insert(uniqueNB);
            // delete 2-hop neighbor
            for (auto t : nb2sets)
                if (t->nbMainAddr == uniqueNB)
                    uncovered.erase(t->nb2hopAddr);
        }
    }

    // greedy select
    while (!uncovered.empty())
    {
        L3Address best;
        int bestCount = -1;

        for (auto nb : nset)
        {
            if (nb->status != 2 || nb->willingness == 0)
                continue;

            int count = 0;
            for (auto t : nb2sets)
                if (t->nbMainAddr == nb->mainAddr && uncovered.count(t->nb2hopAddr))
                    count++;

            if (count > bestCount)
            {
                bestCount = count;
                best = nb->mainAddr;
            }
        }

        if (bestCount <= 0)
            break;

        mpr.insert(best);

        for (auto t : nb2sets)
            if (t->nbMainAddr == best)
                uncovered.erase(t->nb2hopAddr);
    }

    state.SetMprSet(mpr);
}
void OLSR::PopulateMprSelectorSet(
        const L3Address& sender,
        const Ptr<OlsrHello>& hello)
{
    SimTimeType now = simTime();

    int offset = 0;

    int blocks = hello->getLinkCodesArraySize();

    for (int i = 0; i < blocks; ++i)
    {
        uint8_t linkCode = hello->getLinkCodes(i);

        uint8_t neighborType = (linkCode >> 2) & 0x03;

        int count = hello->getLinkCounts(i);

        if (neighborType != 2) // MPR_NEIGH
        {
            offset += count;
            continue;
        }

        for (int j = 0; j < count; ++j)
        {
            L3Address addr = hello->getNbIfaceAddrs(offset + j);

            if (addr == getSelfIPAddress())
            {
                state.addMprSelector(
                        sender,
                        now + neighborHoldTime);
            }
        }

        offset += count;
    }
}

void OLSR::processTc(const Ptr<OlsrTcGroup>& tc, const L3Address& sender)
{
    SimTimeType now = simTime();

    // 1 sender must be symmetric neighbor
    auto nb = state.findNeighbor(sender);
    if (!nb || nb->status != 2)
        return;

    int M = tc->getOriginatorsArraySize();
    int addrIndex = 0;

    for (int i = 0; i < M; i++)
    {
        L3Address origin = tc->getOriginators(i);
        uint16_t ansn = tc->getAnsns(i);
        uint16_t seq = tc->getSeqNums(i);
        uint8_t hop = tc->getHopCounts(i);
        int advCount = tc->getAdvCounts(i);

        // -----------------------------
        // 2 check if newer topology exists
        // -----------------------------

        auto& tset = state.getTopologySet();   //  &

        bool newerExists = false;

        for (auto *t : tset)
        {
            if (t->lastAddr == origin && t->seq > ansn)
            {
                newerExists = true;
                break;
            }
        }

        if (newerExists)
        {
            addrIndex += advCount;
            continue;
        }


        // -----------------------------
        // 3 delete older topology
        // -----------------------------

        for (auto it = tset.begin(); it != tset.end(); )
        {
            OlsrTopologyTuple *t = *it;

            if (t->lastAddr == origin && t->seq < ansn)
            {
                delete t;                  // delete tuple
                it = tset.erase(it);       // TopologySet
            }
            else
            {
                ++it;
            }
        }

        // -----------------------------
        // 4 update topology
        // -----------------------------

        for (int j = 0; j < advCount; j++)
        {
            L3Address dest = tc->getAdvertisedNeighbors(addrIndex++);

            state.addOrUpdateTopology(
                    dest,
                    origin,
                    ansn,
                    now + topologyHoldTime);
        }

        // -----------------------------
        // 5 duplicate detection
        // -----------------------------

        bool duplicate = state.findDup(origin, seq);

        if (!duplicate)
            state.addDup(origin, seq, now + DupHoldTime);

        // -----------------------------
        // 6 forwarding decision
        // -----------------------------

        if (!duplicate && state.isMprSelector(sender))
        {
            CachedTcEntry e;

            e.originator = origin;
            e.ansn = ansn;
            e.msgSeq = seq;
            e.hopCount = hop + 1;
            e.expiry = now + tcInterval;   // cache lifetime

            for (int j = 0; j < advCount; j++)
            {
                e.advertisedNeighbors.push_back(
                    tc->getAdvertisedNeighbors(
                        addrIndex - advCount + j));
            }

            tcCache.push_back(e);
        }
    }

}

std::vector<L3Address> OLSR::getAdvertisedNeighborList()
{
    std::vector<L3Address> res;
    std::vector<OlsrMprSelectorTuple*> mset = state.getmprSelSet();
    for (auto t : mset)
        res.push_back(t->mainAddr);

    return res;
}
void OLSR::constructLocalTc()
{
    CachedTcEntry e;

    e.originator = getSelfIPAddress();
    e.ansn = ++ansn;
    e.msgSeq = ++tcMsgSeq;
    e.hopCount = 0;
    e.hopLimit = 255;
    e.advertisedNeighbors = getAdvertisedNeighborList();
    e.expiry = simTime() + topologyHoldTime;

    tcCache.push_back(e);
}
void OLSR::sendPacket(Packet *packet, /*const L3Address& destAddr, no need*/ int destPort,double delay){
    Enter_Method("sendPacket");
    if (!packet) {
        EV_ERROR << "sendPacketTo: called with null packet\n";
        return;
    }
    take(packet);
    int port = (destPort == 0 ? (int)olsrUDPPort : destPort);

    L3Address dest = multicastGroup;

    EV_INFO << "sendPacket: sending to multicast " << dest << " port=" << port
            << " originalPacket='" << packet->getName() << "'\n";

    // HELLO -> TTL=1
    // TC    -> TTL=255
    int hopLimit = 255;
    std::string name = packet->getName();
    if (name.find("HELLO") != std::string::npos ||
        name.find("Hello") != std::string::npos)
        hopLimit = 1;

    // Interface
    packet->addTagIfAbsent<InterfaceReq>()->setInterfaceId(ifId);

    // IP header fields
    packet->addTagIfAbsent<HopLimitReq>()->setHopLimit(hopLimit);
    packet->addTagIfAbsent<L3AddressReq>()->setDestAddress(dest);

    // UDP
    packet->addTagIfAbsent<L4PortReq>()->setDestPort(port);

    if (delay == 0) {
        /*if (packet) {
            emit(ctrlPacketsSentSignal, packet);
        }*///data collect
        socket.sendTo(packet, dest, port);
        EV_INFO << "sendPacketTo: sent immediately\n";
    } else {
        auto *timer = new OlsrPacketHolderMessage("olsr-delay-send", KIND_DELAYEDSEND);
        timer->setOwnedPacket(packet);
        scheduleAfter(delay, timer);
        EV_INFO << "sendPacketTo: scheduled for delayed send after " << delay << "\n";
    }
}

void OLSR::sendHelloPacket()
{
    Enter_Method("sendHelloPacket");

    simtime_t now = simTime();

    // create chunk
    auto hello = makeShared<OlsrHello>();
    // packetType: use your enum; _msg generated setter takes numeric (uint8)
    hello->setPacketType(OlsrPktType::Hello);

    // header fields
    hello->setOriginator(getSelfIPAddress());
    hello->setHTime(helloInterval);
    hello->setWillingness(willingness);
    hello->setMsgSeq(++helloMsgSeq);

    // collect link messages from state
    std::vector<uint8_t> linkCodes;
    std::vector<int> linkCounts;
    std::vector<inet::L3Address> flatNbAddrs;

    const auto &links = state.getLinkSet();
    for (auto lt : links) {
        // include only link tuples for our local iface and not expired
        if (lt->localIfaceAddr != getSelfIPAddress()) continue;
        if (lt->time < now) continue;

        // compute linkType: 0=LOST,1=ASYM,2=SYM
        uint8_t linkType = 0;
        if (lt->symTime >= now)
            linkType = 2;//2=SYM
        else if (lt->asymTime >= now)
            linkType = 1;//1=ASYM,
        else
            linkType = 3;//0=LOST

        // neighborType: 0=NOT_NEIGH/UNKNOWN,2=MPR_NEIGH,1=SYM_NEIGH
        uint8_t neighborType = 0;
        std::set<L3Address> mset = state.GetMprSet();
        auto it = mset.find(lt->neighborIfaceAddr);
        if (it != mset.end()){
            neighborType = 2; // MPR_NEIGH
        }
        else{
            auto nb = state.findNeighbor(lt->neighborIfaceAddr);
            if (nb && nb->status == 2)//0=LOST,1=ASYM,2=SYM
                neighborType = 1; // SYM_NEIGH
        }

        uint8_t linkCode = (linkType & 0x03) | ((neighborType << 2) & 0x0c);
        linkCodes.push_back(linkCode);

        // In single-interface simplified model, report only neighbor's iface address
        linkCounts.push_back(1);
        flatNbAddrs.push_back(lt->neighborIfaceAddr);
        /*if(getContainingNode(this)->getId() == 10 || getContainingNode(this)->getId() == 11){
            EV_ERROR<<"flatNbAddrs: "<<lt->neighborIfaceAddr<<"\n";
        }*/
        // If you want to report additional neighbor interface addresses,
        // query state for neighbor's interfaces and append them, incrementing linkCounts.back().
    }

    // fill chunk arrays
    int N = (int)linkCodes.size();
    hello->setLinkCodesArraySize(N);
    hello->setLinkCountsArraySize(N);
    int flatSize = (int)flatNbAddrs.size();
    hello->setNbIfaceAddrsArraySize(flatSize);

    for (int i = 0; i < N; ++i) {
        hello->setLinkCodes(i, linkCodes[i]);
        hello->setLinkCounts(i, linkCounts[i]);
    }
    for (int i = 0; i < flatSize; ++i) {
        hello->setNbIfaceAddrs(i, flatNbAddrs[i]);
    }

    //debug
    if(getContainingNode(this)->getId() == 10 || getContainingNode(this)->getId() == 11){//TODO
        EV_INFO << "\n==== SEND HELLO ====\n";
        EV_INFO << "originator=" << getSelfIPAddress() << endl;

        for(int i=0;i<N;i++)
        {
            uint8_t code = linkCodes[i];

            uint8_t linkType = code & 0x03;
            uint8_t neighborType = (code >> 2) & 0x03;

            EV_INFO << "block "
                    << i
                    << " linkType=" << (int)linkType
                    << " neighborType=" << (int)neighborType
                    << " addr=" << flatNbAddrs[i]
                    << endl;
        }
    }
    // optional: keep total count
    //hello->setCount(flatSize);

    // conservative chunk length estimate (for simulation bookkeeping)
    int addrBytes = (getSelfIPAddress().getType() == L3Address::IPv6 ? 16 : 4);
    int bytes = 0;
    bytes += addrBytes;            // originator
    bytes += 8;                    // hTime (double)
    bytes += 1;                    // willingness
    bytes += 2;                    // msgSeq (uint16)
    bytes += N * 1;                // linkCodes
    bytes += N * 4;                // linkCounts (int)
    bytes += flatSize * addrBytes; // addresses
    hello->setChunkLength(B(bytes));

    // build packet
    Packet *packet = new Packet("OLSR-HELLO", hello);

    // jitter per RFC: <= helloInterval/4
    double maxJit = std::min(helloInterval / 4.0, maxPeriodicJitter);
    double jitter = uniform(0, maxJit);

    // send via unified interface
    sendPacket(packet, olsrUDPPort, jitter);

    EV_INFO << "OLSR: sent HELLO seq=" << helloMsgSeq << " links=" << N << " jitter=" << jitter << "\n";
}

void OLSR::sendTcPacket()
{
    SimTime now = simTime();

    tcCache.erase(
        std::remove_if(
            tcCache.begin(),
            tcCache.end(),
            [now](const CachedTcEntry& e)
            {
                return e.expiry < now;
            }),
        tcCache.end());

    if (tcCache.empty())
        return;

    auto group = makeShared<OlsrTcGroup>();

    group->setPacketType(OlsrPktType::Tc);

    group->setGenTime(simTime());

    int M = tcCache.size();

    group->setOriginatorsArraySize(M);
    group->setAnsnsArraySize(M);
    group->setSeqNumsArraySize(M);
    group->setHopCountsArraySize(M);
    //group->setHopLimitsArraySize(M);
    group->setAdvCountsArraySize(M);

    std::vector<L3Address> flat;

    for (int i = 0; i < M; i++)
    {
        auto &e = tcCache[i];

        group->setOriginators(i, e.originator);
        group->setAnsns(i, e.ansn);
        group->setSeqNums(i, e.msgSeq);
        group->setHopCounts(i, e.hopCount);
        //group->setHopLimits(i, e.hopLimit);

        group->setAdvCounts(i,
            e.advertisedNeighbors.size());

        for (auto &a : e.advertisedNeighbors)
            flat.push_back(a);
    }

    group->setAdvertisedNeighborsArraySize(flat.size());

    for (int i = 0; i < (int)flat.size(); i++)
        group->setAdvertisedNeighbors(i, flat[i]);

    int addrBytes = (getSelfIPAddress().getType() == L3Address::IPv6 ? 16 : 4);

    int bytes = 0;

    bytes += 8;

    bytes += M * addrBytes;
    bytes += M * 2;
    bytes += M * 2;
    bytes += M * 1;
    bytes += M * 4;

    int flatSize = flat.size();
    bytes += flatSize * addrBytes;

    group->setChunkLength(B(bytes));

    Packet *pkt = new Packet("OLSR-TC", group);

    double maxJ = std::min(tcCacheInterval / 4.0, maxPeriodicJitter );
    double jitter = uniform(0, maxJ);

    sendPacket(pkt, olsrUDPPort, jitter);

    tcCache.clear();
}


INetfilter::IHook::Result OLSR::ensureRouteForDatagram(Packet *datagram){
    //delete expireTuples
    state.expireTuples(simTime());
    RoutingTableComputation();

    return INetfilter::IHook::ACCEPT;
}
void OLSR::RoutingTableComputation()
{
    clearOlsrRoutes();

    std::map<L3Address, TempRoute> routes;
    L3Address self = getSelfIPAddress();

    // 1-hop
    std::vector<OlsrNeighborTuple*> nset = state.getNbSet();
    for (auto nb : nset)
    {
        if (nb->status != 2)
            continue;

        TempRoute r;
        r.dest = nb->mainAddr;
        r.nextHop = nb->mainAddr;
        r.distance = 1;

        routes[r.dest] = r;
    }

    // 2-hop
    std::vector<OlsrTwoHopTuple*> nb2sets = state.getnb2HopSet();
    for (auto t : nb2sets)
    {
        if (routes.count(t->nb2hopAddr))
            continue;

        auto it = routes.find(t->nbMainAddr);

        if (it == routes.end())
            continue;

        TempRoute r;

        r.dest = t->nb2hopAddr;
        r.nextHop = it->second.nextHop;
        r.distance = 2;

        routes[r.dest] = r;
    }

    // topology expansion
    int h = 2;

    while (true)
    {
        bool added = false;
        std::vector<OlsrTopologyTuple*> tset = state.getTopologySet();
        for (auto t : tset)
        {
            if (routes.count(t->destination))
                continue;

            auto it = routes.find(t->lastAddr);

            if (it == routes.end())
                continue;

            if (it->second.distance != h)
                continue;

            TempRoute r;

            r.dest = t->destination;
            r.nextHop = it->second.nextHop;
            r.distance = h + 1;

            routes[r.dest] = r;

            added = true;
        }

        if (!added)
            break;

        h++;
    }

    // write to inet routing table
    for (auto &kv : routes)
    {
        auto &r = kv.second;

        createRoute(r.dest, r.nextHop, r.distance);
    }
}
void OLSR::clearOlsrRoutes()
{
    std::vector<IRoute*> toDelete;

    for (int i = 0; i < routingTable->getNumRoutes(); i++)
    {
        IRoute *r = routingTable->getRoute(i);

        if (r->getSourceType() == IRoute::OLSR)
            toDelete.push_back(r);
    }

    for (auto r : toDelete)
        routingTable->deleteRoute(r);
}

IRoute *OLSR::createRoute(const L3Address& destAddr,const L3Address& nextHop,unsigned int distance){
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
    newRoute->setSourceType(IRoute::OLSR);//Registration in inet is needed,location is inet/networklayer/contract/IRoute.h in SourceType
    newRoute->setSource(this);

    // adding protocol-specific fields
    OLSR_RouteData *newProtocolData = new OLSR_RouteData();
    newProtocolData->setDistance(distance);

    newRoute->setProtocolData(newProtocolData);

    EV_DETAIL << "Adding new route " << newRoute << endl;

    routingTable->addRoute(newRoute);

    return newRoute;
}
void OLSR::printRouteTable(){
    EV_ERROR <<"\n=== routingTable print ===\n";
    for (int i = 0; i < routingTable->getNumRoutes(); ++i) {
        IRoute *rt = routingTable->getRoute(i);
        EV_ERROR << "dest: "<< rt->getDestinationAsGeneric() << "/" << "next hop"
           << " -> " << rt->getNextHopAsGeneric() << endl;
    }
}

OLSR::~OLSR(){
    if(helloMsgTimer)cancelAndDelete(helloMsgTimer);
    if(tcCacheTimer)cancelAndDelete(tcCacheTimer);
    if(tcGroupTimer)cancelAndDelete(tcGroupTimer);
}

}
}
