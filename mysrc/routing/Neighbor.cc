/*
 * Neighbor.cc
 *
 *  Created on: Aug 31, 2025
 *      Author: yychen
 */

#include "Neighbor.h"
#include "MDTRouting.h" // we assume MDTRouting expose sendControlPacket(...)
#include "MdtRelayTag_m.h" // generated from MDTRelayTag.msg, adjust include path
#include "inet/networklayer/common/L3AddressTag_m.h"
#include "inet/networklayer/contract/ipv4/Ipv4Address.h" // for IPv4Address::ALLONES_ADDRESS
#include "inet/common/packet/ChunkQueue.h"
#include "inet/common/packet/Packet.h"

#include "MDTControlPackets_m.h"
#include <cmath>
#include <algorithm>
#include <iterator>
#include <vector>
#include <limits>

namespace mysrc {
namespace routing {

Define_Module(Neighbor);

static double coordDistanceSquared(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.empty() || b.empty()) return std::numeric_limits<double>::infinity();
    size_t n = std::min(a.size(), b.size());
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double d = a[i] - b[i];
        s += d * d;
    }
    return s;
}
// dest between point and line
static double pointSegmentDistanceSquared(const std::vector<double>& p,
                                          const std::vector<double>& a,
                                          const std::vector<double>& b)
{
    size_t n = std::min({p.size(), a.size(), b.size()});
    if (n == 0) return std::numeric_limits<double>::infinity();

    // ab = b - a
    std::vector<double> ab(n);
    std::vector<double> ap(n);
    for (size_t i = 0; i < n; ++i) {
        ab[i] = b[i] - a[i];
        ap[i] = p[i] - a[i];
    }

    // denom = |ab|^2
    double denom = 0.0;
    for (size_t i = 0; i < n; ++i) denom += ab[i] * ab[i];

    if (denom <= 0.0) {
        // same
        return coordDistanceSquared(p, a);
    }

    // t = dot(ap,ab) / denom
    double num = 0.0;
    for (size_t i = 0; i < n; ++i) num += ap[i] * ab[i];
    double t = num / denom;

    //
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;

    //
    std::vector<double> proj(n);
    for (size_t i = 0; i < n; ++i) proj[i] = a[i] + t * ab[i];

    return coordDistanceSquared(p, proj);
}

void Neighbor::initialize(int stage)
{
    if (stage == INITSTAGE_LOCAL) {
        EV_INFO << "Neighbor::initialize(stage=" << stage << ")\n";
        // owner is sibling named "core" in MDT compound module
        cModule *m = getModuleByPath("^.core");
        if (!m) {
            EV_ERROR << "Neighbor::initialize: failed to find core module at '^.core'\n";
        } else {
            owner = check_and_cast<MDTRouting*>(m);
        }

        //rng
        rng = owner->rng;
        //debug
        if(getContainingNode(owner)->getId() == 7){
            debug = true;
            EV_INFO<<"debug node 7\n";
        }

        // read parameters (NED/ini)
        keepAliveInterval = simtime_t(par("keepAliveInterval").doubleValue());
        KANInterval = simtime_t(par("KANInterval").doubleValue());
        weaklinkKANInterval = simtime_t(par("weaklinkKANInterval").doubleValue());
        presentKANInterval = KANInterval;

        neighborStaleTimeout = simtime_t(par("neighborStaleTimeout").doubleValue());
        knownNodeStaleTimeout = simtime_t(par("knownNodeStaleTimeout").doubleValue());
        tupleTimeout = simtime_t(par("tupleTimeout").doubleValue());//not used

        EV_INFO << "Neighbor params: keepAliveInterval=" << keepAliveInterval
                << ", neighborStaleTimeout=" << neighborStaleTimeout
                << ", tupleTimeout=" << tupleTimeout << "\n";

        keepAliveTimer = new cMessage("neighbor-keepalive-timer");
        KANTimer = new cMessage("neighbor-KAN-timer");
        initiailizationTimer = new cMessage("neighbor-init-timer");
        initInterval = par("initInterval").doubleValue();

        dtTimer = new cMessage("neighbor-dt-timer");
        dtInterval = par("dtInterval").doubleValue();
        scheduleAt(simTime() + uniform(1, dtInterval), dtTimer);

        // schedule first keepalive with jitter
        double jitter = uniform(0, keepAliveInterval.dbl());
        scheduleAt(simTime() + simtime_t(jitter), keepAliveTimer);
        EV_INFO << "Neighbor: scheduled first keepalive after " << jitter << "s\n";

        double jit = uniform(0, presentKANInterval.dbl());
        scheduleAt(simTime() + simtime_t(jit), KANTimer);


        // create DT manager using references to this->members
        dtManager = std::make_unique<NeighborDT>(knownNodeCoords,
                                                 dt,
                                                 vhToAddrs,
                                                 addrToVh,
                                                 pendingRebuildAddrs,
                                                 dtNeighbors,
                                                 minSimplexNodes);
    }
}

void Neighbor::handleMessage(cMessage *msg)
{
    if (msg == keepAliveTimer) {
        EV_INFO << "Neighbor::handleMessage - keepAliveTimer triggered at " << simTime() << "\n";
        // broadcast keepalive on the subnet (IPv4 all-ones)
        L3Address broadcastAddr = L3Address(inet::Ipv4Address::ALLONES_ADDRESS);
        sendKeepAliveTo(broadcastAddr);
        //delete time out node ,but keep node in known nodes

        //rpc
        Coord pos = owner->mobility->getCurrentPosition();
        Coord vel = owner->mobility->getCurrentVelocity();
        std::vector<double> posVec = {pos.x,pos.y,pos.z};
        std::vector<double> velVec = {vel.x,vel.y,vel.z};
        logPosition(this,posVec,velVec,0);
        logNeighborList(this,0,neighbors);
        logLocalKnownNodes(this,0,knownNodeCoords);
        logDtNeighborSet(this,0,dtNeighbors);

        // schedule next
        scheduleAt(simTime() + keepAliveInterval, keepAliveTimer);
        EV_INFO << "Neighbor::handleMessage - next keepalive scheduled at " << simTime() + keepAliveInterval << "\n";
    }
    else if(msg == KANTimer){
        sendKeepAliveNotificationTo();
        refreshEntrys();

        //rpc
        Coord pos = owner->mobility->getCurrentPosition();
        Coord vel = owner->mobility->getCurrentVelocity();
        std::vector<double> posVec = {pos.x,pos.y,pos.z};
        std::vector<double> velVec = {vel.x,vel.y,vel.z};
        logPosition(this,posVec,velVec,0);
        logNeighborList(this,0,neighbors);
        logLocalKnownNodes(this,0,knownNodeCoords);
        logDtNeighborSet(this,0,dtNeighbors);

        std::vector<L3Address> toDelete;
        for(auto &[v,value] : weaklinkNeighbors){
            if(value.timeout <= simTime()){
                toDelete.push_back(v);
            }
        }
        for(L3Address addr:toDelete){
            weaklinkNeighbors.erase(addr);
            neighbors.erase(addr);
        }
        if(weaklinkNeighbors.empty()){
            presentKANInterval = KANInterval;
        }
        else{
            presentKANInterval = weaklinkKANInterval;
        }
        scheduleAt(simTime() + presentKANInterval, KANTimer);
    }
    else if(msg == dtTimer){
        EV_INFO << "Neighbor::handleMessage - dtTimer triggered, recomputing DT neighbors\n";
        updateDTNeighbors();
        if(debug){
            printDT();
            printDTNeighbors();
            printKnownCoords();
        }
        scheduleAt(simTime() + dtInterval, dtTimer);
    }
    else if(msg == initiailizationTimer){
        if(owner->getNodeattachstate() == true && owner->getNodeinitializingstate() == true){
            cancelEvent(initiailizationTimer);
            delete initiailizationTimer;
            initiailizationTimer = nullptr;
        }
        else{
            owner->setNodeinitializingstate(false);
        }
        RPC_EV("INFO", "node connect, node=");
    }
    else {
        EV_WARN << "Neighbor::handleMessage - unexpected message: " << msg->getName() << "\n";
        delete msg;
    }
}

Neighbor::~Neighbor()
{
    if (keepAliveTimer) cancelAndDelete(keepAliveTimer);
    if (KANTimer) cancelAndDelete(KANTimer);
    if (dtTimer) cancelAndDelete(dtTimer);
    if (initiailizationTimer) cancelAndDelete(initiailizationTimer);
}

bool Neighbor::existsNeighborInDiameterSphere(const NeighborEntry& entry, L3Address& nxt) const
{
    std::vector<double> selfCoord;
    if (!owner->getSelfCoordinate(selfCoord)) {
        EV_WARN << "existsNeighborInDiameterSphere: cannot get self coordinate\n";
        return false;
    }
    if (entry.coords.empty()) {
        EV_WARN << "existsNeighborInDiameterSphere: entry has no coords\n";
        return false;
    }
    size_t dim = std::min({ selfCoord.size(), entry.coords.size() });
    if (dim == 0) return false;

    // compute core and radius
    std::vector<double> center(dim);
    for (size_t i = 0; i < dim; ++i)
        center[i] = 0.5 * (selfCoord[i] + entry.coords[i]);
    // (7/8)^2 * 0.25 * dist^2
    double radius2 = 0.25 * (49.0 / 64.0) * coordDistanceSquared(selfCoord, entry.coords);
    if (!std::isfinite(radius2)) return false;

    const double EPS = 1e-9;

    // save (addr, neighborEntry, distSelf2, distEntry2, distToSegment2)
    struct Candidate {
        L3Address addr;
        const NeighborEntry* ne;
        double distSelf2;
        double distEntry2;
        double distToSegment2;
    };
    std::vector<Candidate> preferred; // dist(self) < dist(entry)
    std::vector<Candidate> others;    //

    for (const auto &kv : neighbors) {
        const L3Address &nbrAddr = kv.first;
        const NeighborEntry &ne = kv.second;

        // skip entry and self
        if (nbrAddr == entry.addr) continue;
        if (nbrAddr == owner->getSelfIPAddress()) continue;

        if (ne.coords.empty()) continue;
        if (ne.coords.size() < dim) continue;

        // if in the ball
        double d2_center = 0.0;
        for (size_t i = 0; i < dim; ++i) {
            double d = ne.coords[i] - center[i];
            d2_center += d * d;
        }
        if (d2_center > radius2 + EPS) continue; //

        //  dist^2 between self and entry
        double distSelf2 = 0.0;
        double distEntry2 = 0.0;
        for (size_t i = 0; i < dim; ++i) {
            double ds = ne.coords[i] - selfCoord[i];
            double de = ne.coords[i] - entry.coords[i];
            distSelf2 += ds * ds;
            distEntry2 += de * de;
        }

        // dest from point to line : self-entry
        double distToSeg2 = pointSegmentDistanceSquared(ne.coords, selfCoord, entry.coords);

        Candidate c{ nbrAddr, &ne, distSelf2, distEntry2, distToSeg2 };
        if (distSelf2 + EPS < distEntry2) { //
            preferred.push_back(c);
        } else {
            others.push_back(c);
        }
    }

    auto chooseBest = [&](const std::vector<Candidate>& vec)->bool {
        if (vec.empty()) return false;
        //
        const Candidate* best = &vec[0];
        for (size_t i = 1; i < vec.size(); ++i) {
            const Candidate* cur = &vec[i];
            // close to the line better
            if (cur->distToSegment2 + EPS < best->distToSegment2) {
                best = cur;
            } else if (std::abs(cur->distToSegment2 - best->distToSegment2) <= EPS) {
                // close to self better
                if (cur->distSelf2 + EPS < best->distSelf2) {
                    best = cur;
                }
            }
        }
        // write nxt
        nxt = best->addr;
        return true;
    };

    // preferred first
    if (chooseBest(preferred)) return true;
    if (chooseBest(others)) return true;

    return false;
}

void Neighbor::processKeepAliveChunk(const Ptr<KeepAlive>& ka, const L3Address& from)
{
    Enter_Method("processKeepAliveChunk");
    if (!ka) {
        EV_ERROR << "processKeepAliveChunk: null chunk\n";
        return;
    }

    L3Address src = from;
    simtime_t ts = simTime();
    bool attached = ka->getAttached();
    EV_INFO << "Neighbor::processKeepAliveChunk src=" << from << " this node=" << owner->getSelfAddress()<< " ts=" << ts << "\n";

    std::vector<NodeEntry> nodeset = ka->getNeighbors();

    bool other_node_hear_me = false;

    auto &entry = neighbors[src];
    entry.addr = src;
    entry.timeout = simTime() + neighborStaleTimeout;
    entry.attachedToDT = attached;

    // handle coords if present
    int dim = ka->getCoordsArraySize();
    entry.dim = dim;
    entry.coords.clear();
    for (int i = 0; i < dim; ++i) {
        double c = ka->getCoords(i);
        entry.coords.push_back(c);
        EV_INFO << "Neighbor coords "<< i<< ":"<< c <<"\n";
    }

    auto knEntry = entry;
    knEntry.timeout = simTime() + knownNodeStaleTimeout;
    knownNodeCoords[src] = knEntry;

    for (const auto& nodeEntry : nodeset) {
        L3Address nodeAddr = nodeEntry.addr;
        if(nodeAddr == owner->getSelfIPAddress()){
            EV_ERROR<<"the node can hear me and add it to my pneighbors\n";
            other_node_hear_me = true;
        }
        if(knownNodeCoords.find(nodeAddr) == knownNodeCoords.end()){
            EV_ERROR<<"ka message carry new node and add it to knownNodeCoords,new addr: "<<nodeAddr<<"\n";
        }
        NeighborEntry& kncEntry = knownNodeCoords[nodeAddr];

        kncEntry.addr = nodeAddr;
        kncEntry.coords = nodeEntry.coords;
        kncEntry.attachedToDT = nodeEntry.attachedToDT;

        kncEntry.timeout = simTime() + knownNodeStaleTimeout;

        kncEntry.dim = kncEntry.coords.size();
        //route
        //Routetuple tuple("physics",owner->getSelfAddress(),L3Address(),src,nodeAddr,-1,-1,simTime() + owner->activeRouteTimeout);
        //owner->ensureRoute(nodeAddr, src, tuple);
    }

    if(weaklinkNeighbors.count(src)){
        weaklinkNeighbors.erase(src);
    }
    L3Address nxt = L3Address();
    if(existsNeighborInDiameterSphere(entry,nxt)){
        Routetuple tuple("physics",src,L3Address(),nxt,owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);
        owner->ensureRoute(src, nxt, tuple);
    }
    else{
        Routetuple tuple("physics",src,L3Address(),L3Address(),owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);
        owner->ensureRoute(src, src, tuple);
    }
    // route table
    //Routetuple tuple("physics",src,L3Address(),L3Address(),owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);

    //Routetuple tuple_reverse("physics",owner->getSelfAddress(),L3Address(),L3Address(),src,-1,-1,simTime() + owner->activeRouteTimeout);

    //owner->ensureRoute(src, src, tuple);
    //owner->ensureRoute(owner->getSelfAddress(),src,  tuple_reverse);

    //join
    EV_INFO<<"nodestate:"<<owner->getNodeattachstate()<<owner->getNodeinitializingstate()<<"\n";
    if(owner->getNodeattachstate() == false && owner->getNodeinitializingstate() == false){
        if(ka->getAttached()){
            EV_INFO << "Neighbor::processKeepAliveChunk begin join src=" << from << " this node=" << owner->getSelfAddress()<< "\n";
            owner->sendJoinRequest(from);//
            owner->setNodeinitializingstate(true);
            scheduleAt(simTime() + initInterval, initiailizationTimer);
        }
    }
}
void Neighbor::processKeepAliveWeakLinkChunk(const Ptr<KeepAlive>& ka, const L3Address& from){
    Enter_Method("processKeepAliveWeakLinkChunk");
    if (!ka) {
        EV_ERROR << "processKeepAliveWeakLinkChunk: null chunk\n";
        return;
    }
    L3Address src = from;
    simtime_t ts = simTime();
    bool attached = ka->getAttached();
    EV_INFO << "Neighbor::processKeepAliveWeakLinkChunk src=" << from << " this node=" << owner->getSelfAddress()<< " ts=" << ts << "\n";

    std::vector<NodeEntry> nodeset = ka->getNeighbors();

    auto &entry = knownNodeCoords[src];
    entry.addr = src;
    entry.timeout = simTime() + knownNodeStaleTimeout;
    entry.attachedToDT = attached;

    // handle coords if present
    int dim = ka->getCoordsArraySize();
    entry.dim = dim;
    entry.coords.clear();
    for (int i = 0; i < dim; ++i) {
        double c = ka->getCoords(i);
        entry.coords.push_back(c);
        EV_INFO << "Neighbor coords "<< i<< ":"<< c <<"\n";
    }

    for (const auto& nodeEntry : nodeset) {
        L3Address nodeAddr = nodeEntry.addr;
        if(knownNodeCoords.find(nodeAddr) == knownNodeCoords.end()){
            EV_ERROR<<"processKeepAliveWeakLinkChunk: ka message carry new node and add it to knownNodeCoords,new addr: "<<nodeAddr<<"\n";
        }
        NeighborEntry& kncEntry = knownNodeCoords[nodeAddr];

        kncEntry.addr = nodeAddr;
        kncEntry.coords = nodeEntry.coords;
        kncEntry.attachedToDT = nodeEntry.attachedToDT;

        kncEntry.timeout = simTime() + knownNodeStaleTimeout;

        kncEntry.dim = kncEntry.coords.size();
    }
    L3Address nxt = L3Address();
    if(existsNeighborInDiameterSphere(entry,nxt)){
        Routetuple tuple("physics",src,L3Address(),nxt,owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);
        owner->ensureRoute(src, nxt, tuple);
    }
    else{
        if(weaklinkNeighbors.empty()){
            if (KANTimer->isScheduled())
                cancelEvent(KANTimer);
            scheduleAt(simTime() + weaklinkKANInterval, KANTimer);
        }
        auto weEntry = entry;
        weEntry.timeout = simTime() + weaklinkKANInterval * 2;
        weaklinkNeighbors[src] = weEntry;

        auto pEntry = entry;
        pEntry.timeout = simTime() + weaklinkKANInterval * 2;
        neighbors[src] = pEntry;

        Routetuple tuple("physics",src,L3Address(),L3Address(),owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);
        owner->ensureRoute(src, src, tuple);
    }
}
void Neighbor::processKeepAliveNotificationChunk(const Ptr<KeepAliveNotification>& ka, const L3Address& from,bool weaklink){
    Enter_Method("processKeepAliveNotificationChunk");
    if (!ka) {
        EV_ERROR << "processKeepAliveNotificationChunk: null chunk\n";
        return;
    }

    L3Address src = from;
    simtime_t ts = simTime();
    bool attached = ka->getAttached();
    EV_INFO << "Neighbor::processKeepAliveChunk src=" << from << " this node=" << owner->getSelfAddress()<< " ts=" << ts << "\n";

    auto &entry = neighbors[src];
    entry.addr = src;
    entry.timeout = simTime() + neighborStaleTimeout;
    entry.attachedToDT = attached;

    // handle coords if present
    int dim = ka->getCoordsArraySize();
    entry.dim = dim;
    entry.coords.clear();
    for (int i = 0; i < dim; ++i) {
        double c = ka->getCoords(i);
        entry.coords.push_back(c);
        EV_INFO << "Neighbor coords "<< i<< ":"<< c <<"\n";
    }

    auto knEntry = entry;
    knEntry.timeout = simTime() + knownNodeStaleTimeout;
    knownNodeCoords[src] = knEntry;

    //Routetuple tuple("physics",src,L3Address(),L3Address(),owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);

    //owner->ensureRoute(src, src, tuple);
    if(!weaklink && weaklinkNeighbors.count(src)){
        weaklinkNeighbors.erase(src);
    }

    L3Address nxt = L3Address();
    if(existsNeighborInDiameterSphere(entry,nxt)){
        Routetuple tuple("physics",src,L3Address(),nxt,owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);
        owner->ensureRoute(src, nxt, tuple);
        if(weaklink)neighbors.erase(src);
    }
    else{
        if(weaklink){
            if(weaklinkNeighbors.empty()){
                if (KANTimer->isScheduled())
                    cancelEvent(KANTimer);
                scheduleAt(simTime() + weaklinkKANInterval, KANTimer);
            }
            auto weEntry = entry;
            weEntry.timeout = simTime() + weaklinkKANInterval * 2;
            weaklinkNeighbors[src] = weEntry;

            Routetuple tuple("physics",src,L3Address(),L3Address(),owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);
            owner->ensureRoute(src, src, tuple);
        }
        else{
            Routetuple tuple("physics",src,L3Address(),L3Address(),owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);
            owner->ensureRoute(src, src, tuple);
        }
    }
}
void Neighbor::sendKeepAliveTo(const L3Address& dest)
{
    if (!owner) {
        EV_ERROR << "Neighbor::sendKeepAliveTo: owner is null, dropping keepalive to " << dest << "\n";
        return;
    }

    auto ka = makeShared<KeepAlive>(); // Ptr<KeepAlive>


    ka->setPacketType(MDTControlPacketType::KA);
    //ka->setSourceAddr(owner->getSelfIPAddress());
    //ka->setTimestamp(simTime());//???
    //int dim = 3;
    //ka->setDim(dim);
    ka->setAttached(owner->getNodeattachstate());

    Coord pos = owner->mobility->getCurrentPosition();
    ka->setCoordsArraySize(3);
    ka->setCoords(0, pos.x);
    ka->setCoords(1, pos.y);
    ka->setCoords(2, pos.z);

    std::vector<NodeEntry> entrys;
    for (const auto& pair : neighbors) {
        const NeighborEntry& neighbor = pair.second;
        NodeEntry nodeEntry;

        nodeEntry.addr = neighbor.addr;
        nodeEntry.attachedToDT = neighbor.attachedToDT;

        nodeEntry.coords.reserve(3);
        if (neighbor.coords.size() >= 3) {
            nodeEntry.coords = {neighbor.coords[0], neighbor.coords[1], neighbor.coords[2]};
        } else {
            nodeEntry.coords = neighbor.coords;
            while (nodeEntry.coords.size() < 3) {
                nodeEntry.coords.push_back(0.0);
            }
        }
        entrys.push_back(nodeEntry);
    }
    ka->setNeighbors(entrys);
    // === IMPORTANT: compute and set chunk length ===
    // Must match exactly the bytes written by your serializer.
    // In our serializer we write:
    // 1 byte packetType + 4 bytes src IPv4 + 4 bytes timestamp (ms, uint32)
    // + 4 bytes dim (uint32) + 8*dim bytes for coords (double as uint64)
    // + 1 byte attached flag
    int bytes = 0;
    bytes += 1;                // packet type
    //bytes += 4;                // src IPv4
    //bytes += 4;                // timestamp (uint32 ms)
    //bytes += 4;                // dim (uint32)
    int dim = 3;
    bytes += 8 * dim;          // coords[] (each double -> 8 bytes)
    bytes += 1;                // attached

    for (const auto &ne : entrys) {
        bytes += 4;            // addr IPv4
        bytes += 1;            // attached flag (1 byte)
        int ndim = ne.coords.size();
        if (ndim < 0) ndim = 0;
        //bytes += 4;            // dim (uint32)
        bytes += 8 * ndim;     // coords: ndim * 8 bytes (double as uint64)
        //bytes += 8;            // timeout (uint64 ms)
    }

    ka->setChunkLength(inet::B(bytes)); // set chunk length explicitly

    Packet *p = new Packet("mdt-keepalive", ka);

    EV_INFO << "Neighbor::sendKeepAliveTo -> dest=" << dest
            << " pktType=" << ka->getPacketType()
            //<< " src=" << ka->getSourceAddr()
            << " attached=" << ka->getAttached()
            << " coords=" << ka->getCoords(0) << "," << ka->getCoords(1) << "," << ka->getCoords(2)
            << endl;


    EV_INFO << "Neighbor::sendKeepAliveTo sending KeepAlive to " << dest << " at " << simTime() << "\n";
    owner->sendPacketTo(p, dest, 0);
}
void Neighbor::sendKeepAliveNotificationTo(){
    if (!owner) {
        EV_ERROR << "Neighbor::sendKeepAliveNotificationTo: owner is null, dropping KAN " "\n";
        return;
    }
    auto ka = makeShared<KeepAliveNotification>(); // Ptr<KeepAlive>
    ka->setPacketType(MDTControlPacketType::KAN);
    //ka->setSourceAddr(owner->getSelfIPAddress());
    //ka->setTimestamp(simTime());//???
    //int dim = 3;
    //ka->setDim(dim);
    ka->setAttached(owner->getNodeattachstate());

    Coord pos = owner->mobility->getCurrentPosition();
    ka->setCoordsArraySize(3);
    ka->setCoords(0, pos.x);
    ka->setCoords(1, pos.y);
    ka->setCoords(2, pos.z);

    int bytes = 0;
    bytes += 1;                // packet type
    //bytes += 4;                // src IPv4
    //bytes += 4;                // timestamp (uint32 ms)
    //bytes += 4;                // dim (uint32)
    int dim = 3;
    bytes += 8 * dim;          // coords[] (each double -> 8 bytes)
    bytes += 1;                // attached
    ka->setChunkLength(inet::B(bytes)); // set chunk length explicitly

    Packet *p = new Packet("mdt-keepaliveNotification", ka);
    L3Address broadcastAddr = L3Address(inet::Ipv4Address::ALLONES_ADDRESS);
    owner->sendPacketTo(p, broadcastAddr, 0);
}
void Neighbor::refreshEntrys()
{
    std::vector<L3Address> timeoutnodes;
    for (auto it = neighbors.begin(); it != neighbors.end(); ) {
        if (it->second.timeout < simTime()) {
            L3Address removedKey = it->first;      // keep it
            it = neighbors.erase(it);              // then delete
            EV_ERROR << "erase timeout node in pneighbors: " << removedKey << "\n";
            timeoutnodes.push_back(removedKey);    // add to set
        } else {
            ++it;
        }
    }
    //if(!timeoutnodes.empty()){
    EV_ERROR<<"delete expired pneighbors in route table\n";
    owner->deleteexpiredneighbors(timeoutnodes);
    //}
    /*for (auto it = dtNeighbors.begin(); it != dtNeighbors.end(); ) {
        if (it->second.timeout < simTime()) {
            it = dtNeighbors.erase(it);
            EV_INFO<<"erase timeout node in dtneighbors:"<<it->first<<"\n";
        } else {
            ++it;
        }
    }*/
    EV_INFO<<"finish refresh entrys\n";

}

bool Neighbor::isPhysicalNeighbor(const L3Address& addr) const
{
    auto it = neighbors.find(addr);
    //bool is = it != neighbors.end() && (simTime() - it->second.lastHeard) < neighborStaleTimeout;
    bool is = it != neighbors.end() && it->second.timeout > simTime();                    //TODO:time state neighborStaleTimeout or owner->activeRouteTimeout
    EV_DEBUG << "Neighbor::isPhysicalNeighbor(" << addr << ") => " << is << "\n";
    return is;
}


bool Neighbor::getNeighborCoords(const L3Address& addr, std::vector<double>& outCoords) const
{
    auto it = neighbors.find(addr);
    if (it == neighbors.end()) return false;
    const NeighborEntry &e = it->second;
    if (e.dim <= 0) return false;
    outCoords.clear();
    for (int i = 0; i < e.dim; ++i)
        outCoords.push_back(e.coords[i]);
    return true;
}

bool Neighbor::addOrUpdateKnownNode(const std::vector<NeighborEntry>& entrys)
{
    //inet::L3Address selfAddr = owner->getSelfAddress();
    simtime_t now = simTime();
    bool changed = false; // track whether we should return true

    for (const auto &neNew : entrys) {
        const L3Address &addr = neNew.addr;
        if(addr.isUnspecified()){
            continue;
        }
        // If the new entry is expired, ensure local copy is removed and skip
        if (neNew.timeout != SIMTIME_ZERO && neNew.timeout <= now) {
            auto it = knownNodeCoords.find(addr);
            if (it != knownNodeCoords.end()) {
                EV_INFO << "addOrUpdateKnownNode: received expired entry for " << addr
                        << "\n";
                if(it->second.timeout < neNew.timeout){
                    NeighborEntry ne = neNew;
                    ne.addr = addr;
                    knownNodeCoords[addr] = ne;
                    EV_INFO <<"but the entry is newer than the local one ,replace the older one\n";
                    changed = true;
                }
                //knownNodeCoords.erase(it);
            } else {
                knownNodeCoords[addr] = neNew;
                changed = true;
                EV_DETAIL << "addOrUpdateKnownNode: received expired entry for " << addr << " (no local copy)\n";
            }
            continue;
        }

        // Otherwise the new entry is fresh (timeout == SIMTIME_ZERO means 'no timeout' or indefinite)
        auto it = knownNodeCoords.find(addr);
        if (it == knownNodeCoords.end()) {
            // insert new entry
            NeighborEntry ne = neNew; // copy
            // ensure addr field set
            ne.addr = addr;
            // if dim not set but coords present, set dim
            if (ne.dim <= 0 && !ne.coords.empty()) ne.dim = (int)ne.coords.size();
            // If no timeout provided, leave as SIMTIME_ZERO
            knownNodeCoords[addr] = ne;
            EV_INFO << "addOrUpdateKnownNode: inserted new known node " << addr
                    << " attached=" << (ne.attachedToDT ? "true" : "false")
                    << " timeout=" << (ne.timeout==SIMTIME_ZERO ? std::string("SIMTIME_ZERO") : ne.timeout.str())
                    << "\n";

            // per your requirement: insertion of a node that was not present -> return true
            changed = true;
        }
        else {
            // update existing entry: be conservative with timeout, update coords if provided
            NeighborEntry &neLocal = it->second;

            // Save previous attached state to detect false->true transition
            bool prevAttached = neLocal.attachedToDT;

            // timeout: choose the later expiration (max)
            if (neLocal.timeout == SIMTIME_ZERO) {
                // local had no timeout, keep as is unless incoming provides one (we prefer longer)
                if (neNew.timeout != SIMTIME_ZERO)
                    neLocal.timeout = neNew.timeout;
            } else if (neNew.timeout == SIMTIME_ZERO) {
                // incoming says indefinite, adopt it
                //neLocal.timeout = SIMTIME_ZERO;
                EV_ERROR<<"the timeout of new node info is SIMTIME_ZERO";
            } else {
                // both set -> take later
                neLocal.timeout = std::max(neLocal.timeout, neNew.timeout);
            }

            // attached state: take newest value (incoming is considered more recent)
            if(neNew.timeout > neLocal.timeout){
                neLocal.attachedToDT = neNew.attachedToDT;
                // coords/dim: if incoming provides coords, update; otherwise keep existing
                if (!neNew.coords.empty()) {
                    neLocal.coords = neNew.coords;
                    neLocal.dim = neNew.dim > 0 ? neNew.dim : (int)neNew.coords.size();
                }
                changed = true;
            }
            if (!prevAttached && neLocal.attachedToDT) {
                changed = true;
            }

            EV_INFO << "addOrUpdateKnownNode: updated known node " << addr
                      << " attached=" << (neLocal.attachedToDT ? "true" : "false")
                      << " timeout=" << (neLocal.timeout==SIMTIME_ZERO ? std::string("SIMTIME_ZERO") : neLocal.timeout.str())
                      << "\n";
        }
    }

    // Optionally: remove local entries that have expired (in case caller didn't pass them)
    // This keeps knownNodeCoords clean; comment out if you prefer separate garbage collection.
    /*for (auto it = knownNodeCoords.begin(); it != knownNodeCoords.end(); ) {
        if (it->second.timeout != SIMTIME_ZERO && it->second.timeout <= now) {
            EV_INFO << "addOrUpdateKnownNode: local entry expired for " << it->first << " -> erasing\n";
            it = knownNodeCoords.erase(it);
        } else {
            ++it;
        }
    }*/

    return changed;
}
bool Neighbor::addOrUpdateDTNeighbors(const std::vector<NeighborEntry>& entrys){
    //inet::L3Address selfAddr = owner->getSelfAddress();
    simtime_t now = simTime();
    bool changed = false; // track whether we should return true

    for (const auto &neNew : entrys) {
        const L3Address &addr = neNew.addr;
        if(addr.isUnspecified()){
            continue;
        }
        if (neNew.timeout != SIMTIME_ZERO && neNew.timeout <= now) {
            auto it = dtNeighbors.find(addr);
            if (it != dtNeighbors.end()) {
                EV_INFO << "addOrUpdateDTNeighbors: received expired entry for " << addr
                        << "\n";
                //dtNeighbors.erase(it);
            } else {
                EV_INFO << "addOrUpdateKnownNode: received expired entry for " << addr << " (no local copy)\n";
            }
            continue;
        }

        auto it = dtNeighbors.find(addr);
        if (it == dtNeighbors.end()) {
            // insert new entry
            NeighborEntry ne = neNew; // copy
            // ensure addr field set
            ne.addr = addr;
            // if dim not set but coords present, set dim
            if (ne.dim <= 0 && !ne.coords.empty()) ne.dim = (int)ne.coords.size();
            // If no timeout provided, leave as SIMTIME_ZERO
            dtNeighbors[addr] = ne;
            EV_INFO << "addOrUpdatedtNeighbors: inserted new known node " << addr
                    << " attached=" << (ne.attachedToDT ? "true" : "false")
                    << " timeout=" << (ne.timeout==SIMTIME_ZERO ? std::string("SIMTIME_ZERO") : ne.timeout.str())
                    << "\n";
            changed = true;
        }
        else {
            NeighborEntry &neLocal = it->second;

            // Save previous attached state to detect false->true transition
            bool prevAttached = neLocal.attachedToDT;

            if (neLocal.timeout == SIMTIME_ZERO) {
                if (neNew.timeout != SIMTIME_ZERO)
                    neLocal.timeout = neNew.timeout;
            } else if (neNew.timeout == SIMTIME_ZERO) {
                //neLocal.timeout = SIMTIME_ZERO;
                EV_ERROR<<"the timeout of new node info is SIMTIME_ZERO";
            } else {
                neLocal.timeout = std::max(neLocal.timeout, neNew.timeout);
            }

            if(neNew.timeout > neLocal.timeout){
                neLocal.attachedToDT = neNew.attachedToDT;
                if (!neNew.coords.empty()) {
                    neLocal.coords = neNew.coords;
                    neLocal.dim = neNew.dim > 0 ? neNew.dim : (int)neNew.coords.size();
                }
                changed = true;
            }

            if (!prevAttached && neLocal.attachedToDT) {
                 changed = true;
             }

            EV_INFO << "addOrUpdatedtNeighbors: updated known node " << addr
                      << " attached=" << (neLocal.attachedToDT ? "true" : "false")
                      << " timeout=" << (neLocal.timeout==SIMTIME_ZERO ? std::string("SIMTIME_ZERO") : neLocal.timeout.str())
                      << "\n";
        }
    }

    /*for (auto it = dtNeighbors.begin(); it != dtNeighbors.end(); ) {
        if (it->second.timeout != SIMTIME_ZERO && it->second.timeout <= now) {
            EV_INFO << "addOrUpdatedtNeighbors: local entry expired for " << it->first << " -> erasing\n";
            it = dtNeighbors.erase(it);
        } else {
            ++it;
        }
    }*/
    return changed;
}

bool Neighbor::getDTNeighborCoords(const L3Address& addr, std::vector<double>& outCoords) const
{
    auto it = dtNeighbors.find(addr);
    if (it == dtNeighbors.end()) return false;
    outCoords = it->second.coords;
    return true;
}
bool Neighbor::getKnownNodeCoords(const L3Address& addr, std::vector<double>& outCoords) const
{
    auto it = knownNodeCoords.find(addr);
    if (it == knownNodeCoords.end()) return false;
    outCoords = it->second.coords;
    return true;
}

Point3 Neighbor::makePoint(const std::vector<double>& c) {
    return Point3(c.size()>0?c[0]:0.0,
                  c.size()>1?c[1]:0.0,
                  c.size()>2?c[2]:0.0);
}

void Neighbor::insertNodeToDT(const inet::L3Address& addr, const Point3& p) {
    dtManager->insertNodeToDT(addr,p);
}

//
// Remove an address from DT:
// - Remove addr -> vh mapping
// - Remove addr from vhToAddrs[vh]
// - If that vh no longer has any addresses, remove vh from DT and erase mapping
//
void Neighbor::removeNodeFromDT(const inet::L3Address& addr) {
    dtManager->removeNodeFromDT(addr,pendingRebuildThreshold);
}
// --- 2) rebuild function: full, safe dt rebuild from knownNodeCoords ---
void Neighbor::rebuildDTFromKnownNodes() {
    // When you want to process pending rebuilds but obey minRebuildInterval and simTime(),
    // call dtManager->rebuildDTFromKnownNodes with a predicate to only include non-expired & attached nodes:
    auto includePred = [this](const NeighborEntry &ne) -> bool {
        // include only attachedToDT and not expired
        if (!ne.attachedToDT) return false;
        if (ne.timeout != SIMTIME_ZERO && ne.timeout <= simTime()) return false;
        return true;
    };
    dtManager->rebuildDTFromKnownNodes(includePred);
}

//
// Rebuild the DT neighbor cache for self:
// - enumerate incident cells around self vertex (3D safe using incident_cells circulator)
// - collect unique neighbor vertex handles, expand them to addresses (vhToAddrs)
// - fill dtNeighbors with knownNodeCoords entries only for attachedToDT and not expired
//
void Neighbor::rebuildDTNeighborCache() {
    inet::L3Address selfAddr = owner->getSelfAddress();
    auto includePred = [this](const NeighborEntry &ne) -> bool {
        // include only attachedToDT and not expired
        if (!ne.attachedToDT) return false;
        if (ne.timeout != SIMTIME_ZERO && ne.timeout <= simTime()) return false;
        return true;
    };
    dtManager->rebuildDTNeighborCache(selfAddr, includePred);
}

//
// Full initial construction of DT from knownNodeCoords (used at first run).
//
void Neighbor::computeDTNeighbors_CGAL() {
    // ensure dtManager exists
    if (!dtManager) createDTManager();

    // Ensure self coordinate is present (store into knownNodeCoords if missing)
    inet::L3Address self = owner->getSelfAddress();
    if (knownNodeCoords.find(self) == knownNodeCoords.end()) {
        std::vector<double> c;
        if (owner->getSelfCoordinate(c)) {
            NeighborEntry &me = knownNodeCoords[self];
            me.coords = c;
            me.addr = self;
            me.dim = (int)c.size();
            me.attachedToDT = true;
            me.timeout = SIMTIME_ZERO; // local self has no timeout (was bug: simTime())
            me.type = "self";
            EV_ERROR << "computeDTNeighbors_CGAL: stored self coords and marked attachedToDT\n";
        } else {
            EV_WARN << "computeDTNeighbors_CGAL: cannot obtain self coordinate\n";
        }
    } else {
        // ensure flags are correct if self already present
        auto &me = knownNodeCoords[self];
        me.attachedToDT = true;
        me.addr = self;
        me.timeout = SIMTIME_ZERO;
    }

    // create include predicate: attachedToDT and not expired (or self)
    auto includePred = [this, self](const NeighborEntry &ne) -> bool {
        if (!ne.attachedToDT) return false;
        if (ne.timeout != SIMTIME_ZERO && ne.timeout <= simTime()) return false;
        return true;
    };

    // delegate full rebuild to dtManager (this will repopulate dt and addrToVh from knownNodeCoords)
    dtManager->rebuildDTFromKnownNodes(includePred);

    // After rebuild, ensure self is mapped in addrToVh (defensive)
    if (addrToVh.find(self) == addrToVh.end()) {
        EV_ERROR << "computeDTNeighbors_CGAL: self not present in addrToVh after rebuild -> inserting self directly\n";
        auto it = knownNodeCoords.find(self);
        if (it != knownNodeCoords.end()) {
            Point3 p = makePoint(it->second.coords);
            dtManager->insertNodeToDT(self, p);
        } else {
            EV_ERROR << "computeDTNeighbors_CGAL: self coords not found in knownNodeCoords when trying direct insert\n";
        }
    }

    // rebuild neighbor cache around self; dtManager expects self present in addrToVh
    dtManager->rebuildDTNeighborCache(self, includePred);

    isInitialized = true;
    EV_INFO << "computeDTNeighbors_CGAL: initial dt built, vertices=" << addrToVh.size() << "\n";
}
//
// Incremental update (called periodically or on disturbance).
// If disturbance==true, reset timer (cancelAndDelete previous).
//
// --- 3) updateDTNeighbors: call rebuild handler at start; avoid immediate dt.remove ---
void Neighbor::updateDTNeighbors(bool disturbance) {
    // reset timer as you had
    Enter_Method("updateDTNeighbors");
    if (disturbance) {
        if (dtTimer != nullptr) {
            cancelAndDelete(dtTimer);
            dtTimer = nullptr;
        }
        dtTimer = new cMessage("neighbor-dt-timer");
        scheduleAt(simTime() + dtInterval, dtTimer);
        EV_ERROR<<"update dt owning to disturbance\n";
    }

    // Ensure dtManager exists
    if (!dtManager) createDTManager();

    EV_ERROR<<"handle pending rebuilds (deferred cleanups)\n";
    if (!pendingRebuildAddrs.empty()) {
        // Option A: attempt rebuild immediately (preferred)
        if (simTime() - lastRebuildTime >= minRebuildInterval) {
            // rebuild via dtManager with predicate
            auto includePred = [this](const NeighborEntry &ne) -> bool {
                if (!ne.attachedToDT) return false;
                if (ne.timeout != SIMTIME_ZERO && ne.timeout <= simTime()) return false;
                return true;
            };
            dtManager->rebuildDTFromKnownNodes(includePred);
            EV_ERROR<<"rebuild neighbor cache for self\n";
            dtManager->rebuildDTNeighborCache(owner->getSelfAddress(), includePred);
            pendingRebuildAddrs.clear();
            lastRebuildTime = simTime();
            ++rebuildCount;
            EV_ERROR << "updateDTNeighbors: performed deferred rebuild\n";
        } else {
            EV_DETAIL << "updateDTNeighbors: pending rebuild deferred (min interval)\n";
        }
    }

    if (!isInitialized) {
        EV_ERROR <<"not isInitialized compute dt directly\n";
        computeDTNeighbors_CGAL();
        return;
    }

    // collect expired addresses (we won't call dt.remove for these)
    std::vector<inet::L3Address> expired;
    for (auto &kv : knownNodeCoords) {
        if (kv.second.timeout != SIMTIME_ZERO && kv.second.timeout + 10 <= simTime())//TODO
            expired.push_back(kv.first);
        if(kv.second.addr == owner->getSelfAddress()){
            kv.second.timeout = simTime() + 10;
        }
    }

    // For expired addresses: only remove from addrToVh/vhToAddrs and mark for rebuild
    for (const auto &a : expired) {
        // use the safe removeNodeFromDT above (doesn't call dt.remove)
        removeNodeFromDT(a);
    }

    // create include predicate once
    auto includePred = [this](const NeighborEntry &ne) -> bool {
        if (!ne.attachedToDT) return false;
        if (ne.timeout != SIMTIME_ZERO && ne.timeout <= simTime()) return false;
        return true;
    };

    EV_ERROR <<"update/insert remaining nodes\n";
    for (auto &kv : knownNodeCoords) {
        const inet::L3Address &addr = kv.first;
        NeighborEntry &ne = kv.second;
        if (ne.timeout != SIMTIME_ZERO && ne.timeout <= simTime()) continue;
        if (!ne.attachedToDT) continue;

        Point3 p = makePoint(ne.coords);

        // use dtManager's addr->vh map
        auto &addrToVhRef = dtManager->addrToVhRef();
        if (addrToVhRef.find(addr) != addrToVhRef.end()) {
            Vertex_handle vh = addrToVhRef[addr];
            // if supported, try dt.move; if fails, defer to rebuild
    #if CGAL_VERSION_NR >= 1050500000
            try {
                if (vh->point() != p) dtManager->dtRef().move(vh, p);
            } catch (...) {
                EV_ERROR << "updateDTNeighbors: dt.move failed for " << addr << ", deferring to rebuild\n";
                pendingRebuildAddrs.insert(addr);
            }
    #else
            if (vh->point() != p) pendingRebuildAddrs.insert(addr);
    #endif
        } else {
            // new node -> insert via dtManager
            dtManager->insertNodeToDT(addr, p);
        }
    }

    EV_ERROR <<"finally rebuild neighbor cache for self\n";
    dtManager->rebuildDTNeighborCache(owner->getSelfAddress(), includePred);

    // if pending set big enough and min interval passed, rebuild now
    if (!pendingRebuildAddrs.empty() && (int)pendingRebuildAddrs.size() >= pendingRebuildThreshold) {
        if (simTime() - lastRebuildTime >= minRebuildInterval) {
            EV_WARN << "updateDTNeighbors: pending size reached -> immediate rebuild\n";
            // use dtManager rebuild
            dtManager->rebuildDTFromKnownNodes(includePred);
            dtManager->rebuildDTNeighborCache(owner->getSelfAddress(), includePred);
            pendingRebuildAddrs.clear();
            lastRebuildTime = simTime();
            ++rebuildCount;
        }
    }

    std::vector<L3Address> dtUniqueKeys;
    for (const auto& dtPair : dtNeighbors) {
        const L3Address& key = dtPair.first;
        if (neighbors.find(key) == neighbors.end()) {
            dtUniqueKeys.push_back(key);
        }
    }
    if(!dtUniqueKeys.empty()){
        EV_ERROR<<"find new dt neighbor and route discover\n";
        owner->findroutefornewdt(dtUniqueKeys);
    }

}
// isDTNeighborOfTarget: ensure target is added to local DT and determine whether self is a DT-neighbor of target.
// Return true if: (a) self is in the DT neighbor set of target (same simplex), OR (b) there are no other valid neighbors
// that satisfy attached & not-expired (i.e., effectively only self and target exist).
bool Neighbor::isDTNeighborOfTarget(const std::vector<double>& targetCoord, const L3Address& addr)
{
    inet::L3Address selfAddr = owner->getSelfAddress();
    auto includePred = [this](const NeighborEntry &ne) -> bool {
        // include only attachedToDT and not expired
        if (!ne.attachedToDT) return false;
        if (ne.timeout != SIMTIME_ZERO && ne.timeout <= simTime()) return false;
        return true;
    };
    if (addrToVh.find(selfAddr) == addrToVh.end()) {
        auto kself = knownNodeCoords.find(selfAddr);
        if (kself == knownNodeCoords.end()) {
            std::vector<double> sc;
            if (owner->getSelfCoordinate(sc)) {
                NeighborEntry ne;
                ne.addr = selfAddr;
                ne.coords = sc;
                ne.dim = (int)sc.size();
                ne.attachedToDT = true;
                ne.timeout = SIMTIME_ZERO;
                knownNodeCoords[selfAddr] = ne;
                dtManager->insertNodeToDT(selfAddr, dtManager->makePoint(sc));
            } else {
                EV_ERROR<<"can not get self coords\n";
            }
        } else {
            knownNodeCoords[selfAddr].attachedToDT = true;
            dtManager->insertNodeToDT(selfAddr, dtManager->makePoint(kself->second.coords));
        }
    }
    return dtManager->isDTNeighborOfTarget(targetCoord,addr,selfAddr,includePred);
}

std::vector<NeighborEntry> Neighbor::computeDTNeighborsForCoord(const std::vector<double>& targetCoord, const L3Address& addr)
{
    auto includePred = [this](const NeighborEntry &ne) -> bool {
        // include only attachedToDT and not expired
        if (!ne.attachedToDT) return false;
        if (ne.timeout != SIMTIME_ZERO && ne.timeout <= simTime()) return false;
        return true;
    };
    const L3Address selfAddr = owner->getSelfAddress();
    if (addrToVh.find(selfAddr) == addrToVh.end()) {
        auto kself = knownNodeCoords.find(selfAddr);
        if (kself == knownNodeCoords.end()) {
            std::vector<double> sc;
            if (owner->getSelfCoordinate(sc)) {
                NeighborEntry ne;
                ne.addr = selfAddr;
                ne.coords = sc;
                ne.dim = (int)sc.size();
                ne.attachedToDT = true;
                ne.timeout = SIMTIME_ZERO;
                knownNodeCoords[selfAddr] = ne;
                dtManager->insertNodeToDT(selfAddr, dtManager->makePoint(sc));
            } else {
                EV_ERROR<<"can not get self coords\n";
            }
        } else {
            knownNodeCoords[selfAddr].attachedToDT = true;
            dtManager->insertNodeToDT(selfAddr, dtManager->makePoint(kself->second.coords));
        }
    }
    return dtManager->computeDTNeighborsForCoord(targetCoord,addr,selfAddr,includePred);

}

void Neighbor::computeMinimalNeighborSubset_CGAL_3D() {
    const L3Address selfAddr = owner->getSelfAddress();
    auto includePred = [this](const NeighborEntry &ne) -> bool {
        // include only attachedToDT and not expired
        if (!ne.attachedToDT) return false;
        if (ne.timeout != SIMTIME_ZERO && ne.timeout <= simTime()) return false;
        return true;
    };
    dtManager->computeMinimalNeighborSubset_CGAL_3D(selfAddr, includePred);
}
std::map<L3Address, NeighborEntry> Neighbor::getMinSimplexNodes() {
    this->computeMinimalNeighborSubset_CGAL_3D();
    return minSimplexNodes;
}
void Neighbor::printDT() const
{
    EV_INFO << "=====  DT content (@ " << simTime() << ")  =====\n";
    EV_INFO << "DT vertices=" << dt.number_of_vertices()
            << " cells=" << dt.number_of_finite_cells() << "\n";
    for (const auto& kv : addrToVh) {
        const L3Address& addr = kv.first;
        Vertex_handle vh      = kv.second;
        if (dt.is_infinite(vh)) continue;
        Point3 p = dt.point(vh);
        auto k   = knownNodeCoords.find(addr);
        bool att = (k != knownNodeCoords.end()) ? k->second.attachedToDT : false;
        EV_INFO << "  addr=" << addr
                << " pos=(" << p.x() << "," << p.y() << "," << p.z()
                << ") attached=" << att << "\n";
    }
    EV_INFO << "=====  DT end  =====\n";
}

void Neighbor::printDTNeighbors() const
{
    EV_INFO << "=====  My DT neighbors (@ " << simTime() << ")  =====\n";
    L3Address self = owner->getSelfAddress();
    for (const auto& kv : dtNeighbors) {
        const NeighborEntry& ne = kv.second;
        EV_INFO << "  nbr=" << ne.addr
                << " pos=(" << ne.coords[0] << "," << ne.coords[1] << "," << ne.coords[2]
                << ") timeout=" << ne.timeout
                << (ne.timeout < simTime() ? "  **EXPIRED**" : "") << "\n";
    }
    EV_INFO << "=====  total " << dtNeighbors.size() << "  =====\n";
}

void Neighbor::printKnownCoords() const
{
    EV_INFO << "=====  knownNodeCoords (@ " << simTime() << ")  =====\n";
    for (const auto& kv : knownNodeCoords) {
        const NeighborEntry& ne = kv.second;
        EV_INFO << "  addr=" << ne.addr
                << " pos=(" << ne.coords[0] << "," << ne.coords[1] << "," << ne.coords[2]
                << ") timeout=" << ne.timeout
                << " attached=" << ne.attachedToDT
                << (ne.timeout < simTime() ? "  **EXPIRED**" : "") << "\n";
    }
    EV_INFO << "=====  total " << knownNodeCoords.size() << "  =====\n";
}

void Neighbor::finish()
{
    EV_ERROR << "Neighbor::finish() called. neighbors.size()=" << neighbors.size() << "\n";
    EV_ERROR << "dtneighbor num:"<<dtNeighbors.size()<<"\n";
    for(auto &[key,value]:neighbors){
        EV_ERROR <<"Pneighbor: "<<key<<"attached:"<<value.attachedToDT<<"timeout: "<<value.timeout<<"\n";
    }
    for(auto &[key,value]:dtNeighbors){
        EV_ERROR <<"dtneighbor: "<<key<<"attached:"<<value.attachedToDT<<"timeout: "<<value.timeout<<"\n";
    }
    for(auto &[key,value]:knownNodeCoords){
        EV_ERROR <<"knownnode: "<<key<<"attached:"<<value.attachedToDT<<"timeout: "<<value.timeout<<"\n";
    }
    for(auto &[key,value]:minSimplexNodes){
        EV_ERROR <<"minSimplexNodes: "<<key<<"attached:"<<value.attachedToDT<<"timeout: "<<value.timeout<<"\n";
    }
}



} // namespace routing
} // namespace mysrc


