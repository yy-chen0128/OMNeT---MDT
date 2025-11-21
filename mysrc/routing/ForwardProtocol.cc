/*
 * Forward_Protocol.cc
 *
 *  Created on: Aug 31, 2025
 *      Author: yychen
 */
#include "ForwardProtocol.h"

#include "MDTRouting.h"
#include "ForwardProtocol.h"

#include "Neighbor.h"
#include "MdtRelayTag_m.h"
#include "inet/common/TagBase_m.h"
#include "inet/networklayer/common/L3AddressTag_m.h"
#include "inet/networklayer/ipv4/Ipv4Header_m.h"
#include "inet/networklayer/common/L3Tools.h"

namespace mysrc {
namespace routing {
using namespace inet;

Define_Module(ForwardProtocol);//TODO:short cut

void ForwardProtocol::initialize(int stage)
{
    if (stage == INITSTAGE_LOCAL){
        core = check_and_cast<MDTRouting*>(getModuleByPath("^.core"));
    }
}

void ForwardProtocol::handleMessage(cMessage *msg)
{
    // no timers here
    delete msg;
}
/**
 * Helper: compute Euclidean distance between two coordinate vectors (may be dim N)
 */
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

/**
 * Main forwarding decision per the 6-step table.
 */
ForwardDecision ForwardProtocol::decide(Packet *packet)///not use now
{
    ForwardDecision dec;
    if (!packet || !core) {
        dec.action = ForwardAction::NONE;
        EV_ERROR<<"null packet or core\n";
        dec.reason = "null packet or core";
        return dec;
    }
    const auto& networkHeader = getNetworkProtocolHeader(packet);
    const L3Address& dest = networkHeader->getDestinationAddress();
    const L3Address& src = networkHeader->getSourceAddress();

    if (dest.isUnspecified()) {
        dec.action = ForwardAction::NONE;
        dec.reason = "dest unknown";
        EV_ERROR<<"dest unknown\n";
        return dec;
    }

    L3Address self = core->getSelfAddress();
    if (self == dest) {
        dec.action = ForwardAction::NONE;
        dec.reason = "I am destination";
        EV_ERROR<<"self is the destination\n";//not cover
        return dec;
    }

    // Rule 2: physical neighbor
    auto P = core->getPNeighborNodes();
    for (auto &[key,value] : P) {
        if (key == dest) {
            dec.action = ForwardAction::FORWARD_TO_NEXT_HOP;
            dec.nextHop = key;
            dec.useMdtRelay = false;
            dec.reason = "dest is physical neighbor (rule 2)";
            EV_ERROR<<"Pneighbor is dest:"<<dest<<"\n";
            return dec;
        }
    }

    // Rule 3: check MDT relay tag
    L3Address relayAddr;
    bool relayPresent = false;                          //?????  how this tag be set?
    if (auto rtag = packet->findTag<MdtRelayTag>()) {
        // adapt to your tag API -- example:
        if (rtag->getHasRelay()) {
            EV_ERROR<<"try rule 3\n";
            relayAddr = rtag->getRelayAddr();
            relayPresent = !relayAddr.isUnspecified();
        }
    }

    if (relayPresent && relayAddr != self) {//TODO
        Routetuple t;
        EV_ERROR<<"rule three:MDT forward\n";//not cover
        if (core->findTupleByDestandUpdate(relayAddr, t)) {
            dec.action = ForwardAction::FORWARD_TO_NEXT_HOP;
            dec.nextHop = t.succ;      // next hop is tuple successor
            dec.useMdtRelay = true;
            dec.relayAddr = relayAddr;
            dec.reason = "rule 3: use MDT tuple to relay";
            EV_ERROR<<"rule 3 and relay = "<<relayAddr<<"\n";
            core->updateRoutetuple(src,L3Address(),t.succ,dest);
            return dec;
        }
    }

    // Rules 4 & 5: greedy steps (require coordinates)
    std::vector<double> destCoord;
    bool haveDestCoord = core->getKnownNodeCoordinate(dest, destCoord);//TODO how to find dest coord
    if(!haveDestCoord){
        dec.action = ForwardAction::WAIT_ROUTE_DISCOVERY;
        dec.nextHop = L3Address();
        dec.useMdtRelay = false;
        dec.reason = "wait route discovery";
        return dec;
    }

    EV_INFO<<"haveDestCoord = "<<haveDestCoord<<"\n";//not cover
    if (haveDestCoord) {
        // Rule 4: greedy on physical neighbors
        std::vector<double> selfCoord;
        core->getSelfCoordinate(selfCoord);
        double bestDist = coordDistance(selfCoord, destCoord);
        EV_ERROR<<"dest between self and dest is: "<<bestDist<<"\n";
        L3Address best = L3Address();
        for (auto &[v,value] : P) {
            std::vector<double> vcoord;
            if (!core->getNeighborNodeCoordinate(v, vcoord)) continue;
            double d = coordDistance(vcoord, destCoord);
            if (d < bestDist) {
                bestDist = d; best = v;
                EV_ERROR<<"better dest = "<<best<<"nearer dest = "<<bestDist<<"\n";
            }
        }
        //TODO which strategy should I take?
        /*std::vector<double> selfCoord;
        core->getSelfCoordinate(selfCoord);
        const double selfDist = coordDistance(selfCoord, destCoord);   // self -> dest
        EV_ERROR << "self->dest = " << selfDist << "\n";

        L3Address best = L3Address();
        double bestSelfDist = std::numeric_limits<double>::max();      //

        for (auto &[v, value] : P) {
            std::vector<double> vcoord;
            if (!core->getNeighborNodeCoordinate(v, vcoord)) continue;

            double vDestDist = coordDistance(vcoord, destCoord);       // v -> dest
            if (vDestDist >= selfDist) continue;                       //

            double vSelfDist = coordDistance(vcoord, selfCoord);       // v -> self
            if (vSelfDist < bestSelfDist) {                            //
                bestSelfDist = vSelfDist;
                best = v;
                EV_ERROR << "closer to dest & nearer to self: " << best
                         << " (selfDist=" << vSelfDist << ")\n";
            }
        }*/

        if (!best.isUnspecified() && best != self) {//local minimum
            dec.action = ForwardAction::FORWARD_TO_NEXT_HOP;
            dec.nextHop = best;
            dec.useMdtRelay = false;
            dec.reason = "rule 4 greedy step1";
            core->updateRoutetuple(src,L3Address(),best,dest);
            return dec;
        }

        // Rule 5: DT neighbors / tuples (N_u)
        auto N = core->getDTNeighbors();
        double bestDist2 = 1e300;
        L3Address best2 = L3Address();
        for (auto &v : N) {
            std::vector<double> vcoord;
            if (!core->getDTNeighborNodeCoordinate(v.first, vcoord) || v.first == core->getSelfIPAddress()) continue;
            double d = coordDistance(vcoord, destCoord);
            if (d < bestDist2) {
                bestDist2 = d;
                best2 = v.first;
                EV_ERROR<<"in dt neighbors find better dest = "<<best2<<"nearer dest = "<<bestDist2<<"\n";
            }
        }
        if (!best2.isUnspecified() && best2 != self) {
            Routetuple t;
            if (core->findTupleByDestandUpdate(best2, t)) {
                dec.action = ForwardAction::FORWARD_TO_NEXT_HOP;
                dec.nextHop = t.succ;
                dec.useMdtRelay = true;
                dec.relayAddr = best2;
                dec.reason = "rule 5 greedy step2 via DT tuple";
                EV_ERROR<<"rule 5 and relay = "<<best2<<"\n";
                core->updateRoutetuple(src,L3Address(),t.succ,dest);
                return dec;
            }
            /*else{//TODO
                L3Address nxt = core->findNearestPneighborForTarget(destCoord);
                Routetuple tuple("physics",core->getSelfIPAddress(),L3Address(),nxt,dest,-1,-1,simTime() + core->activeRouteTimeout / 2);
                core->ensureRoute(dest,nxt,tuple);
                dec.action = ForwardAction::FORWARD_TO_NEXT_HOP;
                dec.nextHop = nxt;
                dec.useMdtRelay = true;
                dec.relayAddr = nxt;
                dec.reason = "rule 5 greedy step2 via DT tuple";
                EV_ERROR<<"rule 5 and relay = "<<nxt<<"\n";
                //core->updateRoutetuple(src,L3Address(),t.succ,dest);
                return dec;
            }*/
        }
    }

    // Rule 6: no forwarding needed
    dec.action = ForwardAction::NONE;
    dec.reason = "no forwarding (closest)";
    EV_ERROR<<"error rule 6 is provoked!";
    return dec;
}
bool ForwardProtocol::shortcutforMDTrelay(const L3Address& dest,L3Address& relay,L3Address ban){
    L3Address newrelay = L3Address();
    std::vector<double> destCoord;
    std::vector<double> relayCoord;
    std::vector<double> selfCoord;
    core->getSelfCoordinate(selfCoord);
    L3Address self = core->getSelfIPAddress();
    std::map<L3Address, NeighborEntry> DT = core->getDTNeighbors();
    std::map<L3Address, NeighborEntry> K = core->getKnownNodes();
    for(auto &[v,value] : K){
        if(v == dest){
            if(value.timeout <= simTime()){
                return false;
            }
        }
        if(v == relay){
            if(value.timeout <= simTime()){
                return false;
            }
        }
    }
    if(core->getKnownNodeCoordinate(dest,destCoord) && core->getKnownNodeCoordinate(relay,relayCoord)){
        double bestDist = coordDistance(relayCoord, destCoord);
        L3Address best = L3Address();
        for (auto &[v,value] : DT){
            std::vector<double> vcoord;
            if(value.timeout <= simTime() || v == self || v == ban || !value.attachedToDT)continue;
            if(!core->getDTNeighborNodeCoordinate(v, vcoord))continue;
            double d = coordDistance(vcoord, destCoord);
            if(d < bestDist){
                newrelay = v;
                bestDist = d;
            }
        }
        if (!newrelay.isUnspecified() && newrelay != self){
            relay = newrelay;
            return true;
        }
        else{
            return false;
        }
    }
    else{
        return false;
    }

}
bool ForwardProtocol::findnexthop(const L3Address& dest,L3Address& res,L3Address ban){
    L3Address nxthop = L3Address();
    std::vector<double> destCoord;
    std::vector<double> selfCoord;
    core->getSelfCoordinate(selfCoord);
    L3Address self = core->getSelfIPAddress();
    if(core->getKnownNodeCoordinate(dest,destCoord)){
        double bestDist = coordDistance(selfCoord, destCoord);
        EV_ERROR<<"findnexthop:dest between self and dest is: "<<bestDist<<"\n";
        L3Address best = L3Address();
        auto P = core->getPNeighborNodes();
        for (auto &[v,value] : P) {
            std::vector<double> vcoord;
            if(!value.attachedToDT)continue;
            if (!core->getNeighborNodeCoordinate(v, vcoord) || v == ban) continue;
            double d = coordDistance(vcoord, destCoord);
            if (d < bestDist) {
                bestDist = d; best = v;
                EV_ERROR<<"better dest = "<<best<<"nearer dest = "<<bestDist<<"\n";
            }
        }
        if (!best.isUnspecified() && best != self){
            nxthop = best;
            res = nxthop;
            RPC_EV("INFO", "find nxthop: "<<nxthop);
            return true;
        }
        else{//local min
            RPC_EV("INFO", "local min and dt: ");
            EV_ERROR<<"findnexthop: local min\n";
            double bestDist2 = 1e300;
            L3Address best2 = L3Address();
            auto N = core->getDTNeighbors();
            for (auto &v : N) {

                std::vector<double> vcoord;
                if(!v.second.attachedToDT)continue;
                if (!core->getDTNeighborNodeCoordinate(v.first, vcoord) || v.first == self || v.first == dest || v.first == ban) continue;

                double d = coordDistance(vcoord, destCoord);

                if (d < bestDist2) {
                    bestDist2 = d;
                    best2 = v.first;
                    EV_ERROR<<"findnexthop: in dt neighbors find better dest = "<<best2<<"nearer dest = "<<bestDist2<<"\n";
                }
            }
            if (!best2.isUnspecified() && best2 != self) {//TODO
                res = best2;
            }
            else{
                EV_ERROR<<"findnexthop: error:can not find route for expired route,\n";
            }
            return false;
        }

    }
    else{
        EV_ERROR<<"findnexthop: can not find the coords of the node "<<dest<<"\n";
    }
    return true;
}

}
}
