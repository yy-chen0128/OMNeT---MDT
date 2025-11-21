/*
 * MDT_Join.cc
 *
 *  Created on: Aug 31, 2025
 *      Author: yychen
 */
#include "JoinProtocol.h"

#include "MDTRouting.h"
#include "inet/networklayer/common/L3AddressTag_m.h"

#include <cmath>
namespace mysrc {
namespace routing {
using namespace inet;
Define_Module(JoinProtocol);

void JoinProtocol::initialize(int stage){
    if (stage == INITSTAGE_LOCAL){
        core = check_and_cast<MDTRouting*>(getModuleByPath("^.core"));
    }
}
void JoinProtocol::handleMessage(cMessage *msg){

}
Packet *JoinProtocol::createJoinRequestTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<double>& coords){
    auto jr = makeShared<JoinRequest>();

    jr->setPacketType(MDTControlPacketType::JRQ);
    jr->setSourceAddr(src);
    int dim = 3;
    jr->setDim(dim);
    jr->setDestAddr(dest);
    jr->setPredAddr(pred);


    if(coords.empty()){
        Coord pos = core->mobility->getCurrentPosition();
        jr->setCoordsArraySize(3);
        jr->setCoords(0, pos.x);
        jr->setCoords(1, pos.y);
        jr->setCoords(2, pos.z);
    }
    else{
        jr->setCoordsArraySize(coords.size());
        assert(coords.size() == 3);                  //TODO
        jr->setCoords(0, coords[0]);
        jr->setCoords(1, coords[1]);
        jr->setCoords(2, coords[2]);
    }

    // === IMPORTANT: compute and set chunk length ===
    // Must match exactly the bytes written by your serializer.
    // In our serializer we write:
    // 1 byte packetType + 4 bytes src IPv4 + 4 bytes timestamp (ms, uint32)
    // + 4 bytes dim (uint32) + 8*dim bytes for coords (double as uint64)
    // + 1 byte attached flag
    int bytes = 0;
    bytes += 1;                // packet type
    bytes += 4;                // src IPv4
    bytes += 4;                 //dest
    bytes += 4;                 //pred

    bytes += 4;                // dim (uint32)
    bytes += 8 * dim;          // coords[] (each double -> 8 bytes)
    jr->setChunkLength(inet::B(bytes)); // set chunk length explicitly

    return new Packet("mdt-joinrequest", jr);
}
bool JoinProtocol::isNearestDTNode(const std::vector<double>& dest,int dim){
    std::map<L3Address, NeighborEntry> knownnodes = core->getKnownNodes();
    std::vector<double> coord;
    core->getSelfCoordinate(coord);

    double selfDist = 0.0;
    for (int i = 0; i < dim; i++) {
        double d = (i < (int)coord.size() ? coord[i] : 0.0) - dest[i];
        selfDist += d * d;
    }
    selfDist = sqrt(selfDist);
    double minDist = selfDist;
    L3Address minNode = core->getSelfAddress();

    for (const auto& kv : knownnodes) {
        const L3Address& addr = kv.first;
        const NeighborEntry& ne = kv.second;

        if (!ne.attachedToDT) continue;//DT node only
        const std::vector<double>& coords = ne.coords;
        if ((int)coords.size() < dim) continue;

        double dist = 0.0;
        for (int i = 0; i < dim; i++) {
            double d = coords[i] - dest[i];
            dist += d * d;
        }
        dist = sqrt(dist);
        if (dist < minDist) {
            minDist = dist;
            minNode = addr;
        }
    }
    return (minNode == core->getSelfAddress());
}
void JoinProtocol::processJoinRequestChunk(const Ptr<JoinRequest>& jrq){
    Enter_Method("processJoinRequestChunk");
    if (!jrq) {
        EV_ERROR << "processJoinRequestChunk: null chunk\n";
        return;
    }

    std::vector<double> destcoords;
    int dim = jrq->getDim();
    for (int i = 0; i < dim; ++i) {
        double c = jrq->getCoords(i);
        destcoords.push_back(c);
    }
    if(jrq->getDestAddr() != jrq->getSourceAddr() && jrq->getDestAddr() == core->getSelfIPAddress()){
        //send to the node who sent the ka message with attached value equaling to true
        EV_INFO<<"the node who sent the ka message with attached value equaling to true:"<<core->getSelfIPAddress()<<"receive jrq message\n";
        std::map<L3Address, NeighborEntry> dts = core->getDTNeighbors();
        for(const auto& [key, value] : dts){
            EV_INFO<<"dtneighbors:"<<key<<"\n";
        }
        if(core->isDTNeighborOfTarget(destcoords,jrq->getSourceAddr())){//shortcut
            EV_ERROR<<"shortcut 01\n";
            std::vector<NeighborEntry> entrys = core->computeDTNeighborsForCoord(destcoords,jrq->getSourceAddr());
            Packet* p = createJoinReplyTo(jrq->getSourceAddr(),core->getSelfIPAddress(),L3Address(),entrys);

            Routetuple rt("physics",core->getSelfAddress(),L3Address(),jrq->getSourceAddr(),jrq->getSourceAddr(),-1,-1,simTime() + core->activeRouteTimeout);
            core->ensureRoute(jrq->getSourceAddr(),jrq->getSourceAddr(),rt);

            core->sendPacketTo(p,jrq->getSourceAddr(),0);
            //Routetuple tuple("physics",src,L3Address(),L3Address(),owner->getSelfAddress(),-1,-1,simTime() + owner->activeRouteTimeout);
        }
        else{
            EV_ERROR<<"the node who sent the ka message with attached value equaling to true forward the message to a more near node\n";
            std::vector<NeighborEntry> entrys = core->computeDTNeighborsForCoord(destcoords,jrq->getSourceAddr());
            //which is the nearer node?
            std::vector<double> selfcoords;
            core->getSelfCoordinate(selfcoords);

            auto nearest = std::min_element(
                entrys.begin(), entrys.end(),
                [&selfcoords](const NeighborEntry& a, const NeighborEntry& b)
                {
                    auto dist2 = [](const std::vector<double>& u,
                                    const std::vector<double>& v)
                    {
                        double d2 = 0.0;
                        for (size_t i = 0; i < u.size(); ++i)
                        {
                            double tmp = u[i] - v[i];
                            d2 += tmp * tmp;
                        }
                        return d2;
                    };
                    return dist2(a.coords, selfcoords) < dist2(b.coords, selfcoords);
                });
            if(nearest != entrys.end()){
                EV_ERROR<<"nexthop: "<<nearest->addr<<"\n";
            }
            else{
                EV_ERROR<<"processJoinRequestChunk:process jrq fail,con not find next hop\n";
            }
            Packet* p = createJoinRequestTo(jrq->getSourceAddr(),jrq->getSourceAddr(),core->getSelfAddress(),destcoords);

            Routetuple rt("physics",core->getSelfAddress(),L3Address(),jrq->getSourceAddr(),jrq->getSourceAddr(),-1,-1,simTime() + core->activeRouteTimeout);
            core->ensureRoute(jrq->getSourceAddr(),jrq->getSourceAddr(),rt);

            core->sendPacketTo(p,jrq->getSourceAddr(),0);
        }
        return;
    }
    else if(jrq->getDestAddr() == jrq->getSourceAddr()){
        if(core->isDTNeighborOfTarget(destcoords,jrq->getSourceAddr())){
            //send reply
            EV_ERROR<<"send join reply to "<<jrq->getSourceAddr()<<"\n";
            EV_ERROR<<"update route pred = "<<jrq->getPredAddr()<<"src = "<<jrq->getSourceAddr()<<"\n";
            std::vector<NeighborEntry> entrys = core->computeDTNeighborsForCoord(destcoords,jrq->getSourceAddr());
            Packet* p = createJoinReplyTo(jrq->getSourceAddr(),core->getSelfIPAddress(),L3Address(),entrys);

            Routetuple rt("physics",core->getSelfAddress(),L3Address(),jrq->getPredAddr(),jrq->getSourceAddr(),-1,-1,simTime() + core->activeRouteTimeout);
            core->ensureRoute(jrq->getSourceAddr(),jrq->getPredAddr(),rt);

            core->sendPacketTo(p,jrq->getSourceAddr(),0);
        }
        else{
            //continue forward and route table update
            EV_ERROR<<"continue forward and route table update "<<jrq->getSourceAddr()<<"\n";
            EV_ERROR<<"update route pred = "<<jrq->getPredAddr()<<"src = "<<jrq->getSourceAddr()<<"\n";
            Packet* p = createJoinRequestTo(jrq->getSourceAddr(),jrq->getSourceAddr(),core->getSelfAddress(),destcoords);

            Routetuple rt("physics",core->getSelfAddress(),L3Address(),jrq->getPredAddr(),jrq->getSourceAddr(),-1,-1,simTime() + core->activeRouteTimeout);
            core->ensureRoute(jrq->getSourceAddr(),jrq->getPredAddr(),rt);

            core->sendPacketTo(p,jrq->getSourceAddr(),0);
        }
        return;
    }
    else{
        EV_ERROR<<"process jrq error,not at stage one or two";
        return;
    }

}
Packet *JoinProtocol::createJoinReplyTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<NeighborEntry>& entrys){
    auto jr = makeShared<JoinReply>();

    jr->setPacketType(MDTControlPacketType::JRP);
    jr->setSourceAddr(src);
    jr->setDestAddr(dest);
    jr->setPredAddr(pred);

    jr->setNeighbors(entrys);
    // === IMPORTANT: compute and set chunk length ===
    // Must match exactly the bytes written by your serializer.
    // In our serializer we write:
    // 1 byte packetType + 4 bytes src IPv4 + 4 bytes timestamp (ms, uint32)
    // + 4 bytes dim (uint32) + 8*dim bytes for coords (double as uint64)
    // + 1 byte attached flag
    int bytes = 0;
    bytes += 1;                // packet type
    bytes += 4;                // src IPv4
    bytes += 4;                 //dest
    bytes += 4;                 //pred
    bytes += 4;                // neighborCount (uint32)

    // for each neighbor: addr(4) + attached(1) + dim(4) + coords(8*dim) + timeout(8)
    for (const auto &ne : entrys) {
        bytes += 4;            // addr IPv4
        bytes += 1;            // attached flag (1 byte)
        int ndim = ne.dim;
        if (ndim < 0) ndim = 0;
        bytes += 4;            // dim (uint32)
        bytes += 8 * ndim;     // coords: ndim * 8 bytes (double as uint64)
        bytes += 8;            // timeout (uint64 ms)
    }

    jr->setChunkLength(inet::B(bytes)); // set chunk length explicitly

    EV_INFO << "JoinProtocol::createJoinReplyTo -> dest=" << dest
            << " pred=" << jr->getPredAddr()
            << " src=" << jr->getSourceAddr()
            << " dest = " << jr->getDestAddr()
            << endl;

    return new Packet("mdt-joinreply", jr);
}

void JoinProtocol::processJoinReplyChunk(const Ptr<JoinReply>& jrp){
    Enter_Method("processJoinReplyChunk");
    if (!jrp) {
        EV_ERROR << "processJoinReplyChunk null chunk\n";
        return;
    }
    EV_INFO<<"receive join reply\n";
    if(jrp->getDestAddr() == core->getSelfAddress()){//target node
        EV_ERROR<<"receive join reply change node state and send NSRQ\n";
        //change node state
        core->setNodeattachstate(true);
        core->setNodeinitializingstate(false);
        //update neighbor table
        std::vector<NeighborEntry> entrys = jrp->getNeighbors();
        core->addOrUpdateKnownNode(entrys);
        //known nodes state?
        std::map<L3Address, NeighborEntry>  known = core->getKnownNodes();
        for(auto [key,value]:known){
            EV_ERROR<<"known nodes: "<<key<<": "<<value.attachedToDT<<"\n";
        }
        for(auto key:entrys){
            EV_ERROR<<"entrys nodes: "<<key.addr<<key.attachedToDT<<"\n";
        }
        core->updateDTNeighbors();
        //begin Neighbor-Set Request
        core->sendNSRQ(entrys,jrp->getSourceAddr());
        //change MaintenanceProtocol stage in MDTRouting
    }
    else if(core->getNodeattachstate()){//forward node
        std::vector<NeighborEntry> entrys = jrp->getNeighbors();
        EV_ERROR<<"receive join reply and forward reply\n";
        Packet* p = createJoinReplyTo(jrp->getDestAddr(),jrp->getSourceAddr(),core->getSelfAddress(),entrys);

        IRoute *rt = core->findBestMatchingRoute(jrp->getDestAddr());
        if(rt == nullptr){
            EV_INFO<<"this node has no route to the new node";
        }
        else{
            MDTData *routingData = check_and_cast<MDTData *>(rt->getProtocolData());
            Routetuple tuple = routingData ->getroutetuple();

            Routetuple tp("physics",jrp->getSourceAddr(),jrp->getPredAddr(),tuple.succ,jrp->getDestAddr(),-1,-1,simTime() + core->activeRouteTimeout);
            EV_INFO<<"Routetuple:src = "<<jrp->getSourceAddr()<<"pred = "<<jrp->getPredAddr()<<"succ = "<<tuple.succ<<"dest = "<<jrp->getDestAddr()<<"\n";
            Routetuple tp_reverse("physics",jrp->getDestAddr(),tuple.succ,jrp->getPredAddr(),jrp->getSourceAddr(),-1,-1,simTime() + core->activeRouteTimeout);
            EV_INFO<<"Routetuple_reverse:"<<jrp->getDestAddr()<<"pred = "<<tuple.succ<<"succ = "<<jrp->getPredAddr()<<"dest = "<<jrp->getSourceAddr()<<"\n";;
            core->ensureRoute(jrp->getDestAddr(),tuple.succ,tp);
            core->ensureRoute(jrp->getSourceAddr(),jrp->getPredAddr(),tp_reverse);

            core->sendPacketTo(p,jrp->getDestAddr(),0);
        }
    }
}





}
}

