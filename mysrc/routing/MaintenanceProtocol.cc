/*
 * MDT_Maintenance.cc
 *
 *  Created on: Aug 31, 2025
 *      Author: yychen
 */
#include "MaintenanceProtocol.h"

#include "MDTRouting.h"
#include <algorithm>
#include <unordered_set>
#include <iterator>



namespace mysrc {
namespace routing {
using namespace inet;
Define_Module(MaintenanceProtocol);

void MaintenanceProtocol::initialize(int stage){
    if (stage == INITSTAGE_LOCAL) {
        EV_INFO << "MaintenanceProtocol::initialize(stage=" << stage << ")\n";
        // owner is sibling named "core" in MDT compound module
        cModule *m = getModuleByPath("^.core");
        if (!m) {
            EV_ERROR << "MaintenanceProtocol::initialize: failed to find core module at '^.core'\n";
        } else {
            core = check_and_cast<MDTRouting*>(m);
        }

        //rng
        rng = core->rng;

        MaintenanceInterval = simtime_t(par("MaintenanceInterval").doubleValue());
        EV_INFO << "MaintenanceProtocol params: MaintenanceInterval=" << MaintenanceInterval<<"\n";
        NSRPTimeout =  simtime_t(par("NSRPTimeout").doubleValue());
        EV_INFO << "NSRPTimeout params: NSRPTimeout=" << NSRPTimeout<<"\n";

        MaintenanceTimer = new cMessage("Maintenance-timer");
        scheduleAt(simTime() + uniform(0, 1), MaintenanceTimer);//Interval will change with stage

        this->stage = MaintenanceState::INACTIVE;
    }
}
void MaintenanceProtocol::handleMessage(cMessage *msg){
    if(msg == MaintenanceTimer){
        if(stage == MaintenanceState::INACTIVE){
            scheduleAt(simTime() + uniform(0, 1), MaintenanceTimer);
        }
        else if(stage == MaintenanceState::JOINSTAGE){//check waitingReplyNodes timeout
            simtime_t now = simTime();
            for (auto it = waitingReplyNodes.begin(); it != waitingReplyNodes.end(); ) {//delete TODO: repair or mark this node Invalid
                if (it->second + NSRPTimeout < now)
                    it = waitingReplyNodes.erase(it);
                else
                    ++it;
            }
            if(waitingReplyNodes.empty()){
                EV_ERROR<<"at join stage wait node reply out of date enter CYCLE_MAINTAIN";
                inquiredNodes.clear();
                stage = MaintenanceState::CYCLE_MAINTAIN;
            }
            scheduleAt(simTime() + uniform(0.1, 2), MaintenanceTimer);//TODO adjust
        }
        else if(stage == MaintenanceState::CYCLE_MAINTAIN){
            stage = MaintenanceState::CYCLE_MAINTAIN_PROCESS;
            scheduleAt(simTime() + MaintenanceInterval, MaintenanceTimer);
        }
        else{//MAINTAIN process and check waitingReplyNodes timeout
            if(inquiredNodes.empty()){//TODO begin CYCLE_MAINTAIN
                beginCycleMaintain();
            }
            simtime_t now = simTime();
            for (auto it = waitingReplyNodes.begin(); it != waitingReplyNodes.end(); ) {//delete TODO: repair or mark this node Invalid
                if (it->second + NSRPTimeout < now)
                    it = waitingReplyNodes.erase(it);
                else
                    ++it;
            }
            if(waitingReplyNodes.empty()){
                EV_ERROR<<"wait node reply out of date enter CYCLE_MAINTAIN";
                inquiredNodes.clear();
                stage = MaintenanceState::CYCLE_MAINTAIN;
            }
            scheduleAt(simTime() + uniform(0.1, 2), MaintenanceTimer);//TODO adjust
        }
    }
    else {
        EV_WARN << "MaintenanceProtocol::handleMessage - unexpected message: " << msg->getName() << "\n";
        delete msg;
    }
}

MaintenanceProtocol::~MaintenanceProtocol(){
    if (MaintenanceTimer) cancelAndDelete(MaintenanceTimer);
}
Packet *MaintenanceProtocol::createNeighborSetNotificationTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<double>& coords){
    auto nsn = makeShared<NeighborSetNotification>();

    nsn->setPacketType(MDTControlPacketType::NSN);

    nsn->setSourceAddr(src);

    int dim = 3;
    nsn->setDim(dim);

    nsn->setDestAddr(dest);
    nsn->setPredAddr(pred);

    if(coords.empty()){
        Coord pos = core->mobility->getCurrentPosition();
        nsn->setCoordsArraySize(3);
        nsn->setCoords(0, pos.x);
        nsn->setCoords(1, pos.y);
        nsn->setCoords(2, pos.z);
    }
    else{
        nsn->setCoordsArraySize(coords.size());
        assert(coords.size() == 3);                  //TODO
        nsn->setCoords(0, coords[0]);
        nsn->setCoords(1, coords[1]);
        nsn->setCoords(2, coords[2]);
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
    nsn->setChunkLength(inet::B(bytes)); // set chunk length explicitly

    return new Packet("mdt-neighborSetNotification", nsn);
}
void MaintenanceProtocol::processNeighborSetNotificationChunk(const Ptr<NeighborSetNotification>& nsn){
    Enter_Method("processNeighborSetNotificationChunk");
    if (!nsn) {
        EV_ERROR << "processNeighborSetNotificationChunk null chunk\n";
        return;
    }
    std::vector<double> destcoords;
    int dim = nsn->getDim();
    for (int i = 0; i < dim; ++i) {
        double c = nsn->getCoords(i);
        destcoords.push_back(c);
    }
    if(nsn->getDestAddr() == core->getSelfAddress()){//target node update neighbor info
        L3Address addr = nsn->getDestAddr();
        NeighborEntry entry = {
                addr,
                simTime() + core->activeRouteTimeout,
                true,
                destcoords,
                3,
                "physics"
        };
        std::vector<NeighborEntry> entrys;
        entrys.push_back(entry);
        core->addOrUpdateDTNeighbors(entrys);
    }
    else{//forward node
        Packet* p = createNeighborSetNotificationTo(nsn->getDestAddr(),nsn->getSourceAddr(),core->getSelfAddress(),destcoords);

        Routetuple rt("physics",core->getSelfAddress(),L3Address(),nsn->getPredAddr(),nsn->getSourceAddr(),-1,-1,simTime() + core->activeRouteTimeout);
        if(core->ensureRoute(nsn->getSourceAddr(),nsn->getPredAddr(),rt) == nullptr){
            EV_ERROR<<"processNeighborSetNotificationChunk: Pred is not pneighbor!\n";
        }

        core->sendPacketTo(p,nsn->getDestAddr(),0);
    }
}
void MaintenanceProtocol::beginCycleMaintain(){
    EV_ERROR<<"beginCycleMaintain:\n";
    std::map<L3Address, NeighborEntry> minSet = core->getMinSimplexNodes();

    std::vector<NeighborEntry> entrys;

    for (const auto& kv : minSet){
        entrys.push_back(kv.second);
    }
    std::map<L3Address, NeighborEntry> DT = core->getDTNeighbors();
    std::map<L3Address, NeighborEntry> N = core->getPNeighborNodes();

    std::set<L3Address> existing;
    for (const auto& kv : N)
        existing.insert(kv.first);

    entrys.erase(
        std::remove_if(entrys.begin(), entrys.end(),
                       [&](const NeighborEntry& e)
                       { return existing.count(e.addr); }),
        entrys.end());

    for (const auto& kv : DT) {//get MDT neighbors
        if (N.find(kv.first) == N.end()) {
            entrys.push_back(kv.second);
        }
    }

    if(entrys.empty()){
        EV_ERROR<<"MinSimplexNodes is empty!\n";
    }
    core->sendNSRQ(entrys,L3Address());
    //in this part we can assume that the nodes in entrys is accessible to this node,so these nodes as dest directly

    //NSN message
    /*std::map<L3Address, NeighborEntry> PNeighborSet = core->getPNeighborNodes();

    for (const auto& kv : PNeighborSet) {
        const L3Address& addr = kv.first;
        if (minSet.count(addr) == 0) {
            EV_ERROR<<"send NSN message to"<<addr<<"\n";
            double delay = uniform(0, 1);
            Packet* pkt = createNeighborSetNotificationTo(addr,core->getSelfAddress(),L3Address());
            //core->sendPacketTo(pkt, addr, 0, delay);
        }
    }*/
}
void MaintenanceProtocol::findroutefornewdt(std::vector<NeighborEntry> newdtset){
    std::vector<NeighborEntry> res;
    for(auto entry:newdtset){
        if(waitingReplyNodes.find(entry.addr) == waitingReplyNodes.end()){
            res.push_back(entry);
        }
    }
    if(res.empty()){
        EV_ERROR<<"new dt neighbor route search has started already\n";
        return;
    }
    else{
        EV_ERROR<<"new dt neighbor route search start\n";
        core->sendNSRQ(res,L3Address());
    }

}
Packet *MaintenanceProtocol::createNeighborSetRequestTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<double>& coords){
    auto nsrq = makeShared<NeighborSetRequest>();

    nsrq->setPacketType(MDTControlPacketType::NSRQ);
    nsrq->setSourceAddr(src);

    int dim = 3;
    nsrq->setDim(dim);

    nsrq->setDestAddr(dest);
    nsrq->setPredAddr(pred);

    if(coords.empty()){
        Coord pos = core->mobility->getCurrentPosition();
        nsrq->setCoordsArraySize(3);
        nsrq->setCoords(0, pos.x);
        nsrq->setCoords(1, pos.y);
        nsrq->setCoords(2, pos.z);
    }
    else{
        nsrq->setCoordsArraySize(coords.size());
        assert(coords.size() == 3);                  //TODO
        nsrq->setCoords(0, coords[0]);
        nsrq->setCoords(1, coords[1]);
        nsrq->setCoords(2, coords[2]);
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
    nsrq->setChunkLength(inet::B(bytes)); // set chunk length explicitly

    return new Packet("mdt-neighborSetRequest", nsrq);
}
void MaintenanceProtocol::processNeighborSetRequestChunk(const Ptr<NeighborSetRequest>& nsrq){
    Enter_Method("processNeighborSetRequestChunk");
    if (!nsrq) {
        EV_ERROR << "processNeighborSetRequestChunk null chunk\n";
        return;
    }
    std::vector<double> destcoords;
    int dim = nsrq->getDim();
    for (int i = 0; i < dim; ++i) {
        double c = nsrq->getCoords(i);
        destcoords.push_back(c);
    }
    if(nsrq->getDestAddr() == core->getSelfAddress()){
        EV_ERROR<<"NSRQ to target node:"<<core->getSelfAddress()<<"\n";
        std::vector<NeighborEntry> entrys = core->computeDTNeighborsForCoord(destcoords,nsrq->getSourceAddr());
        Packet* p = createNeighborSetReplyTo(nsrq->getSourceAddr(),core->getSelfIPAddress(),L3Address(),entrys);

        Routetuple rt("physics",core->getSelfAddress(),L3Address(),nsrq->getPredAddr(),nsrq->getSourceAddr(),-1,-1,simTime() + core->activeRouteTimeout);
        if(core->ensureRoute(nsrq->getSourceAddr(),nsrq->getPredAddr(),rt) == nullptr){//TODO:maybe the similar logic in join protocol can be improved
            EV_ERROR<<"processNeighborSetRequestChunk:Pred is not pneighbor!\n";
        }

        core->sendPacketTo(p,nsrq->getSourceAddr(),0);
    }
    else{//forward
        EV_ERROR<<"NSRQ to forward node:"<<nsrq->getDestAddr()<<"\n";
        Packet* p = createNeighborSetRequestTo(nsrq->getDestAddr(),nsrq->getSourceAddr(),core->getSelfAddress(),destcoords);

        Routetuple rt("physics",core->getSelfAddress(),L3Address(),nsrq->getPredAddr(),nsrq->getSourceAddr(),-1,-1,simTime() + core->activeRouteTimeout);
        if(core->ensureRoute(nsrq->getSourceAddr(),nsrq->getPredAddr(),rt) == nullptr){
            EV_ERROR<<"processNeighborSetRequestChunk forward:Pred is not pneighbor!\n";
        }

        core->sendPacketTo(p,nsrq->getDestAddr(),0);
    }
}

Packet *MaintenanceProtocol::createNeighborSetReplyTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<NeighborEntry>& entrys){
    auto nsrp = makeShared<NeighborSetReply>();

    nsrp->setPacketType(MDTControlPacketType::NSRP);
    nsrp->setSourceAddr(src);
    nsrp->setDestAddr(dest);
    nsrp->setPredAddr(pred);

    nsrp->setNeighbors(entrys);
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

    nsrp->setChunkLength(inet::B(bytes)); // set chunk length explicitly

    return new Packet("mdt-neighborsetreply", nsrp);
}
void MaintenanceProtocol::processNeighborSetReplyChunk(const Ptr<NeighborSetReply>& nsrp){
    Enter_Method("processNeighborSetReplyChunk");
    if (!nsrp) {
        EV_ERROR << "processNeighborSetReplyChunk null chunk\n";
        return;
    }

    if(nsrp->getDestAddr() == core->getSelfAddress()){//target node
        //update neighbor table
        std::vector<NeighborEntry> entrys = nsrp->getNeighbors();
        core->addOrUpdateKnownNode(entrys);
        //update set
        waitingReplyNodes.erase(nsrp->getSourceAddr());
        //remove the node in inquiredNodes
        entrys.erase(
            std::remove_if(entrys.begin(), entrys.end(),
                           [this](const NeighborEntry& e) {
                               return this->inquiredNodes.count(e.addr);
                           }),
            entrys.end());
        if(entrys.empty()){//no new neighbor
            if(waitingReplyNodes.empty()){                                      //change stage
                EV_ERROR<<"no new neighbor and change stage to CYCLE_MAINTAIN\n";
                stage = MaintenanceState::CYCLE_MAINTAIN;
                inquiredNodes.clear();
            }
        }
        else { // update (we will insert C1 logic when in maintenance search stage)
            core->addOrUpdateDTNeighbors(entrys); // keep DT neighbor updates as before

            // decide whether we are in maintenance-search stage or join/init stage
            // replace MaintenanceState::MAINTENANCE_SEARCH with your actual enum value
            EV_ERROR<<"new neighbor and update and check C1 condition\n";
            if (stage == MaintenanceState::CYCLE_MAINTAIN_PROCESS) {
                // only send NSRQ to those new nodes that satisfy C1
                std::vector<NeighborEntry> toRequest;
                toRequest.reserve(entrys.size());
                for (const NeighborEntry &ne : entrys) {
                    const L3Address &cand = ne.addr;
                    // skip if already inquired in this maintenance round
                    if (inquiredNodes.count(cand)) continue;

                    // check C1: does there exist a simplex containing u and cand
                    // that does NOT include any node in inquiredNodes?
                    if (satisfiesC1_CGAL(cand, inquiredNodes, core)) {
                        toRequest.push_back(ne);
                    } else {
                        // not satisfying C1 -> defer/ignore this candidate in this round
                        EV_DETAIL << "C1 false for " << cand << " -> not sending NSRQ in this maintenance round\n";
                    }
                }

                if (!toRequest.empty()) {
                    // send NSRQ only to the filtered nodes, and mark them as inquired
                    EV_ERROR<<"find new Neighbor and keep sending NSRQ\n";
                    core->sendNSRQ(toRequest, nsrp->getSourceAddr());
                    //for (const NeighborEntry &ne : toRequest) {
                      //  inquiredNodes.insert(ne.addr);
                        //waitingReplyNodes.insert(ne.addr); // ensure we wait reply from them
                    //}
                } else {
                    // no new node satisfies C1
                    if (waitingReplyNodes.empty()) {
                        EV_ERROR<<"no new node satisfies C1 and change stage to CYCLE_MAINTAIN\n";
                        stage = MaintenanceState::CYCLE_MAINTAIN;
                        inquiredNodes.clear();
                    }
                    // else keep waiting for outstanding replies
                }
            } else {
                // non-maintenance stage (e.g., join/initialization): behave as before
                EV_ERROR<<"in join stage find new neighbor"<<"\n";
                core->addOrUpdateDTNeighbors(entrys);
                core->sendNSRQ(entrys, nsrp->getSourceAddr());
            }
        }
        //else{//update
          //  core->addOrUpdateDTNeighbors(entrys);
            //send nsrq to new nodes and add them to inquiredNodes
            //core->sendNSRQ(entrys,nsrp->getSourceAddr());
        //}
    }
    else if(core->getNodeattachstate()){//forward node
        std::vector<NeighborEntry> entrys = nsrp->getNeighbors();
        Packet* p = createNeighborSetReplyTo(nsrp->getDestAddr(),nsrp->getSourceAddr(),core->getSelfAddress(),entrys);
        EV_ERROR<<"forward NSRP\n";
        IRoute *rt = core->findBestMatchingRoute(nsrp->getDestAddr());
        if(rt == nullptr){
            EV_INFO<<"this node has no route to the new node";
        }
        else{
            MDTData *routingData = check_and_cast<MDTData *>(rt->getProtocolData());
            Routetuple tuple = routingData ->getroutetuple();

            Routetuple tp("physics",nsrp->getSourceAddr(),nsrp->getPredAddr(),tuple.succ,nsrp->getDestAddr(),-1,-1,simTime() + core->activeRouteTimeout);
            EV_ERROR<<"Routetuple:src = "<<nsrp->getSourceAddr()<<"pred = "<<nsrp->getPredAddr()<<"succ = "<<tuple.succ<<"dest = "<<nsrp->getDestAddr()<<"\n";
            Routetuple tp_reverse("physics",nsrp->getDestAddr(),tuple.succ,nsrp->getPredAddr(),nsrp->getSourceAddr(),-1,-1,simTime() + core->activeRouteTimeout);
            EV_ERROR<<"Routetuple_reverse:"<<nsrp->getDestAddr()<<"pred = "<<tuple.succ<<"succ = "<<nsrp->getPredAddr()<<"dest = "<<nsrp->getSourceAddr()<<"\n";
            if(core->ensureRoute(nsrp->getDestAddr(),tuple.succ,tp) == nullptr){
                EV_ERROR<<"processNeighborSetReplyChunk:tuple.succ is not pneighbor!\n";
            }
            if(core->ensureRoute(nsrp->getSourceAddr(),nsrp->getPredAddr(),tp_reverse) == nullptr){
                EV_ERROR<<"processNeighborSetReplyChunk:Pred is not pneighbor!\n";
            }
            //core->ensureRoute(nsrp->getDestAddr(),tuple.succ,tp);
            //core->ensureRoute(nsrp->getSourceAddr(),nsrp->getPredAddr(),tp_reverse);

            core->sendPacketTo(p,nsrp->getDestAddr(),0);
        }
    }
}
bool MaintenanceProtocol::satisfiesC1_CGAL(const L3Address &vAddr,
                         const std::set<L3Address> &requested, MDTRouting *core)
{
    const double EPS = 1e-9;
    EV_DETAIL << "satisfiesC1_CGAL: begin for vAddr=" << vAddr << "\n";

    // 1) self coord
    std::vector<double> selfCoord;
    if (!core->getSelfCoordinate(selfCoord) || selfCoord.size() != 3) {
        EV_ERROR << "satisfiesC1_CGAL: failed to get self coord\n";
        return false;
    }
    Point3 selfP(selfCoord[0], selfCoord[1], selfCoord[2]);

    // 2) known nodes
    const std::map<L3Address, NeighborEntry> &known = core->getKnownNodes();

    // 3) build local 3D Delaunay triangulation from known nodes + self
    Delaunay3 dt;
    std::map<Delaunay3::Vertex_handle, std::set<L3Address>> vhToAddrs;

    // insert self
    Delaunay3::Vertex_handle selfVh;
    try {
        selfVh = dt.insert(selfP);
    } catch (const std::exception &ex) {
        EV_ERROR << "satisfiesC1_CGAL: CGAL exception inserting self: " << ex.what() << "\n";
        return false;
    } catch (...) {
        EV_ERROR << "satisfiesC1_CGAL: unknown exception inserting self\n";
        return false;
    }
    vhToAddrs[selfVh].insert(core->getSelfAddress());

    // insert known neighbors (only those with 3D coords and attachedToDT)
    for (const auto &kv : known) {
        const L3Address &addr = kv.first;
        const NeighborEntry &ne = kv.second;
        if (!ne.attachedToDT) continue;
        if (ne.coords.size() != 3) continue;
        Point3 p(ne.coords[0], ne.coords[1], ne.coords[2]);
        try {
            Delaunay3::Vertex_handle vh = dt.insert(p);
            vhToAddrs[vh].insert(addr);
        } catch (const std::exception &ex) {
            EV_WARN << "satisfiesC1_CGAL: CGAL exception inserting known node " << addr << ": " << ex.what() << "\n";
            // continue inserting others; we won't abort whole function for a single insert
        } catch (...) {
            EV_WARN << "satisfiesC1_CGAL: unknown exception inserting known node " << addr << "\n";
        }
    }

    // ensure we have vAddr coordinate
    auto itv = known.find(vAddr);
    if (itv == known.end() || itv->second.coords.size() != 3) {
        EV_DETAIL << "satisfiesC1_CGAL: no coords for vAddr=" << vAddr << " -> false\n";
        return false;
    }
    Point3 vP(itv->second.coords[0], itv->second.coords[1], itv->second.coords[2]);

    // find the vertex handle for vP in dt (nearest_vertex, guarded)
    Delaunay3::Vertex_handle vVh;
    try {
        vVh = dt.nearest_vertex(vP);
    } catch (const std::exception &ex) {
        EV_ERROR << "satisfiesC1_CGAL: CGAL exception on nearest_vertex: " << ex.what() << "\n";
        return false;
    } catch (...) {
        EV_ERROR << "satisfiesC1_CGAL: unknown exception on nearest_vertex\n";
        return false;
    }
    if (vVh == Delaunay3::Vertex_handle()) {
        EV_WARN << "satisfiesC1_CGAL: nearest_vertex returned infinite for vAddr\n";
        return false;
    }
    // small coordinate tolerance check
    auto close3 = [&](const Point3 &a, const Point3 &b){
        return (std::fabs(a.x()-b.x()) < EPS && std::fabs(a.y()-b.y()) < EPS && std::fabs(a.z()-b.z()) < EPS);
    };
    Point3 foundP = vVh->point();
    if (!close3(foundP, vP)) {
        // try to locate by vhToAddrs mapping (fallback)
        bool located = false;
        for (const auto &kv : vhToAddrs) {
            if (close3(kv.first->point(), vP)) { vVh = kv.first; located = true; break; }
        }
        if (!located) {
            EV_WARN << "satisfiesC1_CGAL: cannot locate vertex for vAddr by coords -> false\n";
            return false;
        }
    }

    // now check triangulation dimension and branch safely
    int dim = dt.dimension();
    EV_DETAIL << "satisfiesC1_CGAL: dt.dimension()=" << dim
              << ", vertices=" << dt.number_of_vertices() << ", cells=" << dt.number_of_cells() << "\n";

    // --- 3D path: use incident_cells safely ---
    if (dim == 3) {
        if (dt.number_of_vertices() < 4) {
            EV_WARN << "satisfiesC1_CGAL: too few vertices for 3D -> false\n";
            return false;
        }
        if (dt.is_infinite(selfVh)) {
            EV_WARN << "satisfiesC1_CGAL: selfVh infinite -> false\n";
            return false;
        }

        std::vector<Delaunay3::Cell_handle> incidentCells;
        try {
            dt.incident_cells(selfVh, std::back_inserter(incidentCells));
        } catch (const std::exception &ex) {
            EV_ERROR << "satisfiesC1_CGAL: CGAL exception incident_cells: " << ex.what() << "\n";
            return false;
        } catch (...) {
            EV_ERROR << "satisfiesC1_CGAL: unknown exception incident_cells\n";
            return false;
        }

        for (auto ch : incidentCells) {
            if (dt.is_infinite(ch)) continue;
            bool containsV = false;
            for (int i = 0; i < 4; ++i) if (ch->vertex(i) == vVh) { containsV = true; break; }
            if (!containsV) continue;

            // check other verts in this cell (exclude self and v)
            bool anyRequested = false;
            for (int i = 0; i < 4; ++i) {
                Delaunay3::Vertex_handle vh = ch->vertex(i);
                if (vh == selfVh || vh == vVh) continue;
                auto mit = vhToAddrs.find(vh);
                if (mit != vhToAddrs.end()) {
                    for (const auto &a : mit->second) {
                        if (requested.count(a)) { anyRequested = true; break; }
                    }
                    if (anyRequested) break;
                } else {
                    // fallback match by coordinate to known set (tolerant)
                    Point3 pt = vh->point();
                    for (const auto &kv : known) {
                        if (kv.second.coords.size() == 3 &&
                            std::fabs(kv.second.coords[0] - pt.x()) < EPS &&
                            std::fabs(kv.second.coords[1] - pt.y()) < EPS &&
                            std::fabs(kv.second.coords[2] - pt.z()) < EPS) {
                            if (requested.count(kv.first)) { anyRequested = true; break; }
                        }
                    }
                    if (anyRequested) break;
                }
            } // per-vertex
            if (!anyRequested) {
                EV_DETAIL << "satisfiesC1_CGAL: found 3D simplex satisfying C1\n";
                return true;
            }
        } // per cell
        EV_DETAIL << "satisfiesC1_CGAL: no 3D simplex satisfied C1\n";
        return false;
    }

    // --- 2D path: COPLANAR CASE: use 2D triangulation + vertex circulator ---
    if (dim == 2) {
        EV_DETAIL << "satisfiesC1_CGAL: entering 2D fallback\n";

        // collect finite vertices as references and to map to 2D DT
        std::vector<Point3> finitePts;
        std::vector<Delaunay3::Vertex_handle> finiteVhs;
        for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
            finitePts.push_back(vit->point());
            finiteVhs.push_back(vit);
            if (finitePts.size() >= 4) break; // a few refs are enough
        }
        if (finitePts.size() < 3) {
            EV_WARN << "satisfiesC1_CGAL: not enough finite points for 2D projection -> false\n";
            return false;
        }

        // build orthonormal basis (u,v) of plane using first 3 non-collinear finitePts
        auto toArr = [&](const Point3 &p, double out[3]) {
            out[0] = static_cast<double>(p.x()); out[1] = static_cast<double>(p.y()); out[2] = static_cast<double>(p.z());
        };
        double p0[3], p1[3], p2[3];
        toArr(finitePts[0], p0); toArr(finitePts[1], p1); toArr(finitePts[2], p2);
        auto sub = [](const double a[3], const double b[3], double out[3]) { out[0] = a[0]-b[0]; out[1] = a[1]-b[1]; out[2] = a[2]-b[2]; };
        auto cross = [](const double a[3], const double b[3], double out[3]) {
            out[0] = a[1]*b[2] - a[2]*b[1];
            out[1] = a[2]*b[0] - a[0]*b[2];
            out[2] = a[0]*b[1] - a[1]*b[0];
        };
        auto dot = [](const double a[3], const double b[3]) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; };
        auto norm = [&](const double a[3]) { return std::sqrt(dot(a,a)); };
        auto normalize = [&](double a[3]) { double n = norm(a); if (n>0) { a[0]/=n; a[1]/=n; a[2]/=n; } };

        double u[3], tmpv[3], nvec[3], v[3];
        sub(p1,p0,u); normalize(u);
        sub(p2,p0,tmpv);
        cross(u,tmpv,nvec);
        if (norm(nvec) < 1e-12) {
            EV_WARN << "satisfiesC1_CGAL: refs nearly collinear -> cannot build plane -> false\n";
            return false;
        }
        normalize(nvec); cross(nvec,u,v); normalize(v);

        // project function
        auto project_to_2d = [&](const Point3 &pp, double &x, double &y){
            double pv[3] = {static_cast<double>(pp.x()), static_cast<double>(pp.y()), static_cast<double>(pp.z())};
            double w[3]; sub(pv, p0, w);
            x = dot(w, u);
            y = dot(w, v);
        };

        // build 2D Delaunay and mapping from 2D vertex handles to addresses
        typedef CGAL::Delaunay_triangulation_2<K> D2;
        typedef K::Point_2 P2;
        typedef D2::Vertex_handle VH2;
        D2 dt2;
        std::map<VH2, std::set<L3Address>> vh2addrs2;
        std::map<L3Address, VH2> addr2vh2;

        // insert every finite vertex from dt into dt2
        for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
            Point3 p3 = vit->point();
            double x,y; project_to_2d(p3, x, y);
            VH2 vh2 = dt2.insert(P2(x,y));
            // map addresses: vhToAddrs may have multiple addresses per 3D vh
            auto itmap = vhToAddrs.find(vit);
            if (itmap != vhToAddrs.end()) {
                for (const auto &a : itmap->second) {
                    addr2vh2[a] = vh2;
                    vh2addrs2[vh2].insert(a);
                }
            } else {
                // try mapping by coordinate against known set (fallback)
                for (const auto &kv : known) {
                    if (kv.second.coords.size() == 3 &&
                        std::fabs(kv.second.coords[0]-p3.x())<EPS &&
                        std::fabs(kv.second.coords[1]-p3.y())<EPS &&
                        std::fabs(kv.second.coords[2]-p3.z())<EPS) {
                        addr2vh2[kv.first] = vh2;
                        vh2addrs2[vh2].insert(kv.first);
                    }
                }
            }
        }

        // ensure self and vAddr present in addr2vh2
        double sx, sy, vx, vy;
        project_to_2d(selfP, sx, sy);
        project_to_2d(vP, vx, vy);
        if (addr2vh2.find(core->getSelfAddress()) == addr2vh2.end()) {
            VH2 vh2 = dt2.insert(P2(sx,sy));
            addr2vh2[core->getSelfAddress()] = vh2;
            vh2addrs2[vh2].insert(core->getSelfAddress());
        }
        if (addr2vh2.find(vAddr) == addr2vh2.end()) {
            VH2 vh2 = dt2.insert(P2(vx,vy));
            addr2vh2[vAddr] = vh2;
            vh2addrs2[vh2].insert(vAddr);
        }

        auto itSelf2 = addr2vh2.find(core->getSelfAddress());
        auto itV2 = addr2vh2.find(vAddr);
        if (itSelf2 == addr2vh2.end() || itV2 == addr2vh2.end()) {
            EV_WARN << "satisfiesC1_CGAL: missing self or v in 2D mapping -> false\n";
            return false;
        }
        VH2 self2 = itSelf2->second;
        VH2 v2 = itV2->second;

        // use vertex circulator to enumerate neighboring vertices around self2
        typedef D2::Vertex_circulator VC;
        std::vector<VH2> neigh2;
        VC vcirc = dt2.incident_vertices(self2);
        if (vcirc != VC()) {
            VC done = vcirc;
            do {
                if (!dt2.is_infinite(vcirc)) neigh2.push_back(vcirc);
                ++vcirc;
            } while (vcirc != done);
        }

        if (neigh2.empty()) {
            EV_DETAIL << "satisfiesC1_CGAL: no 2D incident neighbors -> false\n";
            return false;
        }

        // each consecutive pair forms triangle (self2, neigh2[i], neigh2[i+1])
        for (size_t i = 0; i < neigh2.size(); ++i) {
            VH2 a = neigh2[i];
            VH2 b = neigh2[(i+1) % neigh2.size()];
            bool aIsV = (a == v2);
            bool bIsV = (b == v2);
            if (!(aIsV || bIsV)) continue;
            VH2 other = aIsV ? b : a;

            bool otherRequested = false;
            auto itOtherAddrs = vh2addrs2.find(other);
            if (itOtherAddrs != vh2addrs2.end()) {
                for (const auto &addr : itOtherAddrs->second) {
                    if (requested.count(addr)) { otherRequested = true; break; }
                }
            }
            if (!otherRequested) {
                EV_DETAIL << "satisfiesC1_CGAL: found 2D triangle satisfying C1\n";
                return true;
            }
        }

        EV_DETAIL << "satisfiesC1_CGAL: no 2D triangle satisfied C1\n";
        return false;
    }

    // --- dim <= 1 fallback: cannot check reliably. conservative decision: return false ---
    EV_WARN << "satisfiesC1_CGAL: dt.dimension() <= 1 -> cannot check -> false\n";
    return false;
}





}
}
