/*
 * Neighbor.h
 *
 *  Created on: Aug 31, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_NEIGHBOR_H_
#define MYSRC_ROUTING_NEIGHBOR_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Point_3.h>

#include <map>
#include <vector>
#include <omnetpp.h>
#include "inet/networklayer/common/L3Address.h"
#include "inet/common/INETDefs.h"
#include "inet/common/packet/Packet.h"
#include "MDTData.h"
#include "MDTControlPackets_m.h"

#include"NeighborEntry.h"
#include "NeighborDT.h"
#include "KeepAlive.h"

//rpc
#include"../RPC/log2vis.h"

namespace mysrc {
namespace routing {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point3;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay3;
using Vertex_handle = Delaunay3::Vertex_handle;
using Cell_handle = Delaunay3::Cell_handle;

using namespace inet;

class MDTRouting; // forward




class INET_API Neighbor : public cSimpleModule
{
  protected:
    MDTRouting *owner = nullptr;
    cRNG *rng = nullptr;
    std::map<L3Address, NeighborEntry> neighbors;
    std::map<L3Address, NeighborEntry> weaklinkNeighbors;
    cMessage *keepAliveTimer = nullptr;
    cMessage *KANTimer = nullptr;
    cMessage *dtTimer = nullptr;
    simtime_t dtInterval = 2;
    cMessage *initiailizationTimer = nullptr;
    simtime_t initInterval = 5;

    std::map<L3Address, NeighborEntry> knownNodeCoords;
    std::map<L3Address, NeighborEntry> dtNeighbors;
    std::map<L3Address, NeighborEntry> minSimplexNodes;

    simtime_t keepAliveInterval = 1;      // default: 5s ; will be set from NED/par
    simtime_t KANInterval = 0.2;
    simtime_t weaklinkKANInterval = 0.15;
    simtime_t presentKANInterval = 0.2;
    simtime_t neighborStaleTimeout = 0.5;
    simtime_t knownNodeStaleTimeout = 15; // default: 15s
    simtime_t tupleTimeout = 10;          // default: 20s

    std::unique_ptr<NeighborDT> dtManager;

    //dt
    Delaunay3 dt;
    std::map<Delaunay3::Vertex_handle, std::set<inet::L3Address>> vhToAddrs;
    std::map<inet::L3Address, Delaunay3::Vertex_handle> addrToVh;
    bool isInitialized = false;
    //safe remove node
    /*std::set<inet::L3Address> pendingRemovals;
    int removeFailures = 0;
    const int removeFailuresThreshold = 3;
    const int pendingRemovalsThreshold = 5;*/

    std::set<inet::L3Address> pendingRebuildAddrs; // addresses that triggered removal/move but deferred
    int rebuildCount = 0;
    int rebuildThreshold = 1; // optional: rebuild every time by default; tune for perf
    int pendingRebuildThreshold = 4; // if pending > this, trigger immediate rebuild
    simtime_t lastRebuildTime = SIMTIME_ZERO;
    simtime_t minRebuildInterval = 0.5; // seconds; don't rebuild more often than this

    void insertNodeToDT(const inet::L3Address& addr, const Point3& p);
    void removeNodeFromDT(const inet::L3Address& addr);
    //void handlePendingRemovalsAndMaybeRebuild();
    void rebuildDTNeighborCache();
    void rebuildDTFromKnownNodes();
    Point3 makePoint(const std::vector<double>& c);

  protected:
    virtual void initialize(int stage) override;
    virtual void handleMessage(cMessage *msg) override;
    virtual void finish() override;

  public:
    Neighbor() {}
    virtual ~Neighbor();

    // Create dtManager from this object's containers (call in initialize() or tests)
    void createDTManager() {
        dtManager = std::make_unique<NeighborDT>(knownNodeCoords,
                                                 dt,
                                                 vhToAddrs,
                                                 addrToVh,
                                                 pendingRebuildAddrs,
                                                 dtNeighbors,
                                                 minSimplexNodes);
    }

    // Accessor for tests / other code
    NeighborDT* getDTManager() { return dtManager.get(); }

    //debug
    bool debug = false;
    void printDT() const;
    void printDTNeighbors() const;
    void printKnownCoords() const;

    /** Called by MDTRouting when a control packet of type NEIGHBOR is received */

    /** Called periodically to broadcast / unicast keepalive messages via owner */
    void sendKeepAliveTo(const L3Address& dest);
    void sendKeepAliveNotificationTo();

    /** Refresh or create tuple entry */
    void refreshEntrys();

    bool isPhysicalNeighbor(const L3Address& addr) const;
    /*std::vector<L3Address> getPhysicalNeighbors() const;*/
    bool existsNeighborInDiameterSphere(const NeighborEntry& entry, L3Address& nxt) const;
    void processKeepAliveChunk(const Ptr<KeepAlive>& keepAliveChunk, const L3Address& from);
    void processKeepAliveNotificationChunk(const Ptr<KeepAliveNotification>& KeepAliveNotificationChunk, const L3Address& from,bool weaklink = false);
    void processKeepAliveWeakLinkChunk(const Ptr<KeepAlive>& keepAliveChunk, const L3Address& from);

    bool getNeighborCoords(const L3Address& addr, std::vector<double>& outCoords) const;
    bool getDTNeighborCoords(const L3Address& addr, std::vector<double>& outCoords) const;

    // add or update a known node's coords (called when receiving KA/JR/JRReply etc)
    bool addOrUpdateKnownNode(const std::vector<NeighborEntry>& entrys);
    bool addOrUpdateDTNeighbors(const std::vector<NeighborEntry>& entrys);

    // compute DT neighbors (Gabriel graph approximation)

    void computeDTNeighbors_CGAL();
    void updateDTNeighbors(bool disturbance = false);
    //bool isDTNeighborOfTarget(const std::vector<double>& targetCoord, int dim);
    bool isDTNeighborOfTarget(const std::vector<double>& targetCoord,const L3Address& addr);
    void computeMinimalNeighborSubset_CGAL_3D();

    std::map<L3Address, NeighborEntry> getMinSimplexNodes();

    std::vector<NeighborEntry> computeDTNeighborsForCoord(const std::vector<double>& targetCoord, const L3Address& addr);

    std::map<L3Address, NeighborEntry> getPNeighborNodes() const{return neighbors;}
    // query DT neighbors
    //std::vector<L3Address> getDTNeighbors() const;
    std::map<L3Address, NeighborEntry> getDTNeighbors() const{return dtNeighbors;}

    // query known nodes list
    std::map<L3Address, NeighborEntry> getKnownNodes() const{return knownNodeCoords;}

    std::map<L3Address, NeighborEntry> getweaklinkNeighbors() const{return weaklinkNeighbors;}
    // get coords for a known node
    bool getKnownNodeCoords(const L3Address& addr, std::vector<double>& outCoords) const;





};

} // namespace routing
} // namespace mysrc



#endif /* MYSRC_ROUTING_NEIGHBOR_H_ */
