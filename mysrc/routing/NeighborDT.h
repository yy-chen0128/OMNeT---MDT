/*
 * NeighborDT.h
 *
 *  Created on: Sep 26, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_NEIGHBORDT_H_
#define MYSRC_ROUTING_NEIGHBORDT_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Point_3.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Point_2.h>
#include <cmath>
#include <sstream>

#include <map>
#include <set>
#include <vector>
#include <functional>
#include <unordered_set>
#include <string>

#include "inet/networklayer/common/L3Address.h"
#include "NeighborEntry.h"

namespace mysrc {
namespace routing {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point3;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay3;
using Vertex_handle = Delaunay3::Vertex_handle;
using Cell_handle = Delaunay3::Cell_handle;

typedef K::Point_2 Point2;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay2;
using Vertex_handle2 = Delaunay2::Vertex_handle;

using namespace inet;

/**
 * NeighborDT - pure C++ class that encapsulates CGAL Delaunay3 operations
 * and the bookkeeping maps used by Neighbor. It DOES NOT depend on OMNeT++
 * time or owner. All collections are *references* provided by the caller.
 *
 * The caller (Neighbor) is responsible for:
 *  - lifetime of the referenced containers
 *  - providing any time-based filtering via includePred (std::function<bool(const NeighborEntry&)>)
 */
class NeighborDT {
public:
    // constructor: takes references to the containers stored in Neighbor
    NeighborDT(std::map<L3Address, NeighborEntry> &knownNodeCoords,
               Delaunay3 &dt,
               std::map<Delaunay3::Vertex_handle, std::set<L3Address>> &vhToAddrs,
               std::map<L3Address, Delaunay3::Vertex_handle> &addrToVh,
               std::set<L3Address> &pendingRebuildAddrs,
               std::map<L3Address, NeighborEntry> &dtNeighbors,
               std::map<L3Address, NeighborEntry> &minSimplexNodes);

    // basic utilities
    Point3 makePoint(const std::vector<double> &c) const;

    // mutating operations (do not perform time checks)
    //void insertNodeToDT(const L3Address &addr, const Point3 &p);
    bool insertNodeToDT(const L3Address &addr, const Point3 &p);
    void removeNodeFromDT(const L3Address &addr, int pendingRebuildThreshold = 1);

    // Rebuild DT completely from knownNodeCoords that satisfy includePred (if provided).
    // includePred returns true for entries that should be included in DT.
    void rebuildDTFromKnownNodes(std::function<bool(const NeighborEntry&)> includePred = nullptr);

    // Rebuild neighbor cache around selfAddr. includePred filters entries.
    void rebuildDTNeighborCache(const L3Address &selfAddr,
                                std::function<bool(const NeighborEntry&)> includePred = nullptr);

    // Compute neighbor entries for a target coordinate/address.
    // selfAddr is used to exclude self, includePred filters which known nodes are considered valid.
    std::vector<NeighborEntry> computeDTNeighborsForCoord(const std::vector<double> &targetCoord,
                                                          const L3Address &addr,
                                                          const L3Address &selfAddr,
                                                          std::function<bool(const NeighborEntry&)> includePred = nullptr);

    // Determine whether selfAddr is a DT neighbor of 'addr' at coordinate targetCoord
    bool isDTNeighborOfTarget(const std::vector<double> &targetCoord,
                              const L3Address &addr,
                              const L3Address &selfAddr,
                              std::function<bool(const NeighborEntry&)> includePred = nullptr);

    // Compute minimal hitting set (minSimplexNodes) for selfAddr
    void computeMinimalNeighborSubset_CGAL_3D(const L3Address &selfAddr,
                                              std::function<bool(const NeighborEntry&)> includePred = nullptr);

    // Accessors to underlying containers (useful for tests)
    std::map<L3Address, Delaunay3::Vertex_handle> &addrToVhRef() { return addrToVhRef_; }
    std::map<Delaunay3::Vertex_handle, std::set<L3Address>> &vhToAddrsRef() { return vhToAddrsRef_; }
    std::set<L3Address> &pendingRebuildRef() { return pendingRebuildRef_; }
    std::map<L3Address, NeighborEntry> &dtNeighborsRef() { return dtNeighborsRef_; }
    std::map<L3Address, NeighborEntry> &minSimplexNodesRef() { return minSimplexNodesRef_; }
    Delaunay3 &dtRef() { return dtRef_; }

private:
    // references to caller-owned containers
    std::map<L3Address, NeighborEntry> &knownNodeCoordsRef_;
    Delaunay3 &dtRef_;
    std::map<Delaunay3::Vertex_handle, std::set<L3Address>> &vhToAddrsRef_;
    std::map<L3Address, Delaunay3::Vertex_handle> &addrToVhRef_;
    std::set<L3Address> &pendingRebuildRef_;
    std::map<L3Address, NeighborEntry> &dtNeighborsRef_;
    std::map<L3Address, NeighborEntry> &minSimplexNodesRef_;
};

} // namespace routing
} // namespace mysrc


#endif /* MYSRC_ROUTING_NEIGHBORDT_H_ */
