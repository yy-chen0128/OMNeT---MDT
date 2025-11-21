/*
 * test_Neighbor_dt.cc
 *
 *  Created on: Sep 26, 2025
 *      Author: yychen
 */

// your header (which now includes createDTManager/getDTManager)

// tests/test_NeighborDT.cpp
// GoogleTest unit tests for NeighborDT (pure-class CGAL/DT manager).
//
// Requirements to compile/run:
//   - gtest available
//   - CGAL available
//   - INET/OMNeT includes for L3Address and simtime_t (we use SIMTIME_ZERO)
//   - Link with CGAL and gtest
//

#include <gtest/gtest.h>
#include <vector>
#include <set>
#include <string>

#include <omnetpp.h> // for SIMTIME_ZERO type
#include "../../inet-4.5.4/src/inet/networklayer/common/L3Address.h"
#include "../../inet-4.5.4/src/inet/networklayer/ipv4/Ipv4Header_m.h"
#include "../../inet-4.5.4/src/inet/networklayer/ipv4/Ipv4Route.h"

#include "../mysrc/routing/Neighbor.h"

#include "../mysrc/routing/NeighborDT.h"
#include "../mysrc/routing/NeighborEntry.h"

using namespace mysrc::routing;
using namespace inet;

// helper constructors
static L3Address makeAddr(const char *ipv4str) {
    return L3Address(Ipv4Address(ipv4str));
}

static NeighborEntry makeNeighborEntry(const L3Address &addr,
                                       const std::vector<double> &coords,
                                       bool attached = true)
{
    NeighborEntry ne;
    ne.addr = addr;
    ne.coords = coords;
    ne.dim = (int)coords.size();
    ne.attachedToDT = attached;
    ne.timeout = SIMTIME_ZERO;
    ne.type = std::string("test");
    return ne;
}

// Test fixture constructs empty containers and a NeighborDT instance
class NeighborDTTest : public ::testing::Test {
protected:
    std::map<L3Address, NeighborEntry> known;
    Delaunay3 dt;
    std::map<Delaunay3::Vertex_handle, std::set<L3Address>> vhToAddrs;
    std::map<L3Address, Delaunay3::Vertex_handle> addrToVh;
    std::set<L3Address> pending;
    std::map<L3Address, NeighborEntry> dtNeighbors;
    std::map<L3Address, NeighborEntry> minSimplex;

    std::unique_ptr<NeighborDT> dman;

    void SetUp() override {
        known.clear();
        vhToAddrs.clear();
        addrToVh.clear();
        pending.clear();
        dtNeighbors.clear();
        minSimplex.clear();

        dman = std::make_unique<NeighborDT>(known, dt, vhToAddrs, addrToVh, pending, dtNeighbors, minSimplex);
    }

    void TearDown() override {
        dman.reset();
    }
};

///////////////////////////////////////////////////////////////////////////////
// TEST 1
//
// Method under test: makePoint(const std::vector<double>&)
// Expected behavior:
//   - Project input 1/2/3-d coordinates into a CGAL Point_3 with missing components = 0.0.
// Test scenarios:
//   - 1D input -> (x,0,0)
//   - 2D input -> (x,y,0)
//   - 3D input -> (x,y,z)
///////////////////////////////////////////////////////////////////////////////
TEST_F(NeighborDTTest, MakePointProject1D2D3D) {
    std::vector<double> v1 = {5.0};
    Point3 p1 = dman->makePoint(v1);
    EXPECT_DOUBLE_EQ(p1.x(), 5.0);
    EXPECT_DOUBLE_EQ(p1.y(), 0.0);
    EXPECT_DOUBLE_EQ(p1.z(), 0.0);

    std::vector<double> v2 = {1.1, 2.2};
    Point3 p2 = dman->makePoint(v2);
    EXPECT_DOUBLE_EQ(p2.x(), 1.1);
    EXPECT_DOUBLE_EQ(p2.y(), 2.2);
    EXPECT_DOUBLE_EQ(p2.z(), 0.0);

    std::vector<double> v3 = {9.9, 8.8, 7.7};
    Point3 p3 = dman->makePoint(v3);
    EXPECT_DOUBLE_EQ(p3.x(), 9.9);
    EXPECT_DOUBLE_EQ(p3.y(), 8.8);
    EXPECT_DOUBLE_EQ(p3.z(), 7.7);
}

///////////////////////////////////////////////////////////////////////////////
// TEST 2
//
// Method under test: insertNodeToDT(const L3Address&, const Point3&)
// Expected behavior:
//   - If two different addresses map to same coordinates they should reuse the same
//     Delaunay vertex (detection by dt.nearest_vertex equality + point equality).
//   - addrToVh maps both addresses to the same Vertex_handle, and vhToAddrs[vh]
//     contains both addresses.
// Test scenario:
//   - Insert two addresses with identical coordinates and verify reuse.
///////////////////////////////////////////////////////////////////////////////
TEST_F(NeighborDTTest, Insert_ReusesVertexForSameCoordinates) {
    L3Address a1 = makeAddr("10.0.0.1");
    L3Address a2 = makeAddr("10.0.0.2");
    std::vector<double> coords = {0.0, 0.0, 0.0};

    // populate known (not strictly necessary for insert, but realistic)
    known[a1] = makeNeighborEntry(a1, coords, true);
    known[a2] = makeNeighborEntry(a2, coords, true);

    dman->insertNodeToDT(a1, dman->makePoint(coords));
    dman->insertNodeToDT(a2, dman->makePoint(coords));

    auto &addrToVhRef = dman->addrToVhRef();
    auto &vhToAddrsRef = dman->vhToAddrsRef();

    ASSERT_TRUE(addrToVhRef.find(a1) != addrToVhRef.end());
    ASSERT_TRUE(addrToVhRef.find(a2) != addrToVhRef.end());

    auto vh1 = addrToVhRef[a1];
    auto vh2 = addrToVhRef[a2];
    EXPECT_EQ(vh1, vh2);

    auto it = vhToAddrsRef.find(vh1);
    ASSERT_TRUE(it != vhToAddrsRef.end());
    EXPECT_EQ(it->second.size(), 2u);
    EXPECT_TRUE(it->second.count(a1));
    EXPECT_TRUE(it->second.count(a2));
}

///////////////////////////////////////////////////////////////////////////////
// TEST 3
//
// Method under test: removeNodeFromDT(const L3Address&, int threshold = ...)
// Expected behavior:
//   - Removing an address that is NOT the last ref of a vertex should simply
//     erase the addr->vh mapping and leave the vh->addrs set containing other addrs.
//   - Removing the last address referencing a vertex should erase vhToAddrs entry
//     and add that address to pendingRebuildAddrs (deferred cleanup).
//   - If pendingRebuildAddrs size reaches threshold, the method will trigger a
//     rebuild (rebuildDTFromKnownNodes) which clears pendingRebuildAddrs.
// Test scenario:
//   - Insert two addresses sharing a vertex, remove them sequentially and assert behavior,
//     then remove a third to reach threshold and cause automatic rebuild.
///////////////////////////////////////////////////////////////////////////////
TEST_F(NeighborDTTest, RemoveNode_DeferredCleanupAndRebuildTrigger) {
    L3Address a1 = makeAddr("10.0.0.1");
    L3Address a2 = makeAddr("10.0.0.2");
    L3Address a3 = makeAddr("10.0.0.3");

    // all three share coords except a3 different
    std::vector<double> c01 = {1.0, 0.0, 0.0};
    std::vector<double> c3  = {2.0, 0.0, 0.0};

    known[a1] = makeNeighborEntry(a1, c01, true);
    known[a2] = makeNeighborEntry(a2, c01, true);
    known[a3] = makeNeighborEntry(a3, c3, true);

    // insert a1 and a2 -> same vertex
    dman->insertNodeToDT(a1, dman->makePoint(c01));
    dman->insertNodeToDT(a2, dman->makePoint(c01));

    auto &addrToVhRef = dman->addrToVhRef();
    auto &vhToAddrsRef = dman->vhToAddrsRef();

    ASSERT_TRUE(addrToVhRef.count(a1));
    ASSERT_TRUE(addrToVhRef.count(a2));
    auto vh_shared = addrToVhRef[a1];
    EXPECT_EQ(vh_shared, addrToVhRef[a2]);

    // remove a1 -> vertex still referenced by a2
    dman->removeNodeFromDT(a1, /*threshold*/ 1000); // big threshold so no auto rebuild
    EXPECT_FALSE(addrToVhRef.count(a1));
    ASSERT_TRUE(vhToAddrsRef.find(vh_shared) != vhToAddrsRef.end());
    EXPECT_TRUE(vhToAddrsRef[vh_shared].count(a2));

    // remove a2 -> last reference removed -> vh removed and a2 in pending
    dman->removeNodeFromDT(a2, /*threshold*/ 1000);
    EXPECT_FALSE(addrToVhRef.count(a2));
    EXPECT_TRUE(vhToAddrsRef.find(vh_shared) == vhToAddrsRef.end());
    EXPECT_TRUE(dman->pendingRebuildRef().count(a2));

    // now insert a3 and remove it using threshold = 2 to force rebuild inside removeNodeFromDT
    dman->insertNodeToDT(a3, dman->makePoint(c3));
    // remove a3 with threshold 2 -> pending will contain (a2,a3) and trigger rebuild
    dman->removeNodeFromDT(a3, /*threshold*/ 2);

    // after automatic rebuild, pending should be cleared
    EXPECT_TRUE(dman->pendingRebuildRef().empty());
    // and addrToVh should have entries for known nodes that survive rebuild (depending on include predicate default)
    // At least the known map contains a1,a2,a3; verify addrToVh is consistent (no dangling entries)
    for (const auto &kv : addrToVhRef) {
        ASSERT_TRUE(kv.second != Delaunay3::Vertex_handle());
    }
}

///////////////////////////////////////////////////////////////////////////////
// TEST 4
//
// Methods under test:
//   - rebuildDTFromKnownNodes(includePred)
//   - rebuildDTNeighborCache(selfAddr, includePred)
// Expected behavior:
//   - Given 4 non-coplanar points (a tetrahedron) inserted via knownNodeCoords with
//     attachedToDT=true, a full rebuild should produce a 3D triangulation (dt.dimension()==3).
//   - rebuildDTNeighborCache(self) should populate dtNeighbors for the self vertex with
//     the three other vertices (neighbors).
// Test scenario:
//   - Build tetrahedron points A(self), B, C, D; run rebuild & neighbor cache; assert dimension==3 and
//     dtNeighbors for self has size 3 and contains B,C,D.
///////////////////////////////////////////////////////////////////////////////
TEST_F(NeighborDTTest, RebuildFromKnownNodesAndNeighborCache_Tetrahedron) {
    // coordinates of a simple tetrahedron (non-coplanar)
    // choose points spaced to avoid degeneracy
    L3Address A = makeAddr("10.0.0.1"); // self
    L3Address B = makeAddr("10.0.0.2");
    L3Address C = makeAddr("10.0.0.3");
    L3Address D = makeAddr("10.0.0.4");

    known[A] = makeNeighborEntry(A, {0.0, 0.0, 0.0}, true);
    known[B] = makeNeighborEntry(B, {1.0, 0.0, 0.0}, true);
    known[C] = makeNeighborEntry(C, {0.0, 1.0, 0.0}, true);
    known[D] = makeNeighborEntry(D, {0.0, 0.0, 1.0}, true);

    // full rebuild (no includePred -> include all known)
    dman->rebuildDTFromKnownNodes(nullptr);

    // After rebuild, triangulation should be 3D
    EXPECT_EQ(dman->dtRef().dimension(), 3);

    // build neighbor cache for A (self)
    auto includePred = [](const NeighborEntry &ne)->bool {
        // include only attached; (we already set attached = true)
        return ne.attachedToDT;
    };
    dman->rebuildDTNeighborCache(A, includePred);

    // dtNeighbors should contain B, C, D
    auto &dtN = dman->dtNeighborsRef();
    EXPECT_EQ(dtN.size(), 3u);
    EXPECT_TRUE(dtN.find(B) != dtN.end());
    EXPECT_TRUE(dtN.find(C) != dtN.end());
    EXPECT_TRUE(dtN.find(D) != dtN.end());
}

///////////////////////////////////////////////////////////////////////////////
// TEST 5
//
// Methods under test:
//   - computeDTNeighborsForCoord(targetCoord, addr, selfAddr, includePred)
//   - isDTNeighborOfTarget(targetCoord, addr, selfAddr, includePred)
// Expected behavior:
//   - For a tetrahedron, computeDTNeighborsForCoord for one vertex should return the three others.
//   - isDTNeighborOfTarget should return true when self and target share a simplex (they do in tetrahedron).
// Test scenario:
//   - Reuse tetrahedron constructed above; ask neighbors for B and check that A is among them
//     and that isDTNeighborOfTarget returns true for (self=A, target=B).
///////////////////////////////////////////////////////////////////////////////
TEST_F(NeighborDTTest, ComputeNeighborsForCoordAndIsDTNeighbor_Tetrahedron) {
    L3Address A = makeAddr("10.0.0.1"); // self
    L3Address B = makeAddr("10.0.0.2");
    L3Address C = makeAddr("10.0.0.3");
    L3Address D = makeAddr("10.0.0.4");

    known[A] = makeNeighborEntry(A, {0.0,0.0,0.0}, true);
    known[B] = makeNeighborEntry(B, {1.0,0.0,0.0}, true);
    known[C] = makeNeighborEntry(C, {0.0,1.0,0.0}, true);
    known[D] = makeNeighborEntry(D, {0.0,0.0,1.0}, true);

    dman->rebuildDTFromKnownNodes(nullptr);

    // compute neighbors for B
    std::vector<NeighborEntry> nbList = dman->computeDTNeighborsForCoord(known[B].coords, B, A, nullptr);

    // Expect neighbors to include A,C,D (order unspecified)
    std::set<std::string> s;
    for (const auto &ne : nbList) s.insert(ne.addr.str());
    EXPECT_TRUE(s.count(A.str()));
    EXPECT_TRUE(s.count(C.str()));
    EXPECT_TRUE(s.count(D.str()));
    EXPECT_EQ(s.size(), 3u);

    // isDTNeighborOfTarget: check whether A is a DT-neighbor of B
    bool isNeighbor = dman->isDTNeighborOfTarget(known[B].coords, B, A, nullptr);
    EXPECT_TRUE(isNeighbor);
}

///////////////////////////////////////////////////////////////////////////////
// TEST 6
//
// Method under test:
//   computeMinimalNeighborSubset_CGAL_3D(selfAddr, includePred)
// Expected behavior:
//   - For self=A in tetrahedron, incident_cells form one simplex containing B,C,D.
//   - Minimal hitting set algorithm should select a subset of neighbors that cover all simplices;
//     In tetrahedron with only one simplex incident to A, any single neighbor from that simplex could cover it.
//   - We assert that the returned set is non-empty and contains only known nodes among {B,C,D}.
// Test scenario:
//   - Use tetrahedron and call computeMinimalNeighborSubset_CGAL_3D(A)
///////////////////////////////////////////////////////////////////////////////
TEST_F(NeighborDTTest, ComputeMinimalNeighborSubset_CGAL_3D_Tetrahedron) {
    L3Address A = makeAddr("10.0.0.1"); // self
    L3Address B = makeAddr("10.0.0.2");
    L3Address C = makeAddr("10.0.0.3");
    L3Address D = makeAddr("10.0.0.4");

    known[A] = makeNeighborEntry(A, {0.0,0.0,0.0}, true);
    known[B] = makeNeighborEntry(B, {1.0,0.0,0.0}, true);
    known[C] = makeNeighborEntry(C, {0.0,1.0,0.0}, true);
    known[D] = makeNeighborEntry(D, {0.0,0.0,1.0}, true);

    dman->rebuildDTFromKnownNodes(nullptr);

    // compute minimal subset
    dman->computeMinimalNeighborSubset_CGAL_3D(A, nullptr);

    auto &minSet = dman->minSimplexNodesRef();
    // Should contain at least one of B,C,D
    EXPECT_FALSE(minSet.empty());
    for (const auto &kv : minSet) {
        EXPECT_TRUE(kv.first == B || kv.first == C || kv.first == D);
    }
}



