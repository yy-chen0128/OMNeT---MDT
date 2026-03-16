#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "../include/olsrv2_routing.h"
#include "../include/olsrv2_state.h"

// Mock INET (needed because olsrv2_state.h includes olsrv2_types.h which includes L3Address.h)
#include "mock/inet/networklayer/common/L3Address.h"

using namespace mysrc::olsrv2;

namespace {
int g_failures = 0;

void expectTrue(bool cond, const std::string& msg) {
    if (!cond) {
        ++g_failures;
        std::cerr << "[FAIL] " << msg << "\n";
    }
}

void testDiamondTopology() {
    //      S
    //    /   \
    //   A     B
    //    \   /
    //      D
    // All links cost 1. S->D should have cost 2. Next hop could be A or B.
    
    MainAddress S("10.0.0.1");
    MainAddress A("10.0.0.2");
    MainAddress B("10.0.0.3");
    MainAddress D("10.0.0.4");
    
    Olsrv2State state;
    
    // Neighbors of S
    NeighborTuple nA; nA.main_addr = A; nA.is_symmetric = true; nA.validity.expires_at = 100.0;
    state.upsertNeighbor(nA);
    
    NeighborTuple nB; nB.main_addr = B; nB.is_symmetric = true; nB.validity.expires_at = 100.0;
    state.upsertNeighbor(nB);
    
    // Topology: A->D
    TopologyTuple tAD; tAD.last_addr = A; tAD.dest_addr = D; tAD.ansn = 1; tAD.validity.expires_at = 100.0;
    state.upsertTopology(tAD);
    
    // Topology: B->D
    TopologyTuple tBD; tBD.last_addr = B; tBD.dest_addr = D; tBD.ansn = 1; tBD.validity.expires_at = 100.0;
    state.upsertTopology(tBD);
    
    Olsrv2RoutingEngine engine;
    auto routes = engine.computeRoutes(S, state);
    
    bool foundD = false;
    for (const auto& r : routes) {
        if (r.destination == D) {
            foundD = true;
            expectTrue(r.metric == 2, "Metric to D should be 2");
            expectTrue(r.hop_count == 2, "Hop count to D should be 2");
            expectTrue(r.next_hop == A || r.next_hop == B, "Next hop to D should be A or B");
        }
    }
    expectTrue(foundD, "Route to D must be found");
}

void testBrokenLink() {
    // S -> A -> B
    // Link S->A breaks (neighbor removed or asymmetric)
    
    MainAddress S("10.0.0.1");
    MainAddress A("10.0.0.2");
    MainAddress B("10.0.0.3");
    
    Olsrv2State state;
    
    // Initial state: A is neighbor
    NeighborTuple nA; nA.main_addr = A; nA.is_symmetric = true; nA.validity.expires_at = 100.0;
    state.upsertNeighbor(nA);
    
    // Topology A->B
    TopologyTuple tAB; tAB.last_addr = A; tAB.dest_addr = B; tAB.ansn = 1; tAB.validity.expires_at = 100.0;
    state.upsertTopology(tAB);
    
    Olsrv2RoutingEngine engine;
    auto routes = engine.computeRoutes(S, state);
    expectTrue(routes.size() == 2, "Should have routes to A and B");
    
    // Break link S->A (make it asymmetric or remove)
    state.removeNeighbor(A); // Or set asymmetric
    
    routes = engine.computeRoutes(S, state);
    expectTrue(routes.empty(), "Should have no routes after S->A break");
}

} // namespace

int main() {
    std::cout << "Running test_routing...\n";
    testDiamondTopology();
    testBrokenLink();

    if (g_failures == 0) {
        std::cout << "[PASS] test_routing\n";
        return 0;
    }

    std::cerr << "[FAIL] test_routing failures=" << g_failures << "\n";
    return 1;
}
