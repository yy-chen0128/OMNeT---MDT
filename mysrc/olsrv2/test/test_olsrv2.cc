#include <iostream>
#include <string>
#include <vector>

#include "../include/olsrv2.h"

// Mock INET
#include "mock/inet/common/packet/Packet.h"
#include "mock/inet/networklayer/common/L3Address.h"
#include "mock/message/OLSRv2Packet_m.h"

using namespace mysrc::olsrv2;

namespace {
int g_failures = 0;

void expectTrue(bool cond, const std::string& msg) {
    if (!cond) {
        ++g_failures;
        std::cerr << "[FAIL] " << msg << "\n";
    }
}

const RouteTuple* findRoute(const std::vector<RouteTuple>& routes, const MainAddress& dst) {
    for (const auto& r : routes) {
        if (r.destination == dst)
            return &r;
    }
    return nullptr;
}

void testCoreRecomputeRoutes() {
    const MainAddress origin("10.0.0.1");
    const MainAddress n1("10.0.0.2");
    const MainAddress n2("10.0.0.3");

    Olsrv2Core core(origin);

    NeighborTuple oneHop;
    oneHop.main_addr = n1;
    oneHop.is_symmetric = true;
    oneHop.validity.expires_at = 100.0;
    core.state().upsertNeighbor(oneHop);

    TopologyTuple topo;
    topo.last_addr = n1;
    topo.dest_addr = n2;
    topo.ansn = 7;
    topo.validity.expires_at = 100.0;
    core.state().upsertTopology(topo);

    core.recomputeRoutes();
    const auto& routes = core.state().getRoutes();

    expectTrue(routes.size() == 2, "must have routes to n1 and n2");

    const RouteTuple* r1 = findRoute(routes, n1);
    const RouteTuple* r2 = findRoute(routes, n2);

    expectTrue(r1 != nullptr, "route to one-hop neighbor n1 must exist");
    expectTrue(r2 != nullptr, "route to two-hop node n2 must exist");

    if (r1) {
        expectTrue(r1->next_hop == n1, "n1 next_hop must be n1");
        expectTrue(r1->hop_count == 1, "n1 hop_count must be 1");
    }
    if (r2) {
        expectTrue(r2->next_hop == n1, "n2 next_hop must be n1");
        expectTrue(r2->hop_count == 2, "n2 hop_count must be 2");
    }
}

void testProcessTcAndRecompute() {
    const MainAddress origin("10.0.0.1");
    const MainAddress n1("10.0.0.2"); // My neighbor
    const MainAddress n2("10.0.0.3"); // n1's neighbor

    Olsrv2Core core(origin);

    // Setup: n1 is a neighbor
    NeighborTuple oneHop;
    oneHop.main_addr = n1;
    oneHop.is_symmetric = true;
    oneHop.validity.expires_at = 200.0;
    core.state().upsertNeighbor(oneHop);

    // Create TC from n1 advertising n2
    auto tc = std::make_shared<inet::Olsrv2TcGroup>();
    tc->setOriginatorsArraySize(1);
    tc->setOriginators(0, inet::L3Address("10.0.0.2"));
    tc->setAnsnsArraySize(1);
    tc->setAnsns(0, 10);
    
    tc->setAdvCountsArraySize(1);
    tc->setAdvCounts(0, 1);
    
    tc->setAdvertisedNeighborsArraySize(1);
    tc->setAdvertisedNeighbors(0, inet::L3Address("10.0.0.3"));
    
    // Process TC
    core.processTcAndRecompute(tc.get(), 100.0);
    
    // Verify Routes
    const auto& routes = core.state().getRoutes();
    expectTrue(routes.size() == 2, "Should have route to n2 via n1");
    
    const RouteTuple* r2 = findRoute(routes, n2);
    expectTrue(r2 != nullptr, "Route to n2 must exist");
    if (r2) {
        expectTrue(r2->next_hop == n1, "n2 via n1");
    }
}

void testGlobalTunables() {
    expectTrue(Olsrv2::get_tc_interval() > 0, "tc interval must be positive");
    expectTrue(Olsrv2::get_tc_validity() > 0, "tc validity must be positive");
}
} // namespace

int main() {
    std::cout << "Running test_olsrv2...\n";
    testCoreRecomputeRoutes();
    testProcessTcAndRecompute();
    testGlobalTunables();

    if (g_failures == 0) {
        std::cout << "[PASS] test_olsrv2\n";
        return 0;
    }

    std::cerr << "[FAIL] test_olsrv2 failures=" << g_failures << "\n";
    return 1;
}