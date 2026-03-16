#include <iostream>
#include <string>

#include "../include/nhdp/nhdp.h"
#include "../include/olsrv2_state.h"

#include "inet/common/packet/Packet.h"
#include "inet/networklayer/common/L3Address.h"
#include "../message/OLSRv2Packet_m.h"

using namespace mysrc::olsrv2;

namespace {
int g_failures = 0;

void expectTrue(bool cond, const std::string& msg) {
    if (!cond) {
        ++g_failures;
        std::cerr << "[FAIL] " << msg << "\n";
    }
}

void testNhdpHelloProcessing() {
    Nhdp nhdp;
    nhdp.setOriginator(inet::L3Address("10.0.0.1"));

    // Create a Hello Packet from neighbor 10.0.0.2
    auto hello = std::make_shared<inet::Olsrv2HelloPacket>();
    hello->setOriginator(inet::L3Address("10.0.0.2"));
    hello->setHTime(2.0);
    
    // Process Hello
    double now = 100.0;
    nhdp.processHello(hello.get(), inet::L3Address("10.0.0.2"), now);

    // Verify neighbor state in DB
    auto neighbor = nhdp.getDb().findNeighbor(inet::Ipv4Address("10.0.0.2"));
    expectTrue(neighbor.has_value(), "Neighbor must be added after Hello");
    
    // Verify link state (although implicit in current dummy implementation)
    // We can't access links easily via findNeighbor in current NhdpDb API without iterating.
    // But we can check if neighbor exists.
}

void testNeighborAndLinkLifecycle() {
    // Keep original test for Olsrv2State
    Olsrv2State state;

    const MainAddress local = MainAddress("10.0.0.1");
    const MainAddress n1 = MainAddress("10.0.0.2");
    const MainAddress n2 = MainAddress("10.0.0.3");

    LinkTuple link1;
    link1.local_iface_addr = local;
    link1.neighbor_iface_addr = n1;
    link1.sym_time.expires_at = 100.0;
    link1.heard_time.expires_at = 100.0;

    LinkTuple link2;
    link2.local_iface_addr = local;
    link2.neighbor_iface_addr = n2;
    link2.sym_time.expires_at = 5.0;
    link2.heard_time.expires_at = 5.0;

    state.upsertLink(link1);
    state.upsertLink(link2);
    expectTrue(state.getLinks().size() == 2, "2 links must be inserted");

    NeighborTuple neigh1;
    neigh1.main_addr = n1;
    neigh1.is_symmetric = true;
    neigh1.validity.expires_at = 100.0;
    state.upsertNeighbor(neigh1);

    NeighborTuple neigh2;
    neigh2.main_addr = n2;
    neigh2.is_symmetric = false;
    neigh2.validity.expires_at = 5.0;
    state.upsertNeighbor(neigh2);

    expectTrue(state.getNeighbors().size() == 2, "2 neighbors must be inserted");
    expectTrue(state.findNeighbor(n1).has_value(), "n1 must be present before purge");

    state.purgeExpired(10.0);

    expectTrue(state.getLinks().size() == 1, "expired link must be purged");
    expectTrue(state.getNeighbors().size() == 1, "expired neighbor must be purged");
    expectTrue(state.findNeighbor(n1).has_value(), "n1 must remain after purge");
    expectTrue(!state.findNeighbor(n2).has_value(), "n2 must be removed after purge");

    state.removeNeighbor(n1);
    expectTrue(state.getNeighbors().empty(), "removeNeighbor must erase n1");
}

void testMprSelection() {
    // Topology:
    // Me --(sym)--> N1 --(sym)--> N2
    // Me --(sym)--> N3 --(sym)--> N4
    // N1 covers N2. N3 covers N4.
    // MPR set should be {N1, N3}.

    Nhdp nhdp;
    nhdp.setOriginator(inet::L3Address("10.0.0.1"));
    double now = 100.0;

    // 1. Receive Hello from N1 advertising N2
    auto helloN1 = std::make_shared<inet::Olsrv2HelloPacket>();
    helloN1->setOriginator(inet::L3Address("10.0.0.2")); // N1
    helloN1->setHTime(2.0);
    helloN1->setNeighAddrsArraySize(2);
    helloN1->setNeighAddrs(0, inet::L3Address("10.0.0.1")); // Me (for symmetry)
    helloN1->setNeighAddrs(1, inet::L3Address("10.0.0.3")); // N2
    nhdp.processHello(helloN1.get(), inet::L3Address("10.0.0.2"), now);

    // 2. Receive Hello from N3 advertising N4
    auto helloN3 = std::make_shared<inet::Olsrv2HelloPacket>();
    helloN3->setOriginator(inet::L3Address("10.0.0.4")); // N3
    helloN3->setHTime(2.0);
    helloN3->setNeighAddrsArraySize(2);
    helloN3->setNeighAddrs(0, inet::L3Address("10.0.0.1")); // Me (for symmetry)
    helloN3->setNeighAddrs(1, inet::L3Address("10.0.0.5")); // N4
    nhdp.processHello(helloN3.get(), inet::L3Address("10.0.0.4"), now);

    // Calculate MPR
    auto mprSet = nhdp.calculateMprSet();

    expectTrue(mprSet.size() == 2, "MPR set size should be 2");
    expectTrue(mprSet.count(inet::Ipv4Address("10.0.0.2")), "N1 must be MPR");
    expectTrue(mprSet.count(inet::Ipv4Address("10.0.0.4")), "N3 must be MPR");
}
} // namespace

#include "test_runner.h"

// ... (existing code)

int run_nhdp_tests() {
    omnetpp::SimTime::setScaleExp(-9);
    std::cout << "Running test_nhdp...\n";
    testNhdpHelloProcessing();
    testNeighborAndLinkLifecycle();
    testMprSelection();

    if (g_failures == 0) {
        std::cout << "[PASS] test_nhdp\n";
        return 0;
    }

    std::cerr << "[FAIL] test_nhdp failures=" << g_failures << "\n";
    return 1;
}