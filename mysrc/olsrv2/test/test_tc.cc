#include <iostream>
#include <string>

#include "../include/olsrv2_tc.h"

using namespace mysrc::olsrv2;

namespace {
int g_failures = 0;

void expectTrue(bool cond, const std::string& msg) {
    if (!cond) {
        ++g_failures;
        std::cerr << "[FAIL] " << msg << "\n";
    }
}

void testTcProcessorAppliesTopology() {
    Olsrv2State state;
    Olsrv2TcProcessor tcProcessor;

    TcMessageView tc;
    tc.originator = MainAddress("10.0.0.10");
    tc.ansn = 42;
    tc.validity_duration = 12.5;
    tc.advertised_neighbors.push_back({MainAddress("10.0.0.20")});
    tc.advertised_neighbors.push_back({MainAddress("10.0.0.30")});

    const double now = 100.0;
    tcProcessor.applyTc(tc, now, state);

    const auto topo = state.getTopologyTuples();
    expectTrue(topo.size() == 2, "TC must create 2 topology tuples");

    for (const auto& t : topo) {
        expectTrue(t.last_addr == tc.originator, "last_addr must equal TC originator");
        expectTrue(t.ansn == tc.ansn, "ansn must match TC ansn");
        expectTrue(t.validity.expires_at == now + tc.validity_duration, "expiry must be now + validity_duration");
    }

    state.purgeExpired(now + 13.0);
    expectTrue(state.getTopologyTuples().empty(), "topology tuples must expire after validity window");
}
} // namespace

int main() {
    std::cout << "Running test_message...\n";
    testTcProcessorAppliesTopology();

    if (g_failures == 0) {
        std::cout << "[PASS] test_message\n";
        return 0;
    }

    std::cerr << "[FAIL] test_message failures=" << g_failures << "\n";
    return 1;
}