#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "inet/common/INETDefs.h"
#include "omnetpp/cnullenvir.h"
#include "omnetpp/csimulation.h"
#include "test_runner.h"

int main() {
    omnetpp::SimTime::setScaleExp(-9);

    int argc = 0;
    char **argv = nullptr;
    auto *envir = new omnetpp::cNullEnvir(argc, argv, nullptr);
    auto sim = std::make_unique<omnetpp::cSimulation>("olsrv2-test", envir);
    omnetpp::cSimulation::setActiveSimulation(sim.get());
    
    std::cout << "Starting Test Suite...\n";
    int failures = 0;

    std::cout << "--------------------------------\n";
    std::cout << "Running NHDP Tests:\n";
    if (run_nhdp_tests() != 0) {
        failures++;
        std::cout << "[FAIL] NHDP Tests Failed\n";
    } else {
        std::cout << "[PASS] NHDP Tests Passed\n";
    }

    std::cout << "--------------------------------\n";
    std::cout << "Running OLSRv2 Core Tests:\n";
    if (run_olsrv2_tests() != 0) {
        failures++;
        std::cout << "[FAIL] OLSRv2 Core Tests Failed\n";
    } else {
        std::cout << "[PASS] OLSRv2 Core Tests Passed\n";
    }

    std::cout << "--------------------------------\n";
    std::cout << "Running Routing Tests:\n";
    if (run_routing_tests() != 0) {
        failures++;
        std::cout << "[FAIL] Routing Tests Failed\n";
    } else {
        std::cout << "[PASS] Routing Tests Passed\n";
    }

    std::cout << "--------------------------------\n";
    std::cout << "Running RFC5444 Tests:\n";
    if (run_rfc5444_tests() != 0) {
        failures++;
        std::cout << "[FAIL] RFC5444 Tests Failed\n";
    } else {
        std::cout << "[PASS] RFC5444 Tests Passed\n";
    }

    std::cout << "--------------------------------\n";
    std::cout << "Running TC Tests:\n";
    if (run_tc_tests() != 0) {
        failures++;
        std::cout << "[FAIL] TC Tests Failed\n";
    } else {
        std::cout << "[PASS] TC Tests Passed\n";
    }

    std::cout << "--------------------------------\n";
    std::cout << "Running Integration Tests:\n";
    if (run_integration_tests() != 0) {
        failures++;
        std::cout << "[FAIL] Integration Tests Failed\n";
    } else {
        std::cout << "[PASS] Integration Tests Passed\n";
    }

    std::cout << "--------------------------------\n";
    if (failures == 0) {
        std::cout << "All Tests Passed!\n";
        omnetpp::cSimulation::setActiveSimulation(nullptr);
        return 0;
    } else {
        std::cerr << "Some tests failed (" << failures << " suites failed).\n";
        omnetpp::cSimulation::setActiveSimulation(nullptr);
        return 1;
    }
}
