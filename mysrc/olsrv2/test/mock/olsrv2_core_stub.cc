#include "../include/olsrv2.h"
#include "../include/olsrv2_routing.h"

namespace mysrc::olsrv2 {

// Minimal mock implementation of Olsrv2Core for unit tests that link against it
// but don't compile olsrv2.cc (or to avoid complex dependencies).
// However, test_routing links olsrv2_routing.cc which calls generateTc.
// We need to either link olsrv2.cc or provide a stub.

inet::Packet* Olsrv2Core::generateTc(const NhdpDb& nhdpDb, double now) {
    return nullptr; // Stub
}

// Add other necessary stubs if needed, e.g. constructor/destructor if not inline
// Constructor is inline in olsrv2.h
// processTcAndRecompute is called in processOlsrPacket

void Olsrv2Core::processTcAndRecompute(const inet::Olsrv2TcGroup* tc, double now) {
    // Stub
}

} // namespace mysrc::olsrv2
