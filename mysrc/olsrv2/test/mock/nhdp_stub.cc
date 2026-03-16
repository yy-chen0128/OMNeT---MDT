#include "../include/nhdp/nhdp.h"

namespace mysrc::olsrv2 {

// Minimal mock implementation of Nhdp for unit tests that link against it
// but don't compile nhdp.cc

Nhdp::Nhdp() {}

void Nhdp::setOriginator(const inet::L3Address& originator) {
    originator_ = originator;
}

inet::L3Address Nhdp::getOriginator() const {
    return originator_;
}

inet::Packet* Nhdp::generateHello(double now) {
    return nullptr; // Stub
}

void Nhdp::processHello(const inet::Olsrv2HelloPacket *hello, const inet::L3Address& source, double now) {
    // Stub
}

bool Nhdp::floodingSelector(rfc5444_writer *writer, rfc5444_writer_target *target, void *ptr) {
    return true;
}

bool Nhdp::forwardingSelector(rfc5444_writer_target *target, rfc5444_reader_tlvblock_context *context) {
    return true;
}

} // namespace mysrc::olsrv2
