//
// Created by GuoYudi on 2026/3/10.
//

#ifndef MYSRC_OLSRV2_NHDP_NHDP_H_
#define MYSRC_OLSRV2_NHDP_NHDP_H_

#include <cstdint>
#include <string>

#include "inet/networklayer/common/L3Address.h"
#include "nhdp_db.h"
#include "../../message/OLSRv2Packet_m.h"

namespace mysrc::olsrv2 {

// Forward declarations for RFC5444/OONF bridge types.
// Keep NHDP headers implementation-agnostic and lightweight.
struct netaddr;
struct rfc5444_writer;
struct rfc5444_writer_target;
struct rfc5444_reader_tlvblock_context;

/*! subsystem identifier */
inline constexpr const char *OONF_NHDP_SUBSYSTEM = "nhdp";

/*! special name for 'no' metric/mpr */
inline constexpr const char *CFG_DOMAIN_NO_METRIC_MPR = "-";

/*! special name for 'any' metric/mpr */
inline constexpr const char *CFG_DOMAIN_ANY_METRIC_MPR = "*";

/*!
 * NHDP protocol constants used by OLSRv2 integration.
 */
enum NhdpProtocolConstants : uint16_t {
    /*! maximum number of metric domains */
    NHDP_MAXIMUM_DOMAINS = 4,

    /*! message tlv for transporting IPv4 originator in ipv6 messages */
    NHDP_MSGTLV_IPV4ORIGINATOR = 226,

    /*! message tlv for transporting mac address */
    NHDP_MSGTLV_MAC = 227,

    /*! Address TLV for custom link metric data */
    NHDP_ADDRTLV_LQ_CUSTOM = 228,
};

/*!
 * Core NHDP facade for originator management and forwarding/flooding hooks.
 * This header only exposes interfaces; behavior is implemented in .cc.
 */
class Nhdp
{
  public:
    Nhdp();
    ~Nhdp() = default;

    Nhdp(const Nhdp&) = delete;
    Nhdp& operator=(const Nhdp&) = delete;
    Nhdp(Nhdp&&) = delete;
    Nhdp& operator=(Nhdp&&) = delete;

    // Originator management.
    void setOriginator(const inet::L3Address& originator);
    inet::L3Address getOriginator() const;

    // RFC5444 selectors used by NHDP/OLSRv2 message processing.
    bool floodingSelector(
        rfc5444_writer *writer,
        rfc5444_writer_target *target,
        void *ptr);

    bool forwardingSelector(
        rfc5444_writer_target *target,
        rfc5444_reader_tlvblock_context *context);

    // Hello Message Processing
    void processHello(const inet::Olsrv2HelloPacket *hello, const inet::L3Address& source, double now);
    inet::Packet* generateHello(double now);

    std::set<inet::Ipv4Address> calculateMprSet() const;
    bool shouldForwardFrom(const inet::Ipv4Address& neighbor_originator, double now) const;

    // Access to DB
    const NhdpDb& getDb() const { return db_; }
    NhdpDb& getDb() { return db_; }

  private:
    inet::L3Address originator_;
    NhdpDb db_;
    uint16_t helloSeqNum_ = 0;
    uint64_t nextLinkId_ = 1;
};

} // namespace mysrc::olsrv2

#endif // MYSRC_OLSRV2_NHDP_NHDP_H_
