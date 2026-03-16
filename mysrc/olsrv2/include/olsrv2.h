#ifndef MYSRC_OLSRV2_OLSRV2_H
#define MYSRC_OLSRV2_OLSRV2_H


// Public umbrella header for OLSRv2 core logic.
// External users should include this file instead of individual headers.

#include "olsrv2_routing.h"
#include "olsrv2_state.h"
#include "olsrv2_tc.h"
#include "olsrv2_types.h"
#include "nhdp/nhdp_db.h"
#include "../message/OLSRv2Packet_m.h"
#include "inet/common/packet/Packet.h"

namespace mysrc::olsrv2 {

struct netaddr;
struct rfc5444_reader_tlvblock_context;

// Facade for typical OLSRv2 core workflow in simulation modules.
class Olsrv2Core
{
  public:
    explicit Olsrv2Core(MainAddress originator = MainAddress::UNSPECIFIED_ADDRESS)
        : originator_(originator)
    {
    }

    void setOriginator(MainAddress originator) { originator_ = originator; }
    MainAddress getOriginator() const { return originator_; }

    Olsrv2State& state() { return state_; }
    const Olsrv2State& state() const { return state_; }

    // Ingest one TC and recompute route set.
    void processTcAndRecompute(const inet::Olsrv2TcGroup* tc, double now);

    // Recompute route set using current state.
    void recomputeRoutes()
    {
        state_.setRoutes(routing_engine_.computeRoutes(originator_, state_));
    }

    void purgeExpired(double now) { state_.purgeExpired(now); }
    
    // Generate TC packet
    inet::Packet* generateTc(const NhdpDb& nhdpDb, double now);

  private:
    MainAddress originator_;
    Olsrv2State state_;
    Olsrv2TcProcessor tc_processor_;
    Olsrv2RoutingEngine routing_engine_;
    uint16_t tcSeqNum_ = 0;
    uint16_t ansn_ = 0;
};

class Olsrv2
{
public:
    static uint64_t get_tc_interval(void);
    static uint64_t get_tc_validity(void);

    static bool is_nhdp_routable(struct netaddr* addr);
    static bool is_routable(struct netaddr* addr);

    static bool mpr_shall_process(struct rfc5444_reader_tlvblock_context* context, uint64_t vtime);
    static bool mpr_shall_forwarding(
        struct rfc5444_reader_tlvblock_context* context, struct netaddr* source_address, uint64_t vtime);

    static void generate_tcs(bool generate);
};

} // namespace mysrc::olsrv2

#endif // MYSRC_OLSRV2_OLSRV2_H