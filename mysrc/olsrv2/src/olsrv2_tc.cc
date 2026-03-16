//
// Created by GuoYudi on 2026/3/14.
//
#include "../include/olsrv2_tc.h"

namespace mysrc::olsrv2 {

void Olsrv2TcProcessor::applyTc(const TcMessageView& tc, double now, Olsrv2State& state) const
{
    const double expires = now + tc.validity_duration;

    for (const auto& adv : tc.advertised_neighbors) {
        TopologyTuple tuple;
        tuple.dest_addr = adv.neighbor_main_addr;
        tuple.last_addr = tc.originator;
        tuple.ansn = tc.ansn;
        tuple.validity.expires_at = expires;
        state.upsertTopology(tuple);
    }
}

} // namespace omnet::olsrv2