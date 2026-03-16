#pragma once

#include <vector>

#include "olsrv2_state.h"

namespace mysrc::olsrv2 {

struct TcAdvertisedNeighbor {
    MainAddress neighbor_main_addr;
};

struct TcMessageView {
    MainAddress originator;
    SequenceNumber ansn = 0;
    double validity_duration = 0.0;
    std::vector<TcAdvertisedNeighbor> advertised_neighbors;
};

// Applies validated TC information into topology set.
class Olsrv2TcProcessor
{
  public:
    // now is simulation time in seconds.
    void applyTc(const TcMessageView& tc, double now, Olsrv2State& state) const;
};

} // namespace mysrc::olsrv2
