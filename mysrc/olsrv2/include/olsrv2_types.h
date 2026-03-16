//
// Created by GuoYudi on 2026/3/14.
//

#pragma once

#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "inet/networklayer/common/L3Address.h"

namespace mysrc::olsrv2 {

    using MainAddress = inet::Ipv4Address;
    using SequenceNumber = uint16_t;

    // Basic validity timestamps in simulation time seconds represented as double.
    // The caller can map this to simtime_t if needed.
    struct ValidityTime {
        double expires_at = 0.0;
        bool isValid(double now) const { return now < expires_at; }
    };

    struct LinkTuple {
        MainAddress local_iface_addr;
        MainAddress neighbor_iface_addr;
        ValidityTime sym_time;
        ValidityTime heard_time;
    };

    struct NeighborTuple {
        MainAddress main_addr;
        bool is_symmetric = false;
        std::set<MainAddress> iface_addrs;
        ValidityTime validity;
    };

    struct TopologyTuple {
        MainAddress dest_addr;
        MainAddress last_addr;
        SequenceNumber ansn = 0;
        ValidityTime validity;
    };

    struct RouteTuple {
        MainAddress destination;
        MainAddress next_hop;
        uint32_t metric = 0;
        uint8_t hop_count = 0;
    };

} // namespace mysrc::olsrv2
