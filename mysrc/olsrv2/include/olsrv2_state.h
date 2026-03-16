//
// Created by GuoYudi on 2026/3/14.
//

#ifndef MYSRC_OLSRV2_OLSRV2STATE_H
#define MYSRC_OLSRV2_OLSRV2STATE_H

#include <map>
#include <optional>
#include <set>
#include <vector>

#include "olsrv2_types.h"

namespace mysrc::olsrv2 {

    // Owns OLSRv2 RFC7181 core sets: Link/Neighbor/Topology/Routing.
    class Olsrv2State
    {
    public:
        // Link Set
        void upsertLink(const LinkTuple& tuple);
        void removeLink(const MainAddress& local_iface, const MainAddress& neighbor_iface);
        std::vector<LinkTuple> getLinks() const;

        // Neighbor Set
        void upsertNeighbor(const NeighborTuple& tuple);
        void removeNeighbor(const MainAddress& main_addr);
        std::optional<NeighborTuple> findNeighbor(const MainAddress& main_addr) const;
        std::vector<NeighborTuple> getNeighbors() const;

        // Topology Set
        void upsertTopology(const TopologyTuple& tuple);
        void removeTopologyByDestAndLast(const MainAddress& dest, const MainAddress& last);
        std::vector<TopologyTuple> getTopologyTuples() const;

        // Route Set
        void setRoutes(const std::vector<RouteTuple>& routes);
        const std::vector<RouteTuple>& getRoutes() const;

        // Expire outdated tuples by current time.
        void purgeExpired(double now);

    private:
        struct LinkKey {
            MainAddress local_iface;
            MainAddress neighbor_iface;
            bool operator<(const LinkKey& other) const
            {
                if (local_iface != other.local_iface)
                    return local_iface < other.local_iface;
                return neighbor_iface < other.neighbor_iface;
            }
        };

        struct TopologyKey {
            MainAddress dest;
            MainAddress last;
            bool operator<(const TopologyKey& other) const
            {
                if (dest != other.dest)
                    return dest < other.dest;
                return last < other.last;
            }
        };

        std::map<LinkKey, LinkTuple> link_set_;
        std::map<MainAddress, NeighborTuple> neighbor_set_;
        std::map<TopologyKey, TopologyTuple> topology_set_;
        std::vector<RouteTuple> routing_set_;
    };

} // namespace mysrc::olsrv2

#endif // MYSRC_OLSRV2_OLSRV2STATE_H