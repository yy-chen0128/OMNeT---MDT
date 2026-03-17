//
// Created by GuoYudi on 2026/3/14.
//
#include "../include/olsrv2_state.h"

#include <algorithm>

#include "inet/common/INETDefs.h"

using namespace omnetpp;

namespace mysrc::olsrv2 {

void Olsrv2State::upsertLink(const LinkTuple& tuple)
{
    link_set_[LinkKey{tuple.local_iface_addr, tuple.neighbor_iface_addr}] = tuple;
}

void Olsrv2State::removeLink(const MainAddress& local_iface, const MainAddress& neighbor_iface)
{
    link_set_.erase(LinkKey{local_iface, neighbor_iface});
}

std::vector<LinkTuple> Olsrv2State::getLinks() const
{
    std::vector<LinkTuple> out;
    out.reserve(link_set_.size());
    for (const auto& kv : link_set_)
        out.push_back(kv.second);
    return out;
}

void Olsrv2State::upsertNeighbor(const NeighborTuple& tuple)
{
    neighbor_set_[tuple.main_addr] = tuple;
}

void Olsrv2State::removeNeighbor(const MainAddress& main_addr)
{
    neighbor_set_.erase(main_addr);
}

std::optional<NeighborTuple> Olsrv2State::findNeighbor(const MainAddress& main_addr) const
{
    auto it = neighbor_set_.find(main_addr);
    if (it == neighbor_set_.end())
        return std::nullopt;
    return it->second;
}

std::vector<NeighborTuple> Olsrv2State::getNeighbors() const
{
    std::vector<NeighborTuple> out;
    out.reserve(neighbor_set_.size());
    for (const auto& kv : neighbor_set_)
        out.push_back(kv.second);
    return out;
}

void Olsrv2State::upsertTopology(const TopologyTuple& tuple)
{
    topology_set_[TopologyKey{tuple.dest_addr, tuple.last_addr}] = tuple;
}

void Olsrv2State::removeTopologyByDestAndLast(const MainAddress& dest, const MainAddress& last)
{
    topology_set_.erase(TopologyKey{dest, last});
}

std::vector<TopologyTuple> Olsrv2State::getTopologyTuples() const
{
    std::vector<TopologyTuple> out;
    out.reserve(topology_set_.size());
    for (const auto& kv : topology_set_)
        out.push_back(kv.second);
    return out;
}

void Olsrv2State::setRoutes(const std::vector<RouteTuple>& routes)
{
    routing_set_ = routes;
}

const std::vector<RouteTuple>& Olsrv2State::getRoutes() const
{
    return routing_set_;
}

void Olsrv2State::purgeExpired(double now)
{
    for (auto it = link_set_.begin(); it != link_set_.end();) {
        if (!it->second.sym_time.isValid(now) && !it->second.heard_time.isValid(now))
            it = link_set_.erase(it);
        else
            ++it;
    }

    for (auto it = neighbor_set_.begin(); it != neighbor_set_.end();) {
        if (!it->second.validity.isValid(now))
            it = neighbor_set_.erase(it);
        else
            ++it;
    }

    for (auto it = topology_set_.begin(); it != topology_set_.end();) {
        if (!it->second.validity.isValid(now))
            it = topology_set_.erase(it);
        else
            ++it;
    }
}

void Olsrv2State::printLinkSet() const
{
    EV_INFO << "\n=== LINK SET ===\n";
    for (const auto& kv : link_set_) {
        const auto& t = kv.second;
        EV_INFO << "Neighbor=" << t.neighbor_iface_addr
                << " Local=" << t.local_iface_addr
                << " SymExpire=" << t.sym_time.expires_at
                << " HeardExpire=" << t.heard_time.expires_at
                << omnetpp::endl;
    }
    EV_INFO << "LinkSet size=" << link_set_.size() << "\n";
}

void Olsrv2State::printNeighborSet() const
{
    EV_INFO << "\n=== NEIGHBOR SET ===\n";
    for (const auto& kv : neighbor_set_) {
        const auto& t = kv.second;
        EV_INFO << "NeighborMainAddr=" << t.main_addr
                << " Symmetric=" << (t.is_symmetric ? 1 : 0)
                << " Ifaces=" << t.iface_addrs.size()
                << " Expire=" << t.validity.expires_at
                << omnetpp::endl;
    }
    EV_INFO << "NeighborSet size=" << neighbor_set_.size() << "\n";
}

void Olsrv2State::printTopologySet() const
{
    EV_INFO << "\n=== TOPOLOGY SET ===\n";
    for (const auto& kv : topology_set_) {
        const auto& t = kv.second;
        EV_INFO << "Dest=" << t.dest_addr
                << " LastHop=" << t.last_addr
                << " Ansn=" << t.ansn
                << " Expire=" << t.validity.expires_at
                << omnetpp::endl;
    }
    EV_INFO << "TopologySet size=" << topology_set_.size() << "\n";
}

void Olsrv2State::printRouteSet() const
{
    EV_INFO << "\n=== ROUTE SET ===\n";
    for (const auto& t : routing_set_) {
        EV_INFO << "Dest=" << t.destination
                << " NextHop=" << t.next_hop
                << " Metric=" << t.metric
                << " HopCount=" << static_cast<int>(t.hop_count)
                << omnetpp::endl;
    }
    EV_INFO << "RouteSet size=" << routing_set_.size() << "\n";
}

void Olsrv2State::printOlsrv2State() const
{
    EV_INFO << "\n\n================ OLSRV2 STATE ================\n";
    printLinkSet();
    printNeighborSet();
    printTopologySet();
    printRouteSet();
    EV_INFO << "==============================================\n\n";
}

} // namespace mysrc::olsrv2
