#include "../../include/nhdp/nhdp_db.h"

#include <algorithm>

namespace mysrc::olsrv2 {

void NhdpDb::init()
{
    cleanup();
}

void NhdpDb::cleanup()
{
    neighbors_.clear();
    neighbor_addresses_.clear();
    links_.clear();
    link_addresses_.clear();
    twohop_.clear();
}

void NhdpDb::addOrUpdateNeighbor(const NhdpNeighbor& neighbor)
{
    neighbors_[neighbor.originator] = neighbor;
}

void NhdpDb::removeNeighbor(const NhdpAddress& originator)
{
    neighbors_.erase(originator);
}

std::optional<NhdpNeighbor> NhdpDb::findNeighbor(const NhdpAddress& originator) const
{
    auto it = neighbors_.find(originator);
    if (it == neighbors_.end())
        return std::nullopt;
    return it->second;
}

std::vector<NhdpNeighbor> NhdpDb::getNeighbors() const
{
    std::vector<NhdpNeighbor> result;
    result.reserve(neighbors_.size());
    for (const auto& kv : neighbors_)
        result.push_back(kv.second);
    return result;
}

void NhdpDb::addOrUpdateNeighborAddress(const NhdpNeighborAddress& naddr)
{
    neighbor_addresses_[naddr.neigh_addr] = naddr;
}

void NhdpDb::removeNeighborAddress(const NhdpAddress& neigh_addr)
{
    neighbor_addresses_.erase(neigh_addr);
}

std::optional<NhdpNeighborAddress> NhdpDb::findNeighborAddress(const NhdpAddress& neigh_addr) const
{
    auto it = neighbor_addresses_.find(neigh_addr);
    if (it == neighbor_addresses_.end())
        return std::nullopt;
    return it->second;
}

std::vector<NhdpNeighborAddress> NhdpDb::getNeighborAddresses() const
{
    std::vector<NhdpNeighborAddress> result;
    result.reserve(neighbor_addresses_.size());
    for (const auto& kv : neighbor_addresses_)
        result.push_back(kv.second);
    return result;
}

void NhdpDb::addOrUpdateLink(const NhdpLink& link)
{
    links_[link.id] = link;
}

void NhdpDb::removeLink(uint64_t link_id)
{
    links_.erase(link_id);

    for (auto it = link_addresses_.begin(); it != link_addresses_.end();) {
        if (it->second.link_id == link_id)
            it = link_addresses_.erase(it);
        else
            ++it;
    }

    for (auto it = twohop_.begin(); it != twohop_.end();) {
        if (it->first.link_id == link_id)
            it = twohop_.erase(it);
        else
            ++it;
    }
}

std::optional<NhdpLink> NhdpDb::findLink(uint64_t link_id) const
{
    auto it = links_.find(link_id);
    if (it == links_.end())
        return std::nullopt;
    return it->second;
}

std::vector<NhdpLink> NhdpDb::getLinks() const
{
    std::vector<NhdpLink> result;
    result.reserve(links_.size());
    for (const auto& kv : links_)
        result.push_back(kv.second);
    return result;
}

void NhdpDb::addOrUpdateLinkAddress(const NhdpLinkAddress& laddr)
{
    link_addresses_[laddr.address] = laddr;
}

void NhdpDb::removeLinkAddress(const NhdpAddress& addr)
{
    link_addresses_.erase(addr);
}

std::optional<NhdpLinkAddress> NhdpDb::findLinkAddress(const NhdpAddress& addr) const
{
    auto it = link_addresses_.find(addr);
    if (it == link_addresses_.end())
        return std::nullopt;
    return it->second;
}

void NhdpDb::addOrUpdateTwoHop(const NhdpTwoHop& l2hop)
{
    TwoHopKey key{l2hop.twohop_addr, l2hop.link_id};
    twohop_[key] = l2hop;
}

void NhdpDb::removeTwoHop(const NhdpAddress& addr, uint64_t link_id)
{
    twohop_.erase(TwoHopKey{addr, link_id});
}

std::optional<NhdpTwoHop> NhdpDb::findTwoHop(const NhdpAddress& addr, uint64_t link_id) const
{
    auto it = twohop_.find(TwoHopKey{addr, link_id});
    if (it == twohop_.end())
        return std::nullopt;
    return it->second;
}

std::vector<NhdpTwoHop> NhdpDb::getTwoHops() const
{
    std::vector<NhdpTwoHop> result;
    result.reserve(twohop_.size());
    for (const auto& kv : twohop_)
        result.push_back(kv.second);
    return result;
}

void NhdpDb::setLinkVtime(uint64_t link_id, NhdpTimeMs vtime_ms)
{
    auto it = links_.find(link_id);
    if (it != links_.end())
        it->second.expire_ms = vtime_ms;
}

void NhdpDb::setLinkHeardTime(uint64_t link_id, NhdpTimeMs htime_ms)
{
    auto it = links_.find(link_id);
    if (it != links_.end())
        it->second.heard_expire_ms = htime_ms;
}

void NhdpDb::setLinkSymTime(uint64_t link_id, NhdpTimeMs stime_ms)
{
    auto it = links_.find(link_id);
    if (it != links_.end())
        it->second.sym_expire_ms = stime_ms;
}

void NhdpDb::setTwoHopVtime(const NhdpAddress& addr, uint64_t link_id, NhdpTimeMs vtime_ms)
{
    auto it = twohop_.find(TwoHopKey{addr, link_id});
    if (it != twohop_.end())
        it->second.expire_ms = vtime_ms;
}

void NhdpDb::setNeighborAddressLost(const NhdpAddress& neigh_addr, NhdpTimeMs vtime_ms)
{
    auto it = neighbor_addresses_.find(neigh_addr);
    if (it != neighbor_addresses_.end())
        it->second.lost_expire_ms = vtime_ms;
}

void NhdpDb::clearNeighborAddressLost(const NhdpAddress& neigh_addr)
{
    auto it = neighbor_addresses_.find(neigh_addr);
    if (it != neighbor_addresses_.end())
        it->second.lost_expire_ms = 0;
}

bool NhdpDb::isNeighborAddressLost(const NhdpAddress& neigh_addr, NhdpTimeMs now_ms) const
{
    auto it = neighbor_addresses_.find(neigh_addr);
    if (it == neighbor_addresses_.end())
        return false;
    return it->second.lost_expire_ms > 0 && now_ms >= it->second.lost_expire_ms;
}

bool NhdpDb::isTwoHopLost(const NhdpAddress& addr, uint64_t link_id, NhdpTimeMs now_ms) const
{
    auto it = twohop_.find(TwoHopKey{addr, link_id});
    if (it == twohop_.end())
        return false;
    return it->second.expire_ms > 0 && now_ms >= it->second.expire_ms;
}

void NhdpDb::purgeExpired(NhdpTimeMs now_ms)
{
    for (auto it = twohop_.begin(); it != twohop_.end();) {
        if (it->second.expire_ms > 0 && now_ms >= it->second.expire_ms)
            it = twohop_.erase(it);
        else
            ++it;
    }

    for (auto it = neighbor_addresses_.begin(); it != neighbor_addresses_.end();) {
        if (it->second.lost_expire_ms > 0 && now_ms >= it->second.lost_expire_ms)
            it = neighbor_addresses_.erase(it);
        else
            ++it;
    }

    std::vector<uint64_t> expired_links;
    for (const auto& kv : links_) {
        if (kv.second.expire_ms > 0 && now_ms >= kv.second.expire_ms)
            expired_links.push_back(kv.first);
    }
    for (uint64_t link_id : expired_links)
        removeLink(link_id);
}

} // namespace mysrc::olsrv2

