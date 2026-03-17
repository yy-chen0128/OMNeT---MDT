//
// Created by GuoYudi on 2026/3/10.
//

#ifndef MYSRC_OLSRV2_NHDP_DB_H_
#define MYSRC_OLSRV2_NHDP_DB_H_

#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "inet/networklayer/common/L3Address.h"

namespace mysrc::olsrv2 {

using NhdpAddress = inet::Ipv4Address;
using NhdpDomainId = uint8_t;
using NhdpTimeMs = uint64_t;

inline constexpr uint8_t kNhdpMaximumDomains = 4;

// Link status aligned with RFC6130 semantic ordering.
enum class NhdpLinkStatus : int8_t {
    Pending = -1,
    Lost = 0,
    Symmetric = 1,
    Heard = 2,
};

struct NhdpMetric {
    uint32_t incoming = 0;
    uint32_t outgoing = 0;
};

struct NhdpLinkDomainData {
    NhdpMetric metric;
    NhdpTimeMs last_metric_change_ms = 0;
};

struct NhdpNeighborDomainData {
    NhdpMetric metric;
    NhdpAddress best_out_next_hop = NhdpAddress::UNSPECIFIED_ADDRESS;
    uint32_t best_out_link_metric = 0;
    int best_link_ifindex = -1;
    bool local_is_mpr = false;
    bool neigh_is_mpr = false;
    bool neigh_was_mpr = false;
    uint8_t willingness = 0;
};

struct NhdpTwoHopDomainData {
    NhdpMetric metric;
    uint32_t last_used_outgoing_metric = 0;
};

struct NhdpLink {
    uint64_t id = 0;
    int local_ifindex = -1;
    NhdpAddress neighbor_iface_addr = NhdpAddress::UNSPECIFIED_ADDRESS;
    NhdpTimeMs vtime_ms = 0;
    NhdpTimeMs itime_ms = 0;
    NhdpTimeMs sym_expire_ms = 0;
    NhdpTimeMs heard_expire_ms = 0;
    NhdpTimeMs expire_ms = 0;
    NhdpLinkStatus last_status = NhdpLinkStatus::Pending;
    NhdpLinkStatus status = NhdpLinkStatus::Pending;
    NhdpAddress neighbor_originator = NhdpAddress::UNSPECIFIED_ADDRESS;
    bool local_is_flooding_mpr = false;
    bool neigh_is_flooding_mpr = false;
    uint8_t flooding_willingness = 0;
    bool neigh_was_flooding_mpr = false;
    int process_count = 0;
    std::map<NhdpDomainId, NhdpLinkDomainData> domain_data;
    std::set<NhdpAddress> link_addresses;
    std::set<NhdpAddress> twohop_addresses;
};

struct NhdpLinkAddress {
    NhdpAddress address = NhdpAddress::UNSPECIFIED_ADDRESS;
    uint64_t link_id = 0;
    bool might_be_removed = false;
};

struct NhdpTwoHop {
    NhdpAddress twohop_addr = NhdpAddress::UNSPECIFIED_ADDRESS;
    uint64_t link_id = 0;
    bool same_interface = false;
    NhdpTimeMs expire_ms = 0;
    std::map<NhdpDomainId, NhdpTwoHopDomainData> domain_data;
};

struct NhdpNeighbor {
    NhdpAddress originator = NhdpAddress::UNSPECIFIED_ADDRESS;
    int symmetric_link_count = 0;
    int process_count = 0;
    bool selection_is_mpr = false;
    NhdpAddress old_originator = NhdpAddress::UNSPECIFIED_ADDRESS;
    std::set<uint64_t> link_ids;
    std::set<NhdpAddress> neighbor_addresses;
    std::set<NhdpAddress> link_addresses;
    std::map<NhdpDomainId, NhdpNeighborDomainData> domain_data;
};

struct NhdpNeighborAddress {
    NhdpAddress neigh_addr = NhdpAddress::UNSPECIFIED_ADDRESS;
    NhdpAddress neighbor_originator = NhdpAddress::UNSPECIFIED_ADDRESS;
    int link_addr_count = 0;
    NhdpTimeMs lost_expire_ms = 0;
    bool this_if = false;
    bool might_be_removed = false;
};

class NhdpDb
{
  public:
    void init();
    void cleanup();

    // Neighbor CRUD.
    void addOrUpdateNeighbor(const NhdpNeighbor& neighbor);
    void removeNeighbor(const NhdpAddress& originator);
    std::optional<NhdpNeighbor> findNeighbor(const NhdpAddress& originator) const;
    std::vector<NhdpNeighbor> getNeighbors() const;

    // Neighbor address CRUD.
    void addOrUpdateNeighborAddress(const NhdpNeighborAddress& naddr);
    void removeNeighborAddress(const NhdpAddress& neigh_addr);
    std::optional<NhdpNeighborAddress> findNeighborAddress(const NhdpAddress& neigh_addr) const;
    std::vector<NhdpNeighborAddress> getNeighborAddresses() const;

    // Link CRUD.
    void addOrUpdateLink(const NhdpLink& link);
    void removeLink(uint64_t link_id);
    std::optional<NhdpLink> findLink(uint64_t link_id) const;
    std::vector<NhdpLink> getLinks() const;

    // Link address CRUD.
    void addOrUpdateLinkAddress(const NhdpLinkAddress& laddr);
    void removeLinkAddress(const NhdpAddress& addr);
    std::optional<NhdpLinkAddress> findLinkAddress(const NhdpAddress& addr) const;

    // Two-hop CRUD.
    void addOrUpdateTwoHop(const NhdpTwoHop& l2hop);
    void removeTwoHop(const NhdpAddress& addr, uint64_t link_id);
    std::optional<NhdpTwoHop> findTwoHop(const NhdpAddress& addr, uint64_t link_id) const;
    std::vector<NhdpTwoHop> getTwoHops() const;

    // Timer/state helpers.
    void setLinkVtime(uint64_t link_id, NhdpTimeMs vtime_ms);
    void setLinkHeardTime(uint64_t link_id, NhdpTimeMs htime_ms);
    void setLinkSymTime(uint64_t link_id, NhdpTimeMs stime_ms);
    void setTwoHopVtime(const NhdpAddress& addr, uint64_t link_id, NhdpTimeMs vtime_ms);
    void setNeighborAddressLost(const NhdpAddress& neigh_addr, NhdpTimeMs vtime_ms);
    void clearNeighborAddressLost(const NhdpAddress& neigh_addr);
    bool isNeighborAddressLost(const NhdpAddress& neigh_addr, NhdpTimeMs now_ms) const;
    bool isTwoHopLost(const NhdpAddress& addr, uint64_t link_id, NhdpTimeMs now_ms) const;

    // Periodic cleanup.
    void purgeExpired(NhdpTimeMs now_ms);

    void addOrUpdateMprSelector(const NhdpAddress& selector_originator, NhdpTimeMs expire_ms);
    void removeMprSelector(const NhdpAddress& selector_originator);
    bool isMprSelector(const NhdpAddress& selector_originator, NhdpTimeMs now_ms) const;
    std::vector<NhdpAddress> getMprSelectors() const;

  private:
    struct TwoHopKey {
        NhdpAddress addr;
        uint64_t link_id = 0;
        bool operator<(const TwoHopKey& other) const
        {
            if (addr != other.addr)
                return addr < other.addr;
            return link_id < other.link_id;
        }
    };

    std::map<NhdpAddress, NhdpNeighbor> neighbors_;
    std::map<NhdpAddress, NhdpNeighborAddress> neighbor_addresses_;
    std::map<uint64_t, NhdpLink> links_;
    std::map<NhdpAddress, NhdpLinkAddress> link_addresses_;
    std::map<TwoHopKey, NhdpTwoHop> twohop_;
    std::map<NhdpAddress, NhdpTimeMs> mpr_selectors_;
};

} // namespace mysrc::olsrv2

#endif // MYSRC_OLSRV2_NHDP_DB_H_
