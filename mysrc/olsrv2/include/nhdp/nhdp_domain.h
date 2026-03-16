//
// Created by GuoYudi on 2026/3/10.
//

#ifndef MYSRC_OLSRV2_NHDP_DOMIN_H_
#define MYSRC_OLSRV2_NHDP_DOMIN_H_

#pragma once

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

inline constexpr NhdpDomainId kNhdpDefaultDomainId = 0;
inline constexpr NhdpDomainId kNhdpMaxDomainCount = 4;
inline constexpr const char* kNhdpNoMetricMpr = "-";
inline constexpr const char* kNhdpAnyMetricMpr = "*";

enum class NhdpMetricType : uint8_t {
    HopCount = 0,
    Etx = 1,
    Custom = 2,
};

enum class NhdpMprType : uint8_t {
    Flooding = 0,
    Routing = 1,
    Disabled = 2,
};

struct NhdpDomainConfig {
    NhdpDomainId id = kNhdpDefaultDomainId;
    std::string name;
    NhdpMetricType metric_type = NhdpMetricType::HopCount;
    NhdpMprType mpr_type = NhdpMprType::Flooding;
    uint32_t metric_infinity = UINT32_MAX;
    uint32_t metric_minimum = 1;
    uint32_t metric_maximum = UINT32_MAX;
    bool use_flooding = true;
    bool use_routing = true;
};

struct NhdpDomainLinkMetric {
    NhdpAddress neighbor_iface = NhdpAddress::UNSPECIFIED_ADDRESS;
    uint32_t in_metric = 0;
    uint32_t out_metric = 0;
    bool metric_is_valid = false;
};

struct NhdpDomainNeighborState {
    NhdpAddress neighbor_main = NhdpAddress::UNSPECIFIED_ADDRESS;
    bool local_is_mpr = false;
    bool neigh_is_mpr = false;
    uint8_t willingness = 0;
    uint32_t best_out_metric = UINT32_MAX;
    std::set<NhdpAddress> contributing_ifaces;
};

class NhdpDomain
{
  public:
    explicit NhdpDomain(const NhdpDomainConfig& config);

    const NhdpDomainConfig& getConfig() const;
    NhdpDomainId getId() const;

    // Link-metric view.
    void setLinkMetric(const NhdpDomainLinkMetric& metric);
    std::optional<NhdpDomainLinkMetric> getLinkMetric(const NhdpAddress& neighbor_iface) const;
    void removeLinkMetric(const NhdpAddress& neighbor_iface);
    std::vector<NhdpDomainLinkMetric> getAllLinkMetrics() const;

    // Neighbor per-domain state.
    void setNeighborState(const NhdpDomainNeighborState& state);
    std::optional<NhdpDomainNeighborState> getNeighborState(const NhdpAddress& neighbor_main) const;
    void removeNeighborState(const NhdpAddress& neighbor_main);
    std::vector<NhdpDomainNeighborState> getAllNeighborStates() const;

  private:
    NhdpDomainConfig config_;
    std::map<NhdpAddress, NhdpDomainLinkMetric> link_metrics_;
    std::map<NhdpAddress, NhdpDomainNeighborState> neighbor_states_;
};

class NhdpDomainManager
{
  public:
    bool addDomain(const NhdpDomainConfig& config);
    bool removeDomain(NhdpDomainId domain_id);

    std::optional<NhdpDomainConfig> getDomainConfig(NhdpDomainId domain_id) const;
    std::vector<NhdpDomainConfig> getAllDomainConfigs() const;

    NhdpDomain* findDomain(NhdpDomainId domain_id);
    const NhdpDomain* findDomain(NhdpDomainId domain_id) const;

  private:
    std::map<NhdpDomainId, NhdpDomain> domains_;
};

} // namespace mysrc::olsrv2

#endif //MYSRC_OLSRV2_NHDP_DOMIN_H_