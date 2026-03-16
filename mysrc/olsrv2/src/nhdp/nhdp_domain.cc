#include "../../include/nhdp/nhdp_domain.h"

namespace mysrc::olsrv2 {

NhdpDomain::NhdpDomain(const NhdpDomainConfig& config) : config_(config) {}

const NhdpDomainConfig& NhdpDomain::getConfig() const
{
    return config_;
}

NhdpDomainId NhdpDomain::getId() const
{
    return config_.id;
}

void NhdpDomain::setLinkMetric(const NhdpDomainLinkMetric& metric)
{
    link_metrics_[metric.neighbor_iface] = metric;
}

std::optional<NhdpDomainLinkMetric> NhdpDomain::getLinkMetric(const NhdpAddress& neighbor_iface) const
{
    auto it = link_metrics_.find(neighbor_iface);
    if (it == link_metrics_.end())
        return std::nullopt;
    return it->second;
}

void NhdpDomain::removeLinkMetric(const NhdpAddress& neighbor_iface)
{
    link_metrics_.erase(neighbor_iface);
}

std::vector<NhdpDomainLinkMetric> NhdpDomain::getAllLinkMetrics() const
{
    std::vector<NhdpDomainLinkMetric> result;
    result.reserve(link_metrics_.size());
    for (const auto& kv : link_metrics_)
        result.push_back(kv.second);
    return result;
}

void NhdpDomain::setNeighborState(const NhdpDomainNeighborState& state)
{
    neighbor_states_[state.neighbor_main] = state;
}

std::optional<NhdpDomainNeighborState> NhdpDomain::getNeighborState(const NhdpAddress& neighbor_main) const
{
    auto it = neighbor_states_.find(neighbor_main);
    if (it == neighbor_states_.end())
        return std::nullopt;
    return it->second;
}

void NhdpDomain::removeNeighborState(const NhdpAddress& neighbor_main)
{
    neighbor_states_.erase(neighbor_main);
}

std::vector<NhdpDomainNeighborState> NhdpDomain::getAllNeighborStates() const
{
    std::vector<NhdpDomainNeighborState> result;
    result.reserve(neighbor_states_.size());
    for (const auto& kv : neighbor_states_)
        result.push_back(kv.second);
    return result;
}

bool NhdpDomainManager::addDomain(const NhdpDomainConfig& config)
{
    if (config.id >= kNhdpMaxDomainCount)
        return false;
    auto it = domains_.find(config.id);
    if (it != domains_.end())
        return false;
    domains_.emplace(config.id, NhdpDomain(config));
    return true;
}

bool NhdpDomainManager::removeDomain(NhdpDomainId domain_id)
{
    return domains_.erase(domain_id) > 0;
}

std::optional<NhdpDomainConfig> NhdpDomainManager::getDomainConfig(NhdpDomainId domain_id) const
{
    auto it = domains_.find(domain_id);
    if (it == domains_.end())
        return std::nullopt;
    return it->second.getConfig();
}

std::vector<NhdpDomainConfig> NhdpDomainManager::getAllDomainConfigs() const
{
    std::vector<NhdpDomainConfig> result;
    result.reserve(domains_.size());
    for (const auto& kv : domains_)
        result.push_back(kv.second.getConfig());
    return result;
}

NhdpDomain* NhdpDomainManager::findDomain(NhdpDomainId domain_id)
{
    auto it = domains_.find(domain_id);
    if (it == domains_.end())
        return nullptr;
    return &it->second;
}

const NhdpDomain* NhdpDomainManager::findDomain(NhdpDomainId domain_id) const
{
    auto it = domains_.find(domain_id);
    if (it == domains_.end())
        return nullptr;
    return &it->second;
}

} // namespace mysrc::olsrv2

