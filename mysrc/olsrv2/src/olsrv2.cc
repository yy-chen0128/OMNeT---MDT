// mysrc/olsrv2/src/olsrv2.cc
//
// Created by GuoYudi on 2026/3/10.
//

#include "../include/olsrv2.h"

#include <algorithm>
#include <cstdint>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

#include "inet/networklayer/common/L3Address.h"

namespace mysrc::olsrv2 {

namespace {

// Runtime tunables (milliseconds, aligned with old OONF-style API semantics).
uint64_t g_tc_interval_ms = 5000;
uint64_t g_tc_validity_ms = 300000;
uint64_t g_overwrite_tc_interval_ms = 0;
uint64_t g_overwrite_tc_validity_ms = 0;
bool g_generate_tcs = true;

// Very small duplicate cache for process/forward decision in simulation.
struct DupKey {
    MainAddress originator;
    uint8_t msg_type = 0;
    uint16_t seqno = 0;

    bool operator<(const DupKey& other) const
    {
        if (originator != other.originator)
            return originator < other.originator;
        if (msg_type != other.msg_type)
            return msg_type < other.msg_type;
        return seqno < other.seqno;
    }
};

std::set<DupKey> g_processed_set;
std::set<DupKey> g_forwarded_set;

} // namespace

uint64_t Olsrv2::get_tc_interval(void)
{
    if (g_overwrite_tc_interval_ms != 0)
        return g_overwrite_tc_interval_ms;
    return g_tc_interval_ms;
}

uint64_t Olsrv2::get_tc_validity(void)
{
    if (g_overwrite_tc_validity_ms != 0)
        return g_overwrite_tc_validity_ms;
    return g_tc_validity_ms;
}

bool Olsrv2::is_nhdp_routable(struct netaddr* addr)
{
    return is_routable(addr);
}

bool Olsrv2::is_routable(struct netaddr* addr)
{
    return addr != nullptr;
}

bool Olsrv2::mpr_shall_process(struct rfc5444_reader_tlvblock_context* context, uint64_t /*vtime*/)
{
    if (context == nullptr)
        return false;
    DupKey key{MainAddress::UNSPECIFIED_ADDRESS, 0, 0};
    auto [_, inserted] = g_processed_set.insert(key);
    return inserted;
}

bool Olsrv2::mpr_shall_forwarding(
    struct rfc5444_reader_tlvblock_context* context, struct netaddr* /*source_address*/, uint64_t /*vtime*/)
{
    if (context == nullptr)
        return false;
    return true;
}

void Olsrv2::generate_tcs(bool generate)
{
    g_generate_tcs = generate;
}

// Olsrv2Core Implementation

void Olsrv2Core::processTcAndRecompute(const inet::Olsrv2TcGroup* tc, double now)
{
    if (!tc) return;

    TcMessageView view;
    // Iterate over originators in the TC Group (assuming flattened structure allows multiple TCs?)
    // The .msg defines arrays: originators[], ansns[], seqNums[], hopCounts[], advCounts[], advertisedNeighbors[]
    // This looks like a batch of TCs.
    
    size_t offset = 0;
    for (size_t i = 0; i < tc->getOriginatorsArraySize(); ++i) {
        view.originator = tc->getOriginators(i).toIpv4();
        view.ansn = tc->getAnsns(i);
        view.validity_duration = 300.0; // Default validity (should be in packet)
        
        view.advertised_neighbors.clear();
        int count = tc->getAdvCounts(i);
        for (int j = 0; j < count; ++j) {
            TcAdvertisedNeighbor adv;
            adv.neighbor_main_addr = tc->getAdvertisedNeighbors(offset + j).toIpv4();
            view.advertised_neighbors.push_back(adv);
        }
        offset += count;

        tc_processor_.applyTc(view, now, state_);
    }

    state_.setRoutes(routing_engine_.computeRoutes(originator_, state_));
}

inet::Packet* Olsrv2Core::generateTc(const NhdpDb& nhdpDb, double now)
{
    // Check if we need to generate TC (MPR selection)
    bool am_i_mpr = false;
    // TODO: Check NhdpDb for MPR selector status
    // For now, assume always generate if we have neighbors
    // In strict RFC 7181, check if I am selected as MPR by anyone.
    
    // Iterate over links to find if local_is_flooding_mpr is true
    // nhdpDb doesn't expose links directly in public API efficiently maybe?
    // nhdpDb has getNeighbors().
    // But local_is_flooding_mpr is on Link.
    
    // We assume we advertise all symmetric neighbors for now (Full topology).
    // RFC 7181 allows full topology.
    
    std::vector<MainAddress> advertised_neighbors;
    // Get all symmetric neighbors
    // Note: NhdpDb uses internal structures. We need to access them.
    // Assuming NhdpDb has a way to get neighbors.
    // Let's use getNeighbors() and assume it returns symmetric ones or check status.
    // NhdpNeighbor doesn't have status, NhdpLink has.
    
    // This part requires better NhdpDb access. 
    // For now, dummy implementation to return a packet.
    
    auto pkt = inet::makeShared<inet::Olsrv2TcGroup>();
    pkt->setChunkLength(inet::B(100));
    pkt->setPacketType(inet::Olsrv2Tc);
    pkt->setGenTime(now);

    pkt->setOriginatorsArraySize(1);
    pkt->setOriginators(0, originator_);
    pkt->setAnsnsArraySize(1);
    pkt->setAnsns(0, ansn_++);
    pkt->setSeqNumsArraySize(1);
    pkt->setSeqNums(0, tcSeqNum_++);
    pkt->setHopCountsArraySize(1);
    pkt->setHopCounts(0, 0);
    
    // Add neighbors
    // Dummy: No neighbors added yet because of NhdpDb access.
    pkt->setAdvCountsArraySize(1);
    pkt->setAdvCounts(0, 0);
    
    auto packet = new inet::Packet("OLSRv2-TC");
    packet->insertAtBack(pkt);
    return packet;
}

} // namespace mysrc::olsrv2
