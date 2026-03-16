#include "../../include/nhdp/nhdp.h"
#include "../../include/nhdp/nhdp_db.h"
#include "../../message/OLSRv2Packet_m.h"
#include "inet/common/packet/Packet.h"
#include "inet/common/packet/chunk/Chunk.h"

namespace mysrc::olsrv2 {

Nhdp::Nhdp()
{
    db_.init();
}

void Nhdp::setOriginator(const inet::L3Address& originator)
{
    originator_ = originator;
}

inet::L3Address Nhdp::getOriginator() const
{
    return originator_;
}

bool Nhdp::floodingSelector(
    rfc5444_writer *writer,
    rfc5444_writer_target *target,
    void *ptr)
{
    (void)writer;
    (void)target;
    (void)ptr;
    return true;
}

bool Nhdp::forwardingSelector(
    rfc5444_writer_target *target,
    rfc5444_reader_tlvblock_context *context)
{
    (void)target;
    (void)context;
    return true;
}

void Nhdp::processHello(const inet::Olsrv2HelloPacket *hello, const inet::L3Address& source, double now)
{
    // Basic NHDP processing
    NhdpNeighbor neigh;
    neigh.originator = hello->getOriginator().toIpv4(); 
    
    // Update Neighbor in DB
    db_.addOrUpdateNeighbor(neigh);
    
    // Link processing
    NhdpLink link;
    bool new_link = true;
    
    // Find existing link by source address
    auto linkAddrOpt = db_.findLinkAddress(source.toIpv4());
    if (linkAddrOpt) {
        auto linkOpt = db_.findLink(linkAddrOpt->link_id);
        if (linkOpt) {
            link = *linkOpt;
            new_link = false;
        }
    }
    
    if (new_link) {
        link.id = nextLinkId_++;
        link.local_ifindex = 0; // Default interface
        link.neighbor_iface_addr = source.toIpv4();
        
        // Add link address mapping
        NhdpLinkAddress la;
        la.address = source.toIpv4();
        la.link_id = link.id;
        db_.addOrUpdateLinkAddress(la);
    }
    
    link.neighbor_originator = hello->getOriginator().toIpv4();

    // Check if my address is in the Hello packet (Symmetric check)
    bool symmetric = false;
    inet::L3Address myOriginator = originator_; 
    
    link.twohop_addresses.clear();

    for (size_t i = 0; i < hello->getNeighAddrsArraySize(); ++i) {
        inet::L3Address neighAddr = hello->getNeighAddrs(i);
        if (neighAddr == myOriginator) {
            symmetric = true; 
        } else {
            // Add to 2-hop neighbors
            link.twohop_addresses.insert(neighAddr.toIpv4());
            
            // Update NhdpTwoHop in DB
            NhdpTwoHop th;
            th.twohop_addr = neighAddr.toIpv4();
            th.link_id = link.id;
            th.expire_ms = static_cast<NhdpTimeMs>((now + hello->getHTime().dbl() * 3.0) * 1000);
            db_.addOrUpdateTwoHop(th);
        }
    }
    
    link.status = symmetric ? NhdpLinkStatus::Symmetric : NhdpLinkStatus::Heard;
    
    // Time is in ms
    link.expire_ms = static_cast<NhdpTimeMs>((now + hello->getHTime().dbl() * 3.0) * 1000);
    
    db_.addOrUpdateLink(link);
}

std::set<inet::Ipv4Address> Nhdp::calculateMprSet() const
{
    std::set<inet::Ipv4Address> mprSet;
    std::set<inet::Ipv4Address> N2;
    std::map<inet::Ipv4Address, std::set<inet::Ipv4Address>> coverage; // neighbor -> {2-hop neighbors}

    auto links = db_.getLinks();
    std::set<inet::Ipv4Address> N1;

    // 1. Identify N1 and N2
    for (const auto& link : links) {
        if (link.status == NhdpLinkStatus::Symmetric) {
            inet::Ipv4Address n1 = link.neighbor_iface_addr; // Use interface address for link
            N1.insert(n1);
            
            for (const auto& n2 : link.twohop_addresses) {
                if (n2 != originator_.toIpv4() && N1.find(n2) == N1.end()) {
                    N2.insert(n2);
                    coverage[n1].insert(n2);
                }
            }
        }
    }
    
    // Clean up N2: remove nodes that are in N1 (direct neighbors)
    // (Already done in loop above: N1.find(n2) == N1.end())
    // But we need to be careful about order.
    // Re-check N2 against N1 fully.
    for (auto it = N2.begin(); it != N2.end(); ) {
        if (N1.count(*it)) {
            it = N2.erase(it);
        } else {
            ++it;
        }
    }

    // Update coverage maps to remove non-N2 nodes
    for (auto& kv : coverage) {
        std::set<inet::Ipv4Address> clean;
        for (const auto& addr : kv.second) {
            if (N2.count(addr)) clean.insert(addr);
        }
        kv.second = clean;
    }

    // 2. Greedy Selection
    // a. Select nodes that provide unique coverage
    // (Simplification: skip unique coverage optimization for now, go straight to max coverage)

    while (!N2.empty()) {
        inet::Ipv4Address best_n1;
        size_t max_cover = 0;
        
        for (const auto& kv : coverage) {
            size_t cover = 0;
            for (const auto& n2 : kv.second) {
                if (N2.count(n2)) cover++;
            }
            
            if (cover > max_cover) {
                max_cover = cover;
                best_n1 = kv.first;
            }
        }
        
        if (max_cover == 0) break; // Should not happen
        
        mprSet.insert(best_n1);
        
        // Remove covered
        for (const auto& n2 : coverage[best_n1]) {
            N2.erase(n2);
        }
    }
    
    return mprSet;
}

inet::Packet* Nhdp::generateHello(double now)
{
    auto pkt = inet::makeShared<inet::Olsrv2HelloPacket>();
    pkt->setChunkLength(inet::B(100)); // Arbitrary size
    pkt->setPacketType(inet::Olsrv2Hello);
    pkt->setOriginator(originator_);
    pkt->setHTime(2.0); // 2 seconds
    pkt->setWillingness(7); // Default
    pkt->setMsgSeq(helloSeqNum_++);

    // Calculate MPRs
    std::set<inet::Ipv4Address> mprSet = calculateMprSet();

    auto links = db_.getLinks();
    pkt->setNeighAddrsArraySize(links.size());
    pkt->setLinkStatusArraySize(links.size());
    pkt->setNeighStatusArraySize(links.size());

    for (size_t i = 0; i < links.size(); ++i) {
        const auto& link = links[i];
        pkt->setNeighAddrs(i, inet::L3Address(link.neighbor_iface_addr));
        
        // Link Status
        if (link.status == NhdpLinkStatus::Symmetric) {
            pkt->setLinkStatus(i, 2); // SYMMETRIC
        } else if (link.status == NhdpLinkStatus::Heard) {
            pkt->setLinkStatus(i, 1); // HEARD
        } else {
            pkt->setLinkStatus(i, 0); // LOST/PENDING
        }
        
        // Neighbor Status (MPR)
        if (mprSet.count(link.neighbor_iface_addr)) {
            pkt->setNeighStatus(i, 1); // MPR
        } else {
            pkt->setNeighStatus(i, 0); // Not MPR
        }
    }

    auto packet = new inet::Packet("OLSRv2-Hello");
    packet->insertAtBack(pkt);
    return packet;
}

} // namespace mysrc::olsrv2
