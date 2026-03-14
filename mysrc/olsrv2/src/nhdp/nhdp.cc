#include "../../include/nhdp/nhdp.h"
#include "../../include/nhdp/nhdp_db.h"
#ifdef UNIT_TEST
#include "message/OLSRv2Packet_m.h"
#else
#include "../../message/OLSRv2Packet_m.h"
#endif
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
    
    for (size_t i = 0; i < hello->getNeighAddrsArraySize(); ++i) {
        if (hello->getNeighAddrs(i) == myOriginator) {
            symmetric = true; 
            break;
        }
    }
    
    link.status = symmetric ? NhdpLinkStatus::Symmetric : NhdpLinkStatus::Heard;
    
    // Time is in ms
    link.expire_ms = static_cast<NhdpTimeMs>((now + hello->getHTime().dbl() * 3.0) * 1000);
    
    db_.addOrUpdateLink(link);
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

    // Add neighbors
    // Iterate over links/neighbors in DB
    // Populate linkStatus, neighStatus, neighAddrs
    
    std::vector<inet::L3Address> neighbors;
    std::vector<uint8_t> statuses;
    
    // Dummy implementation: Add all known neighbors as HEARD/SYMMETRIC
    // In a real implementation, we iterate db_.getLinks()
    
    pkt->setNeighAddrsArraySize(neighbors.size());
    for (size_t i = 0; i < neighbors.size(); ++i) {
        pkt->setNeighAddrs(i, neighbors[i]);
        // Set status...
    }

    auto packet = new inet::Packet("OLSRv2-Hello");
    packet->insertAtBack(pkt);
    return packet;
}

} // namespace mysrc::olsrv2
