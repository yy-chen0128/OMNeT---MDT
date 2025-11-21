/*
 * MDTControlPacketsSerializer.cc
 *
 *  Created on: Aug 10, 2025
 *      Author: yychen (modified)
 */

#include "../../mysrc/routing/MDTControlPacketsSerializer.h"

#include "../../mysrc/routing/MDTControlPackets_m.h"
#include "../../mysrc/routing/NeighborEntry.h"
#include "JoinReply.h"
#include "NeighborSetReply.h"
#include "CoordsDiscoverReply.h"
#include "KeepAlive.h"

#include "inet/common/packet/serializer/ChunkSerializerRegistry.h"
#include "inet/common/packet/serializer/ChunkSerializer.h"
#include "inet/common/packet/chunk/FieldsChunk.h"

#include <cstdint>
#include <cstring>

namespace mysrc {
namespace routing {
using namespace inet;

Register_Serializer(MDTControlPacket, MDTControlPacketsSerializer);
Register_Serializer(JoinRequest, MDTControlPacketsSerializer);
Register_Serializer(JoinReply, MDTControlPacketsSerializer);
Register_Serializer(NeighborSetRequest, MDTControlPacketsSerializer);
Register_Serializer(NeighborSetReply, MDTControlPacketsSerializer);
Register_Serializer(KeepAlive, MDTControlPacketsSerializer);
Register_Serializer(NeighborSetNotification, MDTControlPacketsSerializer);
Register_Serializer(CoordsDiscoverRequest, MDTControlPacketsSerializer);
Register_Serializer(CoordsDiscoverReply, MDTControlPacketsSerializer);

// helper: write double as uint64 big-endian
static void writeDoubleAsUint64Be(MemoryOutputStream& stream, double v)
{
    uint64_t bits;
    static_assert(sizeof(bits) == sizeof(v), "double size mismatch");
    std::memcpy(&bits, &v, sizeof(v));
    stream.writeUint64Be(bits);
}
static double readDoubleFromUint64Be(MemoryInputStream& stream)
{
    uint64_t bits = stream.readUint64Be();
    double v;
    std::memcpy(&v, &bits, sizeof(v));
    return v;
}

void MDTControlPacketsSerializer::serialize(MemoryOutputStream& stream, const Ptr<const Chunk>& chunk) const
{
    EV_INFO << "use Serializer####################" << "\n";
    const auto &ctrl = staticPtrCast<const MDTControlPacket>(chunk);
    MDTControlPacketType pktType = static_cast<MDTControlPacketType>(ctrl->getPacketType());

    // write the packet type first (1 byte)
    stream.writeByte(static_cast<uint8_t>(pktType));

    switch (pktType) {
        case MDTControlPacketType::JRQ: { // JoinRequest
            const auto &jrq = CHK(dynamicPtrCast<const JoinRequest>(chunk));
            // source + dest IPv4
            stream.writeIpv4Address(jrq->getSourceAddr().toIpv4());
            stream.writeIpv4Address(jrq->getDestAddr().toIpv4());
            // dim + coords
            uint32_t dim = (uint32_t) jrq->getDim();
            stream.writeUint32Be(dim);
            for (uint32_t i = 0; i < dim; ++i) {
                writeDoubleAsUint64Be(stream, jrq->getCoords(i));
            }
            break;
        }

        case MDTControlPacketType::JRP: { // JoinReply (with neighbor list)
            const auto &jrp = CHK(dynamicPtrCast<const JoinReply>(chunk));
            stream.writeIpv4Address(jrp->getSourceAddr().toIpv4());
            stream.writeIpv4Address(jrp->getDestAddr().toIpv4());

            // neighbors: use your custom JoinReply::getNeighbors()
            const auto &neighbors = jrp->getNeighbors(); // const std::vector<mysrc::routing::NeighborEntry>&
            uint32_t count = (uint32_t) neighbors.size();
            stream.writeUint32Be(count);
            for (uint32_t k = 0; k < count; ++k) {
                const mysrc::routing::NeighborEntry &ne = neighbors[k];
                // addr (IPv4)
                stream.writeIpv4Address(ne.addr.toIpv4());
                // attached flag (1 byte)
                stream.writeUint8(ne.attachedToDT ? 1 : 0);
                // dim and coords
                uint32_t ndim = (uint32_t) ne.dim;
                stream.writeUint32Be(ndim);
                for (uint32_t i = 0; i < ndim; ++i) {
                    double v = (i < (uint32_t)ne.coords.size()) ? ne.coords[i] : 0.0;
                    writeDoubleAsUint64Be(stream, v);
                }
                // --- serialize timeout as absolute simtime in milliseconds (uint64 BE) ---
                // if timeout not set (SIMTIME_ZERO), we write 0
                uint64_t timeout_ms = 0;
                if (ne.timeout != SIMTIME_ZERO) {
                    timeout_ms = (uint64_t) ne.timeout.inUnit(SIMTIME_MS);
                }
                stream.writeUint64Be(timeout_ms);
            }
            break;
        }

        case MDTControlPacketType::NSRQ: {
            const auto &ns = CHK(dynamicPtrCast<const NeighborSetRequest>(chunk));
            stream.writeIpv4Address(ns->getSourceAddr().toIpv4());
            stream.writeIpv4Address(ns->getDestAddr().toIpv4());
            break;
        }

        case MDTControlPacketType::NSRP: {
            const auto &nsr = CHK(dynamicPtrCast<const NeighborSetReply>(chunk));
            stream.writeIpv4Address(nsr->getSourceAddr().toIpv4());
            stream.writeIpv4Address(nsr->getDestAddr().toIpv4());
            const auto &neighbors = nsr->getNeighbors(); // std::vector<mysrc::routing::NeighborEntry>
            uint32_t count = (uint32_t) neighbors.size();
            stream.writeUint32Be(count);
            for (uint32_t k = 0; k < count; ++k) {
                const mysrc::routing::NeighborEntry &ne = neighbors[k];
                stream.writeIpv4Address(ne.addr.toIpv4());
                stream.writeUint8(ne.attachedToDT ? 1 : 0);
                uint32_t ndim = (uint32_t) ne.dim;
                stream.writeUint32Be(ndim);
                for (uint32_t i = 0; i < ndim; ++i)
                    writeDoubleAsUint64Be(stream, (i < (uint32_t)ne.coords.size()) ? ne.coords[i] : 0.0);
                // timeout (uint64 BE ms)
                uint64_t timeout_ms = 0;
                if (ne.timeout != SIMTIME_ZERO) {
                    timeout_ms = (uint64_t) ne.timeout.inUnit(SIMTIME_MS);
                }
                stream.writeUint64Be(timeout_ms);
            }
            break;
        }

        case MDTControlPacketType::KA: { // KeepAlive: message coords[] + attached + neighbors (NodeEntry)
            const auto &ka = CHK(dynamicPtrCast<const KeepAlive>(chunk));

            // 1) message-level coords[]
            uint32_t msg_coords_len = static_cast<uint32_t>(ka->getCoordsArraySize());
            stream.writeUint32Be(msg_coords_len);
            for (uint32_t i = 0; i < msg_coords_len; ++i) {
                double v = ka->getCoords(static_cast<int>(i));
                writeDoubleAsUint64Be(stream, v);
            }

            // 2) message-level attached (1 byte)
            stream.writeUint8(static_cast<uint8_t>(ka->getAttached() ? 1 : 0));

            // 3) neighbors vector
            const auto &neighbors = ka->getNeighbors(); // std::vector<mysrc::routing::NodeEntry>
            uint32_t ncount = static_cast<uint32_t>(neighbors.size());
            stream.writeUint32Be(ncount);

            for (uint32_t k = 0; k < ncount; ++k) {
                const mysrc::routing::NodeEntry &ne = neighbors[k];

                // a) addr (IPv4)
                stream.writeIpv4Address(ne.addr.toIpv4());

                // b) attachedToDT (1 byte)
                stream.writeUint8(ne.attachedToDT ? 1 : 0);

                // c) neighbor coords
                uint32_t ncoords_len = static_cast<uint32_t>(ne.coords.size());
                stream.writeUint32Be(ncoords_len);
                for (uint32_t j = 0; j < ncoords_len; ++j) {
                    writeDoubleAsUint64Be(stream, ne.coords[j]);
                }
            }
            break;
        }


        case MDTControlPacketType::NSN: { // NeighborSetNotification
            const auto &nsn = CHK(dynamicPtrCast<const NeighborSetNotification>(chunk));
            // write addresses
            stream.writeIpv4Address(nsn->getSourceAddr().toIpv4());
            stream.writeIpv4Address(nsn->getPredAddr().toIpv4());
            stream.writeIpv4Address(nsn->getDestAddr().toIpv4());
           // write dim and coords
            uint32_t dim = (uint32_t) nsn->getDim();
            stream.writeUint32Be(dim);
            for (uint32_t i = 0; i < dim; ++i) {
                writeDoubleAsUint64Be(stream, nsn->getCoords(i));
            }
            break;
            }
        case MDTControlPacketType::CDQ: { // CoordsDiscoverRequest
            const auto &cdq = CHK(dynamicPtrCast<const CoordsDiscoverRequest>(chunk));
            // write dest, source, target addresses
            //stream.writeIpv4Address(cdq->getDestAddr().toIpv4());
            stream.writeIpv4Address(cdq->getSourceAddr().toIpv4());
            stream.writeIpv4Address(cdq->getTargetAddr().toIpv4());
            // hopcount (uint32)
            stream.writeUint32Be((uint32_t)cdq->getHopcount());
            break;
        }
       case MDTControlPacketType::CDP: { // CoordsDiscoverReply (with one NeighborEntry)
            const auto &cdp = CHK(dynamicPtrCast<const CoordsDiscoverReply>(chunk));
            // addresses: dest, source, target
            stream.writeIpv4Address(cdp->getDestAddr().toIpv4());
            stream.writeIpv4Address(cdp->getSourceAddr().toIpv4());
            stream.writeIpv4Address(cdp->getTargetAddr().toIpv4());
            // serialize the single neighbor (use same format as JRP/NSRP elements)
           const mysrc::routing::NeighborEntry &ne = cdp->getNeighbor();
            // addr (IPv4)
            stream.writeIpv4Address(ne.addr.toIpv4());
            // attached flag
            stream.writeUint8(ne.attachedToDT ? 1 : 0);
            // dim and coords
            uint32_t ndim = (uint32_t) ne.dim;
            stream.writeUint32Be(ndim);
            for (uint32_t i = 0; i < ndim; ++i)
                writeDoubleAsUint64Be(stream, (i < (uint32_t)ne.coords.size()) ? ne.coords[i] : 0.0);
            // timeout (uint64 BE ms)
            uint64_t timeout_ms = 0;
            if (ne.timeout != SIMTIME_ZERO) {
               timeout_ms = (uint64_t) ne.timeout.inUnit(SIMTIME_MS);
            }
            stream.writeUint64Be(timeout_ms);
            break;
        }
        default:
            throw cRuntimeError("MDTControlPacketsSerializer::serialize: unsupported packet type %d", (int)pktType);
    }
}

const Ptr<Chunk> MDTControlPacketsSerializer::deserialize(MemoryInputStream& stream) const
{
    EV_INFO << "use Deserializer###################" << "\n";
    // read the packet type first
    MDTControlPacketType pktType = static_cast<MDTControlPacketType>(stream.readByte());

    switch (pktType) {
        case MDTControlPacketType::JRQ: {
            auto jrq = makeShared<JoinRequest>();
            jrq->setPacketType(pktType);
            jrq->setSourceAddr(L3Address(stream.readIpv4Address()));
            jrq->setDestAddr(L3Address(stream.readIpv4Address()));
            uint32_t dim = stream.readUint32Be();
            jrq->setDim((int)dim);
            jrq->setCoordsArraySize(dim);
            for (uint32_t i = 0; i < dim; ++i) {
                double c = readDoubleFromUint64Be(stream);
                jrq->setCoords(i, c);
            }
            return jrq;
        }

        case MDTControlPacketType::JRP: {
            auto jrp = makeShared<JoinReply>();
            jrp->setPacketType(pktType);
            jrp->setSourceAddr(L3Address(stream.readIpv4Address()));
            jrp->setDestAddr(L3Address(stream.readIpv4Address()));

            // read neighbors
            uint32_t count = stream.readUint32Be();
            for (uint32_t k = 0; k < count; ++k) {
                mysrc::routing::NeighborEntry ne;
                ne.addr = L3Address(stream.readIpv4Address());
                uint8_t att = stream.readUint8();
                ne.attachedToDT = (att != 0);
                uint32_t ndim = stream.readUint32Be();
                ne.dim = (int)ndim;
                ne.coords.resize(ndim);
                for (uint32_t i = 0; i < ndim; ++i)
                    ne.coords[i] = readDoubleFromUint64Be(stream);
                // read timeout as uint64 ms
                uint64_t timeout_ms = stream.readUint64Be();
                if (timeout_ms == 0) {
                    ne.timeout = SIMTIME_ZERO;
                } else {
                    ne.timeout = simtime_t(timeout_ms, SIMTIME_MS);
                }
                // add to JoinReply's neighbors vector
                jrp->addNeighbor(ne);
            }
            return jrp;
        }

        case MDTControlPacketType::NSRQ: {
            auto ns = makeShared<NeighborSetRequest>();
            ns->setPacketType(pktType);
            ns->setSourceAddr(L3Address(stream.readIpv4Address()));
            ns->setDestAddr(L3Address(stream.readIpv4Address()));
            return ns;
        }

        case MDTControlPacketType::NSRP: {
            auto nsr = makeShared<NeighborSetReply>();
            nsr->setPacketType(pktType);
            nsr->setSourceAddr(L3Address(stream.readIpv4Address()));
            nsr->setDestAddr(L3Address(stream.readIpv4Address()));
            uint32_t count = stream.readUint32Be();
            for (uint32_t k = 0; k < count; ++k) {
                mysrc::routing::NeighborEntry ne;
                ne.addr = L3Address(stream.readIpv4Address());
                ne.attachedToDT = (stream.readUint8() != 0);
                uint32_t ndim = stream.readUint32Be();
                ne.dim = (int)ndim;
                ne.coords.resize(ndim);
                for (uint32_t i = 0; i < ndim; ++i)
                    ne.coords[i] = readDoubleFromUint64Be(stream);
                // read timeout (uint64 BE ms)
                uint64_t timeout_ms = stream.readUint64Be();
                if (timeout_ms == 0) {
                    ne.timeout = SIMTIME_ZERO;
                } else {
                    ne.timeout = simtime_t(timeout_ms, SIMTIME_MS);
                }
                nsr->addNeighbor(ne);
            }
            return nsr;
        }

        case MDTControlPacketType::KA: {
            auto ka = makeShared<KeepAlive>();
            ka->setPacketType(pktType);

            // 1) message-level coords length + elements
            uint32_t msg_coords_len = stream.readUint32Be();
            const uint32_t MAX_MSG_COORDS = 65536;
            if (msg_coords_len > MAX_MSG_COORDS) {
                throw cRuntimeError("KeepAlive deserialization: msg_coords_len %u exceeds MAX_MSG_COORDS", msg_coords_len);
            }
            ka->setCoordsArraySize(static_cast<int>(msg_coords_len));
            for (uint32_t i = 0; i < msg_coords_len; ++i) {
                double v = readDoubleFromUint64Be(stream);
                ka->setCoords(static_cast<int>(i), v);
            }

            // 2) message-level attached
            uint8_t matt = stream.readUint8();
            ka->setAttached(matt != 0);

            // 3) neighbors: read count, then each NodeEntry
            uint32_t ncount = stream.readUint32Be();
            const uint32_t MAX_NEIGHBORS = 16384;
            if (ncount > MAX_NEIGHBORS) {
                throw cRuntimeError("KeepAlive deserialization: neighbors count %u exceeds MAX_NEIGHBORS", ncount);
            }

            for (uint32_t k = 0; k < ncount; ++k) {
                mysrc::routing::NodeEntry ne;

                // a) addr
                ne.addr = L3Address(stream.readIpv4Address());

                // b) attachedToDT
                uint8_t natt = stream.readUint8();
                ne.attachedToDT = (natt != 0);

                // c) neighbor coords len + elements
                uint32_t ncoords_len = stream.readUint32Be();
                const uint32_t MAX_NEIGHBOR_COORDS = 65536;
                if (ncoords_len > MAX_NEIGHBOR_COORDS) {
                    throw cRuntimeError("KeepAlive deserialization: neighbor coords_len %u exceeds MAX_NEIGHBOR_COORDS", ncoords_len);
                }
                ne.coords.resize(ncoords_len);
                for (uint32_t j = 0; j < ncoords_len; ++j) {
                    ne.coords[j] = readDoubleFromUint64Be(stream);
                }

                // add to KeepAlive's neighbors
                ka->addNeighbor(ne);
            }

            return ka;
        }

        case MDTControlPacketType::NSN: {
            auto nsn = makeShared<NeighborSetNotification>();
            nsn->setPacketType(pktType);
            // read addresses (source, pred, dest)
            nsn->setSourceAddr(L3Address(stream.readIpv4Address()));
            nsn->setPredAddr(L3Address(stream.readIpv4Address()));
            nsn->setDestAddr(L3Address(stream.readIpv4Address()));
            // dim + coords
            uint32_t dim = stream.readUint32Be();
            nsn->setDim((int)dim);
            nsn->setCoordsArraySize(dim);
            for (uint32_t i = 0; i < dim; ++i) {
                double c = readDoubleFromUint64Be(stream);
                nsn->setCoords(i, c);
            }
            return nsn;
        }
        case MDTControlPacketType::CDQ: {
            auto cdq = makeShared<CoordsDiscoverRequest>();
            cdq->setPacketType(pktType);
            //cdq->setDestAddr(L3Address(stream.readIpv4Address()));
            cdq->setSourceAddr(L3Address(stream.readIpv4Address()));
            cdq->setTargetAddr(L3Address(stream.readIpv4Address()));
            uint32_t hop = stream.readUint32Be();
            cdq->setHopcount((int)hop);
            return cdq;
        }
        case MDTControlPacketType::CDP: {
            auto cdp = makeShared<CoordsDiscoverReply>();
            cdp->setPacketType(pktType);
            cdp->setDestAddr(L3Address(stream.readIpv4Address()));
            cdp->setSourceAddr(L3Address(stream.readIpv4Address()));
            cdp->setTargetAddr(L3Address(stream.readIpv4Address()));
            // read single neighbor entry
            mysrc::routing::NeighborEntry ne;
            ne.addr = L3Address(stream.readIpv4Address());
            ne.attachedToDT = (stream.readUint8() != 0);
            uint32_t ndim = stream.readUint32Be();
            ne.dim = (int)ndim;
            ne.coords.resize(ndim);
            for (uint32_t i = 0; i < ndim; ++i)
                ne.coords[i] = readDoubleFromUint64Be(stream);
            uint64_t timeout_ms = stream.readUint64Be();
            if (timeout_ms == 0) ne.timeout = SIMTIME_ZERO;
            else ne.timeout = simtime_t(timeout_ms, SIMTIME_MS);
            cdp->setNeighbor(ne);
            return cdp;
        }
        default:
            throw cRuntimeError("MDTControlPacketsSerializer::deserialize: unsupported packet type %d", (int)pktType);
    }
}

} // namespace routing
} // namespace mysrc



