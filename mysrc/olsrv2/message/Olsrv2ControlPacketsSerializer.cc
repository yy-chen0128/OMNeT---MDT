#include "Olsrv2ControlPacketsSerializer.h"

#include "../message/OLSRv2Packet_m.h"

#include "inet/common/packet/serializer/ChunkSerializerRegistry.h"

namespace mysrc::olsrv2 {

using namespace inet;

Register_Serializer(Olsrv2HelloPacket, Olsrv2ControlPacketsSerializer);
Register_Serializer(Olsrv2TcGroup, Olsrv2ControlPacketsSerializer);

void Olsrv2ControlPacketsSerializer::serialize(MemoryOutputStream& stream, const Ptr<const Chunk>& chunk) const
{
    const auto& pkt = staticPtrCast<const Olsrv2ControlPacket>(chunk);

    switch (pkt->getPacketType()) {
        case Olsrv2Hello: {
            const auto& hello = staticPtrCast<const Olsrv2HelloPacket>(chunk);

            stream.writeByte(static_cast<uint8_t>(hello->getPacketType()));
            stream.writeIpv4Address(hello->getOriginator().toIpv4());
            stream.writeUint64Be(static_cast<uint64_t>(hello->getHTime().inUnit(SIMTIME_MS)));
            stream.writeUint64Be(static_cast<uint64_t>(hello->getValidityTime().inUnit(SIMTIME_MS)));
            stream.writeByte(hello->getWillingness());
            stream.writeUint16Be(hello->getMsgSeq());

            const auto n = static_cast<uint16_t>(hello->getNeighAddrsArraySize());
            stream.writeUint16Be(n);

            for (uint16_t i = 0; i < n; ++i) {
                stream.writeIpv4Address(hello->getNeighAddrs(i).toIpv4());
                stream.writeByte(i < hello->getLinkStatusArraySize() ? hello->getLinkStatus(i) : 0);
                stream.writeByte(i < hello->getNeighStatusArraySize() ? hello->getNeighStatus(i) : 0);
                stream.writeUint16Be(static_cast<uint16_t>(i < hello->getAddrCountsArraySize() ? hello->getAddrCounts(i) : 0));
            }
            break;
        }

        case Olsrv2Tc: {
            const auto& tc = staticPtrCast<const Olsrv2TcGroup>(chunk);

            stream.writeByte(static_cast<uint8_t>(tc->getPacketType()));
            stream.writeUint64Be(static_cast<uint64_t>(tc->getGenTime().inUnit(SIMTIME_MS)));
            stream.writeUint64Be(static_cast<uint64_t>(tc->getValidityTime().inUnit(SIMTIME_MS)));

            const auto m = static_cast<uint16_t>(tc->getOriginatorsArraySize());
            stream.writeUint16Be(m);

            for (uint16_t i = 0; i < m; ++i) {
                stream.writeIpv4Address(tc->getOriginators(i).toIpv4());
                stream.writeUint16Be(tc->getAnsns(i));
                stream.writeUint16Be(tc->getSeqNums(i));
                stream.writeByte(tc->getHopCounts(i));
                stream.writeUint32Be(static_cast<uint32_t>(tc->getAdvCounts(i)));
            }

            const auto n = static_cast<uint32_t>(tc->getAdvertisedNeighborsArraySize());
            stream.writeUint32Be(n);
            for (uint32_t i = 0; i < n; ++i)
                stream.writeIpv4Address(tc->getAdvertisedNeighbors(i).toIpv4());

            break;
        }

        default:
            throw cRuntimeError("Unknown OLSRv2 packet type %d", pkt->getPacketType());
    }
}

const Ptr<Chunk> Olsrv2ControlPacketsSerializer::deserialize(MemoryInputStream& stream) const
{
    uint8_t type = stream.readByte();

    switch (static_cast<Olsrv2PktType>(type)) {
        case Olsrv2Hello: {
            auto hello = makeShared<Olsrv2HelloPacket>();
            hello->setPacketType(Olsrv2Hello);
            hello->setOriginator(L3Address(stream.readIpv4Address()));
            hello->setHTime(SimTime(stream.readUint64Be(), SIMTIME_MS));
            hello->setValidityTime(SimTime(stream.readUint64Be(), SIMTIME_MS));
            hello->setWillingness(stream.readByte());
            hello->setMsgSeq(stream.readUint16Be());

            uint16_t n = stream.readUint16Be();
            hello->setNeighAddrsArraySize(n);
            hello->setLinkStatusArraySize(n);
            hello->setNeighStatusArraySize(n);
            hello->setAddrCountsArraySize(n);

            for (uint16_t i = 0; i < n; ++i) {
                hello->setNeighAddrs(i, L3Address(stream.readIpv4Address()));
                hello->setLinkStatus(i, stream.readByte());
                hello->setNeighStatus(i, stream.readByte());
                hello->setAddrCounts(i, static_cast<int>(stream.readUint16Be()));
            }

            return hello;
        }

        case Olsrv2Tc: {
            auto tc = makeShared<Olsrv2TcGroup>();
            tc->setPacketType(Olsrv2Tc);
            tc->setGenTime(SimTime(stream.readUint64Be(), SIMTIME_MS));
            tc->setValidityTime(SimTime(stream.readUint64Be(), SIMTIME_MS));

            uint16_t m = stream.readUint16Be();
            tc->setOriginatorsArraySize(m);
            tc->setAnsnsArraySize(m);
            tc->setSeqNumsArraySize(m);
            tc->setHopCountsArraySize(m);
            tc->setAdvCountsArraySize(m);

            for (uint16_t i = 0; i < m; ++i) {
                tc->setOriginators(i, L3Address(stream.readIpv4Address()));
                tc->setAnsns(i, stream.readUint16Be());
                tc->setSeqNums(i, stream.readUint16Be());
                tc->setHopCounts(i, stream.readByte());
                tc->setAdvCounts(i, static_cast<int>(stream.readUint32Be()));
            }

            uint32_t n = stream.readUint32Be();
            tc->setAdvertisedNeighborsArraySize(n);
            for (uint32_t i = 0; i < n; ++i)
                tc->setAdvertisedNeighbors(i, L3Address(stream.readIpv4Address()));

            return tc;
        }

        default:
            throw cRuntimeError("Unknown OLSRv2 packet type %d", type);
    }
}

} // namespace mysrc::olsrv2

