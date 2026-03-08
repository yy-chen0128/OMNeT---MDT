/*
 * OlsrControlPacketsSerializer.cc
 *
 *  Created on: Mar 7, 2026
 *      Author: yychen
 */


#include "../olsr/OlsrControlPacketsSerializer.h"

#include "../olsr/OlsrControlPackets_m.h"
#include "inet/common/packet/serializer/ChunkSerializerRegistry.h"

namespace mysrc{
namespace olsr {
using namespace inet;

Register_Serializer(OlsrHello, OlsrControlPacketsSerializer);
Register_Serializer(OlsrTcGroup, OlsrControlPacketsSerializer);

void OlsrControlPacketsSerializer::serialize(
        MemoryOutputStream& stream,
        const Ptr<const Chunk>& chunk) const
{
    const auto& olsrPacket =
        staticPtrCast<const OlsrControlPacket>(chunk);

    switch (olsrPacket->getPacketType()) {

        case Hello: {

            const auto& hello =
                staticPtrCast<const OlsrHello>(chunk);

            stream.writeByte(hello->getPacketType());

            stream.writeIpv4Address(
                hello->getOriginator().toIpv4());

            stream.writeUint64Be(
                hello->getHTime().inUnit(SIMTIME_MS));

            stream.writeByte(hello->getWillingness());

            stream.writeUint16Be(hello->getMsgSeq());

            uint8_t linkMsgCount =
                hello->getLinkCodesArraySize();

            stream.writeByte(linkMsgCount);

            for (uint8_t i=0;i<linkMsgCount;i++) {

                stream.writeByte(hello->getLinkCodes(i));

                stream.writeUint16Be(
                    hello->getLinkCounts(i));
            }

            int addrCount =
                hello->getNbIfaceAddrsArraySize();

            for (int i=0;i<addrCount;i++)
                stream.writeIpv4Address(
                    hello->getNbIfaceAddrs(i).toIpv4());

            break;
        }

        case Tc: {

            const auto& tc =
                staticPtrCast<const OlsrTcGroup>(chunk);

            stream.writeByte(tc->getPacketType());

            stream.writeUint64Be(
                tc->getGenTime().inUnit(SIMTIME_MS));

            int M = tc->getOriginatorsArraySize();

            stream.writeUint16Be(M);

            for (int i=0;i<M;i++) {

                stream.writeIpv4Address(
                    tc->getOriginators(i).toIpv4());

                stream.writeUint16Be(tc->getAnsns(i));

                stream.writeUint16Be(tc->getSeqNums(i));

                stream.writeByte(tc->getHopCounts(i));

                stream.writeUint32Be(tc->getAdvCounts(i));
            }

            int N = tc->getAdvertisedNeighborsArraySize();

            stream.writeUint32Be(N);

            for (int i=0;i<N;i++)
                stream.writeIpv4Address(
                    tc->getAdvertisedNeighbors(i).toIpv4());

            break;
        }

        default:
            throw cRuntimeError(
                "Unknown OLSR packet type %d",
                olsrPacket->getPacketType());
    }
}

const Ptr<Chunk> OlsrControlPacketsSerializer::deserialize(
        MemoryInputStream& stream) const
{
    uint8_t type = stream.readByte();

    switch(type) {

        case Hello: {

            auto hello = makeShared<OlsrHello>();

            hello->setPacketType(Hello);

            hello->setOriginator(
                L3Address(stream.readIpv4Address()));

            hello->setHTime(
                SimTime(stream.readUint64Be(),
                        SIMTIME_MS));

            hello->setWillingness(stream.readByte());

            hello->setMsgSeq(stream.readUint16Be());

            uint8_t linkMsgCount = stream.readByte();

            hello->setLinkCodesArraySize(linkMsgCount);
            hello->setLinkCountsArraySize(linkMsgCount);

            int totalAddr = 0;

            for(int i=0;i<linkMsgCount;i++) {

                hello->setLinkCodes(i,
                    stream.readByte());

                uint16_t count =
                    stream.readUint16Be();

                hello->setLinkCounts(i,count);

                totalAddr += count;
            }

            hello->setNbIfaceAddrsArraySize(totalAddr);

            for(int i=0;i<totalAddr;i++)
                hello->setNbIfaceAddrs(i,
                    L3Address(stream.readIpv4Address()));

            return hello;
        }

        case Tc: {

            auto tc = makeShared<OlsrTcGroup>();

            tc->setPacketType(Tc);

            tc->setGenTime(
                SimTime(stream.readUint64Be(),
                        SIMTIME_MS));

            int M = stream.readUint16Be();

            tc->setOriginatorsArraySize(M);
            tc->setAnsnsArraySize(M);
            tc->setSeqNumsArraySize(M);
            tc->setHopCountsArraySize(M);
            tc->setAdvCountsArraySize(M);

            for(int i=0;i<M;i++) {

                tc->setOriginators(i,
                    L3Address(stream.readIpv4Address()));

                tc->setAnsns(i,
                    stream.readUint16Be());

                tc->setSeqNums(i,
                    stream.readUint16Be());

                tc->setHopCounts(i,
                    stream.readByte());

                tc->setAdvCounts(i,
                    stream.readUint32Be());
            }

            int N = stream.readUint32Be();

            tc->setAdvertisedNeighborsArraySize(N);

            for(int i=0;i<N;i++)
                tc->setAdvertisedNeighbors(i,
                    L3Address(stream.readIpv4Address()));

            return tc;
        }

        default: {

            throw cRuntimeError(
                "Unknown OLSR packet type %d", type);
        }
    }
}

}
}
