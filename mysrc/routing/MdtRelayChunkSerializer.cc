/*
 * MdtRelayChunkSerializer.cc
 *
 *  Created on: Oct 26, 2025
 *      Author: yychen
 */


#include "MdtRelayChunkSerializer.h"
#include "MdtRelayChunk_m.h"
#include "inet/common/packet/serializer/ChunkSerializerRegistry.h"

#include "inet/common/packet/serializer/ChunkSerializer.h"
#include "inet/common/packet/chunk/FieldsChunk.h"

#include <cstdint>
#include <cstring>

namespace mysrc {
namespace routing {

using namespace inet;

Register_Serializer(MdtRelayChunk, MdtRelayChunkSerializer);

void MdtRelayChunkSerializer::serialize(MemoryOutputStream& stream, const Ptr<const Chunk>& chunk) const
{
    //const auto &c = staticPtrCast<const MdtRelayChunk>(chunk);
    const auto &c = CHK(dynamicPtrCast<const MdtRelayChunk>(chunk));
    stream.writeByte(c->getHasRelay() ? 1 : 0);
    //stream.writeByte(c->getHasRelay());

    //IPv4Address relay = c->getRelayAddr().toIpv4();
    stream.writeIpv4Address(c->getRelayAddr().toIpv4());

    //IPv4Address pred = c->getPredAddr().toIpv4();
    stream.writeIpv4Address(c->getPredAddr().toIpv4());
}

const Ptr<Chunk> MdtRelayChunkSerializer::deserialize(MemoryInputStream& stream) const
{
    auto chunk = makeShared<MdtRelayChunk>();

    //bool has = stream.readByte();
    uint8_t has = stream.readByte();
    chunk->setHasRelay(static_cast<bool>(has));
    //chunk->setHasRelay(has);

    chunk->setRelayAddr(L3Address(stream.readIpv4Address()));
    chunk->setPredAddr(L3Address(stream.readIpv4Address()));

    //IPv4Address relay = stream.readIpv4Address();
    //IPv4Address pred  = stream.readIpv4Address();
    //chunk->setRelayAddr(L3Address(relay));
    //chunk->setPredAddr(L3Address(pred));

    chunk->setChunkLength(B(1 + 4 + 4)); // 9 bytes

    return chunk;
}

} // namespace routing
} // namespace mysrc

