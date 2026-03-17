#ifndef MYSRC_OLSRV2_OLSRV2CONTROLPACKETSSERIALIZER_H
#define MYSRC_OLSRV2_OLSRV2CONTROLPACKETSSERIALIZER_H

#include "inet/common/packet/serializer/FieldsChunkSerializer.h"

namespace mysrc::olsrv2 {

using namespace inet;

class INET_API Olsrv2ControlPacketsSerializer : public FieldsChunkSerializer
{
  public:
    using FieldsChunkSerializer::serialize;
    using FieldsChunkSerializer::deserialize;

  protected:
    void serialize(MemoryOutputStream& stream, const Ptr<const Chunk>& chunk) const override;
    const Ptr<Chunk> deserialize(MemoryInputStream& stream) const override;

  public:
    Olsrv2ControlPacketsSerializer() : FieldsChunkSerializer() {}
};

} // namespace mysrc::olsrv2

#endif
