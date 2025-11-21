/*
 * MDTControlPacketsSerializer.h
 *
 *  Created on: Aug 10, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_MDTCONTROLPACKETSSERIALIZER_H_
#define MYSRC_ROUTING_MDTCONTROLPACKETSSERIALIZER_H_

#include "inet/common/packet/serializer/FieldsChunkSerializer.h"

namespace mysrc {
namespace routing {

using namespace inet;

/**
 * Converts between AodvControlPackets and binary (network byte order) Aodv control packets.
 */

class MDTControlPacketsSerializer : public FieldsChunkSerializer
{
  protected:
    virtual void serialize(MemoryOutputStream& stream, const Ptr<const Chunk>& chunk) const override;
    virtual const Ptr<Chunk> deserialize(MemoryInputStream& stream) const override;

  public:
    MDTControlPacketsSerializer() : FieldsChunkSerializer() {}
};

}
}


#endif /* MYSRC_ROUTING_MDTCONTROLPACKETSSERIALIZER_H_ */
