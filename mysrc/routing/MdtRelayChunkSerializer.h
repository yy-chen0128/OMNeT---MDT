/*
 * MdtRelayChunkSerializer.h
 *
 *  Created on: Oct 26, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_MDTRELAYCHUNKSERIALIZER_H_
#define MYSRC_ROUTING_MDTRELAYCHUNKSERIALIZER_H_

#include "inet/common/packet/serializer/FieldsChunkSerializer.h"

namespace mysrc {
namespace routing {

using namespace inet;

class MdtRelayChunkSerializer : public FieldsChunkSerializer
{
  protected:
    virtual void serialize(MemoryOutputStream& stream, const Ptr<const Chunk>& chunk) const override;
    virtual const Ptr<Chunk> deserialize(MemoryInputStream& stream) const override;
  public:
    MdtRelayChunkSerializer() : FieldsChunkSerializer() {}
};

} // namespace routing
} // namespace mysrc



#endif /* MYSRC_ROUTING_MDTRELAYCHUNKSERIALIZER_H_ */
