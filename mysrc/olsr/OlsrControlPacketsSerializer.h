/*
 * OlsrControlPacketsSerializer.h
 *
 *  Created on: Mar 7, 2026
 *      Author: yychen
 */

#ifndef MYSRC_OLSR_OLSRCONTROLPACKETSSERIALIZER_H_
#define MYSRC_OLSR_OLSRCONTROLPACKETSSERIALIZER_H_

#include "inet/common/packet/serializer/FieldsChunkSerializer.h"

namespace mysrc {
namespace olsr {
using namespace inet;

class INET_API OlsrControlPacketsSerializer : public FieldsChunkSerializer
{
  protected:
    virtual void serialize(MemoryOutputStream& stream, const Ptr<const Chunk>& chunk) const override;
    virtual const Ptr<Chunk> deserialize(MemoryInputStream& stream) const override;

  public:
    OlsrControlPacketsSerializer() : FieldsChunkSerializer() {}
};

}//olsr
} // namespace inet



#endif /* MYSRC_OLSR_OLSRCONTROLPACKETSSERIALIZER_H_ */
