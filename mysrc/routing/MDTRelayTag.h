/*
 * MDTRelayTag.h
 *
 *  Created on: Sep 1, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_MDTRELAYTAG_H_
#define MYSRC_ROUTING_MDTRELAYTAG_H_
#include <inet/common/TagBase.h>
#include <inet/networklayer/common/L3Address.h>
namespace inet {
class INET_API MDTRelayTag : public TagBase{
public:
    L3Address relayAddr;
    static const auto& getTagName() { static const std::string name = "MDTRelayTag"; return name; }
    static const auto& getTagKind() { static const int kind = 1; return kind; }

    void parsimPack(cCommBuffer* b) const override { b->pack(relayAddr); }
    void parsimUnpack(cCommBuffer* b) override { b->unpack(relayAddr); }
};

}

#endif /* MYSRC_ROUTING_MDTRELAYTAG_H_ */
