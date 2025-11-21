/*
 * Forward_Protocol.h
 *
 *  Created on: Aug 31, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_FORWARDPROTOCOL_H_
#define MYSRC_ROUTING_FORWARDPROTOCOL_H_

#include <omnetpp.h>
#include "MDTControlPackets_m.h"

#include "MdtRelayTag_m.h"
#include "Neighbor.h"

//using namespace omnetpp;
namespace mysrc {
namespace routing {
using namespace inet;
class MDTRouting;
enum class ForwardAction {
    NONE,
    FORWARD_TO_NEXT_HOP,
    DROP,
    WAIT_ROUTE_DISCOVERY
};
struct ForwardDecision {
    ForwardAction action = ForwardAction::NONE;
    L3Address nextHop;
    bool useMdtRelay = false;
    L3Address relayAddr = L3Address();
    std::string reason;
};
class ForwardProtocol : public cSimpleModule
{

  protected:
    MDTRouting *core = nullptr;

    virtual void initialize(int stage) override;
    virtual void handleMessage(cMessage *msg) override;
  public:
    ForwardProtocol() {}
    virtual ~ForwardProtocol() {}

    /** Called by MDTRouting (or by datagramForwardHook) when this node should decide how to forward a packet.
     *  Returns true if the ForwardProtocol took action and forwarded the packet (caller should not forward).
     *  Returns false if the ForwardProtocol did NOT forward it (caller may accept/drop).
     */
    ForwardDecision decide(Packet *packet);
    bool findnexthop(const L3Address& dest,L3Address& res,L3Address ban = L3Address());//return true for greedy false for MDT
    bool shortcutforMDTrelay(const L3Address& dest,L3Address& relay,L3Address ban = L3Address());
};

}
}
#endif /* MYSRC_ROUTING_FORWARDPROTOCOL_H_ */
