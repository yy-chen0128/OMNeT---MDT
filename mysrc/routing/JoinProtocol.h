/*
 * MDT_Join.h
 *
 *  Created on: Aug 31, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_JOINPROTOCOL_H_
#define MYSRC_ROUTING_JOINPROTOCOL_H_
#include <omnetpp.h>

//#include "MDTDataManager.h"
//#include "Coordinate.h"
#include "inet/common/INETDefs.h"
#include "inet/common/packet/Packet.h"
#include "MDTControlPackets_m.h"
#include "Neighbor.h"
#include "MDTData.h"
#include "JoinReply.h"

namespace mysrc {
namespace routing {
using namespace inet;

class MDTRouting;

class INET_API JoinProtocol : public cSimpleModule{
private:

    MDTRouting *core = nullptr;
    //Coordinate nodeCoord;
    virtual void initialize(int stage) override;
    virtual void handleMessage(cMessage *msg) override;
public:
    JoinProtocol(){}
    virtual ~JoinProtocol(){}

    Packet *createJoinRequestTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<double>& coords = std::vector<double>());
    bool isNearestDTNode(const std::vector<double>& dest,int dim);//TODO this method is no use

    void processJoinRequestChunk(const Ptr<JoinRequest>& joinRequestChunk);

    Packet *createJoinReplyTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<NeighborEntry>& entrys);
    void processJoinReplyChunk(const Ptr<JoinReply>& joinReplyChunk);

};

}
}

#endif /* MYSRC_ROUTING_JOINPROTOCOL_H_ */
