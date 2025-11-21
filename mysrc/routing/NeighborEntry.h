/*
 * NeighborEntry.h
 *
 *  Created on: Sep 10, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_NEIGHBORENTRY_H_
#define MYSRC_ROUTING_NEIGHBORENTRY_H_
#include <omnetpp.h>
#include <vector>
#include "inet/networklayer/common/L3Address.h"

#include "inet/common/INETDefs.h"
namespace mysrc {
namespace routing {
using namespace inet;
struct NeighborEntry {
    L3Address addr;

    simtime_t timeout;
    bool attachedToDT = false;
    std::vector<double> coords;
    int dim;
    std::string type;
    //std::map<std::pair<L3Address,L3Address>, Routetuple> tuples; // key: (src,dst)
};
struct NodeEntry{//for ka message only
    L3Address addr;
    std::vector<double> coords;
    bool attachedToDT = false;
};
struct TxPowerEntry {
    double power_W;
    simtime_t ts;
};
struct TxInfo { int txid; double power_W; L3Address src; simtime_t ts; };
struct RecentRxEntry {
    //int txid = -1;               // transmission id if available
    double power_W = NAN;        // received power in W
    simtime_t endTime = SIMTIME_ZERO;
    //int transmitterId = -1;      // transmitter radio id if available
    //Coord startPosition;         // optional for disambiguation
    //std::string transmitterPath; // module path of transmitter if available
};

}
}


#endif /* MYSRC_ROUTING_NEIGHBORENTRY_H_ */
