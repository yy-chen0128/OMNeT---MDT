/*
 * OLSR_Repositories.h
 *
 *  Created on: Mar 6, 2026
 *      Author: yychen
 */

#ifndef MYSRC_OLSR_OLSR_REPOSITORIES_H_
#define MYSRC_OLSR_OLSR_REPOSITORIES_H_

#include <omnetpp.h>
#include <vector>
#include <set>
#include <map>
#include "inet/common/INETDefs.h"
#include "inet/networklayer/common/L3AddressTag_m.h"
#include "inet/networklayer/common/L3Address.h"

using namespace omnetpp;
using namespace inet;

// Note: keep everything in a simple namespace to avoid clashes
namespace mysrc {
namespace olsr {
// convenience typedefs
typedef inet::L3Address L3A;
typedef std::vector<L3A> AddrVector;
typedef simtime_t SimTimeType;

// ---------- Link tuple (local iface <-> neighbor iface) ----------
struct OlsrLinkTuple : public cObject {
    L3A localIfaceAddr;    // local interface address (our iface)
    L3A neighborIfaceAddr; // neighbor's interface address
    SimTimeType symTime;   // time until considered symmetric
    SimTimeType asymTime;
    SimTimeType time;      // expiration time for this tuple

    OlsrLinkTuple() : asymTime(0),symTime(0), time(0) {}
    OlsrLinkTuple(const L3A& local, const L3A& nb, SimTimeType expiry)
        : localIfaceAddr(local), neighborIfaceAddr(nb), symTime(expiry), time(expiry) {}
    OlsrLinkTuple(const L3A& local, const L3A& nb)
        : localIfaceAddr(local),
          neighborIfaceAddr(nb),
          asymTime(0),
          symTime(0),
          time(0) {}

    virtual OlsrLinkTuple *dup() const { return new OlsrLinkTuple(*this); }
};

// ---------- Neighbor tuple (one-hop neighbor, main addr) ----------
struct OlsrNeighborTuple : public cObject {
    L3A mainAddr;          // neighbor's main address
    uint8_t status = 0;    // neighbor/link status flags
    uint8_t willingness = 3; // 0..7
    SimTimeType time;      // expiry

    OlsrNeighborTuple() : time(0) {}
    virtual OlsrNeighborTuple *dup() const { return new OlsrNeighborTuple(*this); }
};

// ---------- 2-hop neighbor tuple ----------
struct OlsrTwoHopTuple : public cObject {
    L3A nbMainAddr;   // neighbor main addr (the 1-hop neighbor)
    L3A nb2hopAddr;   // 2-hop neighbor main addr (reachable via nbMainAddr)
    SimTimeType time; // expiry

    OlsrTwoHopTuple() : time(0) {}
    virtual OlsrTwoHopTuple *dup() const { return new OlsrTwoHopTuple(*this); }
};

// ---------- MPR selector tuple (nodes that selected us as MPR) ----------
struct OlsrMprSelectorTuple : public cObject {
    L3A mainAddr;
    SimTimeType time;

    OlsrMprSelectorTuple() : time(0) {}
    virtual OlsrMprSelectorTuple *dup() const { return new OlsrMprSelectorTuple(*this); }
};

// ---------- Topology tuple (for TC processing) ----------
struct OlsrTopologyTuple : public cObject {
    L3A destination;     // destination main addr
    L3A lastAddr;        // last hop (the neighbor advertising this)
    uint16_t seq;        // sequence number (from TC)
    SimTimeType time;    // expiry

    OlsrTopologyTuple() : seq(0), time(0) {}
    virtual OlsrTopologyTuple *dup() const { return new OlsrTopologyTuple(*this); }
};

// ---------- Duplicate tuple (for duplicate suppression) ----------
struct OlsrDupTuple : public cObject {
    L3A originator;
    uint16_t seqNum;       // message sequence number
    SimTimeType time;      // expiry
    bool retransmitted = false;

    OlsrDupTuple() : seqNum(0), time(0), retransmitted(false) {}
    virtual OlsrDupTuple *dup() const { return new OlsrDupTuple(*this); }
};


}//olsr
} // namespace mysrc



#endif /* MYSRC_OLSR_OLSR_REPOSITORIES_H_ */
