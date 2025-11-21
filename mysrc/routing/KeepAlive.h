/*
 * KeepAlive.h
 *
 *  Created on: Sep 29, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_KEEPALIVE_H_
#define MYSRC_ROUTING_KEEPALIVE_H_

#include "MDTControlPackets_m.h"
#include "NeighborEntry.h"
namespace inet{
using namespace inet;
// Ensure the base symbol is visible: JoinReply_Base should be declared in MDTControlPackets_m.h
class KeepAlive : public KeepAlive_Base {
private:
    std::vector<mysrc::routing::NodeEntry> neighbors;

public:
    // default ctor forwarding to base
    KeepAlive() : KeepAlive_Base() {}

    // copy ctor / assignment to satisfy duplication semantics
    KeepAlive(const KeepAlive& other) : KeepAlive_Base(other), neighbors(other.neighbors) {}
    KeepAlive& operator=(const KeepAlive& other) {
        if (this == &other) return *this;
        KeepAlive_Base::operator=(other);
        neighbors = other.neighbors;
        return *this;
    }

    // dup() override -- OMNeT++/msg classes often use dup() to clone objects
    virtual KeepAlive *dup() const override { return new KeepAlive(*this); }

    virtual ~KeepAlive() override {}

    // neighbors accessor
    const std::vector<mysrc::routing::NodeEntry>& getNeighbors() const { return neighbors; }
    void setNeighbors(const std::vector<mysrc::routing::NodeEntry>& n) { neighbors = n; }
    void addNeighbor(const mysrc::routing::NodeEntry& e) { neighbors.push_back(e); }

    // optional: convenience to clear
    void clearNeighbors() { neighbors.clear(); }
};
Register_Class(KeepAlive);
}



#endif /* MYSRC_ROUTING_KEEPALIVE_H_ */
