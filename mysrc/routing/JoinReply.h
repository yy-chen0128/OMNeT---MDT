/*
 * JoinReply.h
 *
 *  Created on: Sep 11, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_JOINREPLY_H_
#define MYSRC_ROUTING_JOINREPLY_H_

#include "MDTControlPackets_m.h"
#include "NeighborEntry.h"
namespace inet{
using namespace inet;
// Ensure the base symbol is visible: JoinReply_Base should be declared in MDTControlPackets_m.h
class JoinReply : public JoinReply_Base {
private:
    std::vector<mysrc::routing::NeighborEntry> neighbors;

public:
    // default ctor forwarding to base
    JoinReply() : JoinReply_Base() {}

    // copy ctor / assignment to satisfy duplication semantics
    JoinReply(const JoinReply& other) : JoinReply_Base(other), neighbors(other.neighbors) {}
    JoinReply& operator=(const JoinReply& other) {
        if (this == &other) return *this;
        JoinReply_Base::operator=(other);
        neighbors = other.neighbors;
        return *this;
    }

    // dup() override -- OMNeT++/msg classes often use dup() to clone objects
    virtual JoinReply *dup() const override { return new JoinReply(*this); }

    virtual ~JoinReply() override {}

    // neighbors accessor
    const std::vector<mysrc::routing::NeighborEntry>& getNeighbors() const { return neighbors; }
    void setNeighbors(const std::vector<mysrc::routing::NeighborEntry>& n) { neighbors = n; }
    void addNeighbor(const mysrc::routing::NeighborEntry& e) { neighbors.push_back(e); }

    // optional: convenience to clear
    void clearNeighbors() { neighbors.clear(); }
};


Register_Class(JoinReply);
}
#endif /* MYSRC_ROUTING_JOINREPLY_H_ */
