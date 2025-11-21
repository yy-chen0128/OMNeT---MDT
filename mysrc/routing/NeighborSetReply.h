/*
 * NeighborSetReply.h
 *
 *  Created on: Sep 11, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_NEIGHBORSETREPLY_H_
#define MYSRC_ROUTING_NEIGHBORSETREPLY_H_

#include "MDTControlPackets_m.h"
#include "NeighborEntry.h"
#include <vector>
namespace inet{
class NeighborSetReply : public NeighborSetReply_Base {
private:
    std::vector<mysrc::routing::NeighborEntry> neighbors;

public:
    NeighborSetReply() : NeighborSetReply_Base() {}
    NeighborSetReply(const NeighborSetReply& other) : NeighborSetReply_Base(other), neighbors(other.neighbors) {}
    NeighborSetReply& operator=(const NeighborSetReply& other) {
        if (this == &other) return *this;
        NeighborSetReply_Base::operator=(other);
        neighbors = other.neighbors;
        return *this;
    }
    virtual ~NeighborSetReply() override {}

    virtual NeighborSetReply *dup() const override { return new NeighborSetReply(*this); }

    const std::vector<mysrc::routing::NeighborEntry>& getNeighbors() const { return neighbors; }
    void setNeighbors(const std::vector<mysrc::routing::NeighborEntry>& n) { neighbors = n; }
    void addNeighbor(const mysrc::routing::NeighborEntry& e) { neighbors.push_back(e); }
    void clearNeighbors() { neighbors.clear(); }
};

Register_Class(NeighborSetReply);

}
#endif /* MYSRC_ROUTING_NEIGHBORSETREPLY_H_ */
