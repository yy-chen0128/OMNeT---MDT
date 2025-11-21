/*
 * CoordsDiscoverReply.h
 *
 *  Created on: Sep 19, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_COORDSDISCOVERREPLY_H_
#define MYSRC_ROUTING_COORDSDISCOVERREPLY_H_

#include "MDTControlPackets_m.h"
#include "NeighborEntry.h"

namespace inet{
class CoordsDiscoverReply : public CoordsDiscoverReply_Base{
private:
    mysrc::routing::NeighborEntry neighbor;

public:
    CoordsDiscoverReply() : CoordsDiscoverReply_Base() {}
    CoordsDiscoverReply(const CoordsDiscoverReply& other) : CoordsDiscoverReply_Base(other), neighbor(other.neighbor) {}
    CoordsDiscoverReply& operator=(const CoordsDiscoverReply& other) {
        if (this == &other) return *this;
        CoordsDiscoverReply_Base::operator=(other);
        neighbor = other.neighbor;
        return *this;
    }
    virtual ~CoordsDiscoverReply() override {}

    virtual CoordsDiscoverReply *dup() const override { return new CoordsDiscoverReply(*this); }

    const mysrc::routing::NeighborEntry& getNeighbor() const { return neighbor; }
    void setNeighbor(const mysrc::routing::NeighborEntry& n) { neighbor = n; }
    //void addNeighbor(const mysrc::routing::NeighborEntry& e) { neighbors.push_back(e); }
    //void clearNeighbors() { neighbors.clear(); }
};

}
#endif /* MYSRC_ROUTING_COORDSDISCOVERREPLY_H_ */
