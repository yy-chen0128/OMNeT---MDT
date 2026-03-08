//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/.
//

import inet.common.INETDefs;
import inet.common.packet.chunk.Chunk;
import inet.common.packet.Packet;
import inet.networklayer.common.L3Address;


namespace inet;

enum OlsrPktType
{
    Hello = 1;
    Tc = 2;
}
// Base type for OLSR control messages (individual message types inherit from it)
class OlsrControlPacket extends FieldsChunk
{
    // no common fields here: specific message types will extend
    OlsrPktType packetType = static_cast<OlsrPktType>(-1); 
}

// HELLO message
class OlsrHello extends OlsrControlPacket
{
    // packet-level: packetType field inherited from OLSRControlPacket

    // HELLO header fields
    L3Address originator;   // originator/main address
    simtime_t hTime ;  // hello interval (seconds)
    uint8 willingness;      // 0..7
    uint16 msgSeq;          // per-HELLO message sequence number (for debug/dup detection)

    // Per-link messages flattened representation:
    // linkCodes[i] : 1 byte link code for i-th link message
    // linkCounts[i]: number of neighbor iface addresses for i-th link message
    uint8 linkCodes[];         // each is 1 byte in RFC
    int linkCounts[];          // number of addresses in that link message
    L3Address nbIfaceAddrs[];  // flattened neighbour-interface addresses (sum(linkCounts) entries)
}

// OlsrTcGroup: a group-of-TC logical entries packed into one UDP packet.
// Inherit OLSRControlPacket 
class OlsrTcGroup extends OlsrControlPacket
{
    // Group-level metadata
    simtime_t genTime;      // generation time (simtime_t)
    //double vTime @unit(s);  // validity time (VTime) in seconds - lifetime for advertised info
    //uint16 groupSeq;        // group sequence number (optional, for debug)

    // Per-TC parallel arrays (length = M logical TC entries)
    // Each logical TC entry represents one originator's TC: originator addr + ANSN + hopCount + count of advertised neighbors
    L3Address originators[]; // originator addresses (one per logical TC)
    uint16 ansns[];          // ANSN per originator
    uint16 seqNums[];        // message sequence numbers (optional; helpful for de-dup)
    uint8 hopCounts[];       // hopCount for each logical TC (how many hops from originator to here when created)
    int advCounts[];         // number of advertisedNeighbors for each logical TC

    // Flattened list of all advertised neighbor addresses: sum(advCounts) entries
    L3Address advertisedNeighbors[]; // flattened neighbor lists
}

message OlsrPacketHolderMessage//for delay send
{
    Packet *ownedPacket @owned;
}







