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

enum Olsrv2PktType
{
    Olsrv2Hello = 1;
    Olsrv2Tc = 2;
}

// Base type for OLSRv2 control messages
class Olsrv2ControlPacket extends FieldsChunk
{
    Olsrv2PktType packetType = static_cast<Olsrv2PktType>(-1);
}

// OLSRv2 HELLO message (flattened representation)
class Olsrv2HelloPacket extends Olsrv2ControlPacket
{
    L3Address originator;      // Originator/main address
    simtime_t hTime;           // HELLO interval
    simtime_t validityTime;    // RFC 7181 VALIDITY_TIME
    uint8 willingness;         // RFC 7181 willingness
    uint16 msgSeq;             // Message sequence number

    // Neighbor tuples (flattened arrays)
    uint8 linkStatus[];        // Link status per entry
    uint8 neighStatus[];       // Neighbor status per entry
    int addrCounts[];          // Address count per entry
    L3Address neighAddrs[];    // Flattened neighbor interface addresses
}

// OLSRv2 TC group message (flattened representation)
class Olsrv2TcGroup extends Olsrv2ControlPacket
{
    simtime_t genTime;         // Generation time
    simtime_t validityTime;    // RFC 7181 VALIDITY_TIME

    // Per-originator TC metadata
    L3Address originators[];   // Originator addresses
    uint16 ansns[];            // Advertised Neighbor Sequence Number
    uint16 seqNums[];          // Message sequence numbers
    uint8 hopCounts[];         // Hop count per TC
    int advCounts[];           // Advertised-neighbor count per TC

    // Flattened advertised neighbors
    L3Address advertisedNeighbors[];
}

message Olsrv2PacketHolderMessage
{
    Packet *ownedPacket @owned;
}
