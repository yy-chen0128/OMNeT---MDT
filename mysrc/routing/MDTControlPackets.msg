
import inet.common.INETDefs;
import inet.common.packet.chunk.Chunk;
import inet.common.packet.Packet;
import inet.networklayer.common.L3Address;


namespace inet;
//namespace src::routing;

enum MDTControlPacketType
{
    JRQ = 1;
    JRP= 2;
    NSRQ = 3;
	NSRP = 4;
	NSN = 5;
	IT = 6;
	LN = 7;
	FN = 8;
	KA = 9;
	CDQ = 10;
	CDP = 11;
	KAN = 12;

	//JRQ_IPv6 = 11;
    //JRP_IPv6 = 12;
    //NSRQ_IPv6 = 13;
	//NSRP_IPv6 = 14;
	//NSN_IPv6 = 15;
	//IT_IPv6 = 16;
	//LN_IPv6 = 17;
	//FN_IPv6 = 18;
	//KA_IPv6 = 19;
	//CDQ_IPv6 = 20;
	//CDP_IPv6 = 21;
	//KAN_IPv6 = 22;
};

class MDTControlPacket extends FieldsChunk
{
    MDTControlPacketType packetType = static_cast<MDTControlPacketType>(-1);       // packet type, stored in one byte
}

class KeepAlive extends MDTControlPacket {
	//L3Address sourceAddr;
	//simtime_t timestamp;
	@customize(true);
	double coords[];
	//int dim;
	bool attached;
}
class KeepAliveNotification extends MDTControlPacket {
	double coords[];
	//int dim;
	bool attached;
}
class JoinRequest extends MDTControlPacket{
    double coords[];
    //Coordinate nodeCoord;
    int dim;
    L3Address destAddr;
    L3Address sourceAddr;
    L3Address predAddr;
}
class JoinReply extends MDTControlPacket{
    @customize(true);
    L3Address sourceAddr;
    L3Address predAddr;
    L3Address destAddr;
}
class NeighborSetRequest extends MDTControlPacket{
    double coords[];
    //Coordinate nodeCoord;
    int dim;
    L3Address sourceAddr;
    L3Address destAddr;
    L3Address predAddr;
}
class NeighborSetReply extends MDTControlPacket{
    @customize(true);
    L3Address sourceAddr;
    L3Address predAddr;
    L3Address destAddr;
}
class NeighborSetNotification extends MDTControlPacket{
    double coords[];
    int dim;
    L3Address sourceAddr;
    L3Address predAddr;
    L3Address destAddr;
}
class CoordsDiscoverRequest extends MDTControlPacket{
    //L3Address destAddr;
    //int num;
    L3Address sourceAddr;
    L3Address targetAddr;
    int hopcount;
}
class CoordsDiscoverReply extends MDTControlPacket{
    @customize(true);
    L3Address destAddr;
    L3Address sourceAddr;
    L3Address targetAddr;
}
//class InitializationToken extends MDTControlPacket{
    
///}
//class LeaveNotification extends MDTControlPacket{
    
//}
//class FailureNotification extends MDTControlPacket{
    
//}

//
// Represents a timer for delayed sending
//
message PacketHolderForNSMessage//not use
{
    Packet *ownedPacket @owned;
    L3Address nsrqSrcAddr;
    L3Address nsrqDestAddr;
}
message PacketHolderMessage
{
    Packet *ownedPacket @owned;
    L3Address destAddr;
}
