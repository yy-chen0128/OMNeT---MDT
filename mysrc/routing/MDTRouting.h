
#ifndef MYSRC_ROUTING_H_
#define MYSRC_ROUTING_H_

#include <map>
#include <unordered_map>
#include <cmath>
#include <deque>

#include "../../mysrc/routing/MDTControlPackets_m.h"
#include "../../mysrc/routing/MDTData.h"
#include "inet/common/ModuleRefByPar.h"
#include "inet/networklayer/contract/IInterfaceTable.h"
#include "inet/networklayer/contract/IL3AddressType.h"
#include "inet/networklayer/contract/INetfilter.h"
#include "inet/networklayer/contract/IRoutingTable.h"

#include "inet/routing/base/RoutingProtocolBase.h"
#include "inet/transportlayer/contract/udp/UdpSocket.h"
#include "inet/transportlayer/udp/UdpHeader_m.h"

#include "inet/mobility/contract/IMobility.h"

#include"NeighborEntry.h"
#include"CoordsDiscoverReply.h"

//rpc
#include"../RPC/log2vis.h"

//wlan
#include "inet/physicallayer/wireless/common/signal/WirelessSignal.h"
#include "inet/physicallayer/wireless/common/analogmodel/packetlevel/ScalarReception.h"
#include "inet/common/Simsignals.h"
#include "inet/physicallayer/wireless/common/contract/packetlevel/IRadio.h"
#include "inet/physicallayer/wireless/common/contract/packetlevel/IReceiver.h"


namespace mysrc {
namespace routing {
using namespace inet;
using namespace inet::physicallayer;

class ForwardProtocol;
class MaintenanceProtocol;
class JoinProtocol;
class Neighbor;


class INET_API MDTRouting : public RoutingProtocolBase, public NetfilterBase::HookBase, public UdpSocket::ICallback, public cListener
{
  protected:

    //MDT
    bool attached = false;
    bool initializing = false;

    ForwardProtocol* forwardProtocol = nullptr;
    MaintenanceProtocol* maintenanceProtocol = nullptr;
    JoinProtocol* joinProtocol = nullptr;
    Neighbor* neighbor = nullptr;

    //Timer
    cMessage *routeTimer = nullptr;

    // context
    const IL3AddressType *addressType = nullptr; // to support both Ipv4 and v6 addresses.

    // environment
    cModule *host = nullptr;
    cModule *mobilityModule = nullptr;
    //wlan
    cModule *wlan = nullptr;
    cModule *radioModule = nullptr;
    std::unordered_map<int, TxPowerEntry> txPowerMap;
    std::vector<TxInfo> recentTxs;
    simtime_t txPowerEntryLifetime = 1.0;

    std::deque<RecentRxEntry> recentRxDeque;
    size_t recentRxCapacity = 3;
    simtime_t recentEntryLifetime = 0.01; // seconds
    double matchTimeEps = 1e-7; // seconds, adjust if needed

    double kaThreshold_dBm_default = -84; // fallback: threshold for ~225m (approx)
    double kaThreshold_dBm_strict = NAN;      // the "narrow" threshold used when neighbors > 3
    double radioTxPower_dBm = NAN;
    double radioFreq_Hz = 2.412e9; // default 2.412 GHz if can't read
    double pathlossAlpha = 2.0;    // default (use your .ini value)


    ModuleRefByPar<IRoutingTable> routingTable;
    ModuleRefByPar<IInterfaceTable> interfaceTable;
    ModuleRefByPar<INetfilter> networkProtocol;

    UdpSocket socket;

    bool usingIpv6 = false;
    //above are useful

    //parameters: the following parameters are configurable, see the NED file for more info.

    unsigned int mdtUDPPort = 0;


    //useful

    unsigned int netDiameter = 0;

    simtime_t deletePeriod;
    simtime_t myRouteTimeout;

    // TODO:self messages:Timers

    // lifecycle
    simtime_t rebootTime; // the last time when the node rebooted

    // internal
    std::multimap<L3Address, Packet *> targetAddressToDelayedPackets; // queue for the datagrams we have no route for

  protected:
    void sendmessage();
    void createmessage();
    void handlemessage();

    void sendMDTPacket(const Ptr<MDTControlPacket>& packet);



    /* Netfilter hooks */
    Result ensureRouteForDatagram(Packet *datagram,int state);//when a node is going to send a packet or receiving a packet(the packet just arrives network layer)
    //state 0->datagramLocalOutHook 1->datagramPreRoutingHook
    virtual Result datagramPreRoutingHook(Packet *datagram) override { Enter_Method("datagramPreRoutingHook"); return ensureRouteForDatagram(datagram,1); }
    virtual Result datagramForwardHook(Packet *datagram) override{
        Enter_Method("datagramForwardHook");
        std::vector<int> msginfo = getMsgEventInfo(datagram,"FORWARD");
        long long msgId = datagram->getId();
        std::string pkname = datagram->getName();
        if(pkname.empty())pkname = "unknown";
        nlohmann::json payload = {
            {{"bytes_size", datagram->getTotalLength().get() / 8.0}, {"packetType", pkname}} // or packet->getByteLength()
        };
        logMsgEvent(this, "forward", msgId, /*srcNodeId*/msginfo[0], /*dst*/ msginfo[3], msginfo[1], msginfo[2], payload);
        return ACCEPT;
    }
    virtual Result datagramPostRoutingHook(Packet *datagram) override { return ACCEPT; }//about to leave
    virtual Result datagramLocalInHook(Packet *datagram) override {
        Enter_Method("datagramLocalInHook");
        std::vector<int> msginfo = getMsgEventInfo(datagram,"IN");
        long long msgId = datagram->getId();
        std::string pkname = datagram->getName();
        if(pkname.empty())pkname = "unknown";
        nlohmann::json payload = {
            {{"bytes_size", datagram->getTotalLength().get() / 8.0}, {"packetType", pkname}} // or packet->getByteLength()
        };
        logMsgEvent(this, "in", msgId, /*srcNodeId*/msginfo[0], /*dst*/ msginfo[3], msginfo[1], msginfo[2], payload);
        return ACCEPT;
    }//local process the packet(other layer in charge of it)
    virtual Result datagramLocalOutHook(Packet *datagram) override {
        Enter_Method("datagramLocalOutHook");
        EV_ERROR<<"packet send from src: "<<getSelfAddress()<<"\n";
        Result res = ensureRouteForDatagram(datagram,0);

        std::vector<int> msginfo = getMsgEventInfo(datagram,"OUT");
        long long msgId = datagram->getId();
        std::string pkname = datagram->getName();
        if(pkname.empty())pkname = "unknown";
        nlohmann::json payload = {
            {{"bytes_size", datagram->getTotalLength().get() / 8.0}, {"packetType", pkname}} // or packet->getByteLength()
        };
        logMsgEvent(this, "out", msgId, /*srcNodeId*/msginfo[0], /*dst*/ msginfo[3], msginfo[1], msginfo[2], payload);

        return res;
    }
    void delayDatagram(Packet *datagram);
    //rpc
    std::vector<int> getMsgEventInfo(Packet *datagram,std::string type);
    void routetablelog();

    void handleMessageWhenUp(cMessage *msg) override;
    void initialize(int stage) override;
    virtual void finish() override;

    virtual int numInitStages() const override { return NUM_INIT_STAGES; }

    /* Route Discovery */
    void startRouteDiscovery(const L3Address& target, unsigned int timeToLive = 0);

    void completeRouteDiscovery(const L3Address& target);

    bool hasOngoingRouteDiscovery(const L3Address& target);
    void cancelRouteDiscovery(const L3Address& destAddr);

    /* General functions to handle route errors */

    virtual void receiveSignal(cComponent *source, simsignal_t signalID, cObject *obj, cObject *details) override;

    /* Helper functions */

    void processPacket(Packet *pk);
    void clearState();
    void checkIpVersionAndPacketTypeCompatibility(MDTControlPacketType packetType);

    /* UDP callback interface */
    virtual void socketDataArrived(UdpSocket *socket, Packet *packet) override;
    virtual void socketErrorArrived(UdpSocket *socket, Indication *indication) override;
    virtual void socketClosed(UdpSocket *socket) override;

    /* Lifecycle */
    virtual void handleStartOperation(LifecycleOperation *operation) override;
    virtual void handleStopOperation(LifecycleOperation *operation) override;
    virtual void handleCrashOperation(LifecycleOperation *operation) override;

  public:
    MDTRouting();
    virtual ~MDTRouting();

    inet::IMobility *mobility = nullptr;
    simtime_t activeRouteTimeout;
    //rng
    cRNG *rng = getRNG(0);

    bool getNodeattachstate(){return attached;}
    bool getNodeinitializingstate(){return initializing;}
    void setNodeattachstate(bool state){attached = state;}
    void setNodeinitializingstate(bool state){initializing = state;}

    void sendPacketTo(Packet *packet, const L3Address& destAddr, int destPort = 0,double delay = 0);
    L3Address getSelfIPAddress() const;
    void printPacketTags(Packet *packet, const char *context);

    // get self address
    L3Address getSelfAddress() const;

    // delegate to neighbor module: get list of physical neighbors
    //std::vector<L3Address> getPhysicalNeighbors() const;

    // find a tuple by destination (delegates to Neighbor)
    bool findTupleByDestandUpdate(const L3Address& dest, Routetuple &outTuple) const;


    // get coords for a node (delegates to Neighbor)
    bool getNeighborNodeCoordinate(const L3Address& node, std::vector<double>& outCoords) const;
    bool getDTNeighborNodeCoordinate(const L3Address& node, std::vector<double>& outCoords) const;
    bool getKnownNodeCoordinate(const L3Address& node, std::vector<double>& outCoords) const;

    bool updateValidRouteLifeTime(const L3Address& destAddr, simtime_t lifetime);

    bool getSelfCoordinate(std::vector<double> &coordVec) const;

    //std::vector<L3Address> getDTNeighbors() const;

    void updateRoutingTable(IRoute *route,const L3Address& nextHop,const Routetuple& routetuple);
    IRoute *createRoute(const L3Address& destAddr,const L3Address& nextHop,const Routetuple& routetuple);
    IRoute *ensureRoute(const L3Address& dest,const L3Address& nextHop,const Routetuple& rt);

    void sendJoinRequest(const L3Address& destAddr);
    std::map<L3Address, NeighborEntry> getKnownNodes() const;
    std::map<L3Address, NeighborEntry> getDTNeighbors() const;
    std::map<L3Address, NeighborEntry> getMinSimplexNodes() const;
    std::map<L3Address, NeighborEntry> getPNeighborNodes() const;
    std::map<L3Address, NeighborEntry> getweaklinkNeighbors() const;
    void sendNSRQ(const std::vector<NeighborEntry>& entrys,const L3Address& srcAddr);

    //neighbor
    bool isDTNeighborOfTarget(const std::vector<double>& targetCoord, const L3Address& addr);
    std::vector<NeighborEntry> computeDTNeighborsForCoord(const std::vector<double>& targetCoord,const L3Address& addr);


    bool addOrUpdateKnownNode(const std::vector<NeighborEntry>& entrys);
    bool addOrUpdateDTNeighbors(const std::vector<NeighborEntry>& entrys);
    void updateDTNeighbors();

    //Maintenance
    void changeMaintenanceStage(int x);

    //route table
    IRoute *findBestMatchingRoute(const L3Address& dest);
    void updateRoutetuple(const L3Address& sourceAddr,const L3Address& predAddr,const L3Address& succAddr,const L3Address& destAddr,
                            double cost = -1,double error = -1);
    void refreshRoutetuple();
    void findroutefornewdt(const std::vector<L3Address> newdtset);
    void deleteexpiredneighbors(std::vector<L3Address> nodes);
    //dest discover
    //L3Address findNearestPneighborForTarget(const std::vector<double>& targetCoord);
    void sendCoordsDiscoverRequest(const L3Address& targetAddr);
    void processCoordsDiscoverRequestChunk(const Ptr<CoordsDiscoverRequest>& cd);
    void sendCoordsDiscoverReply(const L3Address& destAddr,const L3Address& targetAddr,const NeighborEntry& entry);
    void processCoordsDiscoverReplyChunk(const Ptr<CoordsDiscoverReply>& cd);
};

} // namespace aodv
} // namespace inet

#endif

