/*
 * OLSR.h
 *
 *  Created on: Feb 9, 2026
 *      Author: yychen
 */

#ifndef MYSRC_OLSR_OLSR_H_
#define MYSRC_OLSR_OLSR_H_

#include "inet/common/ModuleRefByPar.h"
#include "inet/networklayer/contract/IInterfaceTable.h"
#include "inet/networklayer/contract/IL3AddressType.h"
#include "inet/networklayer/contract/INetfilter.h"
#include "inet/networklayer/contract/IRoutingTable.h"

#include "inet/routing/base/RoutingProtocolBase.h"
#include "inet/transportlayer/contract/udp/UdpSocket.h"
#include "inet/transportlayer/udp/UdpHeader_m.h"

#include <algorithm>
#include "../olsr/OlsrControlPackets_m.h"
#include "../olsr/OLSR_Repositories.h"
#include "../olsr/OLSR_RouteData.h"
#include "../olsr/OLSR_State.h"

namespace mysrc {
namespace olsr {
using namespace inet;

// TC cache
struct CachedTcEntry
{
    L3Address originator;
    uint16_t ansn;
    uint16_t msgSeq;
    uint8_t hopCount;
    uint8_t hopLimit;
    std::vector<L3Address> advertisedNeighbors;
    SimTimeType expiry;
};
//route compute
struct TempRoute
{
    L3Address dest;
    L3Address nextHop;
    int distance;
};
class INET_API OLSR : public RoutingProtocolBase,
                     public NetfilterBase::HookBase,
                     public UdpSocket::ICallback,
                     public cListener
{
  protected:

    cModule *host = nullptr;
    // timers
    cMessage *helloMsgTimer = nullptr;
    //cMessage *tcMsgTimer = nullptr;
    // socket
    UdpSocket socket;
    L3Address multicastGroup;
    int ifId = -1;
    int olsrUDPPort = 698; //

    // pointers to INET tables
    ModuleRefByPar<IRoutingTable> routingTable;
    ModuleRefByPar<IInterfaceTable> interfaceTable;
    ModuleRefByPar<INetfilter> networkProtocol;
    //state
    OlsrState state;
    //L3Address selfAddress;

    double helloInterval = 0.75;           // seconds
    double maxPeriodicJitter = 0.25;       // seconds
    uint8_t willingness = 3;               // default
    uint16_t helloMsgSeq = 0;

    // tc
    uint16_t tcMsgSeq = 0;
    std::vector<CachedTcEntry> tcCache;   // buffer of received TC entries waiting to be forwarded/aggregated
    double tcCacheInterval = 0.05;      // configurable param (seconds) - short aggregation interval
    double tcInterval = 1.5;   // the longer periodic interval for own TC (example)
    cMessage *tcGroupTimer = nullptr;     // short timer for buffer flush/forward
    cMessage *tcCacheTimer = nullptr;  // long periodic timer for own TC
    int tcMaxPayloadBytes = 1200;         // conservative payload limit for aggregation (avoid fragmentation)
    uint16_t ansn = 0;                    // own advertised neighbor sequence number (increment when ANSN changes)

    const IL3AddressType *addressType = nullptr; // to support both Ipv4 and v6 addresses.now ipv4 only

    //time
    simtime_t DupHoldTime = 10;
    simtime_t neighborHoldTime = 2 * helloInterval;
    simtime_t twoHopHoldTime = 2 * helloInterval;
    simtime_t topologyHoldTime = 2 * tcInterval;
  public:
    OLSR();
    virtual ~OLSR();

    virtual void initialize(int stage) override;
    virtual void handleMessageWhenUp(cMessage *msg) override;
    L3Address getSelfIPAddress() const;
    // Netfilter hooks
    virtual Result datagramPreRoutingHook(Packet *datagram) override { Enter_Method("datagramPreRoutingHook"); return ensureRouteForDatagram(datagram);}
    virtual Result datagramForwardHook(Packet *datagram) override { return ACCEPT; }
    virtual Result datagramPostRoutingHook(Packet *datagram) override { return ACCEPT; }
    virtual Result datagramLocalInHook(Packet *datagram) override { return ACCEPT; }
    virtual Result datagramLocalOutHook(Packet *datagram) override { Enter_Method("datagramLocalOutHook"); return ensureRouteForDatagram(datagram); }
    // hook implement
    Result ensureRouteForDatagram(Packet *datagram);

    void processPacket(Packet *pk);
    // UdpSocket::ICallback
    virtual void socketDataArrived(UdpSocket *socket, Packet *packet) override;
    virtual void socketErrorArrived(UdpSocket *socket, Indication *indication) override;
    virtual void socketClosed(UdpSocket *socket) override;
    /* Lifecycle */
    virtual void handleStartOperation(LifecycleOperation *operation) override;
    virtual void handleStopOperation(LifecycleOperation *operation) override;
    virtual void handleCrashOperation(LifecycleOperation *operation) override;

    // OLSR specifics
    void sendPacket(Packet *packet, int destPort, double delay);
    void sendHelloPacket();
    void sendTcPacket();
    void clearOlsrRoutes();
    IRoute *createRoute(const L3Address& dest, const L3Address& nextHop, unsigned int hopCount);

    //process message
    void processHello(const Ptr<OlsrHello>& hello);
    void LinkSensing(const L3Address& sender, const Ptr<OlsrHello>& hello);
    void PopulateNeighborSet(const L3Address& sender);
    void PopulateTwoHopNeighborSet(const L3Address& sender, const Ptr<OlsrHello>& hello);
    void MprComputation();
    void PopulateMprSelectorSet(const L3Address& sender, const Ptr<OlsrHello>& hello);

    void processTc(const Ptr<OlsrTcGroup>& tc, const L3Address& sourceAddr);
    std::vector<L3Address> getAdvertisedNeighborList();
    void constructLocalTc();

    void RoutingTableComputation();

    //print
    void printRouteTable();

};

}
}
#endif /* MYSRC_OLSR_OLSR_H_ */
