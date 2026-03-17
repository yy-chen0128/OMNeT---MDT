#ifndef MYSRC_OLSRV2_OLSRV2ROUTING_H
#define MYSRC_OLSRV2_OLSRV2ROUTING_H

#include <vector>
#include <map>
#include <set>
#include <memory>
#include <string>

#include "olsrv2_state.h"

#include "inet/common/INETDefs.h"
#include "inet/networklayer/ipv4/IIpv4RoutingTable.h"
#include "inet/networklayer/contract/IInterfaceTable.h"
#include "inet/networklayer/contract/INetfilter.h"
#include "inet/transportlayer/contract/udp/UdpSocket.h"
#include "inet/mobility/contract/IMobility.h"
#include "inet/routing/base/RoutingProtocolBase.h"
#include "inet/networklayer/common/L3Address.h"
#include "nhdp/nhdp.h"
#include "../message/OLSRv2Packet_m.h"

namespace mysrc::olsrv2 {

// Forward declare Core to break circular dependency
class Olsrv2Core;

// Stateless route calculator from current OLSRv2 sets.
class Olsrv2RoutingEngine
{
  public:
    // Computes shortest-path-like routes from originator to reachable nodes.
    // Metric policy: currently hop count.
    std::vector<RouteTuple> computeRoutes(const MainAddress& originator, const Olsrv2State& state) const;
};

// Main OMNeT++ Module for OLSRv2 Routing
class OLSRv2 : public inet::RoutingProtocolBase, public inet::UdpSocket::ICallback, public inet::INetfilter::IHook
{
  public:
    OLSRv2();
    virtual ~OLSRv2();

  protected:
    virtual int numInitStages() const override { return inet::NUM_INIT_STAGES; }
    virtual void initialize(int stage) override;
    virtual void handleMessageWhenUp(omnetpp::cMessage *msg) override;

    // IMHnetRouting interface (not inherited anymore)
    virtual void handleStartOperation(inet::LifecycleOperation *operation) override {}
    virtual void handleStopOperation(inet::LifecycleOperation *operation) override {}
    virtual void handleCrashOperation(inet::LifecycleOperation *operation) override {}

    // UdpSocket::ICallback interface
    virtual void socketDataArrived(inet::UdpSocket *socket, inet::Packet *packet) override;
    virtual void socketErrorArrived(inet::UdpSocket *socket, inet::Indication *indication) override;
    virtual void socketClosed(inet::UdpSocket *socket) override;

    // INetfilter::IHook interface
    virtual inet::INetfilter::IHook::Result datagramPreRoutingHook(inet::Packet *datagram) override;
    virtual inet::INetfilter::IHook::Result datagramPostRoutingHook(inet::Packet *datagram) override { return inet::INetfilter::IHook::ACCEPT; }
    virtual inet::INetfilter::IHook::Result datagramLocalInHook(inet::Packet *datagram) override { return inet::INetfilter::IHook::ACCEPT; }
    virtual inet::INetfilter::IHook::Result datagramLocalOutHook(inet::Packet *datagram) override { return inet::INetfilter::IHook::ACCEPT; }
    virtual inet::INetfilter::IHook::Result datagramForwardHook(inet::Packet *datagram) override { return inet::INetfilter::IHook::ACCEPT; }

  private:
    void processHelloTimer();
    void processTcTimer();
    void sendHello();
    void sendTc();
    void processOlsrPacket(inet::Packet *packet);
    void updateRoutingTable();
    void printProtocolState(const char *tag) const;
    void printNhdpState() const;
    void printIpv4RoutingTable() const;

    // Module parameters
    double helloInterval_ = 2.0;
    double tcInterval_ = 5.0;
    double startJitter_ = 0.0;
    int udpPort_ = 698;
    inet::L3Address multicastGroup_;
    int multicastIfId_ = -1;
    
    // Core logic
    // Use unique_ptr to avoid full definition of Olsrv2Core in header (Pimpl-like)
    // or just pointer. unique_ptr is safer.
    std::unique_ptr<Olsrv2Core> core_;
    Nhdp nhdp_; 

    // INET interactions
    inet::IIpv4RoutingTable *routingTable_ = nullptr;
    inet::IInterfaceTable *interfaceTable_ = nullptr;
    inet::INetfilter *netfilter_ = nullptr;
    inet::IMobility *mobility_ = nullptr;
    inet::UdpSocket socket_;
    
    // Timers
    omnetpp::cMessage *helloTimer_ = nullptr;
    omnetpp::cMessage *tcTimer_ = nullptr;
};

} // namespace mysrc::olsrv2

#endif // MYSRC_OLSRV2_OLSRV2ROUTING_H
