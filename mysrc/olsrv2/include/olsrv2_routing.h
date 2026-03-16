#ifndef MYSRC_OLSRV2_OLSRV2ROUTING_H
#define MYSRC_OLSRV2_OLSRV2ROUTING_H

#include <vector>
#include <map>
#include <set>
#include <memory>

#include "olsrv2_state.h"

// Check if we are in unit test mode or normal OMNeT++ build
#ifdef UNIT_TEST
  #include "mock/omnetpp.h"
  #include "mock/inet/common/packet/Packet.h"
  #include "mock/message/OLSRv2Packet_m.h"
  #include "mock/inet_interfaces.h" 
  #include "nhdp/nhdp.h" // Need Nhdp definition
  
  // Use mock namespace for cSimpleModule etc.
  using namespace omnetpp;
  namespace inet { 
    class LifecycleOperation;
    class RoutingProtocolBase : public omnetpp::cSimpleModule {
    public:
        virtual void handleMessageWhenUp(omnetpp::cMessage *msg) {}
        // Mock handleMessage to forward to handleMessageWhenUp
        virtual void handleMessage(omnetpp::cMessage *msg) override { handleMessageWhenUp(msg); }
        
        // Lifecycle methods
        virtual void handleStartOperation(LifecycleOperation *operation) {}
        virtual void handleStopOperation(LifecycleOperation *operation) {}
        virtual void handleCrashOperation(LifecycleOperation *operation) {}
    }; 
  }
#else
  #include "inet/common/INETDefs.h"
  #include "inet/networklayer/contract/IRoutingTable.h"
  #include "inet/networklayer/contract/INetfilter.h"
  #include "inet/transportlayer/contract/udp/UdpSocket.h"
  #include "inet/mobility/contract/IMobility.h"
  #include "inet/routing/base/RoutingProtocolBase.h"
  #include "nhdp/nhdp.h"
  #include "../message/OLSRv2Packet_m.h"
  
  // In normal build, cSimpleModule is in omnetpp namespace or inet inherits it
#endif

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
class Olsrv2Routing : public inet::RoutingProtocolBase, public inet::UdpSocket::ICallback, public inet::INetfilter::IHook
{
  public:
    Olsrv2Routing();
    virtual ~Olsrv2Routing();

  protected:
    virtual int numInitStages() const override { return inet::NUM_INIT_STAGES; }
    virtual void initialize(int stage) override;
    virtual void handleMessageWhenUp(omnetpp::cMessage *msg) override;

    // IMHnetRouting interface
    virtual void handleStartOperation(inet::LifecycleOperation *operation) override {}
    virtual void handleStopOperation(inet::LifecycleOperation *operation) override {}
    virtual void handleCrashOperation(inet::LifecycleOperation *operation) override {}

    // UdpSocket::ICallback interface
    virtual void onDataAvailable(inet::UdpSocket *socket) override;
    virtual void onError(inet::UdpSocket *socket) override;

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

    // Module parameters
    double helloInterval_ = 2.0;
    double tcInterval_ = 5.0;
    
    // Core logic
    // Use unique_ptr to avoid full definition of Olsrv2Core in header (Pimpl-like)
    // or just pointer. unique_ptr is safer.
    std::unique_ptr<Olsrv2Core> core_;
    Nhdp nhdp_; 
#ifdef UNIT_TEST
#else
#endif

    // INET interactions
    inet::IRoutingTable *routingTable_ = nullptr;
    inet::INetfilter *networkProtocol_ = nullptr;
    inet::IMobility *mobility_ = nullptr;
    inet::UdpSocket socket_;
    
    // Timers
    omnetpp::cMessage *helloTimer_ = nullptr;
    omnetpp::cMessage *tcTimer_ = nullptr;
};

} // namespace mysrc::olsrv2

#endif // MYSRC_OLSRV2_OLSRV2ROUTING_H