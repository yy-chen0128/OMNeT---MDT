#ifndef MOCK_INET_INTERFACES_H
#define MOCK_INET_INTERFACES_H

#include "inet/common/packet/Packet.h"

namespace inet {

class LifecycleOperation;

class IManetRouting {
public:
    virtual ~IManetRouting() = default;
    virtual void handleStartOperation(LifecycleOperation *operation) = 0;
    virtual void handleStopOperation(LifecycleOperation *operation) = 0;
    virtual void handleCrashOperation(LifecycleOperation *operation) = 0;
};

class UdpSocket {
public:
    class ICallback {
    public:
        virtual ~ICallback() = default;
        virtual void onDataAvailable(UdpSocket *socket) = 0;
        virtual void onError(UdpSocket *socket) = 0;
    };
    
    void setOutputGate(void* gate) {}
    void bind(int port) {}
    void setCallback(ICallback* cb) {}
    void processMessage(void* msg) {}
    Packet* receive() { return nullptr; }
    void sendTo(Packet* packet, const L3Address& dest, int port) {}
};

class INetfilter {
public:
    class IHook {
    public:
        enum Result { ACCEPT, DROP, QUEUE };
        virtual ~IHook() = default;
        virtual Result datagramPreRoutingHook(Packet *datagram) = 0;
        virtual Result datagramPostRoutingHook(Packet *datagram) = 0;
        virtual Result datagramLocalInHook(Packet *datagram) = 0;
        virtual Result datagramLocalOutHook(Packet *datagram) = 0;
        virtual Result datagramForwardHook(Packet *datagram) = 0;
    };
    
    void registerHook(int priority, IHook* hook) {}
};

class IMobility {
public:
    virtual ~IMobility() = default;
};

class IRoutingTable {
public:
    virtual ~IRoutingTable() = default;
    L3Address getRouterId() { return L3Address(Ipv4Address(0)); }
    class NetworkInterface {};
    NetworkInterface* getInterfaceByAddress(const Ipv4Address& addr) { return nullptr; }
    NetworkInterface* getInterface(int index) { return nullptr; }
    void addRoute(void* route) {}
};

class Ipv4Route {
public:
    void setDestination(const Ipv4Address& addr) {}
    void setNetmask(const Ipv4Address& mask) {}
    void setGateway(const Ipv4Address& gw) {}
    void setInterface(void* iface) {}
    void setSourceType(int type) {}
    void setMetric(int metric) {}
};

class IRoute {
public:
    static const int MANET = 1;
};

// Template helper for getModuleFromPar
template<typename T>
T* getModuleFromPar(double par, void* ctx) { return nullptr; }

enum InitStage {
    INITSTAGE_LOCAL = 0,
    INITSTAGE_ROUTING_PROTOCOLS = 1,
    NUM_INIT_STAGES = 2
};

} // namespace inet

#endif // MOCK_INET_INTERFACES_H
