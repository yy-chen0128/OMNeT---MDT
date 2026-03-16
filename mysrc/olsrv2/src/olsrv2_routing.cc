//
// Created by GuoYudi on 2026/3/14.
//
#include "../include/olsrv2_routing.h"
#include "../include/olsrv2.h"

#include <map>
#include <queue>
#include <set>

#include "inet/common/ModuleAccess.h"
#include "inet/common/packet/Packet.h"
#include "inet/networklayer/common/L3AddressTag_m.h"
#include "inet/networklayer/common/L3AddressResolver.h"
#include "inet/networklayer/ipv4/Ipv4RoutingTable.h"
#include "inet/networklayer/ipv4/Ipv4Route.h"

using namespace omnetpp;

namespace mysrc::olsrv2 {

std::vector<RouteTuple> Olsrv2RoutingEngine::computeRoutes(const MainAddress& originator, const Olsrv2State& state) const
{
    std::map<MainAddress, std::set<MainAddress>> graph;

    // Build directed graph from topology tuples: last_addr -> dest_addr
    for (const auto& t : state.getTopologyTuples())
        graph[t.last_addr].insert(t.dest_addr);

    // Also add 1-hop symmetric neighbors as originator edges.
    for (const auto& n : state.getNeighbors()) {
        if (n.is_symmetric)
            graph[originator].insert(n.main_addr);
    }

    std::map<MainAddress, uint32_t> dist;
    std::map<MainAddress, MainAddress> prev;
    std::queue<MainAddress> q;

    dist[originator] = 0;
    q.push(originator);

    while (!q.empty()) {
        MainAddress u = q.front();
        q.pop();

        auto it = graph.find(u);
        if (it == graph.end())
            continue;

        for (const auto& v : it->second) {
            if (dist.find(v) != dist.end())
                continue;
            dist[v] = dist[u] + 1;
            prev[v] = u;
            q.push(v);
        }
    }

    std::vector<RouteTuple> routes;
    for (const auto& kv : dist) {
        const MainAddress& dst = kv.first;
        if (dst == originator)
            continue;

        MainAddress hop = dst;
        while (prev.find(hop) != prev.end() && prev[hop] != originator)
            hop = prev[hop];

        RouteTuple rt;
        rt.destination = dst;
        rt.next_hop = hop;
        rt.metric = kv.second;
        rt.hop_count = static_cast<uint8_t>(kv.second);
        routes.push_back(rt);
    }

    return routes;
}

// Module Implementation
Define_Module(Olsrv2Routing);

Olsrv2Routing::Olsrv2Routing()
{
}

Olsrv2Routing::~Olsrv2Routing()
{
    cancelAndDelete(helloTimer_);
    cancelAndDelete(tcTimer_);
}

void Olsrv2Routing::initialize(int stage)
{
    RoutingProtocolBase::initialize(stage);

    if (stage == inet::INITSTAGE_LOCAL) {
        helloInterval_ = par("helloInterval");
        tcInterval_ = par("tcInterval");
        
        helloTimer_ = new inet::cMessage("helloTimer");
        tcTimer_ = new inet::cMessage("tcTimer");
    }
    else if (stage == inet::INITSTAGE_ROUTING_PROTOCOLS) {
        routingTable_ = inet::getModuleFromPar<inet::IIpv4RoutingTable>(par("routingTableModule"), this);
        interfaceTable_ = inet::getModuleFromPar<inet::IInterfaceTable>(par("interfaceTableModule"), this);
        netfilter_ = inet::getModuleFromPar<inet::INetfilter>(par("networkProtocolModule"), this);
        mobility_ = inet::getModuleFromPar<inet::IMobility>(par("mobilityModule"), this);

        netfilter_->registerHook(0, this);

        socket_.setOutputGate(gate("socketOut"));
        socket_.bind(698);
        socket_.setCallback(this);

        // Initialize Core Logic
        inet::Ipv4Address myIp = routingTable_->getRouterId();
        core_ = std::make_unique<Olsrv2Core>(myIp);
        nhdp_.setOriginator(inet::L3Address(myIp));

        scheduleAfter(helloInterval_, helloTimer_);
        scheduleAfter(tcInterval_, tcTimer_);
    }
}

void Olsrv2Routing::handleMessageWhenUp(omnetpp::cMessage *msg)
{
    if (msg->isSelfMessage()) {
        if (msg == helloTimer_) {
            processHelloTimer();
        } else if (msg == tcTimer_) {
            processTcTimer();
        }
    } else {
        socket_.processMessage(msg);
    }
}

void Olsrv2Routing::processHelloTimer()
{
    sendHello();
    scheduleAfter(helloInterval_, helloTimer_);
}

void Olsrv2Routing::processTcTimer()
{
    sendTc();
    scheduleAfter(tcInterval_, tcTimer_);
}

void Olsrv2Routing::sendHello()
{
    double now = simTime().dbl();
    auto pkt = nhdp_.generateHello(now);
    socket_.sendTo(pkt, inet::L3Address(inet::Ipv4Address::ALLONES_ADDRESS), 698);
}

void Olsrv2Routing::sendTc()
{
    double now = simTime().dbl();
    if (core_) {
        auto pkt = core_->generateTc(nhdp_.getDb(), now);
        socket_.sendTo(pkt, inet::L3Address(inet::Ipv4Address::ALLONES_ADDRESS), 698);
    }
}

void Olsrv2Routing::socketDataArrived(inet::UdpSocket *socket, inet::Packet *packet)
{
    processOlsrPacket(packet);
    delete packet;
}

void Olsrv2Routing::socketErrorArrived(inet::UdpSocket *socket, inet::Indication *indication)
{
    // Handle socket error
    delete indication;
}

void Olsrv2Routing::socketClosed(inet::UdpSocket *socket)
{
    // Handle closed
}

void Olsrv2Routing::processOlsrPacket(inet::Packet *packet)
{
    auto chunk = packet->peekAtFront<inet::Olsrv2ControlPacket>();
    if (!chunk) return;

    inet::L3Address source;
    auto l3Tag = packet->getTag<inet::L3AddressInd>();
    if (l3Tag) source = l3Tag->getSrcAddress();

    if (auto hello = inet::dynamicPtrCast<const inet::Olsrv2HelloPacket>(chunk)) {
        nhdp_.processHello(hello.get(), source, simTime().dbl());
    } else if (auto tc = inet::dynamicPtrCast<const inet::Olsrv2TcGroup>(chunk)) {
        if (core_) {
            core_->processTcAndRecompute(tc.get(), simTime().dbl());
            updateRoutingTable();
        }
    }
}

inet::INetfilter::IHook::Result Olsrv2Routing::datagramPreRoutingHook(inet::Packet *datagram)
{
    return inet::INetfilter::IHook::ACCEPT;
}

void Olsrv2Routing::updateRoutingTable()
{
    if (!core_) return;

    const auto& routes = core_->state().getRoutes();
    inet::IIpv4RoutingTable *rt = routingTable_;
    
    for (const auto& r : routes) {
        inet::Ipv4Route *route = new inet::Ipv4Route();
        route->setDestination(inet::Ipv4Address(r.destination.getInt())); 
        route->setNetmask(inet::Ipv4Address::ALLONES_ADDRESS);
        route->setGateway(inet::Ipv4Address(r.next_hop.getInt()));
        
        inet::NetworkInterface *ie = rt->getInterfaceByAddress(inet::Ipv4Address(r.next_hop.getInt()));
        if (!ie) ie = interfaceTable_->getInterface(0); // Fallback

        route->setInterface(ie);
        route->setSourceType(inet::IRoute::MANET);
        route->setMetric(r.metric);
        rt->addRoute(route);
    }
}

} // namespace mysrc::olsrv2