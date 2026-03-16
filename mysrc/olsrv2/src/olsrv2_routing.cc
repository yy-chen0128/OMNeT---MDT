//
// Created by GuoYudi on 2026/3/14.
//
#include "../include/olsrv2_routing.h"
#include "../include/olsrv2.h"

#include <map>
#include <queue>
#include <set>

#ifndef UNIT_TEST
#include "inet/common/ModuleAccess.h"
#include "inet/common/packet/Packet.h"
#include "inet/networklayer/common/L3AddressTag_m.h"
#include "inet/networklayer/common/L3AddressResolver.h"
#include "inet/networklayer/ipv4/Ipv4RoutingTable.h"
#include "inet/networklayer/ipv4/Ipv4Route.h"
#endif

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
        
#ifdef UNIT_TEST
    helloTimer_ = new omnetpp::cMessage("helloTimer");
    tcTimer_ = new omnetpp::cMessage("tcTimer");
#else
    helloTimer_ = new inet::cMessage("helloTimer");
    tcTimer_ = new inet::cMessage("tcTimer");
#endif
    }
    else if (stage == inet::INITSTAGE_ROUTING_PROTOCOLS) {
#ifndef UNIT_TEST
        routingTable_ = inet::getModuleFromPar<inet::IRoutingTable>(par("routingTableModule"), this);
        networkProtocol_ = inet::getModuleFromPar<inet::INetfilter>(par("networkProtocolModule"), this);
        mobility_ = inet::getModuleFromPar<inet::IMobility>(par("mobilityModule"), this);

        networkProtocol_->registerHook(0, this);

        socket_.setOutputGate(gate("socketOut"));
        socket_.bind(698);
        socket_.setCallback(this);

        // Initialize Core Logic
        inet::L3Address myIp = routingTable_->getRouterId();
        core_ = std::make_unique<Olsrv2Core>(myIp.toIpv4());
        nhdp_.setOriginator(myIp);
#else
        core_ = std::make_unique<Olsrv2Core>(inet::Ipv4Address(0));
#endif

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
#ifndef UNIT_TEST
        socket_.processMessage(msg);
#endif
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
#ifdef UNIT_TEST
    double now = 0.0;
#else
    double now = simTime().dbl();
#endif
    auto pkt = nhdp_.generateHello(now);
#ifndef UNIT_TEST
    socket_.sendTo(pkt, inet::L3Address(inet::Ipv4Address::ALLONES_ADDRESS), 698);
#endif
}

void Olsrv2Routing::sendTc()
{
#ifdef UNIT_TEST
    double now = 0.0;
#else
    double now = simTime().dbl();
#endif
    if (core_) {
        auto pkt = core_->generateTc(nhdp_.getDb(), now);
#ifndef UNIT_TEST
        socket_.sendTo(pkt, inet::L3Address(inet::Ipv4Address::ALLONES_ADDRESS), 698);
#endif
    }
}

void Olsrv2Routing::onDataAvailable(inet::UdpSocket *socket)
{
#ifndef UNIT_TEST
    inet::Packet *packet = socket->receive();
    processOlsrPacket(packet);
    delete packet;
#endif
}

void Olsrv2Routing::onError(inet::UdpSocket *socket)
{
    // Handle socket error
}

void Olsrv2Routing::processOlsrPacket(inet::Packet *packet)
{
#ifndef UNIT_TEST
    auto chunk = packet->peekAtFront<inet::Olsrv2ControlPacket>();
    if (!chunk) return;

    inet::L3Address source;
    auto l3Tag = packet->getTag<inet::L3AddressInd>();
    if (l3Tag) source = l3Tag->getSrcAddress();

    if (auto hello = std::dynamic_pointer_cast<const inet::Olsrv2HelloPacket>(chunk)) {
        nhdp_.processHello(hello.get(), source, simTime().dbl());
    } else if (auto tc = std::dynamic_pointer_cast<const inet::Olsrv2TcGroup>(chunk)) {
        if (core_) {
            core_->processTcAndRecompute(tc.get(), simTime().dbl());
            updateRoutingTable();
        }
    }
#endif
}

#ifdef UNIT_TEST
inet::INetfilter::IHook::Result Olsrv2Routing::datagramPreRoutingHook(inet::Packet *datagram)
#else
inet::IHook::Result Olsrv2Routing::datagramPreRoutingHook(inet::Packet *datagram)
#endif
{
#ifdef UNIT_TEST
    return inet::INetfilter::IHook::ACCEPT;
#else
    return inet::IHook::ACCEPT;
#endif
}

void Olsrv2Routing::updateRoutingTable()
{
#ifndef UNIT_TEST
    if (!core_) return;

    const auto& routes = core_->state().getRoutes();
    inet::IRoutingTable *rt = routingTable_;
    
    for (const auto& r : routes) {
        inet::Ipv4Route *route = new inet::Ipv4Route();
        route->setDestination(inet::Ipv4Address(r.destination.getInt())); 
        route->setNetmask(inet::Ipv4Address::ALLONES_MASK);
        route->setGateway(inet::Ipv4Address(r.next_hop.getInt()));
        
        inet::NetworkInterface *ie = rt->getInterfaceByAddress(inet::Ipv4Address(r.next_hop.getInt()));
        if (!ie) ie = rt->getInterface(0); // Fallback

        route->setInterface(ie);
        route->setSourceType(inet::IRoute::MANET);
        route->setMetric(r.metric);
        rt->addRoute(route);
    }
#endif
}

} // namespace mysrc::olsrv2