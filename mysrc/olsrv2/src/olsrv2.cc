//
// Created by GuoYudi on 2026/3/14.
//
#include "../include/olsrv2_routing.h"
#include "../include/olsrv2.h"

#include <map>
#include <queue>
#include <set>
#include <cstring>

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
Define_Module(OLSRv2);

OLSRv2::OLSRv2()
{
}

OLSRv2::~OLSRv2()
{
    cancelAndDelete(helloTimer_);
    cancelAndDelete(tcTimer_);
}

void OLSRv2::initialize(int stage)
{
    RoutingProtocolBase::initialize(stage);

    if (stage == inet::INITSTAGE_LOCAL) {
        helloInterval_ = par("helloInterval");
        tcInterval_ = par("tcInterval");
        startJitter_ = par("startJitter");
        udpPort_ = par("udpPort");
        
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
        socket_.setCallback(this);
        socket_.bind(inet::L3Address(), udpPort_);
        socket_.setBroadcast(true);

        const char *mcast = par("multicastGroup").stringValue();
        if (!mcast || std::strlen(mcast) == 0)
            mcast = "224.0.0.109";
        multicastGroup_ = inet::L3AddressResolver().resolve(mcast);

        const char *ifname = par("interface").stringValue();
        if (ifname && std::strlen(ifname) > 0) {
            auto *ie = interfaceTable_->findInterfaceByName(ifname);
            if (ie)
                multicastIfId_ = ie->getInterfaceId();
        }
        if (multicastIfId_ < 0 && interfaceTable_->getNumInterfaces() > 0) {
            auto *ie = interfaceTable_->getInterface(0);
            if (ie)
                multicastIfId_ = ie->getInterfaceId();
        }
        if (multicastIfId_ >= 0)
            socket_.joinMulticastGroup(multicastGroup_, multicastIfId_);

        // Initialize Core Logic
        inet::Ipv4Address myIp = routingTable_->getRouterId();
        core_ = std::make_unique<Olsrv2Core>(myIp);
        nhdp_.setOriginator(inet::L3Address(myIp));

        const double hello_jitter_s = (startJitter_ > 0.0) ? uniform(0.0, startJitter_) : 0.0;
        const double tc_jitter_s = (startJitter_ > 0.0) ? uniform(0.0, startJitter_) : 0.0;
        scheduleAfter(helloInterval_ + hello_jitter_s, helloTimer_);
        scheduleAfter(tcInterval_ + tc_jitter_s, tcTimer_);
    }
}

void OLSRv2::handleMessageWhenUp(omnetpp::cMessage *msg)
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

void OLSRv2::processHelloTimer()
{
    sendHello();
    scheduleAfter(helloInterval_, helloTimer_);
}

void OLSRv2::processTcTimer()
{
    sendTc();
    scheduleAfter(tcInterval_, tcTimer_);
}

void OLSRv2::sendHello()
{
    double now = simTime().dbl();
    auto pkt = nhdp_.generateHello(now);
    socket_.sendTo(pkt, multicastGroup_, udpPort_);
}

void OLSRv2::sendTc()
{
    double now = simTime().dbl();
    if (core_) {
        auto pkt = core_->generateTc(nhdp_.getDb(), now);
        if (pkt)
            socket_.sendTo(pkt, multicastGroup_, udpPort_);
    }
}

void OLSRv2::socketDataArrived(inet::UdpSocket *socket, inet::Packet *packet)
{
    processOlsrPacket(packet);
    delete packet;
}

void OLSRv2::socketErrorArrived(inet::UdpSocket *socket, inet::Indication *indication)
{
    // Handle socket error
    delete indication;
}

void OLSRv2::socketClosed(inet::UdpSocket *socket)
{
    // Handle closed
}

void OLSRv2::processOlsrPacket(inet::Packet *packet)
{
    auto chunk = packet->peekAtFront<inet::Olsrv2ControlPacket>();
    if (!chunk) return;

    inet::L3Address source;
    auto l3Tag = packet->getTag<inet::L3AddressInd>();
    if (l3Tag) source = l3Tag->getSrcAddress();

    if (auto hello = inet::dynamicPtrCast<const inet::Olsrv2HelloPacket>(chunk)) {
        nhdp_.processHello(hello.get(), source, simTime().dbl());
        printProtocolState("RX_HELLO");
    } else if (auto tc = inet::dynamicPtrCast<const inet::Olsrv2TcGroup>(chunk)) {
        if (core_) {
            core_->processTcAndRecompute(tc.get(), simTime().dbl());
            updateRoutingTable();
        }
    }
}

inet::INetfilter::IHook::Result OLSRv2::datagramPreRoutingHook(inet::Packet *datagram)
{
    return inet::INetfilter::IHook::ACCEPT;
}

void OLSRv2::updateRoutingTable()
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

    printProtocolState("UPDATE_ROUTING_TABLE");
}

void OLSRv2::printProtocolState(const char *tag) const
{
    EV_INFO << "\n\n================ OLSRV2 DEBUG ================\n";
    EV_INFO << "Tag=" << (tag ? tag : "") << " Time=" << simTime() << omnetpp::endl;

    printNhdpState();

    if (core_) {
        core_->state().printOlsrv2State();
    }

    printIpv4RoutingTable();
    EV_INFO << "==============================================\n\n";
}

void OLSRv2::printNhdpState() const
{
    EV_INFO << "\n=== NHDP DB ===\n";

    const auto links = nhdp_.getDb().getLinks();
    for (const auto& l : links) {
        const char *status = "UNKNOWN";
        switch (l.status) {
            case NhdpLinkStatus::Pending: status = "PENDING"; break;
            case NhdpLinkStatus::Lost: status = "LOST"; break;
            case NhdpLinkStatus::Symmetric: status = "SYMMETRIC"; break;
            case NhdpLinkStatus::Heard: status = "HEARD"; break;
        }

        EV_INFO << "LinkId=" << l.id
                << " Status=" << status
                << " NeighIface=" << l.neighbor_iface_addr
                << " NeighOriginator=" << l.neighbor_originator
                << " TwoHopCount=" << l.twohop_addresses.size()
                << omnetpp::endl;
    }
    EV_INFO << "LinkCount=" << links.size() << "\n";

    const auto mpr = nhdp_.calculateMprSet();
    EV_INFO << "\n=== MPR SET (from NHDP) ===\n";
    for (const auto& a : mpr) {
        EV_INFO << "MPR=" << a << omnetpp::endl;
    }
    EV_INFO << "MPR count=" << mpr.size() << "\n";
}

void OLSRv2::printIpv4RoutingTable() const
{
    EV_INFO << "\n=== IPV4 ROUTING TABLE ===\n";

    auto *rt = dynamic_cast<inet::Ipv4RoutingTable *>(routingTable_);
    if (!rt) {
        EV_INFO << "RoutingTable cast failed\n";
        return;
    }

    for (int i = 0; i < rt->getNumRoutes(); ++i) {
        const auto *r = rt->getRoute(i);
        if (!r)
            continue;

        const auto *ie = r->getInterface();
        EV_INFO << "Dest=" << r->getDestination()
                << " Netmask=" << r->getNetmask()
                << " Gateway=" << r->getGateway()
                << " Iface=" << (ie ? ie->getInterfaceName() : "")
                << " Metric=" << r->getMetric()
                << " SrcType=" << static_cast<int>(r->getSourceType())
                << omnetpp::endl;
    }
    EV_INFO << "RouteCount=" << rt->getNumRoutes() << "\n";
}

} // namespace mysrc::olsrv2
