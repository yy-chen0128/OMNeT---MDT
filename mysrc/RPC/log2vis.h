/*
 * log2vis.h
 *
 *  Created on: Oct 9, 2025
 *      Author: yychen
 */

#ifndef MYSRC_RPC_LOG2VIS_H_
#define MYSRC_RPC_LOG2VIS_H_

#include <omnetpp.h>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "inet/networklayer/contract/IRoutingTable.h"
#include "../routing/NeighborEntry.h"
#include "../routing/MDTData.h"

#include "RpcBridge.h"

using json = nlohmann::json;
using namespace omnetpp;
using namespace inet;
using namespace mysrc::routing;

extern RpcBridge* gRpcBridge;
extern void writeJsonToLocalFile(const nlohmann::json &j);

// ------------------------- simple marco -------------------------
#define RPC_EV(level, msg) do { \
    std::ostringstream _oss; _oss << msg; \
    EV << "[" << level << "] " << _oss.str() << endl; \
    if (gRpcBridge) { \
        json j = {{"type","log"}, {"component", this->getFullPath()}, {"level", level}, {"message", _oss.str()}, {"timestamp", simTime().dbl()}}; \
        gRpcBridge->enqueueOutgoingLog(j); \
    } \
} while(0)

// ------------------------- helper-------------------------

inline int moduleIndexOrMinusOne(cModule* module) {
    if (!module) return -1;
    try { return module->getIndex(); } catch(...) { return -1; }
}
static std::string l3AddressToString(const inet::L3Address &addr) {
    return addr.str();
}

static int ipv4StringToNodeId(const std::string &addrStr) {
    auto pos = addrStr.find_last_of('.');
    if (pos == std::string::npos) return -1;
    try {
        int last = std::stoi(addrStr.substr(pos+1));
        return last - 1; // mapping rule
    } catch (...) {
        return -1;
    }
}

// sender
inline void sendJsonFromModule(cModule* module, json j) {
    if (module) {
        try { j["modulePath"] = module->getFullPath(); } catch(...) {}
        int idx = moduleIndexOrMinusOne(module);
        if (idx >= 0) j["nodeId"] = idx;
    }
    // timestamp
    if (j.find("timestamp") == j.end()) j["timestamp"] = simTime().dbl();
    //if (gRpcBridge) gRpcBridge->enqueueOutgoingLog(j);
    writeJsonToLocalFile(j);
}

// ------------------------- interface -------------------------

// position
// posSeq (optional) speed/vel (optional)
inline void logPosition(cModule* module, const std::vector<double>& posVec,
                        const std::vector<double>& velVec, int posSeq = -1)
{
    json j;
    j["type"] = "position";
    j["pos"] = posVec;
    if (!velVec.empty()) j["vel"] = velVec;
    if (posSeq >= 0) j["posSeq"] = posSeq;
    sendJsonFromModule(module, j);
}

// routingTable
// routes: vector of json entries, each entry contains at least destId or dest, nextHopId/nextHop, metric, iface(optional)

static std::vector<json> buildRoutingJsonList(cModule* module, const IRoutingTable *routingTable) {
    std::vector<json> routesJson;
    int n = routingTable->getNumRoutes();
    routesJson.reserve(n);

    for (int i = 0; i < n; ++i) {
        IRoute *rt = routingTable->getRoute(i);
        if (!rt) continue;
        if(rt->getDestinationAsGeneric().str() == "10.0.0.0" || rt->getDestinationAsGeneric().str() == "127.0.0.0")continue;

        json r;
        // destination and prefix
        r["dest"] = ipv4StringToNodeId(rt->getDestinationAsGeneric().str());
        r["prefixLength"] = rt->getPrefixLength();

        // nextHop
        r["nextHop"] = ipv4StringToNodeId(rt->getNextHopAsGeneric().str());

        // maybe metric (if available)
        try {
            r["metric"] = rt->getMetric(); // if IRoute supports getMetric()
        } catch (...) {}

        // optional interface name
        try {
            if (rt->getInterface());
                //r["iface"] = rt->getInterface()->getName();
        } catch (...) {}

        // protocol-specific data
        MDTData *routeDestData = rt ? dynamic_cast<MDTData *>(rt->getProtocolData()) : nullptr;
        if (routeDestData) {
            Routetuple tuple = routeDestData->getroutetuple();
            json t;
            t["source"] = ipv4StringToNodeId(tuple.source.str()); // or tuple.source as string
            t["dest"] = ipv4StringToNodeId(tuple.dest.str());
            t["pred"] = ipv4StringToNodeId(tuple.pred.str()); // if pid or id
            t["succ"] = ipv4StringToNodeId(tuple.succ.str());
            t["timeout"] = tuple.timeout.dbl(); // double or simtime
            r["routetuple"] = t;
        }

        routesJson.push_back(std::move(r));
    }

    return routesJson;
}
inline void logRoutingTable(cModule* module, int tableSeq, const IRoutingTable *routingTable)
{
    json j;
    j["type"] = "routingTable";
    j["tableSeq"] = tableSeq;

    auto routes = buildRoutingJsonList(module, routingTable);
    j["routes"] = routes;
    sendJsonFromModule(module, j);
}
// neighborList
// neighbors: vector of json entries {neighborId, rssi, linkQuality, neighborPath (optional)}


static std::vector<json> exportNeighborMapAsJson(cModule* module, int neighSeq, const std::map<inet::L3Address, NeighborEntry> &neighbors) {
    std::vector<json> neighJson;
    neighJson.reserve(neighbors.size());

    for (const auto &p : neighbors) {
        const inet::L3Address &addr = p.first;
        const NeighborEntry &entry = p.second;

        json n;
        std::string addrStr = l3AddressToString(addr);
        n["addr"] = addrStr;
        int nodeId = ipv4StringToNodeId(addrStr);
        if (nodeId >= 0) n["neighborId"] = nodeId;
        // timeout: simtime_t -> double seconds
        n["timeout"] = entry.timeout.dbl();
        n["attachedToDT"] = entry.attachedToDT;

        // coords: entry.coords is vector<double>
        if (!entry.coords.empty()) {
            n["coords"] = entry.coords; // json will serialize vector<double>
        }

        // optional type
        if (!entry.type.empty())
            n["type"] = entry.type;

        neighJson.push_back(std::move(n));
    }
    return neighJson;
}
inline void logNeighborList(cModule* module, int neighSeq, const std::map<inet::L3Address, NeighborEntry> &neighbors)
{
    json j;
    j["type"] = "neighborList";
    j["neighSeq"] = neighSeq;
    std::vector<json> neighJson = exportNeighborMapAsJson(module,neighSeq,neighbors);
    j["neighbors"] = neighJson;
    sendJsonFromModule(module, j);
}
inline void logLocalKnownNodes(cModule* module, int seq, const std::map<inet::L3Address, NeighborEntry> &knownNodes) {
    nlohmann::json j;
    j["type"] = "localKnown";            // localKnown
    j["seq"] = seq;
    j["desc"] = "local_known_nodes";
    // reuse export function
    std::vector<nlohmann::json> arr = exportNeighborMapAsJson(module, seq, knownNodes);
    j["knownNodes"] = arr;
    sendJsonFromModule(module, j);
}

inline void logDtNeighborSet(cModule* module, int seq, const std::map<inet::L3Address, NeighborEntry> &dtNeighbors) {
    nlohmann::json j;
    j["type"] = "dtNeighborList";        // dtNeighborList
    j["seq"] = seq;
    j["desc"] = "dt_neighbor_set";
    std::vector<nlohmann::json> arr = exportNeighborMapAsJson(module, seq, dtNeighbors);
    j["neighbors"] = arr;
    sendJsonFromModule(module, j);
}
// msgEvent
// eventType: "create","send","recv","forward","drop","deliver"
// messageId: use cMessage::getId()
// srcNodeId/dstNodeId: application-level source/destination id
// prevHop/nextHop: -1 means unknown
// payload: any json mata data (size, flowId, seq )
inline void logMsgEvent(cModule* module,
                        const std::string& eventType,
                        long long messageId,
                        int srcNodeId = -1,
                        int dstNodeId = -1,
                        int prevHop = -1,
                        int nextHop = -1,
                        const json& payload = json::object())
{
    json j;
    j["type"] = "msgEvent";
    j["eventType"] = eventType;
    j["messageId"] = messageId;
    if (srcNodeId >= 0) j["srcNodeId"] = srcNodeId;
    if (dstNodeId >= 0) j["dstNodeId"] = dstNodeId;
    if (prevHop >= 0) j["prevHop"] = prevHop;
    if (nextHop >= 0) j["nextHop"] = nextHop;
    if (!payload.empty()) j["payloadSummary"] = payload;
    sendJsonFromModule(module, j);
}

// finish
inline void logFinish(cModule* module, const std::string& reason = "simulation_finish")
{
    json j;
    j["type"] = "finish";
    j["reason"] = reason;
    sendJsonFromModule(module, j);
}

// ------------------------- examples -------------------------
/*

// 1
RPC_EV("INFO", "routing updated, nextHop=" << nextHop);

// 2
Coord pos = mobility->getLocation();
logPosition(this, pos.x, pos.y, pos.z, mobility->getSpeed(), posSeq++);

// 3
std::vector<json> routes;
routes.push_back({{"destId", 12}, {"nextHopId", 8}, {"metric", 3}, {"iface","wlan0"}});
logRoutingTable(this, tableSeq++, routes);

// 4
std::vector<json> neighs;
neighs.push_back({{"neighborId", 8}, {"rssi",-70}, {"linkQuality",0.9}});
logNeighborList(this, neighSeq++, neighs);

// 5
logMsgEvent(this, "send", msg->getId(), srcId, dstId, -1, nextHop, {{"size", msg->getByteLength()}, {"flowId", flow}, {"seq", seq}});

*/

#endif /* MYSRC_RPC_LOG2VIS_H_ */
