/*
 * Maintenance_Protocol.h
 *
 *  Created on: Aug 31, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_MAINTENANCEPROTOCOL_H_
#define MYSRC_ROUTING_MAINTENANCEPROTOCOL_H_
#include <omnetpp.h>
#include "MDTControlPackets_m.h"
#include "NeighborEntry.h"
#include "NeighborSetReply.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Point_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>


namespace mysrc {
namespace routing {
using namespace inet;

class MDTRouting;

enum class MaintenanceState{
    INACTIVE = 0,
    JOINSTAGE = 1,
    CYCLE_MAINTAIN = 2,
    CYCLE_MAINTAIN_PROCESS = 3,
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point3;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay3;
// prepare 2D triangulation
typedef CGAL::Delaunay_triangulation_2<K> Delaunay2;
typedef K::Point_2 Point2;
typedef Delaunay2::Vertex_handle Vertex_handle2;
typedef Delaunay2::Vertex_circulator Vertex_circulator2;

class INET_API MaintenanceProtocol: public cSimpleModule{
private:
    MDTRouting *core = nullptr;
    cRNG *rng = nullptr;
    simtime_t MaintenanceInterval;
    simtime_t NSRPTimeout;
    cMessage *MaintenanceTimer = nullptr;
    MaintenanceState stage = MaintenanceState::INACTIVE;
    //stage one
    std::map<L3Address,simtime_t>waitingReplyNodes;
    std::set<L3Address> inquiredNodes;
    //std::vector<L3Address> waitingReplyNodes;


    virtual void initialize(int stage) override;
    virtual void handleMessage(cMessage *msg) override;
public:
    MaintenanceProtocol(){}
    virtual~MaintenanceProtocol();

    void setStage(MaintenanceState x){stage = x;}
    MaintenanceState getStage(){return stage;}

    Packet *createNeighborSetRequestTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<double>& coords = std::vector<double>());
    void processNeighborSetRequestChunk(const Ptr<NeighborSetRequest>& neighborSetRequestChunk);

    Packet *createNeighborSetReplyTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<NeighborEntry>& entrys);
    void processNeighborSetReplyChunk(const Ptr<NeighborSetReply>& neighborSetReplyChunk);

    Packet *createNeighborSetNotificationTo(const L3Address& dest,const L3Address& src,const L3Address& pred,const std::vector<double>& coords = std::vector<double>());
    void processNeighborSetNotificationChunk(const Ptr<NeighborSetNotification>& neighborSetNotificationChunk);

    void beginCycleMaintain();
    void findroutefornewdt(std::vector<NeighborEntry> newdtset);
    bool satisfiesC1_CGAL(const L3Address &vAddr,
                             const std::set<L3Address> &requested,MDTRouting *core);

    void AddwaitingReplyNodes(const L3Address& node,simtime_t time){
        auto [it, inserted] = waitingReplyNodes.insert({node, time});
        if (!inserted) {
            EV_DETAIL << "node " << node << " already in map, ignored\n";
        }
    }
    bool AddinquiredNodes(const L3Address& node){
        auto [it,success] = inquiredNodes.insert(node);
        if(success)EV_INFO << "node: " << node << " add in inquiredNodes\n";
        else
            EV_INFO << "node: " << node << "add fail\n";
        return success;
    }
    bool findNodeinInquiredNodes(const L3Address& addr){
        if (inquiredNodes.count(addr) == 1) {
            return true;
        }
        return false;
    }

};

}
}

#endif /* MYSRC_ROUTING_MAINTENANCEPROTOCOL_H_ */
