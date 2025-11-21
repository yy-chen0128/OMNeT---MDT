/*
 * MDTData.h
 *
 *  Created on: Aug 10, 2025
 *      Author: yychen
 */
#ifndef MYSRC_ROUTING_MDTDATA_H_
#define MYSRC_ROUTING_MDTDATA_H_

#include <set>
#include <tuple>
#include <string>
#include "inet/networklayer/common/L3Address.h"

namespace mysrc {
namespace routing {
using namespace inet;
struct Routetuple{
    std::string type;
    L3Address source;
    L3Address pred;
    L3Address succ;
    L3Address dest;
    double cost;
    double error;
    // soft-state bookkeeping
    //simtime_t lastRefresh = SIMTIME_ZERO;
    simtime_t timeout = SIMTIME_ZERO;
    Routetuple(std::string s,L3Address sr,L3Address pre,L3Address suc,L3Address dt,double c,double e,simtime_t t)
        :type(s),source(sr),pred(pre),succ(suc),dest(dt),cost(c),error(e),timeout(t){}
    Routetuple(){
        type = "physics";
        source = L3Address();
        pred = L3Address();
        succ = L3Address();
        dest = L3Address();
        cost = -1;
        error = -1;
        timeout = simTime() + 10;//TODO
    }
};
class INET_API MDTData : public cObject
{
  protected:
    std::set<L3Address> precursorList;
    bool active;
    bool repariable;
    bool beingRepaired;
    bool validDestNum;
    unsigned int destSeqNum;
    simtime_t lifeTime; // expiration or deletion time of the route

    //MDT
    Routetuple routetuple;
    //type,source, pred, succ, dest, cost, error

  public:

    MDTData():routetuple("NONE",L3Address(),L3Address(),L3Address(),L3Address(),0,0,SIMTIME_ZERO)
    {
        active = true;
        repariable = false;
        beingRepaired = false;
        validDestNum = true;
        lifeTime = SIMTIME_ZERO;
        destSeqNum = 0;

        //MDT
    }

    virtual ~MDTData() {}

    unsigned int getDestSeqNum() const { return destSeqNum; }
    void setDestSeqNum(unsigned int destSeqNum) { this->destSeqNum = destSeqNum; }
    bool hasValidDestNum() const { return validDestNum; }
    void setHasValidDestNum(bool hasValidDestNum) { this->validDestNum = hasValidDestNum; }
    bool isBeingRepaired() const { return beingRepaired; }
    void setIsBeingRepaired(bool isBeingRepaired) { this->beingRepaired = isBeingRepaired; }
    bool isRepariable() const { return repariable; }
    void setIsRepariable(bool isRepariable) { this->repariable = isRepariable; }
    const simtime_t& getLifeTime() const { return lifeTime; }
    void setLifeTime(const simtime_t& lifeTime) { this->lifeTime = lifeTime; }
    bool isActive() const { return active; }
    void setIsActive(bool active) { this->active = active; }
    void addPrecursor(const L3Address& precursorAddr) { precursorList.insert(precursorAddr); }
    const std::set<L3Address>& getPrecursorList() const { return precursorList; }
    virtual std::string str() const;

    //MDT
    Routetuple getroutetuple() const{return routetuple;}
    void setroutetuple(const Routetuple& rt){routetuple = rt;}
    simtime_t gettimeout(){return routetuple.timeout;}
    void settimeout(simtime_t t){routetuple.timeout = t;}
    void setroutetupletype(std::string t){
        if(t == "VPoD"){
            routetuple.type = "VPoD";
        }
        else if(t == "physics"){
            routetuple.type = "physics";
            routetuple.cost = -1;
            routetuple.error = -1;
        }
        else{
            EV_ERROR<<"Wrong routetuple type!";
        }
    }
    void setroutetuplepath(L3Address s,L3Address pred,L3Address succ,L3Address d){
        routetuple.source = s;
        routetuple.pred = pred;
        routetuple.succ = succ;
        routetuple.dest = d;
    }
    void setroutetuplecost(double c,double e){
        ASSERT(routetuple.type == "VPoD");
        routetuple.cost = c;
        routetuple.error = e;
    }

};

} //
} //







#endif /* MYSRC_ROUTING_MDTDATA_H_ */
