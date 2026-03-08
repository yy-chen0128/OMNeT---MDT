/*
 * OLSR_State.h
 *
 *  Created on: Mar 6, 2026
 *      Author: yychen
 */

#ifndef MYSRC_OLSR_OLSR_STATE_H_
#define MYSRC_OLSR_OLSR_STATE_H_

#include <unordered_map>
#include "../olsr/OLSR_Repositories.h"

namespace mysrc {
namespace olsr {
class OlsrState : public cObject {
  protected:
    // sets/containers
    std::vector<OlsrLinkTuple*> linkSet;
    std::vector<OlsrNeighborTuple*> nbSet;
    std::vector<OlsrTwoHopTuple*> nb2HopSet;
    std::set<L3Address> mprSet;
    std::vector<OlsrMprSelectorTuple*> mprSelSet;
    std::vector<OlsrTopologyTuple*> topologySet;
    std::vector<OlsrDupTuple*> dupSet;


    // small helpers
    SimTimeType helloValidity; // e.g., allowedHelloLoss * helloInterval
    SimTimeType tcValidity;    // TC entry lifetime
    int nextTcGroupId = 1;

  public:
    OlsrState() : helloValidity(6.0), tcValidity(60.0) {}
    virtual ~OlsrState() { clearAll(); }

    // ------- link/neighbor management -------
    OlsrLinkTuple* findLink(const L3A& local, const L3A& nb);
    //OlsrLinkTuple* addOrUpdateLink(const L3A& local, const L3A& nb, SimTimeType expiry);
    OlsrLinkTuple* addLink(const L3A& local, const L3A& nb);
    OlsrLinkTuple* getOrCreateLink(const L3A& local, const L3A& nb);
    bool removeLink(const L3A& local, const L3A& nb);

    std::vector<OlsrNeighborTuple*>& getNbSet(){return nbSet;}
    OlsrNeighborTuple* findNeighbor(const L3A& mainAddr);
    OlsrNeighborTuple* addOrUpdateNeighbor(const L3A& mainAddr, uint8_t willingness, SimTimeType expiry);
    bool removeNeighbor(const L3A& mainAddr);

    // ------- 2-hop management -------
    std::vector<OlsrTwoHopTuple*>& getnb2HopSet(){return nb2HopSet;}
    OlsrTwoHopTuple* findTwoHop(const L3A& nbMain, const L3A& nb2hop);
    void addOrUpdateTwoHop(const L3Address& selfMainAddr,const L3A& nbMain, const L3A& nb2hop, SimTimeType expiry);
    bool removeTwoHop(const L3A& nbMain, const L3A& nb2hop);

    // ------- topology -------
    std::vector<OlsrTopologyTuple*>& getTopologySet(){return topologySet;}
    OlsrTopologyTuple* findTopology(const L3A& dest, const L3A& last);
    OlsrTopologyTuple* addOrUpdateTopology(const L3A& dest, const L3A& last, uint16_t seq, SimTimeType expiry);
    bool removeTopology(const L3A& dest, const L3A& last);

    // ------- duplicates -------
    OlsrDupTuple* findDup(const L3A& origin, uint16_t seq);
    void addDup(const L3A& origin, uint16_t seq, SimTimeType expiry);
    void purgeExpiredDup(SimTimeType now);

    //MPR set
    void SetMprSet(std::set<L3Address> mSet)
    {
        mprSet = mSet;
    }

    std::set<L3Address> GetMprSet() const
    {
        return mprSet;
    }

    // ------- MPR selectors -------
    std::vector<OlsrMprSelectorTuple*>& getmprSelSet(){return mprSelSet;}
    void addMprSelector(const L3A& addr, SimTimeType expiry);
    void removeMprSelector(const L3A& addr);
    bool isMprSelector(const L3A& addr) const;

    // ------- housekeeping -------
    void expireTuples(SimTimeType now); // walk all sets and remove expired tuples
    void clearAll(); // free all stored cObjects

    //
    const std::vector<OlsrLinkTuple*>& getLinkSet() const { return linkSet; }

    // utility
    uint32_t nextTcGroupIdAndInc() { return nextTcGroupId++; }

    //print
    void printLinkSet();
    void printNeighborSet();
    void printTwoHopSet();
    void printMprSet();
    void printMprSelectorSet();
    void printTopologySet();
    void printDupSet();
    void printOlsrState();
};
}//olsr
} // namespace mysrc



#endif /* MYSRC_OLSR_OLSR_STATE_H_ */
