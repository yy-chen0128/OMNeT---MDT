/*
 * OLSR_State.cc
 *
 *  Created on: Mar 6, 2026
 *      Author: yychen
 */

#include "../olsr/OLSR_State.h"

namespace mysrc {
namespace olsr {
//
// Helper: deletePointerAndErase - deletes pointer at index i from vector and erases element
//
template <typename T>
static void deletePointerAndErase(std::vector<T*>& vec, size_t i)
{
    if (i < vec.size()) {
        delete vec[i];
        vec.erase(vec.begin() + i);
    }
}

//
// findLink
//
OlsrLinkTuple* OlsrState::findLink(const L3A& local, const L3A& nb)
{
    for (auto *lt : linkSet) {
        if (lt->localIfaceAddr == local && lt->neighborIfaceAddr == nb)
            return lt;
    }
    return nullptr;
}

//
// addOrUpdateLink
//
OlsrLinkTuple* OlsrState::addLink(const L3A& local, const L3A& nb)
{
    OlsrLinkTuple *lt = new OlsrLinkTuple(local, nb);
    linkSet.push_back(lt);

    EV_DETAIL << "OLSR: add link tuple local="
              << local.str()
              << " nb=" << nb.str() << "\n";

    return lt;
}
OlsrLinkTuple* OlsrState::getOrCreateLink(const L3A& local, const L3A& nb)
{
    OlsrLinkTuple *lt = findLink(local, nb);

    if (!lt)
        lt = addLink(local, nb);

    return lt;
}
//
// removeLink
//
bool OlsrState::removeLink(const L3A& local, const L3A& nb)
{
    for (size_t i = 0; i < linkSet.size(); ++i) {
        OlsrLinkTuple *lt = linkSet[i];
        if (lt->localIfaceAddr == local && lt->neighborIfaceAddr == nb) {
            deletePointerAndErase(linkSet, i);
            EV_DETAIL << "OLSR: removed link tuple local=" << local.str() << " nb=" << nb.str() << "\n";
            return true;
        }
    }
    return false;
}

//
// Neighbor: findNeighbor
//
OlsrNeighborTuple* OlsrState::findNeighbor(const L3A& mainAddr)
{
    for (auto *nb : nbSet) {
        if (nb->mainAddr == mainAddr)
            return nb;
    }
    return nullptr;
}

//
// addOrUpdateNeighbor
//
OlsrNeighborTuple* OlsrState::addOrUpdateNeighbor(const L3A& mainAddr, uint8_t willingness, SimTimeType expiry)
{
    OlsrNeighborTuple *nt = findNeighbor(mainAddr);
    if (nt) {
        nt->willingness = willingness;
        nt->time = expiry;
        EV_DETAIL << "OLSR: updated neighbor " << mainAddr.str() << " willingness=" << (int)willingness
                  << " expiry=" << expiry << "\n";
        return nt;
    }
    nt = new OlsrNeighborTuple();
    nt->mainAddr = mainAddr;
    nt->willingness = willingness;
    nt->time = expiry;
    nbSet.push_back(nt);
    EV_DETAIL << "OLSR: added neighbor " << mainAddr.str() << " willingness=" << (int)willingness
              << " expiry=" << expiry << "\n";
    return nt;
}

//
// removeNeighbor
//
bool OlsrState::removeNeighbor(const L3A& mainAddr)
{
    for (size_t i = 0; i < nbSet.size(); ++i) {
        if (nbSet[i]->mainAddr == mainAddr) {
            deletePointerAndErase(nbSet, i);
            EV_DETAIL << "OLSR: removed neighbor " << mainAddr.str() << "\n";
            return true;
        }
    }
    return false;
}

//
// Two-hop: findTwoHop
//
OlsrTwoHopTuple* OlsrState::findTwoHop(const L3A& nbMain, const L3A& nb2hop)
{
    for (auto *t : nb2HopSet) {
        if (t->nbMainAddr == nbMain && t->nb2hopAddr == nb2hop)
            return t;
    }
    return nullptr;
}

//
// addOrUpdateTwoHop
//
void OlsrState::addOrUpdateTwoHop(const L3Address& selfMainAddr ,const L3A& nbMain, const L3A& nb2hop, SimTimeType expiry)
{
    if (nb2hop == selfMainAddr)
        return;
    OlsrTwoHopTuple *t = findTwoHop(nbMain, nb2hop);
    if (t) {
        t->time = expiry;
        EV_INFO << "OLSR: updated 2-hop nbMain=" << nbMain.str() << " nb2hop=" << nb2hop.str()
                  << " expiry=" << expiry << "\n";
        return;
    }
    t = new OlsrTwoHopTuple();
    t->nbMainAddr = nbMain;
    t->nb2hopAddr = nb2hop;
    t->time = expiry;
    nb2HopSet.push_back(t);
    EV_INFO << "OLSR: added 2-hop nbMain=" << nbMain.str() << " nb2hop=" << nb2hop.str()
              << " expiry=" << expiry << "\n";
}

//
// removeTwoHop
//
bool OlsrState::removeTwoHop(const L3A& nbMain, const L3A& nb2hop)
{
    for (size_t i = 0; i < nb2HopSet.size(); ++i) {
        if (nb2HopSet[i]->nbMainAddr == nbMain && nb2HopSet[i]->nb2hopAddr == nb2hop) {
            deletePointerAndErase(nb2HopSet, i);
            EV_DETAIL << "OLSR: removed 2-hop nbMain=" << nbMain.str() << " nb2hop=" << nb2hop.str() << "\n";
            return true;
        }
    }
    return false;
}

//
// Topology: findTopology
//
OlsrTopologyTuple* OlsrState::findTopology(const L3A& dest, const L3A& last)
{
    for (auto *tt : topologySet) {
        if (tt->destination == dest && tt->lastAddr == last)
            return tt;
    }
    return nullptr;
}

//
// addOrUpdateTopology
//
OlsrTopologyTuple* OlsrState::addOrUpdateTopology(const L3A& dest, const L3A& last, uint16_t seq, SimTimeType expiry)
{
    OlsrTopologyTuple *tt = findTopology(dest, last);
    if (tt) {
        // Update seq/time only if seq is newer or time extended.
        // Simple policy: if incoming seq >= stored seq then update.
        if (seq > tt->seq) {
            tt->seq = seq;
            tt->time = expiry;
            EV_DETAIL << "OLSR: updated topology dest=" << dest.str() << " last=" << last.str()
                      << " seq=" << seq << " expiry=" << expiry << "\n";
        }
        return tt;
    }
    tt = new OlsrTopologyTuple();
    tt->destination = dest;
    tt->lastAddr = last;
    tt->seq = seq;
    tt->time = expiry;
    topologySet.push_back(tt);
    EV_DETAIL << "OLSR: added topology dest=" << dest.str() << " last=" << last.str()
              << " seq=" << seq << " expiry=" << expiry << "\n";
    return tt;
}

//
// removeTopology
//
bool OlsrState::removeTopology(const L3A& dest, const L3A& last)
{
    for (size_t i = 0; i < topologySet.size(); ++i) {
        if (topologySet[i]->destination == dest && topologySet[i]->lastAddr == last) {
            deletePointerAndErase(topologySet, i);
            EV_DETAIL << "OLSR: removed topology dest=" << dest.str() << " last=" << last.str() << "\n";
            return true;
        }
    }
    return false;
}

//
// Duplicate set: findDup
//
OlsrDupTuple* OlsrState::findDup(const L3A& origin, uint16_t seq)
{
    for (auto *d : dupSet) {
        if (d->originator == origin && d->seqNum == seq)
            return d;
    }
    return nullptr;
}

//
// addDup
//
void OlsrState::addDup(const L3A& origin, uint16_t seq, SimTimeType expiry)
{
    OlsrDupTuple *d = findDup(origin, seq);
    if (d) {
        d->time = expiry;
        EV_DETAIL << "OLSR: updated dup origin=" << origin.str() << " seq=" << seq << " expiry=" << expiry << "\n";
        return;
    }
    d = new OlsrDupTuple();
    d->originator = origin;
    d->seqNum = seq;
    d->time = expiry;
    dupSet.push_back(d);
    EV_DETAIL << "OLSR: added dup origin=" << origin.str() << " seq=" << seq << " expiry=" << expiry << "\n";
}

//
// purgeExpiredDup
//
void OlsrState::purgeExpiredDup(SimTimeType now)
{
    for (int i = (int)dupSet.size() - 1; i >= 0; --i) {
        if (dupSet[i]->time <= now) {
            EV_DETAIL << "OLSR: purging dup origin=" << dupSet[i]->originator.str()
                      << " seq=" << dupSet[i]->seqNum << "\n";
            deletePointerAndErase(dupSet, i);
        }
    }
}

//
// MPR selectors management
//
void OlsrState::addMprSelector(const L3A& addr, SimTimeType expiry)
{
    for (auto *m : mprSelSet) {
        if (m->mainAddr == addr) {
            m->time = expiry;
            EV_DETAIL << "OLSR: updated mpr selector " << addr.str() << " expiry=" << expiry << "\n";
            return;
        }
    }
    OlsrMprSelectorTuple *ms = new OlsrMprSelectorTuple();
    ms->mainAddr = addr;
    ms->time = expiry;
    mprSelSet.push_back(ms);
    EV_DETAIL << "OLSR: added mpr selector " << addr.str() << " expiry=" << expiry << "\n";
}

void OlsrState::removeMprSelector(const L3A& addr)
{
    for (int i = (int)mprSelSet.size() - 1; i >= 0; --i) {
        if (mprSelSet[i]->mainAddr == addr) {
            deletePointerAndErase(mprSelSet, i);
            EV_DETAIL << "OLSR: removed mpr selector " << addr.str() << "\n";
            return;
        }
    }
}

bool OlsrState::isMprSelector(const L3A& addr) const
{
    for (auto *m : mprSelSet)
        if (m->mainAddr == addr)
            return true;
    return false;
}

//
// expireTuples: remove expired elements from all sets
//
void OlsrState::expireTuples(SimTimeType now)
{
    EV_DETAIL << "OLSR: expireTuples now=" << now << "\n";

    // linkSet
    for (int i = (int)linkSet.size() - 1; i >= 0; --i) {
        if (linkSet[i]->time <= now) {
            EV_DETAIL << "OLSR: expiring link local=" << linkSet[i]->localIfaceAddr.str()
                      << " nb=" << linkSet[i]->neighborIfaceAddr.str() << "\n";
            deletePointerAndErase(linkSet, i);
        }
    }

    // nbSet
    for (int i = (int)nbSet.size() - 1; i >= 0; --i) {
        if (nbSet[i]->time <= now) {
            EV_DETAIL << "OLSR: expiring neighbor " << nbSet[i]->mainAddr.str() << "\n";
            deletePointerAndErase(nbSet, i);
        }
    }

    // nb2HopSet
    for (int i = (int)nb2HopSet.size() - 1; i >= 0; --i) {
        if (nb2HopSet[i]->time <= now) {
            EV_DETAIL << "OLSR: expiring 2hop nbMain=" << nb2HopSet[i]->nbMainAddr.str()
                      << " nb2hop=" << nb2HopSet[i]->nb2hopAddr.str() << "\n";
            deletePointerAndErase(nb2HopSet, i);
        }
    }

    // mprSelSet
    for (int i = (int)mprSelSet.size() - 1; i >= 0; --i) {
        if (mprSelSet[i]->time <= now) {
            EV_DETAIL << "OLSR: expiring mpr selector " << mprSelSet[i]->mainAddr.str() << "\n";
            deletePointerAndErase(mprSelSet, i);
        }
    }

    // topologySet
    for (int i = (int)topologySet.size() - 1; i >= 0; --i) {
        if (topologySet[i]->time <= now) {
            EV_DETAIL << "OLSR: expiring topology dest=" << topologySet[i]->destination.str()
                      << " last=" << topologySet[i]->lastAddr.str() << "\n";
            deletePointerAndErase(topologySet, i);
        }
    }

    // dupSet
    purgeExpiredDup(now);

}

//
// clearAll
//
void OlsrState::clearAll()
{
    for (auto *p : linkSet) delete p;
    linkSet.clear();

    for (auto *p : nbSet) delete p;
    nbSet.clear();

    for (auto *p : nb2HopSet) delete p;
    nb2HopSet.clear();

    for (auto *p : mprSelSet) delete p;
    mprSelSet.clear();

    for (auto *p : topologySet) delete p;
    topologySet.clear();

    for (auto *p : dupSet) delete p;
    dupSet.clear();

}

void OlsrState::printLinkSet()
{
    EV_INFO << "\n=== LINK SET ===\n";

    for (auto *t : linkSet) {
        EV_INFO << "Neighbor=" << t->neighborIfaceAddr
                << " Local=" << t->localIfaceAddr
                << " SymTime=" << t->symTime
                << " AsymTime=" << t->asymTime
                << " Expire=" << t->time
                << endl;
    }

    EV_INFO << "LinkSet size=" << linkSet.size() << "\n";
}
void OlsrState::printNeighborSet()
{
    EV_INFO << "\n=== NEIGHBOR SET ===\n";

    for (auto *t : nbSet) {
        EV_INFO << "NeighborMainAddr=" << t->mainAddr
                << " Status=" << t->status
                << " Willingness=" << (int)t->willingness
                << endl;
    }

    EV_INFO << "NeighborSet size=" << nbSet.size() << "\n";
}
void OlsrState::printTwoHopSet()
{
    EV_INFO << "\n=== TWO HOP SET ===\n";

    for (auto *t : nb2HopSet) {
        EV_INFO << "Neighbor=" << t->nbMainAddr
                << " TwoHop=" << t->nb2hopAddr
                << " Expire=" << t->time
                << endl;
    }
    EV_INFO << "TwoHopSet size=" << nb2HopSet.size() << "\n";
}
void OlsrState::printMprSet()
{
    EV_INFO << "\n=== MPR SET ===\n";

    for (auto &a : mprSet) {
        EV_INFO << "MPR=" << a << endl;
    }

    EV_INFO << "MPR count=" << mprSet.size() << "\n";
}
void OlsrState::printMprSelectorSet()
{
    EV_INFO << "\n=== MPR SELECTOR SET ===\n";

    for (auto *t : mprSelSet) {
        EV_INFO << "Selector=" << t->mainAddr
                << " Expire=" << t->time
                << endl;
    }

    EV_INFO << "MPR Selector count=" << mprSelSet.size() << "\n";
}
void OlsrState::printTopologySet()
{
    EV_INFO << "\n=== TOPOLOGY SET ===\n";

    for (auto *t : topologySet) {
        EV_INFO << "Dest=" << t->destination
                << " LastHop=" << t->lastAddr
                << " Seq=" << t->seq
                << " Expire=" << t->time
                << endl;
    }

    EV_INFO << "TopologySet size=" << topologySet.size() << "\n";
}
void OlsrState::printDupSet()
{
    EV_INFO << "\n=== DUPLICATE SET ===\n";

    for (auto *t : dupSet) {
        EV_INFO << "Originator=" << t->originator
                << " Seq=" << t->seqNum
                << " Retransmitted=" << t->retransmitted
                << " Expire=" << t->time
                << endl;
    }

    EV_INFO << "DupSet size=" << dupSet.size() << "\n";
}
void OlsrState::printOlsrState()
{
    EV_INFO << "\n\n================ OLSR STATE ================\n";

    printLinkSet();
    printNeighborSet();
    printTwoHopSet();
    printMprSet();
    printMprSelectorSet();
    printTopologySet();
    printDupSet();

    EV_INFO << "============================================\n\n";
}

}//olsr
} // namespace mysrc




