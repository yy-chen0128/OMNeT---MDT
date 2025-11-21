/*
 * test_dt_rounds.cc
 *
 *  Created on: Sep 27, 2025
 *      Author: yychen
 */
// tests/test_dt_rounds.cc
#include <gtest/gtest.h>

#include <regex>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <map>
#include <set>
#include <iostream>

#include <omnetpp.h> // for SIMTIME_ZERO type if NeighborEntry uses simtime_t
#include "../../inet-4.5.4/src/inet/networklayer/common/L3Address.h"
#include "../../inet-4.5.4/src/inet/networklayer/ipv4/Ipv4Header_m.h"
#include "../../inet-4.5.4/src/inet/networklayer/ipv4/Ipv4Route.h"

#include "../mysrc/routing/NeighborDT.h"
#include "../mysrc/routing/NeighborEntry.h"

using namespace std;
using namespace inet;
using namespace mysrc::routing;

// ---- parser types ----
struct NodeDesc {
    int idx;
    array<double,3> coord;
    vector<int> neighbors;
};

struct RoundData {
    vector<NodeDesc> initial;
    vector<int> deleted;
    vector<NodeDesc> afterDelete;
    vector<int> added;
    vector<NodeDesc> finalNodes;
};

// helper: trim
static string trim(const string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a==string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a,b-a+1);
}

// parse index list like [1, 2, 3]
static vector<int> parseIndexList(const string &s) {
    vector<int> result;
    size_t l = s.find('[');
    size_t r = s.find(']');
    if (l==string::npos || r==string::npos || r<=l) return result;
    string inside = s.substr(l+1, r-l-1);
    stringstream ss(inside);
    string token;
    while (getline(ss, token, ',')) {
        token = trim(token);
        if (!token.empty()) result.push_back(stoi(token));
    }
    return result;
}

static NodeDesc parseNodeLine(const string &line) {
    // Node 0: [12.8, 2.8, 9.1]  Neighbor#=9  Neighbors={1, 2, 3, 4, 5, 6, 14, 17, 18}
    NodeDesc nd;
    nd.idx = -1;
    nd.coord = {0.0,0.0,0.0};
    nd.neighbors.clear();
    std::smatch m;
    std::regex r_idx(R"(Node\s+(\d+)\s*:\s*\[)");
    if (std::regex_search(line, m, r_idx)) nd.idx = stoi(m[1].str());
    std::regex r_coords(R"(\[\s*([^\]]+)\s*\])");
    if (std::regex_search(line, m, r_coords)) {
        string inside = m[1].str();
        stringstream ss(inside);
        string tok;
        vector<double> vals;
        while (getline(ss, tok, ',')) vals.push_back(stod(trim(tok)));
        for (int i=0;i<3 && i<(int)vals.size();++i) nd.coord[i]=vals[i];
    }
    std::regex r_n(R"(Neighbors=\{([^\}]*)\})");
    if (std::regex_search(line, m, r_n)) {
        string inside = m[1].str();
        stringstream ss(inside);
        string tok;
        while (getline(ss, tok, ',')) {
            tok = trim(tok);
            if (!tok.empty()) nd.neighbors.push_back(stoi(tok));
        }
    }
    return nd;
}

// parse a multi-round text (one or more "=== Round X ===" blocks)
static vector<RoundData> parseRoundsFromText(const string &text) {
    vector<RoundData> rounds;
    stringstream ss(text);
    string line;
    RoundData cur;
    enum State { SEEK_ROUND, IN_INITIAL, AFTER_DELETE, AFTER_ADD } state = SEEK_ROUND;
    bool started = false;
    while (getline(ss, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (line.rfind("=== Round", 0) == 0) {
            if (started) {
                rounds.push_back(cur);
                cur = RoundData();
            }
            started = true;
            state = IN_INITIAL;
            continue;
        }
        if (line.rfind("Deleted nodes:", 0) == 0) {
            cur.deleted = parseIndexList(line);
            state = AFTER_DELETE;
            continue;
        }
        if (line.rfind("Added nodes:", 0) == 0) {
            cur.added = parseIndexList(line);
            state = AFTER_ADD;
            continue;
        }
        // parse node line
        NodeDesc nd = parseNodeLine(line);
        if (nd.idx < 0) continue;
        if (state == IN_INITIAL) cur.initial.push_back(nd);
        else if (state == AFTER_DELETE) cur.afterDelete.push_back(nd);
        else if (state == AFTER_ADD) cur.finalNodes.push_back(nd);
    }
    if (started) rounds.push_back(cur);
    return rounds;
}

// helper: map index -> L3Address
static L3Address idxToAddr(int idx) {
    std::ostringstream o; o << "10.0.0." << (idx+1);
    return L3Address(Ipv4Address(o.str().c_str()));
}
std::string idxToIp(int idx) {
    return "10.0.0." + std::to_string(idx + 1);
}

// Convert expected neighbor indices -> set<string> of address strings
static unordered_set<string> idxVecToAddrStrSet(const vector<int> &v) {
    unordered_set<string> s;
    for (int i : v) s.insert(idxToAddr(i).str());
    return s;
}

// compute neighbor addr set from result vector<NeighborEntry>
static unordered_set<string> resultToAddrStrSet(const vector<NeighborEntry> &res) {
    unordered_set<string> s;
    for (const auto &ne : res) s.insert(ne.addr.str());
    return s;
}

// Utility: find NodeDesc by idx in a vector (returns nullptr if not found)
static const NodeDesc* findNodeDesc(const vector<NodeDesc> &vec, int idx) {
    for (const auto &nd : vec) if (nd.idx == idx) return &nd;
    return nullptr;
}

// debug helpers - paste near other helpers in test_dt_rounds.cc




// Run the standard round test as you requested:
//  - pick a self (first initial not in deleted list)
//  - pick other (first initial not self)
//  - initial -> validate self and other against initial expectations
//  - perform deletions -> rebuild -> validate self and other using afterDelete expectations
//  - perform additions -> rebuild -> validate self and other using final expectations
static void runRoundTest(const RoundData &rd) {
    // containers used by NeighborDT
    map<L3Address, NeighborEntry> known;
    CGAL::Delaunay_triangulation_3<CGAL::Exact_predicates_inexact_constructions_kernel> dt;
    map<Delaunay3::Vertex_handle, set<L3Address>> vhToAddrs;
    map<L3Address, Delaunay3::Vertex_handle> addrToVh;
    set<L3Address> pending;
    map<L3Address, NeighborEntry> dtNeighbors;
    map<L3Address, NeighborEntry> minSimplex;

    NeighborDT dman(known, dt, vhToAddrs, addrToVh, pending, dtNeighbors, minSimplex);

    // helper to populate known from a list of NodeDesc
    auto populateKnown = [&](const vector<NodeDesc> &nodes) {
        known.clear();
        for (const auto &nd : nodes) {
            NeighborEntry ne;
            ne.addr = idxToAddr(nd.idx);
            ne.coords = { nd.coord[0], nd.coord[1], nd.coord[2] };
            ne.dim = 3;
            ne.attachedToDT = true;
            ne.timeout = SIMTIME_ZERO;
            known[ne.addr] = ne;
        }
    };

    unordered_set<int> deletedSet(rd.deleted.begin(), rd.deleted.end());
    // collect candidates for self
    std::vector<int> selfCandidates;
    for (const auto &nd : rd.initial) {
        if (!deletedSet.count(nd.idx)) {
            selfCandidates.push_back(nd.idx);
        }
    }
    ASSERT_FALSE(selfCandidates.empty()) << "No suitable self found (all initial nodes deleted?)";
    int selfIdx = selfCandidates[rand() % selfCandidates.size()];

    // collect candidates for other
    std::vector<int> otherCandidates;
    for (const auto &nd : rd.initial) {
        if (nd.idx != selfIdx && !deletedSet.count(nd.idx)) {
            otherCandidates.push_back(nd.idx);
        }
    }
    ASSERT_FALSE(otherCandidates.empty()) << "No suitable other found (all initial nodes same as self?)";
    int otherIdx = otherCandidates[rand() % otherCandidates.size()];

    // --- INITIAL stage ---
    populateKnown(rd.initial);
    dman.rebuildDTFromKnownNodes();
    // compute self neighbors
    const NodeDesc* selfNode = findNodeDesc(rd.initial, selfIdx);
    ASSERT_NE(selfNode, nullptr);
    L3Address selfAddr = idxToAddr(selfIdx);
    auto selfResInitial = dman.computeDTNeighborsForCoord(vector<double>{selfNode->coord[0], selfNode->coord[1], selfNode->coord[2]},
                                                          selfAddr, selfAddr);
    auto expectSelfInitial = idxVecToAddrStrSet(selfNode->neighbors);
    auto gotSelfInitial = resultToAddrStrSet(selfResInitial);
    EXPECT_EQ(gotSelfInitial, expectSelfInitial) << "Initial: self node " << selfIdx << " neighbors mismatch";


    // compute other neighbors (if exists) and compare against initial expectation for that node
    if (otherIdx != -1) {
        const NodeDesc* otherNode = findNodeDesc(rd.initial, otherIdx);
        ASSERT_NE(otherNode, nullptr);
        L3Address otherAddr = idxToAddr(otherIdx);
        auto otherResInitial = dman.computeDTNeighborsForCoord(vector<double>{otherNode->coord[0], otherNode->coord[1], otherNode->coord[2]},
                                                               otherAddr, selfAddr);
        auto expectOtherInitial = idxVecToAddrStrSet(otherNode->neighbors);
        auto gotOtherInitial = resultToAddrStrSet(otherResInitial);
        EXPECT_EQ(gotOtherInitial, expectOtherInitial) << "Initial: other node " << otherIdx << " neighbors mismatch";

    }

    // --- AFTER DELETE stage ---
    // remove nodes in rd.deleted: call removeNodeFromDT and erase from known
    for (int di : rd.deleted) {
        L3Address a = idxToAddr(di);
        dman.removeNodeFromDT(a,1);
        known.erase(a);
    }
    // rebuild whole DT to reflect deletions
    dman.rebuildDTFromKnownNodes();

    // validate self using afterDelete list (if self present in afterDelete)
    const NodeDesc* selfAfter = findNodeDesc(rd.afterDelete, selfIdx);
    if (selfAfter != nullptr) {
        auto selfResAfter = dman.computeDTNeighborsForCoord(vector<double>{selfAfter->coord[0], selfAfter->coord[1], selfAfter->coord[2]},
                                                             selfAddr, selfAddr);
        auto expectSelfAfter = idxVecToAddrStrSet(selfAfter->neighbors);
        auto gotSelfAfter = resultToAddrStrSet(selfResAfter);
        EXPECT_EQ(gotSelfAfter, expectSelfAfter) << "After-delete: self " << selfIdx << " neighbors mismatch";

    } else {
        // if self not present in afterDelete, it probably was deleted; assert that
        ASSERT_TRUE(deletedSet.count(selfIdx) == 0) << "Self unexpectedly missing from afterDelete but not in deleted list";
    }

    // other after-delete check (if other still present)
    const NodeDesc* otherAfter = findNodeDesc(rd.afterDelete, otherIdx);
    if (otherAfter != nullptr) {
        L3Address otherAddr = idxToAddr(otherIdx);
        auto otherResAfter = dman.computeDTNeighborsForCoord(vector<double>{otherAfter->coord[0], otherAfter->coord[1], otherAfter->coord[2]},
                                                              otherAddr, selfAddr);
        auto expectOtherAfter = idxVecToAddrStrSet(otherAfter->neighbors);
        auto gotOtherAfter = resultToAddrStrSet(otherResAfter);
        EXPECT_EQ(gotOtherAfter, expectOtherAfter) << "After-delete: other " << otherIdx << " neighbors mismatch";

    }

    // --- FINAL (after additions) stage ---
    // Add nodes from rd.finalNodes that are not present yet in known
    for (const auto &nd : rd.finalNodes) {
        L3Address a = idxToAddr(nd.idx);
        if (known.find(a) == known.end()) {
            NeighborEntry ne;
            ne.addr = a;
            ne.coords = { nd.coord[0], nd.coord[1], nd.coord[2] };
            ne.dim = 3;
            ne.attachedToDT = true;
            ne.timeout = SIMTIME_ZERO;
            known[a] = ne;
            // insert into DT incrementally as well
            if(dman.insertNodeToDT(a, dman.makePoint(ne.coords))){

            }
            else{
                GTEST_LOG_(INFO) <<"insert error\n";
            }
        }
    }
    // full rebuild to compute final neighbor cache
    dman.rebuildDTFromKnownNodes();

    // prepare expect/got for self (compute outside so available for diagnostics)
    std::unordered_set<std::string> expectSelfFinal; // may remain empty if selfFinal==nullptr
    std::unordered_set<std::string> gotSelfFinal;
    const NodeDesc* selfFinal = findNodeDesc(rd.finalNodes, selfIdx);
    if (selfFinal != nullptr) {
        auto selfResFinal = dman.computeDTNeighborsForCoord(
            std::vector<double>{ selfFinal->coord[0], selfFinal->coord[1], selfFinal->coord[2] },
            selfAddr, selfAddr);
        expectSelfFinal = idxVecToAddrStrSet(selfFinal->neighbors);
        gotSelfFinal = resultToAddrStrSet(selfResFinal);
        EXPECT_EQ(gotSelfFinal, expectSelfFinal) << "Final: self " << selfIdx << " neighbors mismatch";
    } else {
        // If there's no final description for self, we can't assert; optionally warn
        GTEST_LOG_(INFO) << "No final node description for selfIdx=" << selfIdx << "\n";
    }

    // other final check
    const NodeDesc* otherFinal = findNodeDesc(rd.finalNodes, otherIdx);
    if (otherFinal != nullptr) {
        L3Address otherAddr = idxToAddr(otherIdx);
        auto otherResFinal = dman.computeDTNeighborsForCoord(vector<double>{otherFinal->coord[0], otherFinal->coord[1], otherFinal->coord[2]},
                                                             otherAddr, selfAddr);
        auto expectOtherFinal = idxVecToAddrStrSet(otherFinal->neighbors);
        auto gotOtherFinal = resultToAddrStrSet(otherResFinal);
        EXPECT_EQ(gotOtherFinal, expectOtherFinal) << "Final: other " << otherIdx << " neighbors mismatch";
    }
}

// -------------------- the three rounds original text (embedded) --------------------
static const string DATA_ROUNDS = R"DATA(
=== Round 1 ===
Node 0: [12.8, 2.8, 9.1]  Neighbor#=9  Neighbors={1, 2, 3, 4, 5, 6, 14, 17, 18}
Node 1: [10.4, 1.5, 8.4]  Neighbor#=10  Neighbors={0, 2, 6, 7, 9, 10, 13, 14, 17, 18}
Node 2: [14.9, 7.3, 2.5]  Neighbor#=11  Neighbors={0, 1, 3, 4, 5, 7, 8, 14, 16, 17, 18}
Node 3: [17.9, 7.6, 9.6]  Neighbor#=10  Neighbors={0, 2, 4, 5, 6, 11, 12, 15, 16, 17}
Node 4: [10.0, 12.9, 11.0]  Neighbor#=12  Neighbors={0, 2, 3, 7, 8, 9, 11, 15, 16, 17, 18, 19}
Node 5: [18.6, 6.6, 5.0]  Neighbor#=9  Neighbors={0, 2, 3, 6, 8, 12, 14, 15, 16}
Node 6: [13.9, 1.6, 11.1]  Neighbor#=12  Neighbors={0, 1, 3, 5, 9, 10, 11, 12, 13, 14, 17, 19}
Node 7: [3.8, 11.8, 5.8]  Neighbor#=10  Neighbors={1, 2, 4, 8, 9, 10, 13, 17, 18, 19}
Node 8: [8.7, 18.0, 0.3]  Neighbor#=10  Neighbors={2, 4, 5, 7, 13, 14, 15, 16, 18, 19}
Node 9: [5.4, 9.1, 13.8]  Neighbor#=9  Neighbors={1, 4, 6, 7, 10, 13, 17, 18, 19}
Node 10: [1.4, 6.1, 11.5]  Neighbor#=8  Neighbors={1, 6, 7, 9, 11, 13, 17, 19}
Node 11: [20.0, 5.2, 19.5]  Neighbor#=8  Neighbors={3, 4, 6, 10, 12, 15, 17, 19}
Node 12: [19.9, 3.7, 14.4]  Neighbor#=7  Neighbors={3, 5, 6, 11, 14, 15, 17}
Node 13: [2.5, 2.6, 6.1]  Neighbor#=7  Neighbors={1, 6, 7, 8, 9, 10, 18}
Node 14: [14.7, 6.1, 0.5]  Neighbor#=9  Neighbors={0, 1, 2, 5, 6, 8, 12, 16, 18}
Node 15: [17.4, 18.0, 14.1]  Neighbor#=9  Neighbors={3, 4, 5, 8, 11, 12, 16, 17, 19}
Node 16: [16.7, 12.2, 7.5]  Neighbor#=8  Neighbors={2, 3, 4, 5, 8, 14, 15, 17}
Node 17: [12.3, 5.6, 12.2]  Neighbor#=15  Neighbors={0, 1, 2, 3, 4, 6, 7, 9, 10, 11, 12, 15, 16, 18, 19}
Node 18: [6.1, 4.5, 2.4]  Neighbor#=10  Neighbors={0, 1, 2, 4, 7, 8, 9, 13, 14, 17}
Node 19: [3.4, 12.0, 18.4]  Neighbor#=9  Neighbors={4, 6, 7, 8, 9, 10, 11, 15, 17}

Deleted nodes: [0, 3, 5, 7, 11, 12, 14, 15, 16, 18]
Node 1: [10.4, 1.5, 8.4]  Neighbor#=7  Neighbors={2, 4, 6, 9, 10, 13, 17}
Node 2: [14.9, 7.3, 2.5]  Neighbor#=6  Neighbors={1, 4, 6, 8, 13, 17}
Node 4: [10.0, 12.9, 11.0]  Neighbor#=8  Neighbors={1, 2, 8, 9, 10, 13, 17, 19}
Node 6: [13.9, 1.6, 11.1]  Neighbor#=7  Neighbors={1, 2, 9, 10, 13, 17, 19}
Node 8: [8.7, 18.0, 0.3]  Neighbor#=6  Neighbors={2, 4, 9, 10, 13, 19}
Node 9: [5.4, 9.1, 13.8]  Neighbor#=8  Neighbors={1, 4, 6, 8, 10, 13, 17, 19}
Node 10: [1.4, 6.1, 11.5]  Neighbor#=8  Neighbors={1, 4, 6, 8, 9, 13, 17, 19}
Node 13: [2.5, 2.6, 6.1]  Neighbor#=7  Neighbors={1, 2, 4, 6, 8, 9, 10}
Node 17: [12.3, 5.6, 12.2]  Neighbor#=7  Neighbors={1, 2, 4, 6, 9, 10, 19}
Node 19: [3.4, 12.0, 18.4]  Neighbor#=6  Neighbors={4, 6, 8, 9, 10, 17}

Added nodes: [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
Node 1: [10.4, 1.5, 8.4]  Neighbor#=11  Neighbors={2, 6, 9, 10, 17, 20, 21, 22, 23, 25, 28}
Node 2: [14.9, 7.3, 2.5]  Neighbor#=12  Neighbors={1, 4, 8, 13, 21, 22, 23, 24, 25, 26, 28, 29}
Node 4: [10.0, 12.9, 11.0]  Neighbor#=9  Neighbors={2, 9, 17, 19, 21, 24, 25, 26, 29}
Node 6: [13.9, 1.6, 11.1]  Neighbor#=10  Neighbors={1, 9, 10, 17, 19, 20, 21, 22, 27, 28}
Node 8: [8.7, 18.0, 0.3]  Neighbor#=8  Neighbors={2, 13, 19, 23, 24, 25, 26, 27}
Node 9: [5.4, 9.1, 13.8]  Neighbor#=10  Neighbors={1, 4, 6, 10, 13, 17, 19, 21, 22, 25}
Node 10: [1.4, 6.1, 11.5]  Neighbor#=9  Neighbors={1, 6, 9, 13, 17, 19, 21, 22, 25}
Node 13: [2.5, 2.6, 6.1]  Neighbor#=6  Neighbors={2, 8, 9, 10, 22, 25}
Node 17: [12.3, 5.6, 12.2]  Neighbor#=10  Neighbors={1, 4, 6, 9, 10, 19, 20, 21, 27, 29}
Node 19: [3.4, 12.0, 18.4]  Neighbor#=11  Neighbors={4, 6, 8, 9, 10, 17, 20, 24, 25, 26, 29}
Node 20: [16.6, 0.3, 10.3]  Neighbor#=9  Neighbors={1, 6, 17, 19, 21, 22, 27, 28, 29}
Node 21: [12.1, 4.2, 10.2]  Neighbor#=14  Neighbors={1, 2, 4, 6, 9, 10, 17, 20, 22, 23, 25, 27, 28, 29}
Node 22: [4.9, 1.1, 6.3]  Neighbor#=11  Neighbors={1, 2, 6, 9, 10, 13, 20, 21, 23, 25, 28}
Node 23: [15.7, 6.3, 2.8]  Neighbor#=9  Neighbors={1, 2, 8, 21, 22, 26, 27, 28, 29}
Node 24: [11.2, 18.6, 3.2]  Neighbor#=7  Neighbors={2, 4, 8, 19, 25, 26, 29}
Node 25: [2.6, 9.9, 4.0]  Neighbor#=11  Neighbors={1, 2, 4, 8, 9, 10, 13, 19, 21, 22, 24}
Node 26: [17.7, 18.4, 7.5]  Neighbor#=8  Neighbors={2, 4, 8, 19, 23, 24, 27, 29}
Node 27: [17.6, 4.8, 3.8]  Neighbor#=9  Neighbors={6, 8, 17, 20, 21, 23, 26, 28, 29}
Node 28: [12.4, 0.1, 5.8]  Neighbor#=8  Neighbors={1, 2, 6, 20, 21, 22, 23, 27}
Node 29: [17.5, 13.6, 9.5]  Neighbor#=10  Neighbors={2, 4, 17, 19, 20, 21, 23, 24, 26, 27}

=== Round 2 ===
Node 0: [6.0, 7.3, 9.7]  Neighbor#=14  Neighbors={1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, 16, 18, 19}
Node 1: [0.4, 6.6, 13.1]  Neighbor#=10  Neighbors={0, 2, 3, 5, 8, 11, 13, 16, 17, 18}
Node 2: [3.8, 0.6, 13.7]  Neighbor#=8  Neighbors={0, 1, 3, 7, 8, 9, 10, 11}
Node 3: [7.9, 3.2, 18.2]  Neighbor#=8  Neighbors={0, 1, 2, 9, 10, 11, 18, 19}
Node 4: [19.5, 19.6, 15.7]  Neighbor#=9  Neighbors={5, 6, 7, 9, 10, 12, 14, 16, 17}
Node 5: [1.1, 18.9, 16.2]  Neighbor#=8  Neighbors={1, 4, 6, 11, 14, 16, 17, 18}
Node 6: [7.7, 16.0, 0.7]  Neighbor#=9  Neighbors={0, 4, 5, 7, 8, 13, 14, 16, 17}
Node 7: [14.9, 4.6, 5.0]  Neighbor#=10  Neighbors={0, 2, 4, 6, 8, 9, 10, 12, 13, 14}
Node 8: [2.2, 1.3, 11.2]  Neighbor#=9  Neighbors={0, 1, 2, 6, 7, 10, 11, 13, 17}
Node 9: [17.2, 4.5, 11.5]  Neighbor#=10  Neighbors={0, 2, 3, 4, 7, 10, 12, 14, 15, 19}
Node 10: [17.5, 1.4, 15.1]  Neighbor#=9  Neighbors={2, 3, 4, 7, 8, 9, 12, 15, 19}
Node 11: [0.8, 4.6, 17.1]  Neighbor#=7  Neighbors={1, 2, 3, 5, 8, 16, 18}
Node 12: [16.3, 10.3, 17.7]  Neighbor#=9  Neighbors={0, 4, 7, 9, 10, 14, 15, 18, 19}
Node 13: [4.1, 10.4, 6.7]  Neighbor#=9  Neighbors={0, 1, 6, 7, 8, 14, 16, 17, 18}
Node 14: [16.4, 14.8, 18.0]  Neighbor#=11  Neighbors={0, 4, 5, 6, 7, 9, 12, 13, 16, 18, 19}
Node 15: [15.5, 7.2, 16.7]  Neighbor#=6  Neighbors={0, 9, 10, 12, 18, 19}
Node 16: [1.9, 15.0, 15.3]  Neighbor#=10  Neighbors={0, 1, 4, 5, 6, 11, 13, 14, 17, 18}
Node 17: [0.3, 16.6, 3.1]  Neighbor#=7  Neighbors={1, 4, 5, 6, 8, 13, 16}
Node 18: [7.4, 11.3, 18.8]  Neighbor#=11  Neighbors={0, 1, 3, 5, 11, 12, 13, 14, 15, 16, 19}
Node 19: [11.9, 4.8, 18.6]  Neighbor#=8  Neighbors={0, 3, 9, 10, 12, 14, 15, 18}

Deleted nodes: [0, 3, 4, 7, 9, 10, 11, 16, 17, 18]
Node 1: [0.4, 6.6, 13.1]  Neighbor#=6  Neighbors={2, 5, 6, 8, 13, 19}
Node 2: [3.8, 0.6, 13.7]  Neighbor#=6  Neighbors={1, 5, 8, 13, 15, 19}
Node 5: [1.1, 18.9, 16.2]  Neighbor#=7  Neighbors={1, 2, 6, 12, 13, 14, 19}
Node 6: [7.7, 16.0, 0.7]  Neighbor#=7  Neighbors={1, 5, 8, 12, 13, 14, 15}
Node 8: [2.2, 1.3, 11.2]  Neighbor#=5  Neighbors={1, 2, 6, 13, 15}
Node 12: [16.3, 10.3, 17.7]  Neighbor#=6  Neighbors={5, 6, 13, 14, 15, 19}
Node 13: [4.1, 10.4, 6.7]  Neighbor#=9  Neighbors={1, 2, 5, 6, 8, 12, 14, 15, 19}
Node 14: [16.4, 14.8, 18.0]  Neighbor#=5  Neighbors={5, 6, 12, 13, 19}
Node 15: [15.5, 7.2, 16.7]  Neighbor#=6  Neighbors={2, 6, 8, 12, 13, 19}
Node 19: [11.9, 4.8, 18.6]  Neighbor#=7  Neighbors={1, 2, 5, 12, 13, 14, 15}

Added nodes: [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
Node 1: [0.4, 6.6, 13.1]  Neighbor#=10  Neighbors={2, 5, 8, 13, 20, 21, 23, 26, 27, 28}
Node 2: [3.8, 0.6, 13.7]  Neighbor#=10  Neighbors={1, 8, 19, 20, 21, 22, 25, 26, 28, 29}
Node 5: [1.1, 18.9, 16.2]  Neighbor#=8  Neighbors={1, 13, 14, 21, 23, 24, 26, 27}
Node 6: [7.7, 16.0, 0.7]  Neighbor#=6  Neighbors={8, 13, 22, 23, 24, 25}
Node 8: [2.2, 1.3, 11.2]  Neighbor#=9  Neighbors={1, 2, 6, 13, 20, 22, 23, 25, 26}
Node 12: [16.3, 10.3, 17.7]  Neighbor#=8  Neighbors={14, 15, 19, 21, 22, 24, 28, 29}
Node 13: [4.1, 10.4, 6.7]  Neighbor#=14  Neighbors={1, 5, 6, 8, 14, 20, 21, 22, 23, 24, 25, 27, 28, 29}
Node 14: [16.4, 14.8, 18.0]  Neighbor#=9  Neighbors={5, 12, 13, 19, 21, 23, 24, 27, 28}
Node 15: [15.5, 7.2, 16.7]  Neighbor#=6  Neighbors={12, 19, 21, 22, 24, 29}
Node 19: [11.9, 4.8, 18.6]  Neighbor#=9  Neighbors={2, 12, 14, 15, 21, 22, 26, 28, 29}
Node 20: [6.1, 3.5, 10.6]  Neighbor#=7  Neighbors={1, 2, 8, 13, 25, 28, 29}
Node 21: [4.6, 12.3, 18.2]  Neighbor#=13  Neighbors={1, 2, 5, 12, 13, 14, 15, 19, 23, 26, 27, 28, 29}
Node 22: [11.8, 0.8, 3.4]  Neighbor#=12  Neighbors={2, 6, 8, 12, 13, 15, 19, 24, 25, 26, 28, 29}
Node 23: [3.7, 17.1, 5.1]  Neighbor#=9  Neighbors={1, 5, 6, 8, 13, 14, 21, 24, 27}
Node 24: [14.6, 16.6, 3.6]  Neighbor#=11  Neighbors={5, 6, 12, 13, 14, 15, 22, 23, 27, 28, 29}
Node 25: [6.9, 1.5, 6.6]  Neighbor#=8  Neighbors={2, 6, 8, 13, 20, 22, 28, 29}
Node 26: [1.1, 0.3, 16.4]  Neighbor#=8  Neighbors={1, 2, 5, 8, 19, 21, 22, 27}
Node 27: [1.9, 17.4, 16.6]  Neighbor#=8  Neighbors={1, 5, 13, 14, 21, 23, 24, 26}
Node 28: [6.6, 6.1, 10.8]  Neighbor#=12  Neighbors={1, 2, 12, 13, 14, 19, 20, 21, 22, 24, 25, 29}
Node 29: [10.7, 4.5, 12.0]  Neighbor#=11  Neighbors={2, 12, 13, 15, 19, 20, 21, 22, 24, 25, 28}

=== Round 3 ===
Node 0: [6.8, 1.9, 3.9]  Neighbor#=8  Neighbors={1, 2, 3, 5, 6, 7, 8, 9}
Node 1: [2.0, 0.5, 0.5]  Neighbor#=6  Neighbors={0, 3, 4, 5, 6, 9}
Node 2: [7.6, 6.5, 5.0]  Neighbor#=6  Neighbors={0, 3, 6, 7, 8, 9}
Node 3: [4.4, 0.9, 4.1]  Neighbor#=8  Neighbors={0, 1, 2, 4, 5, 6, 7, 8}
Node 4: [0.6, 7.6, 6.4]  Neighbor#=7  Neighbors={1, 3, 5, 6, 7, 8, 9}
Node 5: [1.7, 2.3, 3.5]  Neighbor#=7  Neighbors={0, 1, 3, 4, 6, 7, 9}
Node 6: [3.2, 2.1, 0.8]  Neighbor#=8  Neighbors={0, 1, 2, 3, 4, 5, 7, 9}
Node 7: [5.1, 7.2, 4.8]  Neighbor#=8  Neighbors={0, 2, 3, 4, 5, 6, 8, 9}
Node 8: [3.9, 9.4, 8.1]  Neighbor#=6  Neighbors={0, 2, 3, 4, 7, 9}
Node 9: [4.7, 9.8, 2.1]  Neighbor#=8  Neighbors={0, 1, 2, 4, 5, 6, 7, 8}

Deleted nodes: [1, 4, 5, 8, 9]
Node 0: [6.8, 1.9, 3.9]  Neighbor#=4  Neighbors={2, 3, 6, 7}
Node 2: [7.6, 6.5, 5.0]  Neighbor#=4  Neighbors={0, 3, 6, 7}
Node 3: [4.4, 0.9, 4.1]  Neighbor#=4  Neighbors={0, 2, 6, 7}
Node 6: [3.2, 2.1, 0.8]  Neighbor#=4  Neighbors={0, 2, 3, 7}
Node 7: [5.1, 7.2, 4.8]  Neighbor#=4  Neighbors={0, 2, 3, 6}

Added nodes: [10, 11, 12, 13, 14]
Node 0: [6.8, 1.9, 3.9]  Neighbor#=6  Neighbors={2, 3, 6, 11, 13, 14}
Node 2: [7.6, 6.5, 5.0]  Neighbor#=6  Neighbors={0, 7, 10, 11, 13, 14}
Node 3: [4.4, 0.9, 4.1]  Neighbor#=8  Neighbors={0, 6, 7, 10, 11, 12, 13, 14}
Node 6: [3.2, 2.1, 0.8]  Neighbor#=7  Neighbors={0, 3, 7, 10, 12, 13, 14}
Node 7: [5.1, 7.2, 4.8]  Neighbor#=7  Neighbors={2, 3, 6, 10, 11, 12, 14}
Node 10: [0.9, 6.3, 0.1]  Neighbor#=8  Neighbors={2, 3, 6, 7, 11, 12, 13, 14}
Node 11: [1.5, 8.3, 8.9]  Neighbor#=8  Neighbors={0, 2, 3, 7, 10, 12, 13, 14}
Node 12: [2.9, 7.4, 4.2]  Neighbor#=6  Neighbors={3, 6, 7, 10, 11, 14}
Node 13: [8.7, 2.7, 3.3]  Neighbor#=7  Neighbors={0, 2, 3, 6, 10, 11, 14}
Node 14: [6.9, 4.1, 4.5]  Neighbor#=9  Neighbors={0, 2, 3, 6, 7, 10, 11, 12, 13}
)DATA";

// instantiate and run parser, then create tests per round
static vector<RoundData> ALL_ROUNDS = parseRoundsFromText(DATA_ROUNDS);

// convenience: create a single NeighborDT for each test run with fresh containers (done in runRoundTest)

// Test cases: one test per round
TEST(DTRounds, Round1) {
    ASSERT_GE((int)ALL_ROUNDS.size(), 1);
    runRoundTest(ALL_ROUNDS[0]);
}

TEST(DTRounds, Round2) {
    ASSERT_GE((int)ALL_ROUNDS.size(), 2);
    runRoundTest(ALL_ROUNDS[1]);
}

TEST(DTRounds, Round3) {
    ASSERT_GE((int)ALL_ROUNDS.size(), 3);
    runRoundTest(ALL_ROUNDS[2]);
}






