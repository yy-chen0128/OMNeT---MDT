/*
 * NeighborDT.cc
 *
 *  Created on: Sep 26, 2025
 *      Author: yychen
 */

#include "NeighborDT.h"
#include <algorithm>
#include <iterator>
#include <unordered_map>

namespace mysrc {
namespace routing {

NeighborDT::NeighborDT(std::map<L3Address, NeighborEntry> &knownNodeCoords,
                       Delaunay3 &dt,
                       std::map<Delaunay3::Vertex_handle, std::set<L3Address>> &vhToAddrs,
                       std::map<L3Address, Delaunay3::Vertex_handle> &addrToVh,
                       std::set<L3Address> &pendingRebuildAddrs,
                       std::map<L3Address, NeighborEntry> &dtNeighbors,
                       std::map<L3Address, NeighborEntry> &minSimplexNodes)
    : knownNodeCoordsRef_(knownNodeCoords),
      dtRef_(dt),
      vhToAddrsRef_(vhToAddrs),
      addrToVhRef_(addrToVh),
      pendingRebuildRef_(pendingRebuildAddrs),
      dtNeighborsRef_(dtNeighbors),
      minSimplexNodesRef_(minSimplexNodes)
{
}

static inline void point3_to_vec(const Point3 &p, double &x, double &y, double &z) {
    x = static_cast<double>(p.x());
    y = static_cast<double>(p.y());
    z = static_cast<double>(p.z());
}

static inline void vec_sub(const double a[3], const double b[3], double out[3]) {
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}
static inline double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static inline void cross3(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}
static inline double norm3(const double a[3]) {
    return std::sqrt(dot3(a,a));
}
static inline void normalize3(double a[3]) {
    double n = norm3(a);
    if (n > 0.0) { a[0]/=n; a[1]/=n; a[2]/=n; }
}

Point3 NeighborDT::makePoint(const std::vector<double> &c) const {
    double x = c.size() > 0 ? c[0] : 0.0;
    double y = c.size() > 1 ? c[1] : 0.0;
    double z = c.size() > 2 ? c[2] : 0.0;
    return Point3(x, y, z);
}

bool NeighborDT::insertNodeToDT(const L3Address &addr, const Point3 &p) {
    // If addr already exists and maps to a vertex, ensure it's up-to-date (but don't move)
    auto itAddr = addrToVhRef_.find(addr);
    if (itAddr != addrToVhRef_.end()) {
        Vertex_handle existingVh = itAddr->second;
        if (existingVh != Vertex_handle() && existingVh->point() == p)
            return false;
        // otherwise caller can trigger a rebuild if move required (pure class won't call move)
    }

    // Try to find an existing vertex at same coordinates
    Vertex_handle nearest = dtRef_.nearest_vertex(p);
    if (nearest != Vertex_handle() && nearest->point() == p) {
        // reuse existing vertex
        addrToVhRef_[addr] = nearest;
        vhToAddrsRef_[nearest].insert(addr);
        return true;
    }

    // Insert new vertex
    Vertex_handle vh = dtRef_.insert(p);
    addrToVhRef_[addr] = vh;
    vhToAddrsRef_[vh].insert(addr);
    return true;
}



void NeighborDT::removeNodeFromDT(const L3Address &addr, int pendingRebuildThreshold) {
    auto itAddr = addrToVhRef_.find(addr);
    if (itAddr == addrToVhRef_.end()) {
        return;
    }
    Vertex_handle vh = itAddr->second;

    auto vhIt = vhToAddrsRef_.find(vh);
    if (vhIt == vhToAddrsRef_.end()) {
        // inconsistent: just erase addr->vh mapping
        addrToVhRef_.erase(itAddr);
        return;
    }

    size_t erased = vhIt->second.erase(addr);
    addrToVhRef_.erase(itAddr);

    if (erased == 0) {
        return;
    }

    if (!vhIt->second.empty()) {
        return;
    }

    // empty: erase mapping, and mark for rebuild (do not call dt.remove here)
    vhToAddrsRef_.erase(vhIt);
    pendingRebuildRef_.insert(addr);

    // Optionally trigger immediate rebuild if threshold reached (no timing in this class)
    if ((int)pendingRebuildRef_.size() >= pendingRebuildThreshold) {
        // immediate rebuild using default includePred (groups caller must decide valid entries)
        rebuildDTFromKnownNodes();
    }
}

void NeighborDT::rebuildDTFromKnownNodes(std::function<bool(const NeighborEntry&)> includePred) {
    // full rebuild
    dtRef_.clear();
    addrToVhRef_.clear();
    vhToAddrsRef_.clear();
    dtNeighborsRef_.clear();
    minSimplexNodesRef_.clear();

    // iterate knownNodeCoordsRef_ and insert those passing includePred (or all if null)
    for (const auto &kv : knownNodeCoordsRef_) {
        const L3Address &addr = kv.first;
        const NeighborEntry &ne = kv.second;
        if (includePred && !includePred(ne)) continue;

        Point3 p = makePoint(ne.coords);
        insertNodeToDT(addr, p);
    }

    // caller can rebuild neighbor cache by calling rebuildDTNeighborCache separately,
    // but doing it here is convenient:
    // If caller wants a self-based cache, call rebuildDTNeighborCache(selfAddr,...)
    pendingRebuildRef_.clear();
}

void NeighborDT::rebuildDTNeighborCache(const L3Address &selfAddr,
                                        std::function<bool(const NeighborEntry&)> includePred) {
    dtNeighborsRef_.clear();

    auto itSelf = addrToVhRef_.find(selfAddr);
    if (itSelf == addrToVhRef_.end()) {
        EV_ERROR<<"self is not at dt can not rebuild dt\n";
        return;
    }
    Vertex_handle selfVh = itSelf->second;
    if (dtRef_.is_infinite(selfVh)) return;

    /*if (dtRef_.dimension() != 3) {
        // fallback: include all attached & included known nodes except self
        for (const auto &kv : knownNodeCoordsRef_) {
            if (kv.first == selfAddr) continue;
            if (includePred && !includePred(kv.second)) continue;
            dtNeighborsRef_[kv.first] = kv.second;
        }
        return;
    }*/
    if (dtRef_.dimension() != 3) {
        int dim = dtRef_.dimension();
        if (dim <= 1) {
            // 1D / 0D: keep conservative fallback (all attached & included nodes except self)
            EV_ERROR << "rebuildDTNeighborCache: dt.dimension()=" << dim << " (1D/0D fallback)\n";
            simtime_t now = simTime();
            for (const auto &kv : knownNodeCoordsRef_) {
                if (kv.first == selfAddr) continue;
                if (includePred && !includePred(kv.second)) continue;
                if (!kv.second.attachedToDT) continue;
                dtNeighborsRef_[kv.first] = kv.second;
            }
            EV_ERROR << "rebuildDTNeighborCache (1D fallback): dtNeighbors.size()=" << dtNeighborsRef_.size() << "\n";
            return;
        }

        // dim == 2 -> attempt 2D projection + 2D Delaunay
        EV_ERROR << "rebuildDTNeighborCache: dt.dimension()=2 -> using 2D projection Delaunay\n";

        // pick three non-collinear reference points among known to define plane
        std::vector<Point3> refs;
        std::vector<inet::L3Address> refAddrs;
        for (const auto &kv : knownNodeCoordsRef_) {
            if (includePred && !includePred(kv.second)) continue;
            if (!kv.second.attachedToDT) continue;
            if (kv.second.coords.size() < 3) continue;
            refs.emplace_back(makePoint(kv.second.coords)); // makePoint available
            refAddrs.push_back(kv.first);
            if (refs.size() >= 3) break;
        }
        if (refs.size() < 3) {
            // cannot form plane reliably -> fallback to conservative
            EV_ERROR << "rebuildDTNeighborCache: not enough 3D refs for plane -> fallback to conservative\n";
            for (const auto &kv : knownNodeCoordsRef_) {
                if (kv.first == selfAddr) continue;
                if (includePred && !includePred(kv.second)) continue;
                if (!kv.second.attachedToDT) continue;
                dtNeighborsRef_[kv.first] = kv.second;
            }
            return;
        }

        // convert refs to double arrays
        double p0v[3], p1v[3], p2v[3];
        point3_to_vec(refs[0], p0v[0], p0v[1], p0v[2]);
        point3_to_vec(refs[1], p1v[0], p1v[1], p1v[2]);
        point3_to_vec(refs[2], p2v[0], p2v[1], p2v[2]);
        double u[3], tmp[3], n[3], v[3];
        vec_sub(p1v, p0v, u); normalize3(u);
        vec_sub(p2v, p0v, tmp);
        cross3(u, tmp, n);
        if (norm3(n) < 1e-12) {
            // nearly collinear -> fallback to conservative
            for (const auto &kv : knownNodeCoordsRef_) {
                if (kv.first == selfAddr) continue;
                if (includePred && !includePred(kv.second)) continue;
                if (!kv.second.attachedToDT) continue;
                dtNeighborsRef_[kv.first] = kv.second;
            }
            return;
        }
        normalize3(n);
        cross3(n, u, v); normalize3(v);

        EV_ERROR<<"Build 2D Delaunay\n";
        Delaunay2 dt2;
        std::map<inet::L3Address, Vertex_handle2> addr2vh;
        std::map<Vertex_handle2, std::set<inet::L3Address>> vh2addrs;
        for (const auto &kv : knownNodeCoordsRef_) {
            const auto &la = kv.first;
            const auto &ne = kv.second;
            if (!ne.attachedToDT) continue;
            if (includePred && !includePred(ne)) continue;
            // need at least 3 coords to project; if less, we still can project using zeros for missing dims
            Point3 p3 = makePoint(ne.coords);
            double pv[3]; point3_to_vec(p3, pv[0], pv[1], pv[2]);
            double w[3]; vec_sub(pv, p0v, w);
            double x = dot3(w, u);
            double y = dot3(w, v);
            Vertex_handle2 vh2 = dt2.insert(Point2(x, y));
            addr2vh[la] = vh2;
            vh2addrs[vh2].insert(la);
        }

        EV_ERROR<<"if selfAddr not present, nothing to do\n";
        auto itSelf2 = addr2vh.find(selfAddr);
        if (itSelf2 == addr2vh.end()) {
            // no self in 2D projection -> fallback conservative
            for (const auto &kv : knownNodeCoordsRef_) {
                if (kv.first == selfAddr) continue;
                if (includePred && !includePred(kv.second)) continue;
                if (!kv.second.attachedToDT) continue;
                dtNeighborsRef_[kv.first] = kv.second;
            }
            return;
        }

        EV_ERROR<<"collect incident vertices in 2D using a circulator (portable across CGAL versions)\n";
        {
            EV_ERROR << "rebuildDTNeighborCache: collecting 2D incident vertices for self\n";

            // use CGAL circulator for incident vertices
            Delaunay2::Vertex_circulator vcirc = dt2.incident_vertices(itSelf2->second);
            if (vcirc == Delaunay2::Vertex_circulator()) {
                // no neighbors (isolated vertex) -> nothing to add
                EV_ERROR << "rebuildDTNeighborCache (2D): no incident vertices found\n";
            } else {
                Delaunay2::Vertex_circulator done = vcirc;
                do {
                    // convert circulator position to a vertex handle
                    Vertex_handle2 vh2 = vcirc;
                    ++vcirc; // advance early to keep consistent semantics if we continue/break

                    if (dt2.is_infinite(vh2)) continue;

                    auto it = vh2addrs.find(vh2);
                    if (it == vh2addrs.end()) continue;

                    EV_ERROR<<"expand each associated address and add to dtNeighborsRef_ if it passes includePred\n";
                    for (const L3Address &a : it->second) {
                        if (a == selfAddr) continue;
                        auto kn = knownNodeCoordsRef_.find(a);
                        if (kn == knownNodeCoordsRef_.end()) continue;
                        if (includePred && !includePred(kn->second)) continue;
                        dtNeighborsRef_[a] = kn->second;
                    }
                } while (vcirc != done);
            }
        }

        EV_ERROR << "rebuildDTNeighborCache (2D): dtNeighbors.size()=" << dtNeighborsRef_.size() << "\n";
        return;
    }

    std::vector<Delaunay3::Cell_handle> cells;
    dtRef_.incident_cells(selfVh, std::back_inserter(cells));

    std::vector<Vertex_handle> neighborVhs;
    neighborVhs.reserve(cells.size() * 3);
    for (const auto &ch : cells) {
        if (dtRef_.is_infinite(ch)) continue;
        for (int i = 0; i < 4; ++i) {
            Vertex_handle vh = ch->vertex(i);
            if (vh == selfVh || dtRef_.is_infinite(vh)) continue;
            bool already = false;
            for (const auto &existing : neighborVhs) {
                if (existing == vh) { already = true; break; }
            }
            if (!already) neighborVhs.push_back(vh);
        }
    }

    for (const Vertex_handle &vh : neighborVhs) {
        auto vhIt = vhToAddrsRef_.find(vh);
        if (vhIt == vhToAddrsRef_.end()) continue;
        for (const L3Address &nbAddr : vhIt->second) {
            if (nbAddr == selfAddr) continue;
            auto kn = knownNodeCoordsRef_.find(nbAddr);
            if (kn == knownNodeCoordsRef_.end()) continue;
            const NeighborEntry &ne = kn->second;
            if (includePred && !includePred(ne)) continue;
            dtNeighborsRef_[nbAddr] = ne;
        }
    }
}

std::vector<NeighborEntry> NeighborDT::computeDTNeighborsForCoord(const std::vector<double> &targetCoord,
                                                                   const L3Address &addr,
                                                                   const L3Address &selfAddr,
                                                                   std::function<bool(const NeighborEntry&)> includePred) {
    std::vector<NeighborEntry> result;

    // ensure self present in addrToVhRef_ is caller responsibility; if not present we still try to insert
    // but since this class cannot call owner getSelfCoordinate, caller should pre-insert if needed.

    auto itAddrVH = addrToVhRef_.find(addr);
    if (itAddrVH == addrToVhRef_.end()) {
        // attempt to create from knownNodeCoordsRef_ or targetCoord
        auto kn = knownNodeCoordsRef_.find(addr);
        if (kn == knownNodeCoordsRef_.end()) {
            // create temporary entry in knownNodeCoordsRef_ so that vertex mapping persists
            NeighborEntry ne;
            ne.addr = addr;
            ne.coords = targetCoord;
            ne.dim = (int)targetCoord.size();
            ne.attachedToDT = true;
            ne.timeout = simTime() + 10; // meaningless here
            knownNodeCoordsRef_[addr] = ne;
            insertNodeToDT(addr, makePoint(targetCoord));
        } else {
            // use stored coords (if empty use targetCoord)
            const std::vector<double> coords = kn->second.coords.empty() ? targetCoord : kn->second.coords;
            knownNodeCoordsRef_[addr].attachedToDT = true;
            insertNodeToDT(addr, makePoint(coords));
        }
    }

    auto itTargetVH = addrToVhRef_.find(addr);
    if (itTargetVH == addrToVhRef_.end()) {
        // fallback conservative: return all included knownNodeCoords
        EV_ERROR<<"computeDTNeighborsForCoord:can not insert the target node\n";
        for (const auto &kv : knownNodeCoordsRef_) {
            //if (kv.first == selfAddr) continue;
            if (includePred && !includePred(kv.second)) continue;
            if (!kv.second.attachedToDT) continue;
            result.push_back(kv.second);
        }
        return result;
    }

    Vertex_handle tvh = itTargetVH->second;

    /*if (dtRef_.dimension() != 3) {
        for (const auto &kv : knownNodeCoordsRef_) {
            //if (kv.first == selfAddr) continue;
            if (includePred && !includePred(kv.second)) continue;
            if (!kv.second.attachedToDT) continue;
            result.push_back(kv.second);
        }
        return result;
    }*/
    if (dtRef_.dimension() == 3){
        std::vector<Delaunay3::Cell_handle> cells;
        dtRef_.incident_cells(tvh, std::back_inserter(cells));

        std::unordered_set<Delaunay3::Vertex_handle> neighVhs;
        for (const auto &ch : cells) {
            if (dtRef_.is_infinite(ch)) continue;
            for (int i = 0; i < 4; ++i) {
                auto vh = ch->vertex(i);
                if (vh != tvh && !dtRef_.is_infinite(vh)) neighVhs.insert(vh);
            }
        }

        std::unordered_set<std::string> seen;
        for ( auto &vh : neighVhs) {
            auto vhIt = vhToAddrsRef_.find(vh);
            if (vhIt == vhToAddrsRef_.end()) continue;
            for (const L3Address &nbAddr : vhIt->second) {
                //if (nbAddr == selfAddr) continue;
                auto kn = knownNodeCoordsRef_.find(nbAddr);
                if (kn == knownNodeCoordsRef_.end()) continue;
                NeighborEntry &ne = kn->second;
                if (includePred && !includePred(ne)) continue;
                if (!ne.attachedToDT && nbAddr != addr) continue;
                std::string key = nbAddr.str();
                if (nbAddr == selfAddr) {
                    ne.timeout == 0 ? ne.timeout += 10:ne.timeout += 1;
                }
                if (seen.insert(key).second) result.push_back(ne);
            }
        }

        return result;
    }

    // --- Non-3D path: differentiate between 2D and 1D/0D ---
    int dim = dtRef_.dimension();
    EV_INFO << "computeDTNeighborsForCoord: dt.dimension()=" << dim << " (non-3D path)\n";

    if (dim <= 1) {
        // conservative fallback: include all attached & included known nodes (except self)
        for (const auto &kv : knownNodeCoordsRef_) {
            if (kv.first == selfAddr) continue;
            if (includePred && !includePred(kv.second)) continue;
            if (!kv.second.attachedToDT) continue;
            result.push_back(kv.second);
        }
        EV_INFO << "computeDTNeighborsForCoord (1D/0D fallback): result.size=" << result.size() << "\n";
        return result;
    }

    // --- dim == 2: do 2D projection + Delaunay_triangulation_2 on-the-fly and use it ---
    // We will:
    //  - pick three non-collinear reference points among knownNodeCoordsRef_
    //  - construct plane basis u,v using those refs
    //  - project all included known nodes into 2D coords and build Delaunay2
    //  - ensure target 'addr' is present in 2D triangulation (insert from known or from targetCoord)
    //  - collect incident vertices around the target 2D vertex and expand to addresses

    EV_INFO << "computeDTNeighborsForCoord: entering 2D projection path\n";

    // typedefs for 2D
    typedef CGAL::Delaunay_triangulation_2<K> Delaunay2;
    typedef K::Point_2 Point2;
    typedef Delaunay2::Vertex_handle Vertex_handle2;

    // gather reference points (must be attached & pass includePred)
    std::vector<Point3> refs;
    std::vector<L3Address> refAddrs;
    for (const auto &kv : knownNodeCoordsRef_) {
        const auto &ne = kv.second;
        if (!ne.attachedToDT) continue;
        if (includePred && !includePred(ne)) continue;
        if (ne.coords.size() < 3) continue;
        refs.push_back(makePoint(ne.coords));
        refAddrs.push_back(kv.first);
        if (refs.size() >= 3) break;
    }
    if (refs.size() < 3) {
        EV_WARN << "computeDTNeighborsForCoord: not enough 3 refs for plane -> fallback conservative\n";
        for (const auto &kv : knownNodeCoordsRef_) {
            if (kv.first == selfAddr) continue;
            if (includePred && !includePred(kv.second)) continue;
            if (!kv.second.attachedToDT) continue;
            result.push_back(kv.second);
        }
        return result;
    }

    // convert three refs to double vectors
    auto to_xyz = [&](const Point3 &p, double out[3]) {
        out[0] = static_cast<double>(p.x());
        out[1] = static_cast<double>(p.y());
        out[2] = static_cast<double>(p.z());
    };
    double p0[3], p1[3], p2[3];
    to_xyz(refs[0], p0);
    to_xyz(refs[1], p1);
    to_xyz(refs[2], p2);

    auto sub = [](const double a[3], const double b[3], double out[3]) {
        out[0] = a[0] - b[0];
        out[1] = a[1] - b[1];
        out[2] = a[2] - b[2];
    };
    auto cross = [](const double a[3], const double b[3], double out[3]) {
        out[0] = a[1]*b[2] - a[2]*b[1];
        out[1] = a[2]*b[0] - a[0]*b[2];
        out[2] = a[0]*b[1] - a[1]*b[0];
    };
    auto dot = [](const double a[3], const double b[3]) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    auto norm = [&](const double a[3]) {
        return std::sqrt(dot(a,a));
    };
    auto normalize = [&](double a[3]) {
        double n = norm(a);
        if (n > 0.0) { a[0]/=n; a[1]/=n; a[2]/=n; }
    };

    double u[3], tmpv[3], nvec[3], v[3];
    sub(p1, p0, u); normalize(u);
    sub(p2, p0, tmpv);
    cross(u, tmpv, nvec);
    double nn = norm(nvec);
    if (nn < 1e-12) {
        EV_WARN << "computeDTNeighborsForCoord: refs nearly collinear -> fallback\n";
        for (const auto &kv : knownNodeCoordsRef_) {
            if (kv.first == selfAddr) continue;
            if (includePred && !includePred(kv.second)) continue;
            if (!kv.second.attachedToDT) continue;
            result.push_back(kv.second);
        }
        return result;
    }
    normalize(nvec);
    cross(nvec, u, v); normalize(v);

    // Build 2D Delaunay and mapping structures
    Delaunay2 dt2;
    std::map<L3Address, Vertex_handle2> addr2vh;
    std::map<Vertex_handle2, std::set<L3Address>> vh2addrs;

    auto project_to_2d = [&](const std::vector<double> &coords, double &x, double &y){
        Point3 p3 = makePoint(coords);
        double pv[3]; pv[0] = static_cast<double>(p3.x()); pv[1] = static_cast<double>(p3.y()); pv[2] = static_cast<double>(p3.z());
        double w[3]; sub(pv, p0, w);
        x = dot(w, u);
        y = dot(w, v);
    };

    // insert all included known nodes into dt2
    for (const auto &kv : knownNodeCoordsRef_) {
        const auto &la = kv.first;
        const auto &ne = kv.second;
        if (!ne.attachedToDT) continue;
        if (includePred && !includePred(ne)) continue;
        double x,y;
        project_to_2d(ne.coords, x, y);
        Vertex_handle2 vh2 = dt2.insert(Point2(x,y));
        addr2vh[la] = vh2;
        vh2addrs[vh2].insert(la);
    }

    // ensure target addr is present in dt2 (use its coords if known, else use targetCoord)
    if (addr2vh.find(addr) == addr2vh.end()) {
        auto kn = knownNodeCoordsRef_.find(addr);
        std::vector<double> coords = (kn == knownNodeCoordsRef_.end() || kn->second.coords.empty()) ? targetCoord : kn->second.coords;
        double x,y;
        project_to_2d(coords, x, y);
        Vertex_handle2 vh2 = dt2.insert(Point2(x,y));
        addr2vh[addr] = vh2;
        vh2addrs[vh2].insert(addr);
    }

    // find target vertex handle in dt2
    auto itTarget2 = addr2vh.find(addr);
    if (itTarget2 == addr2vh.end()) {
        EV_WARN << "computeDTNeighborsForCoord: target not in 2D mapping -> fallback\n";
        for (const auto &kv : knownNodeCoordsRef_) {
            if (kv.first == selfAddr) continue;
            if (includePred && !includePred(kv.second)) continue;
            if (!kv.second.attachedToDT) continue;
            result.push_back(kv.second);
        }
        return result;
    }

    // collect 2D incident vertices around target using a circulator (portable)
    std::vector<Vertex_handle2> neigh2;
    {
        typename Delaunay2::Vertex_circulator vcirc = dt2.incident_vertices(itTarget2->second);
        if (vcirc != 0) {
            typename Delaunay2::Vertex_circulator done = vcirc;
            // The circulator iterates over neighboring vertices; we need to guard against
            // empty/degenerate cases where itTarget2->second may be isolated.
            do {
                if (!dt2.is_infinite(vcirc)) neigh2.push_back(vcirc);
                ++vcirc;
            } while (vcirc != done);
        }
    }

    // expand vertex handles to addresses and build result
    std::unordered_set<std::string> seen;
    for (const auto &vh2 : neigh2) {
        if (dt2.is_infinite(vh2)) continue;
        auto it = vh2addrs.find(vh2);
        if (it == vh2addrs.end()) continue;
        for (const L3Address &nbAddr : it->second) {
            if (nbAddr == selfAddr) continue;
            auto kn = knownNodeCoordsRef_.find(nbAddr);
            if (kn == knownNodeCoordsRef_.end()) continue;
            const NeighborEntry &ne = kn->second;
            if (includePred && !includePred(ne)) continue;
            if (!ne.attachedToDT && nbAddr != addr) continue;
            std::string key = nbAddr.str();
            if (seen.insert(key).second) result.push_back(ne);
        }
    }

    EV_INFO << "computeDTNeighborsForCoord (2D): found " << result.size() << " neighbors\n";
    return result;
}

bool NeighborDT::isDTNeighborOfTarget(const std::vector<double> &targetCoord,
                                      const L3Address &addr,
                                      const L3Address &selfAddr,
                                      std::function<bool(const NeighborEntry&)> includePred) {
    // Ensure target inserted (caller responsibility to insert self as needed)
    if (addrToVhRef_.find(addr) == addrToVhRef_.end()) {
        auto kn = knownNodeCoordsRef_.find(addr);
        if (kn == knownNodeCoordsRef_.end()) {
            NeighborEntry ne;
            ne.addr = addr;
            ne.coords = targetCoord;
            ne.dim = (int)targetCoord.size();
            ne.attachedToDT = true;
            ne.timeout = simTime() + 10;
            knownNodeCoordsRef_[addr] = ne;
            insertNodeToDT(addr, makePoint(targetCoord));
        } else {
            knownNodeCoordsRef_[addr].attachedToDT = true;
            insertNodeToDT(addr, makePoint(kn->second.coords.empty() ? targetCoord : kn->second.coords));
        }
    }

    auto itTargetVH = addrToVhRef_.find(addr);
    if (itTargetVH == addrToVhRef_.end()) {
        // if there are no other valid neighbors except self and addr, return true
        for (const auto &kv : knownNodeCoordsRef_) {
            const auto &a = kv.first;
            const auto &ne = kv.second;
            if (a == selfAddr || a == addr) continue;
            if (!ne.attachedToDT) continue;
            if (includePred && !includePred(ne)) continue;
            return false;
        }
        return true;
    }

    if (dtRef_.dimension() != 3) {
        if (dtRef_.dimension() == 2) {
            EV_INFO<<"isDTNeighborOfTarget:dim 2 compute dt\n";
            auto neighs = computeDTNeighborsForCoord(targetCoord, addr, selfAddr, includePred);
            for (const auto &n : neighs) if (n.addr == selfAddr) return true;
        }

        bool hasOtherValid = false;
        for (const auto &kv : knownNodeCoordsRef_) {
            const auto &a = kv.first;
            const auto &ne = kv.second;
            if (a == selfAddr || a == addr) continue;
            if (!ne.attachedToDT) continue;
            if (includePred && !includePred(ne)) continue;
            hasOtherValid = true;
            break;
        }
        if (!hasOtherValid) return true;
        Point3 tp = makePoint(targetCoord);
        auto itSelf = addrToVhRef_.find(selfAddr);
        if (itSelf == addrToVhRef_.end()) return false;
        return (dtRef_.nearest_vertex(tp) == itSelf->second);
    }

    Vertex_handle tvh = itTargetVH->second;
    std::vector<Delaunay3::Cell_handle> cells;
    dtRef_.incident_cells(tvh, std::back_inserter(cells));

    std::unordered_set<Delaunay3::Vertex_handle> neighVhs;
    for (const auto &ch : cells) {
        if (dtRef_.is_infinite(ch)) continue;
        for (int i = 0; i < 4; ++i) {
            auto vh = ch->vertex(i);
            if (vh != tvh && !dtRef_.is_infinite(vh)) neighVhs.insert(vh);
        }
    }

    auto itSelfVh = addrToVhRef_.find(selfAddr);
    if (itSelfVh == addrToVhRef_.end()) return false;
    Delaunay3::Vertex_handle selfVh = itSelfVh->second;
    bool foundSelf = (neighVhs.find(selfVh) != neighVhs.end());
    if (foundSelf) return true;

    // If no other valid neighbors exist, return true (fallback)
    bool hasOtherValid = false;
    for (const auto &kv : knownNodeCoordsRef_) {
        const auto &a = kv.first;
        const auto &ne = kv.second;
        if (a == selfAddr || a == addr) continue;
        if (!ne.attachedToDT) continue;
        if (includePred && !includePred(ne)) continue;
        hasOtherValid = true;
        break;
    }
    if (!hasOtherValid) return true;
    return false;
}

/*void NeighborDT::computeMinimalNeighborSubset_CGAL_3D(const L3Address &selfAddr,
                                                      std::function<bool(const NeighborEntry&)> includePred) {
    minSimplexNodesRef_.clear();

    auto itSelf = addrToVhRef_.find(selfAddr);
    if (itSelf == addrToVhRef_.end()) return;
    Vertex_handle selfVh = itSelf->second;
    if (dtRef_.is_infinite(selfVh)) return;

    std::vector<Cell_handle> cells;
    dtRef_.incident_cells(selfVh, std::back_inserter(cells));
    if (cells.empty()) return;

    std::vector<std::set<L3Address>> simplices;
    for (const Cell_handle &ch : cells) {
        if (dtRef_.is_infinite(ch)) continue;
        std::set<L3Address> one;
        for (int i = 0; i < 4; ++i) {
            Vertex_handle vh = ch->vertex(i);
            if (vh == selfVh) continue;
            auto aIt = vhToAddrsRef_.find(vh);
            if (aIt == vhToAddrsRef_.end()) continue;
            for (const L3Address &addr : aIt->second) {
                auto kn = knownNodeCoordsRef_.find(addr);
                if (kn == knownNodeCoordsRef_.end()) continue;
                const NeighborEntry &ne = kn->second;
                if (!ne.attachedToDT) continue;
                if (includePred && !includePred(ne)) continue;
                if (addr == selfAddr) continue;
                one.insert(addr);
            }
        }
        if (!one.empty()) simplices.emplace_back(std::move(one));
    }

    if (simplices.empty()) return;

    std::set<int> remaining;
    for (int i = 0; i < (int)simplices.size(); ++i) remaining.insert(i);

    std::unordered_map<std::string, std::set<int>> coverMap;
    std::unordered_map<std::string, L3Address> keyToAddr;
    for (int i = 0; i < (int)simplices.size(); ++i) {
        for (const auto &a : simplices[i]) {
            std::string k = a.str();
            coverMap[k].insert(i);
            keyToAddr[k] = a;
        }
    }

    std::set<L3Address> cover;
    while (!remaining.empty()) {
        std::string bestKey;
        int bestCover = 0;
        for (const auto &kv : coverMap) {
            int c = 0;
            for (int idx : kv.second) if (remaining.count(idx)) ++c;
            if (c > bestCover) { bestCover = c; bestKey = kv.first; }
        }
        if (bestCover == 0) break;
        L3Address chosen = keyToAddr[bestKey];
        cover.insert(chosen);

        std::vector<int> toErase;
        for (int idx : coverMap[bestKey]) if (remaining.count(idx)) toErase.push_back(idx);
        for (int idx : toErase) remaining.erase(idx);
        coverMap.erase(bestKey);
    }

    if (!remaining.empty()) {
        for (int idx : remaining) {
            if (!simplices[idx].empty()) cover.insert(*simplices[idx].begin());
        }
    }

    for (const L3Address &a : cover) {
        auto kn = knownNodeCoordsRef_.find(a);
        if (kn != knownNodeCoordsRef_.end()) minSimplexNodesRef_[a] = kn->second;
    }
}*/
void NeighborDT::computeMinimalNeighborSubset_CGAL_3D(
        const L3Address &selfAddr,
        std::function<bool(const NeighborEntry&)> includePred)
{
    minSimplexNodesRef_.clear();

    // find self vertex handle
    auto itSelf = addrToVhRef_.find(selfAddr);
    if (itSelf == addrToVhRef_.end()) {
        EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: self not in addrToVh\n";
        return;
    }
    Vertex_handle selfVh = itSelf->second;

    if (dtRef_.is_infinite(selfVh)) {
        EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: selfVh is infinite\n";
        return;
    }

    int dim = dtRef_.dimension();
    EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: dt.dimension()=" << dim
              << " number_of_vertices=" << dtRef_.number_of_vertices() << "\n";

    // --- 3D path (original) ---
    if (dim == 3) {
        if (dtRef_.number_of_vertices() < 4) {
            EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: too few vertices for 3D\n";
            return;
        }

        std::vector<Cell_handle> cells;
        try {
            dtRef_.incident_cells(selfVh, std::back_inserter(cells));
        } catch (const std::exception &ex) {
            EV_ERROR << "computeMinimalNeighborSubset_CGAL_3D: CGAL exception incident_cells: " << ex.what() << "\n";
            return;
        } catch (...) {
            EV_ERROR << "computeMinimalNeighborSubset_CGAL_3D: unknown exception incident_cells\n";
            return;
        }

        if (cells.empty()) {
            EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: no incident cells (3D)\n";
            return;
        }

        std::vector<std::set<L3Address>> simplices;
        for (const Cell_handle &ch : cells) {
            if (dtRef_.is_infinite(ch)) continue;
            std::set<L3Address> one;
            for (int i = 0; i < 4; ++i) {
                Vertex_handle vh = ch->vertex(i);
                if (vh == selfVh) continue;
                auto aIt = vhToAddrsRef_.find(vh);
                if (aIt == vhToAddrsRef_.end()) continue;
                for (const L3Address &addr : aIt->second) {
                    auto kn = knownNodeCoordsRef_.find(addr);
                    if (kn == knownNodeCoordsRef_.end()) continue;
                    const NeighborEntry &ne = kn->second;
                    if (!ne.attachedToDT) continue;
                    if (includePred && !includePred(ne)) continue;
                    if (addr == selfAddr) continue;
                    one.insert(addr);
                }
            }
            if (!one.empty()) simplices.emplace_back(std::move(one));
        }

        if (simplices.empty()) {
            EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: simplices empty (3D)\n";
            return;
        }

        // greedy hitting set (as before)
        std::set<int> remaining;
        for (int i = 0; i < (int)simplices.size(); ++i) remaining.insert(i);

        std::unordered_map<std::string, std::set<int>> coverMap;
        std::unordered_map<std::string, L3Address> keyToAddr;
        for (int i = 0; i < (int)simplices.size(); ++i) {
            for (const auto &a : simplices[i]) {
                std::string k = a.str();
                coverMap[k].insert(i);
                keyToAddr[k] = a;
            }
        }

        std::set<L3Address> cover;
        while (!remaining.empty()) {
            std::string bestKey;
            int bestCover = 0;
            for (const auto &kv : coverMap) {
                int c = 0;
                for (int idx : kv.second) if (remaining.count(idx)) ++c;
                if (c > bestCover) { bestCover = c; bestKey = kv.first; }
            }
            if (bestCover == 0) break;
            L3Address chosen = keyToAddr[bestKey];
            cover.insert(chosen);

            std::vector<int> toErase;
            auto it = coverMap.find(bestKey);
            if (it != coverMap.end()) {
                for (int idx : it->second) if (remaining.count(idx)) toErase.push_back(idx);
                for (int idx : toErase) remaining.erase(idx);
                coverMap.erase(it);
            } else {
                break;
            }
        }

        if (!remaining.empty()) {
            for (int idx : remaining) {
                if (!simplices[idx].empty()) cover.insert(*simplices[idx].begin());
            }
        }

        for (const L3Address &a : cover) {
            auto kn = knownNodeCoordsRef_.find(a);
            if (kn != knownNodeCoordsRef_.end()) minSimplexNodesRef_[a] = kn->second;
        }

        EV_INFO << "computeMinimalNeighborSubset_CGAL_3D (3D): selected size=" << minSimplexNodesRef_.size() << "\n";
        return;
    }

    // --- 2D path (coplanar points) ---
    if (dim == 2) {
        EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: entering 2D fallback\n";

        // typedefs for 2D triangulation
        typedef CGAL::Delaunay_triangulation_2<K> Delaunay2;
        typedef K::Point_2 Point2;
        typedef Delaunay2::Vertex_handle Vertex_handle2;
        typedef Delaunay2::Vertex_circulator Vertex_circulator2;

        // gather three reference non-collinear points from finite vertices
        std::vector<Point3> refs3;
        refs3.reserve(3);
        for (auto vit = dtRef_.finite_vertices_begin(); vit != dtRef_.finite_vertices_end() && refs3.size() < 3; ++vit) {
            refs3.push_back(vit->point());
        }
        if (refs3.size() < 3) {
            EV_WARN << "computeMinimalNeighborSubset_CGAL_3D: not enough finite vertices to build plane\n";
            return;
        }
        auto to_arr = [&](const Point3 &p, double out[3]) {
            out[0] = static_cast<double>(p.x());
            out[1] = static_cast<double>(p.y());
            out[2] = static_cast<double>(p.z());
        };
        double p0[3], p1[3], p2[3];
        to_arr(refs3[0], p0); to_arr(refs3[1], p1); to_arr(refs3[2], p2);

        auto sub = [](const double a[3], const double b[3], double out[3]) {
            out[0] = a[0] - b[0]; out[1] = a[1] - b[1]; out[2] = a[2] - b[2];
        };
        auto cross = [](const double a[3], const double b[3], double out[3]) {
            out[0] = a[1]*b[2] - a[2]*b[1];
            out[1] = a[2]*b[0] - a[0]*b[2];
            out[2] = a[0]*b[1] - a[1]*b[0];
        };
        auto dot = [](const double a[3], const double b[3]) {
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
        };
        auto norm = [&](const double a[3]) {
            return std::sqrt(dot(a,a));
        };
        auto normalize = [&](double a[3]) {
            double n = norm(a);
            if (n > 0.0) { a[0]/=n; a[1]/=n; a[2]/=n; }
        };

        double u[3], tmpv[3], nvec[3], v[3];
        sub(p1, p0, u); normalize(u);
        sub(p2, p0, tmpv);
        cross(u, tmpv, nvec);
        if (norm(nvec) < 1e-12) {
            EV_WARN << "computeMinimalNeighborSubset_CGAL_3D: refs nearly collinear -> cannot build plane\n";
            return;
        }
        normalize(nvec);
        cross(nvec, u, v); normalize(v);

        auto project_to_2d = [&](const Point3 &p3, double &x, double &y) {
            double pv[3]; pv[0] = static_cast<double>(p3.x()); pv[1] = static_cast<double>(p3.y()); pv[2] = static_cast<double>(p3.z());
            double w[3]; sub(pv, p0, w);
            x = dot(w, u);
            y = dot(w, v);
        };

        // Build 2D triangulation and mapping
        Delaunay2 dt2;
        std::map<L3Address, Vertex_handle2> addr2vh;
        std::map<Vertex_handle2, std::set<L3Address>> vh2addrs2;

        // insert all finite 3D vertices mapped to 2D
        for (auto vit = dtRef_.finite_vertices_begin(); vit != dtRef_.finite_vertices_end(); ++vit) {
            Point3 p3 = vit->point();
            double x,y; project_to_2d(p3, x, y);
            Vertex_handle2 vh2 = dt2.insert(Point2(x,y));
            // map addresses that correspond to this 3D vh
            auto itmap = vhToAddrsRef_.find(vit);
            if (itmap != vhToAddrsRef_.end()) {
                for (const auto &addr : itmap->second) {
                    addr2vh[addr] = vh2;
                    vh2addrs2[vh2].insert(addr);
                }
            } else {
                // try match by coords against known map if mapping missing
                for (const auto &kv : knownNodeCoordsRef_) {
                    if (kv.second.coords.size() == 3) {
                        if (std::fabs(kv.second.coords[0] - p3.x()) < 1e-9 &&
                            std::fabs(kv.second.coords[1] - p3.y()) < 1e-9 &&
                            std::fabs(kv.second.coords[2] - p3.z()) < 1e-9) {
                            addr2vh[kv.first] = vh2;
                            vh2addrs2[vh2].insert(kv.first);
                        }
                    }
                }
            }
        }

        // ensure self present in addr2vh (if not, insert using self 3D point)
        Point3 selfP = selfVh->point();
        double sx, sy; project_to_2d(selfP, sx, sy);
        if (addr2vh.find(selfAddr) == addr2vh.end()) {
            Vertex_handle2 self2 = dt2.insert(Point2(sx, sy));
            addr2vh[selfAddr] = self2;
            vh2addrs2[self2].insert(selfAddr);
        }
        Vertex_handle2 self2 = addr2vh[selfAddr];

        // collect neighbor vertex handles around self2 using circulator
        std::vector<Vertex_handle2> neigh2;
        {
            Vertex_circulator2 vcirc = dt2.incident_vertices(self2);
            if (vcirc != Vertex_circulator2()) {
                Vertex_circulator2 done = vcirc;
                do {
                    if (!dt2.is_infinite(vcirc)) neigh2.push_back(vcirc);
                    ++vcirc;
                } while (vcirc != done);
            }
        }
        if (neigh2.empty()) {
            EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: no 2D incident vertices\n";
            return;
        }

        // build simplices: for each consecutive pair (v_i, v_{i+1}), triangle (self2, v_i, v_{i+1}),
        // simplex nodes = addresses on v_i union addresses on v_{i+1}
        std::vector<std::set<L3Address>> simplices;
        for (size_t i = 0; i < neigh2.size(); ++i) {
            Vertex_handle2 va = neigh2[i];
            Vertex_handle2 vb = neigh2[(i+1) % neigh2.size()];
            std::set<L3Address> one;
            auto ita = vh2addrs2.find(va);
            if (ita != vh2addrs2.end()) for (const auto &a : ita->second) {
                if (a == selfAddr) continue;
                auto kn = knownNodeCoordsRef_.find(a);
                if (kn == knownNodeCoordsRef_.end()) continue;
                if (!kn->second.attachedToDT) continue;
                if (includePred && !includePred(kn->second)) continue;
                one.insert(a);
            }
            auto itb = vh2addrs2.find(vb);
            if (itb != vh2addrs2.end()) for (const auto &a : itb->second) {
                if (a == selfAddr) continue;
                auto kn = knownNodeCoordsRef_.find(a);
                if (kn == knownNodeCoordsRef_.end()) continue;
                if (!kn->second.attachedToDT) continue;
                if (includePred && !includePred(kn->second)) continue;
                one.insert(a);
            }
            if (!one.empty()) simplices.emplace_back(std::move(one));
        }

        if (simplices.empty()) {
            EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: simplices empty (2D)\n";
            return;
        }

        // greedy hitting set (same as 3D)
        std::set<int> remaining;
        for (int i = 0; i < (int)simplices.size(); ++i) remaining.insert(i);

        std::unordered_map<std::string, std::set<int>> coverMap;
        std::unordered_map<std::string, L3Address> keyToAddr;
        for (int i = 0; i < (int)simplices.size(); ++i) {
            for (const auto &a : simplices[i]) {
                std::string k = a.str();
                coverMap[k].insert(i);
                keyToAddr[k] = a;
            }
        }

        std::set<L3Address> cover;
        while (!remaining.empty()) {
            std::string bestKey;
            int bestCover = 0;
            for (const auto &kv : coverMap) {
                int c = 0;
                for (int idx : kv.second) if (remaining.count(idx)) ++c;
                if (c > bestCover) { bestCover = c; bestKey = kv.first; }
            }
            if (bestCover == 0) break;
            L3Address chosen = keyToAddr[bestKey];
            cover.insert(chosen);

            std::vector<int> toErase;
            auto it = coverMap.find(bestKey);
            if (it != coverMap.end()) {
                for (int idx : it->second) if (remaining.count(idx)) toErase.push_back(idx);
                for (int idx : toErase) remaining.erase(idx);
                coverMap.erase(it);
            } else break;
        }

        if (!remaining.empty()) {
            for (int idx : remaining) {
                if (!simplices[idx].empty()) cover.insert(*simplices[idx].begin());
            }
        }

        for (const L3Address &a : cover) {
            auto kn = knownNodeCoordsRef_.find(a);
            if (kn != knownNodeCoordsRef_.end()) minSimplexNodesRef_[a] = kn->second;
        }

        EV_INFO << "computeMinimalNeighborSubset_CGAL_3D (2D): selected size=" << minSimplexNodesRef_.size() << "\n";
        return;
    }

    // dim <= 1: cannot compute simplices meaningfully
    EV_DETAIL << "computeMinimalNeighborSubset_CGAL_3D: dt dimension <= 1 -> skipping\n";
    return;
}

} // namespace routing
} // namespace mysrc


