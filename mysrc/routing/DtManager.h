/*
 * DtManager.h
 *
 *  Created on: Dec 7, 2025
 *      Author: yychen
 */

#ifndef MYSRC_ROUTING_DTMANAGER_H_
#define MYSRC_ROUTING_DTMANAGER_H_

#include <omnetpp.h>
#include <inet/common/INETDefs.h>
#include <inet/networklayer/common/L3Address.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // for 2D kernel too
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <map>
#include <set>
#include <vector>
#include <unordered_map>
#include <functional>
#include <optional>
#include <cmath>
#include <algorithm>

#include <cstdint>
#include <unordered_set>

#include"NeighborEntry.h"

using namespace omnetpp;
using namespace inet;

class DtManager {
public:
    /*struct NeighborEntry {
        bool attachedToDT = false;
        std::vector<double> coords; // expected {x,y[,z]} ; dim optionally present in coords.size()
        int dim = 3; // logical dimension indicator (usually 3)
        // other fields (timeout, addr, ...) are out of scope for DtManager
    };*/

    enum class UpdateMode { FULL, INCREMENTAL };

    DtManager(UpdateMode mode = UpdateMode::INCREMENTAL,
              double coordinateEps = 1e-9)
        : mode(mode), eps(coordinateEps),
          currentDimension(3), using2D(false) {}


    ~DtManager() = default;

    // --- Public API -------------------------------------------------------

    // Full build from provided map (only nodes with attachedToDT==true considered)
    void buildFromMap(const std::map<L3Address, mysrc::routing::NeighborEntry>& nodes) {
        // clear current structures
        clearAll();

        initialized = true;

        // collect points (only attached)
        std::vector<std::pair<Point3_3D, L3Address>> pts3;
        for (const auto &kv : nodes) {
            if (!kv.second.attachedToDT) continue;
            auto p = makePoint3(kv.second);
            if (!p) continue;
            pts3.emplace_back(*p, kv.first);
        }

        // special handling: <=3 nodes -> we may avoid building CGAL at all
        if ((int)pts3.size() <= 3) {
            // keep dt empty, mark currentDimension accordingly
            currentDimension = (int)pts3.size() - 1; // 0-based heuristic
            using2D = false;
            // store vertex list for neighbor queries via special-case
            smallSetPoints3 = pts3;
            return;
        }

        // detect coplanarity
        bool coplanar = arePointsCoplanar(pts3);

        if (coplanar) {
            // build 2D Delaunay by projection
            build2DFrom3DPoints(pts3);
            using2D = true;
            currentDimension = 2;
        } else {
            // build 3D Delaunay
            build3D(pts3);
            using2D = false;
            currentDimension = 3;
        }
    }

    // Update from map. If mode==FULL we call buildFromMap; if INCREMENTAL do local updates.
    // selfAddr: address of this node. newNeighborsOut: append newly-discovered neighbors of self in this update.
    void updateFromMap(const std::map<L3Address, mysrc::routing::NeighborEntry>& nodes,
                       const L3Address& selfAddr,
                       std::vector<L3Address>& newNeighborsOut)
    {
        if (!initialized) {
                buildFromMap(nodes);
                std::vector<L3Address> cur = getDTNeighbors(selfAddr);
                newNeighborsOut.insert(newNeighborsOut.end(), cur.begin(), cur.end());
                return;
        }

        if (mode == UpdateMode::FULL) {
            // full rebuild and compute new neighbors diff
            std::set<L3Address> prev = getDTNeighborsSetUnsafe(selfAddr);
            buildFromMap(nodes);
            std::set<L3Address> cur = getDTNeighborsSetUnsafe(selfAddr);
            for (const auto &a : cur) if (prev.find(a) == prev.end()) newNeighborsOut.push_back(a);
            return;
        }
        //for fallback
        std::set<L3Address> prevneighbors = getDTNeighborsSetUnsafe(selfAddr);

        // INCREMENTAL: attempt minimal changes
        // 1) gather should-exist set (only attached)
        std::set<L3Address> shouldExist;
        for (const auto &kv : nodes) if (kv.second.attachedToDT) shouldExist.insert(kv.first);

        // 2) special-case sizes (0/1/2/3) -> if nodes small, we will clear CGAL and handle via smallSetPoints3
        size_t shouldSize = shouldExist.size();
        if (shouldSize <= 3) {
            // compute prev self neighbors
            std::set<L3Address> prev = getDTNeighborsSetUnsafe(selfAddr);

            // clear CGAL and fill smallSetPoints3
            clearAll();
            std::vector<std::pair<Point3_3D,L3Address>> pts3;
            for (const auto &addr : shouldExist) {
                auto p = makePoint3(nodes.at(addr));
                if (p) pts3.emplace_back(*p, addr);
            }
            smallSetPoints3 = pts3;
            currentDimension = (int)std::max<int>(0, (int)pts3.size() - 1);
            using2D = false;

            std::set<L3Address> cur = getDTNeighborsSetUnsafe(selfAddr);
            for (const auto &a : cur) if (prev.find(a) == prev.end()) newNeighborsOut.push_back(a);
            return;
        }

        // 3) For >=4 points: we must decide whether coplanar or full 3D
        // collect candidate pts3
        std::vector<std::pair<Point3_3D,L3Address>> pts3;
        pts3.reserve(shouldExist.size());
        for (const auto &addr : shouldExist) {
            auto it = nodes.find(addr);
            if (it == nodes.end()) continue;
            auto p = makePoint3(it->second);
            if (p) pts3.emplace_back(*p, addr);
        }
        if (pts3.size() < 4) {
            // fallback: treat similar to <=3 case
            std::set<L3Address> prev = getDTNeighborsSetUnsafe(selfAddr);
            clearAll();
            smallSetPoints3 = pts3;
            using2D = false;
            currentDimension = (int)std::max<int>(0, (int)pts3.size() - 1);
            std::set<L3Address> cur = getDTNeighborsSetUnsafe(selfAddr);
            for (const auto &a : cur) if (prev.find(a) == prev.end()) newNeighborsOut.push_back(a);
            return;
        }

        bool coplanar = arePointsCoplanar(pts3);

        if (coplanar) {
            // must switch into 2D mode and do a full 2D rebuild
            std::set<L3Address> prev = getDTNeighborsSetUnsafe(selfAddr);
            clearAll(); // clear 3D
            build2DFrom3DPoints(pts3);
            using2D = true;
            currentDimension = 2;
            std::set<L3Address> cur = getDTNeighborsSetUnsafe(selfAddr);
            for (const auto &a : cur) if (prev.find(a) == prev.end()) newNeighborsOut.push_back(a);
            return;
        }

        // 4) Non-coplanar (proper 3D). Attempt incremental updates against existing 3D dt.
        if (using2D) {
            // currently in 2D mode, but incoming is 3D ; rebuild 3D from map (full rebuild)
            std::set<L3Address> prev = getDTNeighborsSetUnsafe(selfAddr);
            clearAll();
            build3D(pts3);
            using2D = false;
            currentDimension = 3;
            std::set<L3Address> cur = getDTNeighborsSetUnsafe(selfAddr);
            for (const auto &a : cur) if (prev.find(a) == prev.end()) newNeighborsOut.push_back(a);
            return;
        }

        // incremental changes on 3D structure:
        // (a) remove those present in dt but not in shouldExist
        std::vector<L3Address> toRemove;
        for (const auto &kv : addrToVh3) {
            if (shouldExist.find(kv.first) == shouldExist.end()) toRemove.push_back(kv.first);
        }
        for (const auto &a : toRemove) {
            try {
                dt3.remove(addrToVh3[a]);
            } catch (...) {
                // fallback to full rebuild if local removals fail
                EV_WARN << "DtManager: remove failed -> full rebuild\n";
                std::set<L3Address> prev = getDTNeighborsSetUnsafe(selfAddr);
                clearAll();
                build3D(pts3);
                using2D = false;
                currentDimension = 3;
                std::set<L3Address> cur = getDTNeighborsSetUnsafe(selfAddr);
                for (const auto &x : cur) if (prev.find(x) == prev.end()) newNeighborsOut.push_back(x);
                return;
            }
            addrToVh3.erase(a);
        }

        // (b) insert or move existing points
        std::set<L3Address> prev = getDTNeighborsSetUnsafe(selfAddr);
        for (const auto &pinfo : pts3) {
            const L3Address &addr = pinfo.second;
            const Point3_3D &p = pinfo.first;
            auto it = addrToVh3.find(addr);
            if (it == addrToVh3.end()) {
                // insert new
                //auto vh = dt3.insert(std::make_pair(p, addr));
                auto vh = insertVertexWithInfo3(dt3, p, addr);
                if (vh != Delaunay3::Vertex_handle()) addrToVh3.emplace(addr, vh);
                else {
                    EV_WARN << "DtManager: 3D insert failed -> full rebuild\n";
                    clearAll();
                    build3D(pts3);
                    using2D = false;
                    currentDimension = 3;
                    std::set<L3Address> cur = getDTNeighborsSetUnsafe(selfAddr);
                    for (const auto &x : cur) if (prev.find(x) == prev.end()) newNeighborsOut.push_back(x);
                    return;
                }
            } else {//not safe enough
                // try move if needed
                // try safe move if needed
                /*auto vh = it->second;
                if (!pointsEqual3(vh->point(), p)) {
                    // safe approach: refresh handle, check collision, else defer to rebuild
                    auto storedVhIt = addrToVh3.find(addr);
                    Delaunay3::Vertex_handle currVh = Delaunay3::Vertex_handle();
                    if (storedVhIt != addrToVh3.end()) {
                        // do not trust stored handle; find live handle by info
                        currVh = findVertexHandleByInfo(addr);
                    }

                    // if we didn't find an up-to-date handle, we cannot safely move/remove -> defer rebuild
                    if (currVh == Delaunay3::Vertex_handle()) {
                        EV_WARN << "DtManager: live vertex handle not found for " << addr << ", deferring full rebuild\n";
                        pendingRebuildAddrs.insert(addr);
                    } else {
                        // if points already equal -> nothing to do
                        if (pointsEqual3(currVh->point(), p)) {
                            // already at correct place
                        } else {
                            // collision check: is another vertex already at target p?
                            auto near = dt3.nearest_vertex(p);
                            if (near != Delaunay3::Vertex_handle() && near != currVh && pointsEqual3(near->point(), p)) {
                                EV_WARN << "DtManager: move target collides for " << addr << " with " << near->info() << " -> defer rebuild\n";
                                pendingRebuildAddrs.insert(addr);
                            } else {
                                // don't call dt3.move here to avoid segfault; defer to rebuild for safety
                                EV_WARN << "DtManager: skipping dt3.move for " << addr << " (safe-mode) -> will rebuild later\n";
                                pendingRebuildAddrs.insert(addr);
                            }
                        }
                    }

                }*/
            }
        }

        // (c) finally compute newly added self neighbors
        std::set<L3Address> cur = getDTNeighborsSetUnsafe(selfAddr);
        for (const auto &a : cur) if (prev.find(a) == prev.end()) newNeighborsOut.push_back(a);
    }

    // Return vector of DT neighbors of an address (works for small sets, 3D dt, or 2D dt)
    std::vector<L3Address> getDTNeighbors(const L3Address& addr) const  {
        std::vector<L3Address> out;
        // small sets
        if (!smallSetPoints3.empty()) {
            // if only 1 point -> no neighbors
            if (smallSetPoints3.size() <= 1) return out;
            for (const auto &p : smallSetPoints3) {
                if (p.second == addr) continue;
                out.push_back(p.second);
            }
            return out;
        }

        if (using2D) {
            auto it = addrToVh2.find(addr);
            if (it == addrToVh2.end()) return out;
            auto vh = it->second;
            // special-case dimension 2/1 in dt2
            if (dt2.dimension() == 1) {
                // incident vertices along line
                auto circ = dt2.incident_vertices(vh);
                auto done = circ;
                if (circ == 0) return out;
                do {
                    if (!dt2.is_infinite(circ)) out.push_back(circ->info());
                    ++circ;
                } while (circ != done);
                return out;
            }
            // general 2D case
            auto vc = dt2.incident_vertices(vh);
            if (vc == 0) return out;
            auto done = vc;
            do {
                if (!dt2.is_infinite(vc)) out.push_back(vc->info());
                ++vc;
            } while (vc != done);
            return out;
        } else {
            // 3D
            auto it = addrToVh3.find(addr);
            if (it == addrToVh3.end()) return out;
            auto vh = it->second;
            if (dt3.number_of_vertices() == 2) {
                for (auto vit = dt3.finite_vertices_begin(); vit != dt3.finite_vertices_end(); ++vit) {
                    if (vit->info() == addr) continue;
                    out.push_back(vit->info());
                }
                return out;
            }
            /*auto vc = dt3.incident_vertices(vh);
            if (vc == 0) return out;
            auto done = vc;
            do {
                if (!dt3.is_infinite(vc)) out.push_back(vc->info());
                ++vc;
            } while (vc != done);
            return out;*/

            // Collect neighbors by scanning finite cells that contain vh.
            // This avoids using incident_vertices which may not be available as const in some CGAL setups.
            std::set<L3Address> seen;
            for (auto cit = dt3.finite_cells_begin(); cit != dt3.finite_cells_end(); ++cit) {
                // skip infinite cells (finite_cells_* usually already ensures this)
                if (dt3.is_infinite(cit)) continue;

                // check whether this cell contains vh
                bool contains = false;
                for (int vi = 0; vi < 4; ++vi) {
                    if (cit->vertex(vi) == vh) { contains = true; break; }
                }
                if (!contains) continue;

                // collect other vertices from this tetrahedron
                for (int vi = 0; vi < 4; ++vi) {
                    auto vvh = cit->vertex(vi);
                    if (vvh != vh && !dt3.is_infinite(vvh)) {
                        const L3Address &a = vvh->info();
                        if (seen.insert(a).second) { // newly inserted -> append
                            out.push_back(a);
                        }
                    }
                }
            }
            return out;
        }
    }

    // getMinSimplexNodes implementation
    std::vector<L3Address> getMinSimplexNodes(const L3Address &selfAddr) {
        // Degenerate / trivial cases: 2D, small set (<=3), or self not present in 3D triangulation
        if (using2D || smallSetPoints3.size() <= 3) {
            return getDTNeighbors(selfAddr);
        }

        auto it = addrToVh3.find(selfAddr);
        if (it == addrToVh3.end()) {
            // not present in 3D DT -> fallback to neighbors
            return getDTNeighbors(selfAddr);
        }
        auto vh = it->second;

        // collect all finite incident cells (tetrahedra) that contain vh
        std::vector<std::vector<L3Address>> simplices;
        {
            // iterate all finite cells and pick those that contain vh
            for (auto cit = dt3.finite_cells_begin(); cit != dt3.finite_cells_end(); ++cit) {
                // skip infinite cells just in case (finite_cells_* should already do that)
                if (dt3.is_infinite(cit)) continue;

                // check whether this cell contains vh as one of its vertices
                bool contains = false;
                for (int vi = 0; vi < 4; ++vi) {
                    if (cit->vertex(vi) == vh) { contains = true; break; }
                }
                if (!contains) continue;

                // collect other vertices (the 3 neighbors) from this tetrahedron
                std::vector<L3Address> neighs;
                neighs.reserve(3);
                for (int vi = 0; vi < 4; ++vi) {
                    auto vvh = cit->vertex(vi);
                    if (vvh != vh && !dt3.is_infinite(vvh)) {
                        neighs.push_back(vvh->info());
                    }
                }
                if (!neighs.empty()) {
                    std::sort(neighs.begin(), neighs.end());
                    neighs.erase(std::unique(neighs.begin(), neighs.end()), neighs.end());
                    simplices.push_back(std::move(neighs));
                }
            }

            // if no simplices were found, fall back to neighbors
            if (simplices.empty()) {
                return getDTNeighbors(selfAddr);
            }
        }


        if (simplices.empty()) {
            return getDTNeighbors(selfAddr);
        }

        // Build universe of neighbor addresses appearing in any simplex
        std::map<L3Address, int> addrToIndex;
        std::vector<L3Address> indexToAddr;
        indexToAddr.reserve(32);
        for (const auto &s : simplices) {
            for (const auto &a : s) {
                if (addrToIndex.find(a) == addrToIndex.end()) {
                    int idx = (int)indexToAddr.size();
                    addrToIndex.emplace(a, idx);
                    indexToAddr.push_back(a);
                }
            }
        }
        int m = (int)indexToAddr.size();
        if (m == 0) return {};

        // If small universe, use uint64_t bitmask exact search else greedy fallback
        const int BIT_LIMIT = 63; // safe bit limit for uint64_t
        if (m <= BIT_LIMIT) {
            std::vector<std::uint64_t> simplexMasks;
            simplexMasks.reserve(simplices.size());
            for (const auto &s : simplices) {
                std::uint64_t mask = 0;
                for (const auto &a : s) {
                    mask |= (std::uint64_t)1 << addrToIndex.at(a);
                }
                simplexMasks.push_back(mask);
            }

            // Try exact search by increasing subset size k
            std::vector<int> solutionIndices;
            bool found = false;

            // recursive DFS with pruning; picks 'remain' items starting at start, maintains curMask
            std::function<void(int,int,std::uint64_t,std::vector<int>&)> dfs =
            [&](int start, int remain, std::uint64_t curMask, std::vector<int> &chosen) {
                if (found) return;
                if (remain == 0) {
                    // check cover: every simplexMask must intersect curMask
                    for (const auto &sm : simplexMasks) {
                        if ((curMask & sm) == 0) return; // this chosen set misses a simplex
                    }
                    solutionIndices = chosen;
                    found = true;
                    return;
                }
                if ((int)indexToAddr.size() - start < remain) return; // not enough left
                for (int i = start; i < m; ++i) {
                    chosen.push_back(i);
                    dfs(i + 1, remain - 1, curMask | ((std::uint64_t)1 << i), chosen);
                    chosen.pop_back();
                    if (found) return;
                }
            };

            for (int k = 1; k <= m && !found; ++k) {
                std::vector<int> chosen;
                dfs(0, k, 0ULL, chosen);
                if (found) break;
            }

            if (found) {
                std::vector<L3Address> res;
                res.reserve(solutionIndices.size());
                for (int idx : solutionIndices) res.push_back(indexToAddr[idx]);
                return res;
            }

            // fallback to greedy cover if exact search didn't find (should be rare)
            std::vector<int> greedyIdx = greedy_set_cover_ints(m, simplexMasks);
            std::vector<L3Address> greedyRes;
            greedyRes.reserve(greedyIdx.size());
            for (int idx : greedyIdx) greedyRes.push_back(indexToAddr[idx]);
            return greedyRes;
        } else {
            // Universe too large for exact bitmask method -> greedy approximation using sets
            std::vector<std::unordered_set<int>> simplexSets;
            simplexSets.reserve(simplices.size());
            for (const auto &s : simplices) {
                std::unordered_set<int> st;
                for (const auto &a : s) st.insert(addrToIndex.at(a));
                simplexSets.push_back(std::move(st));
            }

            std::vector<int> chosenIdx;
            std::vector<char> covered(simplexSets.size(), 0);
            int uncoveredCount = (int)simplexSets.size();

            while (uncoveredCount > 0) {
                int bestIdx = -1;
                int bestCover = -1;
                for (int j = 0; j < m; ++j) {
                    int coverCnt = 0;
                    for (size_t si = 0; si < simplexSets.size(); ++si) {
                        if (covered[si]) continue;
                        if (simplexSets[si].count(j)) ++coverCnt;
                    }
                    if (coverCnt > bestCover) { bestCover = coverCnt; bestIdx = j; }
                }
                if (bestIdx < 0 || bestCover == 0) break;
                chosenIdx.push_back(bestIdx);
                for (size_t si = 0; si < simplexSets.size(); ++si) {
                    if (!covered[si] && simplexSets[si].count(bestIdx)) {
                        covered[si] = 1;
                        --uncoveredCount;
                    }
                }
            }

            std::vector<L3Address> res;
            res.reserve(chosenIdx.size());
            for (int idx : chosenIdx) res.push_back(indexToAddr[idx]);
            return res;
        }
    }


    // utility: clear all internal structures
    void clearAll() {
        // Clear 3D
        dt3.clear();
        addrToVh3.clear();
        // Clear 2D
        dt2.clear();
        addrToVh2.clear();
        // small set
        smallSetPoints3.clear();
        using2D = false;
        currentDimension = 3;
    }

    // set update mode at runtime
    void setUpdateMode(UpdateMode m) { mode = m; }

private:
    // ------------------- CGAL kernels & typedefs -------------------------
    using K3 = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point3_3D = K3::Point_3;

    using Vb3 = CGAL::Triangulation_vertex_base_with_info_3<L3Address, K3>;
    using Tds3 = CGAL::Triangulation_data_structure_3<Vb3>;
    using Delaunay3 = CGAL::Delaunay_triangulation_3<K3, Tds3>;

    // 2D
    using K2 = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point2_2D = K2::Point_2;
    using Vb2 = CGAL::Triangulation_vertex_base_with_info_2<L3Address, K2>;
    using Tds2 = CGAL::Triangulation_data_structure_2<Vb2>;
    using Delaunay2 = CGAL::Delaunay_triangulation_2<K2, Tds2>;

    // ------------------- members ---------------------------------------
    Delaunay3 dt3;
    std::map<L3Address, Delaunay3::Vertex_handle> addrToVh3;

    Delaunay2 dt2;
    std::map<L3Address, Delaunay2::Vertex_handle> addrToVh2;

    // for small sizes (<=3) we avoid CGAL and keep explicit list
    std::vector<std::pair<Point3_3D, L3Address>> smallSetPoints3;


    UpdateMode mode;
    double eps;
    int currentDimension;
    bool using2D;

    bool initialized = false;

    // ------------------- helpers ---------------------------------------
    static double absd(double x) { return x < 0 ? -x : x; }

    std::optional<Point3_3D> makePoint3(const mysrc::routing::NeighborEntry& ne) const {
        if (ne.coords.size() < 2) return std::nullopt;
        double x = ne.coords[0];
        double y = ne.coords[1];
        double z = 0.0;
        if (ne.coords.size() >= 3) z = ne.coords[2];
        return Point3_3D(x, y, z);
    }

    bool pointsEqual3(const Point3_3D &a, const Point3_3D &b) const {
        return (absd(a.x() - b.x()) <= eps) && (absd(a.y() - b.y()) <= eps) && (absd(a.z() - b.z()) <= eps);
    }

    // quick coplanar test for a set of >=4 points (returns true if all points coplanar within eps)
    bool arePointsCoplanar(const std::vector<std::pair<Point3_3D, L3Address>>& pts) const {
        if (pts.size() < 4) return true; // vacuously
        // find first 3 non-collinear points
        size_t n = pts.size();
        bool found = false;
        Point3_3D p0, p1, p2;
        for (size_t i = 0; i < n && !found; ++i) {
            for (size_t j = i+1; j < n && !found; ++j) {
                for (size_t k = j+1; k < n && !found; ++k) {
                    auto v1x = pts[j].first.x() - pts[i].first.x();
                    auto v1y = pts[j].first.y() - pts[i].first.y();
                    auto v1z = pts[j].first.z() - pts[i].first.z();
                    auto v2x = pts[k].first.x() - pts[i].first.x();
                    auto v2y = pts[k].first.y() - pts[i].first.y();
                    auto v2z = pts[k].first.z() - pts[i].first.z();
                    // cross product
                    double cx = v1y * v2z - v1z * v2y;
                    double cy = v1z * v2x - v1x * v2z;
                    double cz = v1x * v2y - v1y * v2x;
                    double norm = std::sqrt(cx*cx + cy*cy + cz*cz);
                    if (norm > eps) {
                        p0 = pts[i].first;
                        p1 = pts[j].first;
                        p2 = pts[k].first;
                        found = true;
                    }
                }
            }
        }
        if (!found) {
            // All points are collinear -> treat as coplanar (degenerate)
            return true;
        }
        // plane normal
        double nx = (p1.y()-p0.y())*(p2.z()-p0.z()) - (p1.z()-p0.z())*(p2.y()-p0.y());
        double ny = (p1.z()-p0.z())*(p2.x()-p0.x()) - (p1.x()-p0.x())*(p2.z()-p0.z());
        double nz = (p1.x()-p0.x())*(p2.y()-p0.y()) - (p1.y()-p0.y())*(p2.x()-p0.x());
        double nlen = std::sqrt(nx*nx + ny*ny + nz*nz);
        if (nlen <= eps) return true; // degeneracy -> treat as coplanar

        // normalize
        nx /= nlen; ny /= nlen; nz /= nlen;

        // check each point's distance to plane
        for (const auto &pp : pts) {
            const auto &q = pp.first;
            double vx = q.x() - p0.x();
            double vy = q.y() - p0.y();
            double vz = q.z() - p0.z();
            double dist = absd(nx*vx + ny*vy + nz*vz);
            if (dist > 1e-6) return false; // not coplanar
        }
        return true;
    }

    // build 3D dt from list of (Point3,L3Address)
    void build3D(const std::vector<std::pair<Point3_3D,L3Address>>& pts) {
        dt3.clear(); addrToVh3.clear();
        for (const auto &pinfo : pts) {
            auto vh = insertVertexWithInfo3(dt3, pinfo.first, pinfo.second);
            if (vh != Delaunay3::Vertex_handle()) addrToVh3.emplace(pinfo.second, vh);
        }
    }

    // build 2D dt via projection of 3D points onto plane determined by first three non-collinear pts
    void build2DFrom3DPoints(const std::vector<std::pair<Point3_3D,L3Address>>& pts) {
        dt2.clear(); addrToVh2.clear();

        // find three non-collinear (or fallback)
        Point3_3D p0 = pts[0].first, p1 = pts[1].first, p2 = pts[2].first;
        bool found = false;
        for (size_t i = 0; i < pts.size() && !found; ++i) {
            for (size_t j = i+1; j < pts.size() && !found; ++j) {
                for (size_t k = j+1; k < pts.size() && !found; ++k) {
                    auto v1 = pts[j].first - pts[i].first;
                    auto v2 = pts[k].first - pts[i].first;
                    // cross
                    double cx = v1.y()*v2.z() - v1.z()*v2.y();
                    double cy = v1.z()*v2.x() - v1.x()*v2.z();
                    double cz = v1.x()*v2.y() - v1.y()*v2.x();
                    double norm = std::sqrt(cx*cx + cy*cy + cz*cz);
                    if (norm > eps) {
                        p0 = pts[i].first; p1 = pts[j].first; p2 = pts[k].first;
                        found = true;
                    }
                }
            }
        }
        if (!found) {
            // collinear fallback: choose arbitrary orthonormal basis
            // choose p0 = pts[0]
            p0 = pts[0].first;
            // pick p1 = any != p0
            p1 = pts[1].first;
            // choose p2 = p1 + orthogonal small vector
            p2 = Point3_3D(p0.x(), p0.y(), p0.z() + 1.0);
        }

        // build orthonormal basis (u,v) on plane
        // u = normalize(p1 - p0)
        double ux = p1.x() - p0.x();
        double uy = p1.y() - p0.y();
        double uz = p1.z() - p0.z();
        double ulen = std::sqrt(ux*ux + uy*uy + uz*uz);
        if (ulen <= eps) {
            // degenerate fallback
            ux = 1; uy = 0; uz = 0; ulen = 1;
        }
        ux/=ulen; uy/=ulen; uz/=ulen;

        // n = cross(p1-p0, p2-p0)
        double vx = p2.x() - p0.x();
        double vy = p2.y() - p0.y();
        double vz = p2.z() - p0.z();
        double nx = uy * vz - uz * vy;
        double ny = uz * vx - ux * vz;
        double nz = ux * vy - uy * vx;
        double nlen = std::sqrt(nx*nx + ny*ny + nz*nz);
        if (nlen <= eps) {
            // fallback normal
            nx = 0; ny = 0; nz = 1; nlen = 1;
        }
        nx/=nlen; ny/=nlen; nz/=nlen;

        // v = n x u
        double vx2 = ny * uz - nz * uy;
        double vy2 = nz * ux - nx * uz;
        double vz2 = nx * uy - ny * ux;
        // no need to normalize v2 because u and n both normalized and orthogonal -> v2 normalized

        // now project each point q to coords (dot(q-p0,u), dot(q-p0,v))
        std::vector<std::pair<Point2_2D, L3Address>> pts2;
        pts2.reserve(pts.size());
        for (const auto &pp : pts) {
            double qx = pp.first.x() - p0.x();
            double qy = pp.first.y() - p0.y();
            double qz = pp.first.z() - p0.z();
            double coord_u = qx*ux + qy*uy + qz*uz;
            double coord_v = qx*vx2 + qy*vy2 + qz*vz2;
            pts2.emplace_back(Point2_2D(coord_u, coord_v), pp.second);
        }

        // Build dt2
        for (const auto &pp : pts2) {
            auto vh = insertVertexWithInfo2(dt2, pp.first, pp.second);
            if (vh != Delaunay2::Vertex_handle()) addrToVh2.emplace(pp.second, vh);
        }

        // if dt2 ended up dimension < 2 (collinear), fallback to pairwise neighbors handled by smallSetPoints3
        if (dt2.dimension() < 2) {
            smallSetPoints3.clear();
            for (const auto &pp : pts) smallSetPoints3.emplace_back(pp);
            // clear dt2 and mapping to keep deterministic behavior
            dt2.clear();
            addrToVh2.clear();
            using2D = false;
            currentDimension = 1;
        }
    }

    // get DT neighbors set for self (unsafe means returns set even if using smallset)
    std::set<L3Address> getDTNeighborsSetUnsafe(const L3Address &selfAddr) const {
        std::set<L3Address> out;
        auto vec = getDTNeighbors(selfAddr);
        for (const auto &a : vec) out.insert(a);
        return out;
    }

    static std::vector<int> greedy_set_cover_ints(int universeSize, const std::vector<std::uint64_t> &simplexMasks) {
        std::vector<int> chosen;
        std::vector<char> covered(simplexMasks.size(), 0);
        int remaining = (int)simplexMasks.size();

        while (remaining > 0) {
            int bestIdx = -1;
            int bestCover = -1;
            for (int j = 0; j < universeSize; ++j) {
                int cnt = 0;
                std::uint64_t bit = (std::uint64_t)1 << j;
                for (size_t si = 0; si < simplexMasks.size(); ++si) {
                    if (!covered[si] && (simplexMasks[si] & bit)) ++cnt;
                }
                if (cnt > bestCover) { bestCover = cnt; bestIdx = j; }
            }
            if (bestIdx < 0 || bestCover == 0) break; // nothing more to cover
            chosen.push_back(bestIdx);
            std::uint64_t chosenBit = (std::uint64_t)1 << bestIdx;
            for (size_t si = 0; si < simplexMasks.size(); ++si) {
                if (!covered[si] && (simplexMasks[si] & chosenBit)) {
                    covered[si] = 1;
                    --remaining;
                }
            }
        }
        return chosen;
    }

    // Attempt to safely move the vertex 'vh' (which corresponds to 'addr') to point p.
    // Returns true if operation succeeded and addrToVh3 has been updated correctly.
    // Returns false if the caller should give up (e.g., mark for rebuild).
    bool safeMoveVertex3(const L3Address &addr,
                                    Delaunay3::Vertex_handle vh,
                                    const Point3_3D &p)
    {
        // basic sanity
        if (vh == Delaunay3::Vertex_handle()) {
            EV_WARN << "safeMoveVertex3: invalid vertex handle for " << addr << "\n";
            return false;
        }
        if (dt3.is_infinite(vh)) {
            EV_WARN << "safeMoveVertex3: vertex is infinite for " << addr << " -> cannot move\n";
            return false;
        }

        // no-op if already equal (within eps)
        if (pointsEqual3(vh->point(), p)) {
            return true;
        }

        // check for collision: is there another vertex at target location?
        // If nearest vertex at p is not vh and is within eps -> collision, avoid move.
        {
            auto nearest = dt3.nearest_vertex(p);
            if (nearest != Delaunay3::Vertex_handle() && nearest != vh) {
                if (pointsEqual3(nearest->point(), p)) {
                    EV_WARN << "safeMoveVertex3: target location collides with existing vertex (addr="
                            << nearest->info() << ") for addr=" << addr << " -> defer rebuild\n";
                    return false;
                }
            }
        }

    #if CGAL_VERSION_NR >= 1050300000
        // Try move under try/catch (some CGAL may throw, but note: segfault won't be caught)
        try {
            dt3.move(vh, p);

            // After move, update mapping: find vertex whose info()==addr.
            // Ideally vh remains valid and vh->info()==addr, but to be safe refresh.
            if (vh->info() == addr) {
                addrToVh3[addr] = vh;
                return true;
            }

            // Otherwise search nearby vertices for one with matching info.
            bool found = false;
            for (auto vit = dt3.finite_vertices_begin(); vit != dt3.finite_vertices_end(); ++vit) {
                if (!dt3.is_infinite(vit) && vit->info() == addr) {
                    addrToVh3[addr] = vit;
                    found = true;
                    break;
                }
            }
            if (!found) {
                EV_WARN << "safeMoveVertex3: moved vertex but handle with info()==addr not found -> request rebuild\n";
                return false;
            }
            return true;
        } catch (const std::exception &ex) {
            EV_WARN << "safeMoveVertex3: dt3.move threw exception for " << addr << " : " << ex.what() << " -> fallback\n";
            // fall through to remove+insert fallback
        } catch (...) {
            EV_WARN << "safeMoveVertex3: dt3.move unknown exception for " << addr << " -> fallback\n";
        }
    #endif

        // Fallback: try remove + insert, but check validity first
        try {
            // As a precaution, verify vh still appears in finite_vertices (avoid removing an invalid handle)
            bool stillPresent = false;
            for (auto vit = dt3.finite_vertices_begin(); vit != dt3.finite_vertices_end(); ++vit) {
                // compare by pointer to the underlying vertex object to avoid relying on operator==
                if (&*vit == &*vh) { stillPresent = true; break; }

                // alternatively (safer if info() is unique), you can use:
                // if (vit->info() == vh->info() && pointsEqual3(vit->point(), vh->point())) { stillPresent = true; break; }
            }
            if (!stillPresent) {
                EV_WARN << "safeMoveVertex3: vertex handle not present anymore for " << addr << " -> rebuild\n";
                return false;
            }

            // remove vh
            dt3.remove(vh);
            addrToVh3.erase(addr);

            // insert new vertex
            auto newVh = dt3.insert(p);
            if (newVh == Delaunay3::Vertex_handle()) {
                EV_WARN << "safeMoveVertex3: fallback insert failed for " << addr << " -> rebuild\n";
                return false;
            }
            newVh->info() = addr;
            addrToVh3[addr] = newVh;
            return true;
        } catch (const std::exception &ex) {
            EV_WARN << "safeMoveVertex3: remove/insert fallback threw for " << addr << " : " << ex.what() << " -> rebuild\n";
            return false;
        } catch (...) {
            EV_WARN << "safeMoveVertex3: remove/insert fallback unknown error for " << addr << " -> rebuild\n";
            return false;
        }

        // unreachable
        return false;
    }


    // helper: insert point into 3D triangulation and attach info
    inline Delaunay3::Vertex_handle insertVertexWithInfo3(Delaunay3 &dt, const Point3_3D &p, const L3Address &addr) {
        Delaunay3::Vertex_handle vh = dt.insert(p);
        if (vh != Delaunay3::Vertex_handle()) {
            // attach info (overwrite existing info if any)
            vh->info() = addr;
        } else {
            // insertion returned null handle (very rare) - nothing to do
        }
        return vh;
    }

    // helper: insert point into 2D triangulation and attach info (symmetry/consistency)
    inline Delaunay2::Vertex_handle insertVertexWithInfo2(Delaunay2 &dt, const Point2_2D &p, const L3Address &addr) {
        Delaunay2::Vertex_handle vh = dt.insert(p);
        if (vh != Delaunay2::Vertex_handle()) vh->info() = addr;
        return vh;
    }


    // find a current vertex handle in dt3 with info()==addr, or return null handle
    Delaunay3::Vertex_handle findVertexHandleByInfo(const L3Address &addr) const {
        for (auto vit = dt3.finite_vertices_begin(); vit != dt3.finite_vertices_end(); ++vit) {
            if (!dt3.is_infinite(vit) && vit->info() == addr) {
                return vit;
            }
        }
        return Delaunay3::Vertex_handle();
    }

};



#endif /* MYSRC_ROUTING_DTMANAGER_H_ */
