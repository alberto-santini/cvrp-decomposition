//
// Created by Alberto Santini on 04/06/2021.
//

#include "ArcBasedDecomposition.h"
#include "../../PALNS-CVRP/CVRPInstance.h"
#include "../../PALNS-CVRP/CVRPSolution.h"
#include "../../PALNS-CVRP/PALNS-CVRP.h"
#include "../../SRP-Utils/ParFileParser.h"
#include "../SubProblem.h"
#include <algorithm>
#include <cassert>
#include <set>

namespace {
    bool isSubPath(const ArcBasedDecomposition::Path& path, const ArcBasedDecomposition::Arc& arc) {
        if(path.size() < 2ul) {
            return false;
        }

        for(auto pos = 0u; pos < path.size() - 1u; ++pos) {
            // Arc is already part of this path.
            if(arc.o == path[pos] && arc.d == path[pos + 1u]) {
                return true;
            }

            // The reverse arc is already part of this path.
            if(arc.d == path[pos] && arc.o == path[pos + 1u]) {
                return true;
            }
        }

        return false;
    }

    bool isSubPathOfAny(const std::vector<ArcBasedDecomposition::Path>& paths, const ArcBasedDecomposition::Arc& arc) {
        return std::any_of(paths.begin(), paths.end(), [&arc](const auto& path) -> bool { return isSubPath(path, arc); });
    }

    bool attachesToPath(ArcBasedDecomposition::Path& path, const ArcBasedDecomposition::Arc& arc, ArcBasedDecomposition::CustomerSet& custs) {
        if(arc.d == path.front()) {
            // I have a path v1, v2, v3, ...
            // and now I consider an arc (v0, v1)
            // The path needs to become v0, v1, v2, ...
            path.insert(path.begin(), arc.o);
            custs.push_back(arc.o);
            return true;
        } else if(arc.o == path.back()) {
            // I have a path ..., v98, v99, v100
            // and now I consider an arc (v100, 101)
            // The path needs to become ..., v99, v100, v101
            path.push_back(arc.d);
            custs.push_back(arc.d);
            return true;
        }

        return false;
    }

    bool attachesToAnyPath(std::vector<ArcBasedDecomposition::Path>& paths, const ArcBasedDecomposition::Arc& arc, ArcBasedDecomposition::CustomerSet& custs) {
        return std::any_of(paths.begin(), paths.end(), [&arc, &custs](auto& path) -> bool { return attachesToPath(path, arc, custs); });
    }

    void mergePaths(std::vector<ArcBasedDecomposition::Path>& paths, ArcBasedDecomposition::CustomerSet& custs) {
        // Check if there are any paths to merge.
        // For example, if before I had two paths
        // P1: v0, v1, v2
        // P2: v3, v4, v5
        // And the new arc is (v2, v3), I would now have, e.g.:
        // P1: v0, v1, v2, v3
        // P2: v3, v4, v5
        // This needs to be merged into:
        // P1: v0, v1, v2, v3, v4, v5

        for(auto p1 = 0ul; p1 < paths.size(); ++p1) {
            if(paths[p1].empty()) {
                continue;
            }
            bool merged = false;

            for(auto p2 = p1 + 1ul; p2 < paths.size(); ++p2) {
                if(paths[p2].empty()) {
                    continue;
                }

                if(paths[p1].front() == paths[p2].back()) {
                    const auto link = paths[p1].front();

                    paths[p2].insert(paths[p2].end(), paths[p1].begin() + 1, paths[p1].end());
                    paths[p1].clear();

                    // The "linking" customer was counted twice.
                    if(std::count(custs.begin(), custs.end(), link) > 1u) {
                        custs.erase(std::find(custs.begin(), custs.end(), link));
                    }

                    merged = true;
                    break;
                }

                if(paths[p1].back() == paths[p2].front()) {
                    const auto link = paths[p1].back();

                    paths[p1].insert(paths[p1].end(), paths[p2].begin() + 1, paths[p2].end());
                    paths[p2].clear();

                    // The "linking" customer was counted twice.
                    if(std::count(custs.begin(), custs.end(), link) > 1u) {
                        custs.erase(std::find(custs.begin(), custs.end(), link));
                    }

                    merged = true;
                    break;
                }
            }

            if(merged) {
                break;
            }
        }
    }
} // namespace

void ArcBasedDecomposition::addArcToPaths(const Arc& arc, std::vector<ArcBasedDecomposition::Path>& paths,
                                          ArcBasedDecomposition::CustomerSet& custsInPaths) const {
    if(isSubPathOfAny(paths, arc)) {
        return;
    }

    if(attachesToAnyPath(paths, arc, custsInPaths)) {
        mergePaths(paths, custsInPaths);
    } else {
        paths.push_back({arc.o, arc.d});
        custsInPaths.push_back(arc.o);
        custsInPaths.push_back(arc.d);
    }

    paths.erase(std::remove_if(paths.begin(), paths.end(), [](const auto& p) { return p.empty(); }), paths.end());
}

std::pair<ArcBasedDecomposition::CustomerSet, std::vector<ArcBasedDecomposition::Path>> ArcBasedDecomposition::getPaths(const vector<Arc>& arcs) const {
    CustomerSet custs;
    std::vector<Path> paths;

    for(const auto& arc : arcs) {
        // If the arc is already in some path, skip it.
        if(isSubPathOfAny(paths, arc)) {
            continue;
        }

        if(attachesToAnyPath(paths, arc, custs)) {
            // Check if there are any paths to merge.
            mergePaths(paths, custs);
        } else {
            // Start a new path.
            paths.push_back({arc.o, arc.d});
            custs.push_back(arc.o);
            custs.push_back(arc.d);
        }

        // Check there are no duplicates.
        assert(std::set(custs.begin(), custs.end()).size() == custs.size());

        // Clean up empty paths.
        paths.erase(std::remove_if(paths.begin(), paths.end(), [](const auto& p) { return p.empty(); }), paths.end());
    }

    return std::make_pair(custs, paths);
}

CVRPSolution ArcBasedDecomposition::decompose(CVRPSolution sol, double elapsedSec) const {
#ifndef BATCH_MODE
    if(i.getN() <= par.m_iDecompositionSize) {
        std::cerr << par.outputPrefix << "[Decomposition] Skipping arc-based decomposition because of too few customers: " << i.getN() << " vs. "
                  << par.m_iDecompositionSize << " target subproblem size\n";

        return sol;
    }
#endif

    if(par.m_iDecompositionLevel > 1) {
#ifndef BATCH_MODE
        std::cerr << par.outputPrefix << "[Decomposition] Too nested (arc-based decomposition only goes down one level): skipping\n";
        return sol;
#endif
    }

    auto arcs = getArcsToFix(sol);

    std::sort(arcs.begin(), arcs.end());
    arcs.erase(std::unique(arcs.begin(), arcs.end()), arcs.end());
    assert(std::set(arcs.begin(), arcs.end()).size() == arcs.size());

    CustomerSet custsInPaths;
    std::vector<Path> paths;
    std::size_t nArcsUsed = par.m_iDecompositionSize / 2;

    while(paths.size() + (i.getN() - custsInPaths.size()) > par.m_iDecompositionSize && nArcsUsed < arcs.size()) {
        const auto& newArc = arcs[nArcsUsed];
        addArcToPaths(newArc, paths, custsInPaths);
        ++nArcsUsed;
    }

    CustomerSet custsNotInPaths;
    custsNotInPaths.reserve(i.getN() - custsInPaths.size());
    for(auto id = 1; id <= i.getN(); ++id) {
        if(std::find(custsInPaths.begin(), custsInPaths.end(), id) == custsInPaths.end()) {
            custsNotInPaths.push_back(id);
        }
    }

    assert(custsInPaths.size() + custsNotInPaths.size() == i.getN());

    std::map<int, CustomerSet> spToMp;
    std::map<int, int> mpToSp;

    spToMp[0] = {0};
    mpToSp[0] = 0;

    int custNumber = 1;
    int addedCustomers = 0;

    for(const auto& c : custsNotInPaths) {
        spToMp[custNumber] = {c};
        mpToSp[c] = custNumber;
        ++custNumber;
        ++addedCustomers;
    }

    static_assert(std::is_same_v<CustomerSet, Path>);
    for(const auto& path : paths) {
        assert(!path.empty());

        spToMp[custNumber] = path;

        for(const auto& mpCust : path) {
            mpToSp[mpCust] = custNumber;
        }

        ++custNumber;
        addedCustomers += path.size();
    }

    // Make sure all master problem customers are placed somewhere.
    assert(addedCustomers == i.getN());

    spToMp[custNumber] = {i.getN() + 1};
    mpToSp[i.getN() + 1] = {custNumber};

    int sp_n_vertices = spToMp.size() - 2;
    std::vector<Coordinate<double>> sp_coord(sp_n_vertices + 2);
    std::vector<int> sp_demand(sp_n_vertices + 2, 0);
    std::vector<int> sp_service_time(sp_n_vertices + 2, 0);
    const int sp_n_vehicles = i.getNVehicles();
    const int sp_capacity = i.getMaxLoad();
    const int sp_max_duration = i.getMaxDuration();

    sp_coord[0] = i.getCoord(0);
    sp_demand[0] = i.getDemand(0);
    sp_service_time[0] = i.getDemand(0);

    for(const auto& [sp_vertex_id, mp_vertex_sequence] : spToMp) {
        double cX = 0, cY = 0;

        for(auto j = 0; j < mp_vertex_sequence.size(); ++j) {
            const auto& mp_coord = i.getCoord(mp_vertex_sequence[j]);
            cX += mp_coord.getX();
            cY += mp_coord.getY();
            sp_demand[sp_vertex_id] += i.getDemand(mp_vertex_sequence[j]);
            sp_service_time[sp_vertex_id] += i.getServiceTime(mp_vertex_sequence[j]);

            if(j < mp_vertex_sequence.size() - 1) {
                sp_service_time[sp_vertex_id] += i.getDist(mp_vertex_sequence[j], mp_vertex_sequence[j + 1]);
            }
        }

        sp_coord[sp_vertex_id].setX(cX / mp_vertex_sequence.size());
        sp_coord[sp_vertex_id].setY(cY / mp_vertex_sequence.size());
    }

    sp_coord[sp_n_vertices + 1] = i.getCoord(i.getN() + 1);
    sp_demand[sp_n_vertices + 1] = i.getDemand(i.getN() + 1);
    sp_service_time[sp_n_vertices + 1] = i.getServiceTime(i.getN() + 1);

    Matrix<double> sp_dist_matrix{sp_n_vertices + 2};

    for(auto j = 0; j <= sp_n_vertices + 1; ++j) {
        for(auto k = 0; k <= sp_n_vertices + 1; ++k) {
            sp_dist_matrix.setElement(j, k, i.getDist(spToMp[j].back(), spToMp[k].front()));
        }
    }

#ifndef BATCH_MODE
    std::cout << par.outputPrefix << "[Decomposition] Arc-based decomposition sub-problem has " << sp_n_vertices << " customers\n";
#endif

    CVRPInstance sp_instance{};

    sp_instance.setProblem(sp_n_vertices, sp_coord, sp_demand, sp_service_time, sp_n_vehicles, false, sp_capacity, sp_max_duration);
    sp_instance.setDistMatrix(sp_dist_matrix);

    CVRPSolution sp_sol{sp_instance};
    ParFileParser pfp;

    pfp.loadParFile(SubProblem::cvrp_specific_params);
    pfp.setParameter("MaxIter", "10000");
    pfp.setParameter("DecompositionMethod", std::to_string(par.m_iDecompositionMethod));
    pfp.setParameter("DecompositionIters", std::to_string(par.m_iDecompositionIters));
    pfp.setParameter("DecompositionSize", std::to_string(par.m_iDecompositionSize));
    pfp.setParameter("DecompositionLevel", std::to_string(par.m_iDecompositionLevel));
    pfp.setParameter("OutputPrefix", par.outputPrefix);

    double timeoutSec;
    pfp.getDoubleParameter("TimeoutSec", timeoutSec);
    timeoutSec -= elapsedSec;
    pfp.setParameter("TimeoutSec", std::to_string(timeoutSec));

    PALNS_CVRP solver;
    solver.go(sp_sol, 1 /*thread*/, 1 /*retry*/, false /*print header*/, pfp);

    CVRPSolution new_mp_sol{i};
    new_mp_sol.clear(sp_sol.getNRoutes());

    for(auto j = 0; j < sp_sol.getNRoutes(); ++j) {
        const auto& sp_route = sp_sol.getRoute(j);
        std::vector<int> mpRouteNodes;

        for(auto sp_node : sp_route.getNodes()) {
            for(auto mp_node : spToMp[sp_node]) {
                mpRouteNodes.push_back(mp_node);
            }
        }

        new_mp_sol.setRouteNodes(j, mpRouteNodes);
    }

    return new_mp_sol;
}
