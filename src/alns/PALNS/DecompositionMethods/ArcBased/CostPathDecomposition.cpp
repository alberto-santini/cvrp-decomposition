//
// Created by Alberto Santini on 08/06/2021.
//

#include "CostPathDecomposition.h"
#include <cassert>
#include <utility>

std::vector<ArcBasedDecomposition::Arc> CostPathDecomposition::getArcsToFix(const CVRPSolution& sol) const {
    using CostPath = std::pair<double, Path>;
    std::vector<CostPath> paths;

    for(auto r = 0; r < sol.getNRoutes(); ++r) {
        const auto& route = sol.getRoute(r);
        const auto& nodes = route.getNodes();

        if(route.size() < par.m_iPathLength + 2) {
            continue;
        }
        assert(nodes.front() == 0);
        assert(nodes.back() == i.getN() + 1);

        for(auto it = nodes.begin() + 1; it <= nodes.end() - par.m_iPathLength - 1; ++it) {
            double cost = 0;

            for(auto it2 = it; it2 < it + par.m_iPathLength - 1; ++it2) {
                cost += i.getDist(*it, *it2);
            }

            paths.emplace_back(cost, std::vector<int>(it, it + par.m_iPathLength));
        }
    }

    std::sort(paths.begin(), paths.end());

    std::vector<Arc> arcs;

    for(const auto& cpath : paths) {
        const auto& path = cpath.second;

        for(auto id = 0u; id < path.size() - 2; ++id) {
            const Arc arc = {path[id], path[id + 1], i.getDist(path[id], path[id + 1])};

            if(std::find(arcs.begin(), arcs.end(), arc) == arcs.end()) {
                arcs.push_back(arc);
            }
        }
    }

    return arcs;
}
