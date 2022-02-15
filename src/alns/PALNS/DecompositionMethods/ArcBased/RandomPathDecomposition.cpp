#include "RandomPathDecomposition.h"
#include <cassert>
#include <random>

vector<ArcBasedDecomposition::Arc> RandomPathDecomposition::getArcsToFix(const CVRPSolution& sol) const {
    std::vector<Path> paths;

    for(auto r = 0; r < sol.getNRoutes(); ++r) {
        const auto& route = sol.getRoute(r);
        const auto& nodes = route.getNodes();

        if(route.size() < par.m_iPathLength + 2) {
            continue;
        }
        assert(nodes.front() == 0);
        assert(nodes.back() == i.getN() + 1);

        for(auto it = nodes.begin() + 1; it <= nodes.end() - par.m_iPathLength - 1; ++it) {
            paths.emplace_back(it, it + par.m_iPathLength);
        }
    }

    std::random_shuffle(paths.begin(), paths.end());

    std::vector<Arc> arcs;

    for(const auto& path : paths) {
        for(auto id = 0u; id < path.size() - 2; ++id) {
            const Arc arc = {path[id], path[id + 1], i.getDist(path[id], path[id + 1])};

            if(std::find(arcs.begin(), arcs.end(), arc) == arcs.end()) {
                arcs.push_back(arc);
            }
        }
    }

    return arcs;
}