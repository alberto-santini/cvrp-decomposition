//
// Created by Alberto Santini on 08/06/2021.
//

#include "CostArcDecomposition.h"

std::vector<ArcBasedDecomposition::Arc> CostArcDecomposition::getArcsToFix(const CVRPSolution& sol) const {
    std::vector<Arc> arcs;

    for(auto r = 0; r < sol.getNRoutes(); ++r) {
        const auto& route = sol.getRoute(r);
        const auto& nodes = route.getNodes();

        for(auto j = 0; j < nodes.size() - 1; ++j) {
            if(nodes[j] != 0 && nodes[j + 1] != i.getN() + 1) {
                arcs.push_back({nodes[j], nodes[j + 1], i.getDist(nodes[j], nodes[j + 1])});
            }
        }
    }

    std::sort(arcs.begin(), arcs.end(), [](const Arc& a1, const Arc& a2) -> bool { return a1.cost < a2.cost; });

    return arcs;
}
