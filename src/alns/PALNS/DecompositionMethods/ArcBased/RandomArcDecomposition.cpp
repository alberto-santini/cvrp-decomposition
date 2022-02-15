//
// Created by Alberto Santini on 05/06/2021.
//

#include "RandomArcDecomposition.h"

vector<ArcBasedDecomposition::Arc> RandomArcDecomposition::getArcsToFix(const CVRPSolution& sol) const {
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

    std::random_shuffle(arcs.begin(), arcs.end());

    return arcs;
}
