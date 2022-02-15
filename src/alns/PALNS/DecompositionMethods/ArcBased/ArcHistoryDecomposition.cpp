//
// Created by Alberto Santini on 08/06/2021.
//

#include "ArcHistoryDecomposition.h"

std::vector<ArcBasedDecomposition::Arc> ArcHistoryDecomposition::getArcsToFix(const CVRPSolution& sol) const {
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

    std::sort(arcs.begin(), arcs.end(), [this](const Arc& a1, const Arc& a2) -> bool { return afreq[a1.o - 1][a1.d - 1] > afreq[a2.o - 1][a2.d - 1]; });

    return arcs;
}
