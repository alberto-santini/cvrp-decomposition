//
// Created by Alberto Santini on 03/06/2021.
//

#include "BarycentreSwipeDecomposition.h"
#include <cmath>
#include <utility>

RouteSequenceDecomposition::SPList BarycentreSwipeDecomposition::getSubproblems(const CVRPSolution& mpSol) const {
#ifndef BATCH_MODE
    std::cout << par.outputPrefix << "[Decomposition] Barycentre swipe decomposition: starting (recursion level: " << par.m_iDecompositionLevel << ")\n";
#endif

    auto mpRouteAngleIndex = std::vector<std::pair<double, int>>(mpSol.getNRoutes());

    for(auto i = 0; i < mpSol.getNRoutes(); ++i) {
        double totX = 0, totY = 0;
        const auto& route = mpSol.getRoute(i);

        for(auto node : route.getNodes()) {
            if(node == 0 || node == mpSol.getInstance().getN() + 1) {
                continue;
            }
            const auto coords = mpSol.getInstance().getCoord(node);

            totX += coords.getX();
            totY += coords.getY();
        }

        totX /= route.size();
        totY /= route.size();

        mpRouteAngleIndex[i] = {std::atan2(totY, totX), i};
    }

    // Sort routes by barycentre angle
    std::sort(mpRouteAngleIndex.begin(), mpRouteAngleIndex.end());

    // Select a random starting point and rotate.
    const auto start = std::rand() % mpRouteAngleIndex.size();
    std::rotate(mpRouteAngleIndex.begin(), mpRouteAngleIndex.begin() + start, mpRouteAngleIndex.end());

    auto mpRouteIndices = std::vector<int>(mpRouteAngleIndex.size());
    for(auto i = 0u; i < mpRouteIndices.size(); ++i) {
        mpRouteIndices[i] = mpRouteAngleIndex[i].second;
    }

    return getSubproblemsFromRouteIndices(mpSol, mpRouteIndices);
}
