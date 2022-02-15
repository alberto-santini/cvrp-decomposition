//
// Created by Alberto Santini on 04/06/2021.
//

#include "BarycentreQuadrantDecomposition.h"
#include <array>
#include <iostream>
#include <vector>

RouteSequenceDecomposition::SPList BarycentreQuadrantDecomposition::getSubproblems(const CVRPSolution& mpSol) const {
#ifndef BATCH_MODE
    std::cout << par.outputPrefix << "[Decomposition] Barycentre swipe decomposition: starting (recursion level: " << par.m_iDecompositionLevel << ")\n";
#endif

    if(par.m_iDecompositionLevel > 1) {
#ifndef BATCH_MODE
        std::cout << par.outputPrefix << "[Decomposition] Skipping because quadrant decomposition can only be applied once\n";
#endif

        return {};
    }

    SPList sps{};

    std::array<std::vector<int>, 4> routesByQuadrant;

    for(auto i = 0; i < mpSol.getNRoutes(); ++i) {
        const auto& route = mpSol.getRoute(i);

        if(route.empty()) {
            continue;
        }

        double totX = 0, totY = 0;

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

        if(totX >= 0) {
            if(totY >= 0) {
                routesByQuadrant[0].push_back(i);
            } else {
                routesByQuadrant[1].push_back(i);
            }
        } else {
            if(totY >= 0) {
                routesByQuadrant[3].push_back(i);
            } else {
                routesByQuadrant[2].push_back(i);
            }
        }
    }

    for(const auto& qRoutes : routesByQuadrant) {
        if(!qRoutes.empty()) {
            sps.push_back(std::make_unique<SubProblem>(par, i, mpSol, qRoutes, elapsedSec));
        }
    }

#ifndef BATCH_MODE
    std::cout << par.outputPrefix << "[Decomposition] Found " << sps.size() << " non-empty quadrants\n";
#endif

    return sps;
}
