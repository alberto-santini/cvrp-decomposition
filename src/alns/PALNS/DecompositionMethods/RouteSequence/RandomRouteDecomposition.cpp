//
// Created by Alberto Santini on 02/06/2021.
//

#include "RandomRouteDecomposition.h"
#include <numeric>
#include <random>

RouteSequenceDecomposition::SPList RandomRouteDecomposition::getSubproblems(const CVRPSolution& mpSol) const {
#ifndef BATCH_MODE
    std::cout << par.outputPrefix << "[Decomposition] Random route decomposition: starting (recursion level: " << par.m_iDecompositionLevel << ")\n";
#endif

    auto mpRouteIndices = std::vector<int>(mpSol.getNRoutes());
    std::iota(mpRouteIndices.begin(), mpRouteIndices.end(), 0);
    std::random_shuffle(mpRouteIndices.begin(), mpRouteIndices.end());

    return getSubproblemsFromRouteIndices(mpSol, mpRouteIndices);
}
