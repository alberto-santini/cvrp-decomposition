//
// Created by Alberto Santini on 02/06/2021.
//

#include "RouteSequenceDecomposition.h"

RouteSequenceDecomposition::SPList RouteSequenceDecomposition::getSubproblemsFromRouteIndices(const CVRPSolution& mpSol,
                                                                                              const vector<int>& routeIndices) const {
    auto sps = SPList{};
    auto currentCustomersInSP = 0;
    auto currentRoutesInSP = std::vector<int>{};

    for(const auto id : routeIndices) {
        const auto& mpRoute = mpSol.getRoute(id);

        // Use -2 to remove the depot and the returning depot
        if(currentCustomersInSP + mpRoute.size() - 2 > par.m_iDecompositionSize) {
            sps.push_back(std::make_unique<SubProblem>(par, i, mpSol, currentRoutesInSP, this->elapsedSec));
            currentCustomersInSP = 0;
            currentRoutesInSP.clear();
        }

        currentCustomersInSP += mpRoute.size() - 2;
        currentRoutesInSP.push_back(id);
    }

    if(!currentRoutesInSP.empty()) {
        sps.push_back(std::make_unique<SubProblem>(par, i, mpSol, currentRoutesInSP, this->elapsedSec));
    }

    return sps;
}

CVRPSolution RouteSequenceDecomposition::decompose(CVRPSolution sol, double elapsedSec) const {
    if(i.getN() - 1 <= par.m_iDecompositionSize) {
#ifndef BATCH_MODE
        std::cerr << par.outputPrefix << "[Decomposition] Problem has too few customers: " << i.getN() - 1 << " vs. " << par.m_iDecompositionSize
                  << " min required\n";
#endif

        return sol;
    }

    if(par.m_iDecompositionLevel > 5) {
#ifndef BATCH_MODE
        std::cerr << par.outputPrefix << "[Decomposition] Too nested (> 5 levels): skipping\n";
#endif

        return sol;
    }

    this->elapsedSec = elapsedSec;

    auto subproblems = getSubproblems(sol);
#ifndef BATCH_MODE
    std::cerr << par.outputPrefix << "[Decomposition] Created " << subproblems.size() << " subproblems\n";
#endif

    if(subproblems.size() < 1u) {
#ifndef BATCH_MODE
        std::cerr << par.outputPrefix << "[Decomposition] Too few subproblems (" << subproblems.size() << "), skipping\n";
#endif
        return sol;
    }

    CVRPSolution new_mp_sol{i};
    new_mp_sol.clear(sol.getNRoutes());

    for(const auto& sp : subproblems) {
        const std::optional<CVRPSolution> sp_sol = sp->solve();

        if(sp_sol) {
            mergeMpPartialSolutionToMpSolution(*sp_sol, new_mp_sol);
        }
    }

    return new_mp_sol;
}
