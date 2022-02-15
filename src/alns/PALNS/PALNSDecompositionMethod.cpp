//
// Created by Alberto Santini on 02/06/2021.
//

#include "PALNSDecompositionMethod.h"

#include "DecompositionMethods/ArcBased/ArcHistoryDecomposition.h"
#include "DecompositionMethods/ArcBased/CostArcDecomposition.h"
#include "DecompositionMethods/ArcBased/CostPathDecomposition.h"
#include "DecompositionMethods/ArcBased/PathHistoryDecomposition.h"
#include "DecompositionMethods/ArcBased/RandomArcDecomposition.h"
#include "DecompositionMethods/ArcBased/RandomPathDecomposition.h"
#include "DecompositionMethods/RouteSequence/BarycentreClusteringDecomposition.h"
#include "DecompositionMethods/RouteSequence/BarycentreQuadrantDecomposition.h"
#include "DecompositionMethods/RouteSequence/BarycentreSwipeDecomposition.h"
#include "DecompositionMethods/RouteSequence/RandomRouteDecomposition.h"
#include "DecompositionMethods/RouteSequence/RouteHistoryDecomposition.h"

void DecompositionMethod::mergeMpPartialSolutionToMpSolution(const CVRPSolution& mpParSol, CVRPSolution& mpSol) const {
    for(auto i = 0; i < mpParSol.getNRoutes(); ++i) {
        const auto& route = mpParSol.getRoute(i);
        bool inserted = false;

        for(auto j = 0; j < mpSol.getNRoutes(); ++j) {
            auto& m_route = mpSol.getRoute(j);

            if(m_route.empty()) {
                mpSol.setRouteNodes(j, route.getNodes());
                inserted = true;
                break;
            }
        }

        if(!inserted) {
            std::cerr << par.outputPrefix << "Cannot find an empty MP route to insert the current route from the SP, creating a new one...\n";
            mpSol.createRoute();
            mpSol.setRouteNodes(mpSol.getNRoutes() - 1, route.getNodes());
        }
    }
}

std::unique_ptr<DecompositionMethod> get_decomposition_method(const PALNS& p) {
    const PALNSParameters& m_par = p.m_par;
    const CVRPInstance& m_instance = p.m_instance;

    switch(m_par.m_iDecompositionMethod) {
    case 0: return nullptr;

    case 1: return std::unique_ptr<DecompositionMethod>(new RandomRouteDecomposition(m_par, m_instance));

    case 2: return std::unique_ptr<DecompositionMethod>(new BarycentreSwipeDecomposition(m_par, m_instance));

    case 3: return std::unique_ptr<DecompositionMethod>(new BarycentreQuadrantDecomposition(m_par, m_instance));

    case 4: return std::unique_ptr<DecompositionMethod>(new BarycentreClusteringDecomposition(m_par, m_instance));

    case 5: return std::unique_ptr<DecompositionMethod>(new RouteHistoryDecomposition(m_par, m_instance, p.vfreq));

    case 6: return std::unique_ptr<DecompositionMethod>(new RandomArcDecomposition(m_par, m_instance));

    case 7: return std::unique_ptr<DecompositionMethod>(new RandomPathDecomposition(m_par, m_instance));

    case 8: return std::unique_ptr<DecompositionMethod>(new CostArcDecomposition(m_par, m_instance));

    case 9: return std::unique_ptr<DecompositionMethod>(new CostPathDecomposition(m_par, m_instance));

    case 10: return std::unique_ptr<DecompositionMethod>(new ArcHistoryDecomposition(m_par, m_instance, p.afreq));

    case 11: return std::unique_ptr<DecompositionMethod>(new PathHistoryDecomposition(m_par, m_instance, p.pfreq));

    default: throw std::runtime_error("Unknown decomposition method!");
    }
}