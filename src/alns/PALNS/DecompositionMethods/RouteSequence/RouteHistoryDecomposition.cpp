//
// Created by Alberto Santini on 04/06/2021.
//

#include "RouteHistoryDecomposition.h"
#include "../../hgs-TV/Decomposition/ClusteringUtils/KMedoids.h"
#include <limits>

RouteSequenceDecomposition::SPList RouteHistoryDecomposition::getSubproblems(const CVRPSolution& mpSol) const {
    using namespace genvrp::kmedoids;

    const std::size_t k = std::ceil(i.getN() / par.m_iDecompositionSize);

#ifndef BATCH_MODE
    std::cout << par.outputPrefix << "[Decomposition] Route history clustering decomposition: starting (recursion level: " << par.m_iDecompositionLevel
              << ")\n";
    std::cout << par.outputPrefix << "[Decomposition] Number of clusters: k = " << k << "\n";
#endif

    if(k <= 1) {
#ifndef BATCH_MODE
        std::cout << par.outputPrefix << "[Decomposition] k too small: skipping!\n";
#endif

        return {};
    }

    std::vector<int> emptyRoutes;
    for(auto r = 0; r < mpSol.getNRoutes(); ++r) {
        if(mpSol.getRoute(r).empty()) {
            emptyRoutes.push_back(r);
        }
    }

    std::vector<std::vector<double>> routeDistances(mpSol.getNRoutes(), std::vector<double>(mpSol.getNRoutes(), 0.0));
    for(auto r1 = 0; r1 < mpSol.getNRoutes(); ++r1) {
        const auto& route1 = mpSol.getRoute(r1);

        for(auto r2 = r1 + 1; r2 < mpSol.getNRoutes(); ++r2) {
            const auto& route2 = mpSol.getRoute(r2);
            double score = 0;

            for(auto n1 : route1.getNodes()) {
                if(n1 == 0 || n1 == i.getN() + 1) {
                    continue;
                }

                for(auto n2 : route2.getNodes()) {
                    if(n2 == 0 || n2 == i.getN() + 1) {
                        continue;
                    }

                    score += vfreq[n1 - 1][n2 - 1];
                }
            }

            if(score == 0.0) {
                routeDistances[r1][r2] = std::numeric_limits<double>::max();
            } else {
                routeDistances[r1][r2] = 1.0 / score;
            }
        }
    }

    const auto clustering = kMedoids(k, routeDistances, emptyRoutes);

#ifndef BATCH_MODE
    std::cout << par.outputPrefix << "[Decomposition] Cluster sizes: ";
    for(const auto& c : clustering) {
        std::cout << c.size() << " ";
    }
    std::cout << "\n";
#endif

    std::vector<std::unique_ptr<SubProblem>> sps;

    for(const auto& cluster : clustering) {
        if(!cluster.empty()) {
            sps.push_back(std::make_unique<SubProblem>(par, i, mpSol, cluster, elapsedSec));
        }
    }

    return sps;
}
