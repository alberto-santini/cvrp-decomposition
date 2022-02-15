//
// Created by Alberto Santini on 04/06/2021.
//

#include "BarycentreClusteringDecomposition.h"
#include "../../hgs-TV/Decomposition/ClusteringUtils/KMeans.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

RouteSequenceDecomposition::SPList BarycentreClusteringDecomposition::getSubproblems(const CVRPSolution& mpSol) const {
    using namespace genvrp::kmeans;
    struct Barycentre {
        double x;
        double y;
    };

    const std::size_t k = std::ceil(i.getN() / par.m_iDecompositionSize);

#ifndef BATCH_MODE
    std::cout << par.outputPrefix << "[Decomposition] Barycentre clustering decomposition: starting (recursion level: " << par.m_iDecompositionLevel << ")\n";
    std::cout << par.outputPrefix << "[Decomposition] Number of clusters: k = " << k << "\n";
#endif

    if(k <= 1) {
#ifndef BATCH_MODE
        std::cout << par.outputPrefix << "[Decomposition] k too small: skipping!\n";
#endif

        return {};
    }

    std::vector<int> emptyRoutes;
    std::vector<Barycentre> barycentres;
    for(auto r = 0; r < mpSol.getNRoutes(); ++r) {
        if(mpSol.getRoute(r).empty()) {
            emptyRoutes.push_back(r);
            barycentres.push_back({i.getCoord(0).getX(), i.getCoord(0).getX()});
        } else {
            double sumX = 0, sumY = 0;
            const auto& route = mpSol.getRoute(r);
            for(auto node : route.getNodes()) {
                sumX += i.getCoord(node).getX();
                sumY += i.getCoord(node).getY();
            }
            barycentres.push_back({sumX / route.size(), sumY / route.size()});
        }
    }

    const auto clustering = kMeans(k, barycentres, emptyRoutes);

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
