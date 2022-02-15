//
// Created by Alberto Santini on 04/06/2021.
//

#ifndef DECOMPOSITION_ROUTEHISTORYDECOMPOSITION_H
#define DECOMPOSITION_ROUTEHISTORYDECOMPOSITION_H

#include "../RouteSequenceDecomposition.h"
#include "../../../PALNS-CVRP/CVRPInstance.h"
#include "../../../PALNS-CVRP/CVRPSolution.h"
#include <vector>

struct RouteHistoryDecomposition : public RouteSequenceDecomposition {
  const std::vector<std::vector<double>>& vfreq;

  RouteHistoryDecomposition(const PALNSParameters& par, const CVRPInstance& i, const std::vector<std::vector<double>>& vfreq) : RouteSequenceDecomposition{par, i}, vfreq{vfreq} {}
  RouteSequenceDecomposition::SPList getSubproblems(const CVRPSolution& mpSol) const override;
};


#endif //DECOMPOSITION_ROUTEHISTORYDECOMPOSITION_H
