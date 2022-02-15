//
// Created by Alberto Santini on 05/06/2021.
//

#ifndef DECOMPOSITION_RANDOMARCDECOMPOSITION_H
#define DECOMPOSITION_RANDOMARCDECOMPOSITION_H

#include "../ArcBasedDecomposition.h"
#include "../../PALNSParameters.h"
#include "../../../PALNS-CVRP/CVRPInstance.h"
#include "../../../PALNS-CVRP/CVRPSolution.h"
#include <vector>

struct RandomArcDecomposition : public ArcBasedDecomposition {
  RandomArcDecomposition(const PALNSParameters& par, const CVRPInstance& i) : ArcBasedDecomposition{par, i} {}
  std::vector<Arc> getArcsToFix(const CVRPSolution& sol) const override;
};

#endif //DECOMPOSITION_RANDOMARCDECOMPOSITION_H
