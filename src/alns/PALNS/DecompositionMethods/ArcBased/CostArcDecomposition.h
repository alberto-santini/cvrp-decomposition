//
// Created by Alberto Santini on 08/06/2021.
//

#ifndef DECOMPOSITION_COSTARCDECOMPOSITION_H
#define DECOMPOSITION_COSTARCDECOMPOSITION_H

#include "../ArcBasedDecomposition.h"
#include "../../PALNSParameters.h"
#include "../../../PALNS-CVRP/CVRPInstance.h"
#include "../../../PALNS-CVRP/CVRPSolution.h"
#include <vector>

struct CostArcDecomposition : public ArcBasedDecomposition {
  CostArcDecomposition(const PALNSParameters& par, const CVRPInstance& i) : ArcBasedDecomposition{par, i} {}
  std::vector<Arc> getArcsToFix(const CVRPSolution& sol) const override;
};


#endif //DECOMPOSITION_COSTARCDECOMPOSITION_H
