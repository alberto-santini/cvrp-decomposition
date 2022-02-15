//
// Created by Alberto Santini on 08/06/2021.
//

#ifndef DECOMPOSITION_COSTPATHDECOMPOSITION_H
#define DECOMPOSITION_COSTPATHDECOMPOSITION_H

#include "../ArcBasedDecomposition.h"
#include "../../PALNSParameters.h"
#include "../../../PALNS-CVRP/CVRPInstance.h"
#include "../../../PALNS-CVRP/CVRPSolution.h"
#include <vector>

struct CostPathDecomposition : public ArcBasedDecomposition {
  CostPathDecomposition(const PALNSParameters& par, const CVRPInstance& i) : ArcBasedDecomposition{par, i} {}
  std::vector<Arc> getArcsToFix(const CVRPSolution& sol) const override;
};


#endif //DECOMPOSITION_COSTPATHDECOMPOSITION_H
