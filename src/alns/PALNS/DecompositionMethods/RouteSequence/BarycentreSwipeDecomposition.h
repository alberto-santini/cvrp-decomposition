//
// Created by Alberto Santini on 03/06/2021.
//

#ifndef DECOMPOSITION_BARYCENTRESWIPEDECOMPOSITION_H
#define DECOMPOSITION_BARYCENTRESWIPEDECOMPOSITION_H

#include "../RouteSequenceDecomposition.h"
#include "../../../PALNS-CVRP/CVRPInstance.h"
#include "../../../PALNS-CVRP/CVRPSolution.h"

struct BarycentreSwipeDecomposition : public RouteSequenceDecomposition {
  BarycentreSwipeDecomposition(const PALNSParameters& par, const CVRPInstance& i) : RouteSequenceDecomposition{par, i} {}

  RouteSequenceDecomposition::SPList getSubproblems(const CVRPSolution& mpSol) const override;
};


#endif //DECOMPOSITION_BARYCENTRESWIPEDECOMPOSITION_H
