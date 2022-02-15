//
// Created by Alberto Santini on 29/05/2021.
//

#ifndef DECOMPOSITION_RANDOMROUTEDECOMPOSITION_H
#define DECOMPOSITION_RANDOMROUTEDECOMPOSITION_H

#include "../RouteSequenceDecomposition.h"
#include "../../../PALNS-CVRP/CVRPInstance.h"
#include "../../../PALNS-CVRP/CVRPSolution.h"

struct RandomRouteDecomposition : public RouteSequenceDecomposition {
  RandomRouteDecomposition(const PALNSParameters& par, const CVRPInstance& i) : RouteSequenceDecomposition{par, i} {}

  RouteSequenceDecomposition::SPList getSubproblems(const CVRPSolution& mpSol) const override;
};

#endif //DECOMPOSITION_RANDOMROUTEDECOMPOSITION_H
