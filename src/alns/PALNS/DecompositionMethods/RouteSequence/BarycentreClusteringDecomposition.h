//
// Created by Alberto Santini on 04/06/2021.
//

#ifndef DECOMPOSITION_BARYCENTRECLUSTERINGDECOMPOSITION_H
#define DECOMPOSITION_BARYCENTRECLUSTERINGDECOMPOSITION_H

#include "../RouteSequenceDecomposition.h"
#include "../../../PALNS-CVRP/CVRPInstance.h"
#include "../../../PALNS-CVRP/CVRPSolution.h"

struct BarycentreClusteringDecomposition : public RouteSequenceDecomposition {
  BarycentreClusteringDecomposition(const PALNSParameters& par, const CVRPInstance& i) : RouteSequenceDecomposition{par, i} {}
  RouteSequenceDecomposition::SPList getSubproblems(const CVRPSolution& mpSol) const override;
};


#endif //DECOMPOSITION_BARYCENTRECLUSTERINGDECOMPOSITION_H
