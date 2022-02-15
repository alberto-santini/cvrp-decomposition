//
// Created by Alberto Santini on 09/06/2021.
//

#ifndef DECOMPOSITION_PATHHISTORYDECOMPOSITION_H
#define DECOMPOSITION_PATHHISTORYDECOMPOSITION_H

#include "../ArcBasedDecomposition.h"
#include "../../PALNS.h"
#include "../../PALNSParameters.h"
#include "../../../PALNS-CVRP/CVRPInstance.h"
#include "../../../PALNS-CVRP/CVRPSolution.h"
#include <vector>

struct PathHistoryDecomposition : public ArcBasedDecomposition {
  const PALNS::PathFrequencyMatrix& pfreq;

  PathHistoryDecomposition(const PALNSParameters& par, const CVRPInstance& i, const PALNS::PathFrequencyMatrix& pfreq) : ArcBasedDecomposition{par, i}, pfreq{pfreq} {}
  std::vector<Arc> getArcsToFix(const CVRPSolution& sol) const override;
};

#endif //DECOMPOSITION_PATHHISTORYDECOMPOSITION_H
