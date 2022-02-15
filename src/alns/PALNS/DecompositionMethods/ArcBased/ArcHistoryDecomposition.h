//
// Created by Alberto Santini on 08/06/2021.
//

#ifndef DECOMPOSITION_ARCHISTORYDECOMPOSITION_H
#define DECOMPOSITION_ARCHISTORYDECOMPOSITION_H

#include "../ArcBasedDecomposition.h"
#include "../../PALNSParameters.h"
#include "../../../PALNS-CVRP/CVRPInstance.h"
#include "../../../PALNS-CVRP/CVRPSolution.h"
#include <vector>

struct ArcHistoryDecomposition : public ArcBasedDecomposition {
  const std::vector<std::vector<double>>& afreq;

  ArcHistoryDecomposition(const PALNSParameters& par, const CVRPInstance& i, const std::vector<std::vector<double>>& afreq) : ArcBasedDecomposition{par, i}, afreq{afreq} {}
  std::vector<Arc> getArcsToFix(const CVRPSolution& sol) const override;
};


#endif //DECOMPOSITION_ARCHISTORYDECOMPOSITION_H
