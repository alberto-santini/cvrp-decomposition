//
// Created by Alberto Santini on 29/05/2021.
//

#ifndef DECOMPOSITION_ROUTESEQUENCEDECOMPOSITION_H
#define DECOMPOSITION_ROUTESEQUENCEDECOMPOSITION_H

#include "../PALNSDecompositionMethod.h"
#include "../SubProblem.h"
#include <memory>

struct RouteSequenceDecomposition : public DecompositionMethod {
  using SPList = std::vector<std::unique_ptr<SubProblem>>;

  RouteSequenceDecomposition(const PALNSParameters& par, const CVRPInstance& i) : DecompositionMethod{par, i} {}
  virtual ~RouteSequenceDecomposition() = default;
  virtual SPList getSubproblems(const CVRPSolution& mpSol) const = 0;

  SPList getSubproblemsFromRouteIndices(const CVRPSolution& mpSol, const std::vector<int>& routeIndices) const;
  CVRPSolution decompose(CVRPSolution sol, double elapsedSec) const override;
};

#endif //DECOMPOSITION_ROUTESEQUENCEDECOMPOSITION_H
