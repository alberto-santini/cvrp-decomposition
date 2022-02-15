//
// Created by Alberto Santini on 29/05/2021.
//

#ifndef DECOMPOSITION_PALNSDECOMPOSITIONMETHOD_H
#define DECOMPOSITION_PALNSDECOMPOSITIONMETHOD_H

#include "../PALNS-CVRP/CVRPSolution.h"
#include "../PALNS-CVRP/CVRPInstance.h"
#include "PALNSParameters.h"
#include "SubProblem.h"
#include <memory>

struct PALNS;

struct DecompositionMethod {
  const PALNSParameters& par;
  const CVRPInstance& i;
  mutable double elapsedSec;

  DecompositionMethod(const PALNSParameters& par, const CVRPInstance& i) : par{par}, i{i} {}

  virtual CVRPSolution decompose(CVRPSolution sol, double elapsedSec) const = 0;
  virtual ~DecompositionMethod() = default;

protected:
  void mergeMpPartialSolutionToMpSolution(const CVRPSolution& mpParSol, CVRPSolution& mpSol) const;
};

std::unique_ptr<DecompositionMethod> get_decomposition_method(const PALNS& p);

#endif //DECOMPOSITION_PALNSDECOMPOSITIONMETHOD_H
