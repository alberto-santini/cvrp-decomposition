//
// Created by Alberto Santini on 04/06/2021.
//

#ifndef DECOMPOSITION_ARCBASEDDECOMPOSITION_H
#define DECOMPOSITION_ARCBASEDDECOMPOSITION_H

#include "../PALNSDecompositionMethod.h"
#include "../PALNSParameters.h"
#include "../../PALNS-CVRP/CVRPInstance.h"
#include <tuple>
#include <vector>

struct ArcBasedDecomposition : public DecompositionMethod {
  ArcBasedDecomposition(const PALNSParameters& par, const CVRPInstance& i) : DecompositionMethod{par, i} {}
  ~ArcBasedDecomposition() override = default;

  CVRPSolution decompose(CVRPSolution sol, double elapsedSec) const override;

  struct Arc {
    int o; int d; double cost;
    bool operator==(const Arc& other) const { return o == other.o && d == other.d; }
    bool operator<(const Arc& other) const { return std::tie(o, d) < std::tie(other.o, other.d); }
  };
  using Path = std::vector<int>;
  using CustomerSet = std::vector<int>;

  std::pair<CustomerSet, std::vector<Path>> getPaths(const std::vector<Arc>& arcs) const;
  void addArcToPaths(const Arc& arc, std::vector<Path>& paths, CustomerSet& custsInPaths) const;
  virtual std::vector<Arc> getArcsToFix(const CVRPSolution& mpSol) const = 0;
};


#endif //DECOMPOSITION_ARCBASEDDECOMPOSITION_H
