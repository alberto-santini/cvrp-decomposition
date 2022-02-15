//
// Created by Alberto Santini on 29/05/2021.
//

#ifndef DECOMPOSITION_SUBPROBLEM_H
#define DECOMPOSITION_SUBPROBLEM_H

#include "PALNSParameters.h"
#include "../PALNS-CVRP/CVRPInstance.h"
#include "../PALNS-CVRP/CVRPSolution.h"
#include <map>
#include <vector>
#include <optional>
#include <string>

struct SubProblem {
  static const std::string cvrp_specific_params;

  const PALNSParameters& mp_params;
  const CVRPInstance& mp_instance;
  const CVRPSolution& mp_sol;

  CVRPInstance sp_instance;
  std::map<int, int> mp_to_sp_customers;
  std::map<int, int> sp_to_mp_customers;
  std::string outputPrefix = "";
  double elapsedSec;

  SubProblem(const PALNSParameters& mp_params, const CVRPInstance& mp_instance, const CVRPSolution& mp_sol, const std::vector<int>& mp_sol_routes, double elapsedSec);
  std::optional<CVRPSolution> solve();
};

#endif //DECOMPOSITION_SUBPROBLEM_H
