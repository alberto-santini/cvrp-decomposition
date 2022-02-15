//
// Created by Alberto Santini on 01/06/2021.
//

#include "SubProblem.h"
#include "../PALNS-CVRP/PALNS-CVRP.h"
#include "PALNS.h"
#include <iostream>
#include <optional>

const std::string SubProblem::cvrp_specific_params = "../src/alns/ALNS_SubProblem_Params.txt";

SubProblem::SubProblem(const PALNSParameters& mp_params, const CVRPInstance& mp_instance, const CVRPSolution& mp_sol, const std::vector<int>& mp_sol_routes,
                       double elapsedSec)
    : mp_params{mp_params}, mp_instance{mp_instance}, mp_sol{mp_sol}, elapsedSec{elapsedSec} {
    sp_instance = CVRPInstance();

    int curr_vertex = 0;
    int sp_n_vertices = 0;
    std::vector<Coordinate<double>> sp_coords;
    std::vector<int> sp_demand;
    std::vector<int> sp_service_time;
    const int sp_n_vehicles = mp_sol_routes.size();
    const int sp_capacity = mp_instance.getMaxLoad();
    const int sp_max_duration = mp_instance.getMaxDuration();

    // Map the depot to the depot
    sp_coords.emplace_back(mp_instance.getCoord(0).getX(), mp_instance.getCoord(0).getY());
    mp_to_sp_customers[0] = 0;
    sp_to_mp_customers[0] = 0;
    sp_demand.push_back(mp_instance.getDemand(0));
    sp_service_time.push_back(mp_instance.getServiceTime(0));

    for(auto route_id : mp_sol_routes) {
        const auto& route = mp_sol.getRoute(route_id);
        sp_n_vertices += route.size() - 2;

        for(auto node : route.getNodes()) {
            if(node == 0) {
                continue;
            }
            if(node == mp_instance.getN() + 1) {
                continue;
            }

            ++curr_vertex;

            sp_coords.emplace_back(mp_instance.getCoord(node).getX(), mp_instance.getCoord(node).getY());
            mp_to_sp_customers[node] = curr_vertex;
            sp_to_mp_customers[curr_vertex] = node;
            sp_demand.push_back(mp_instance.getDemand(node));
            sp_service_time.push_back(mp_instance.getServiceTime(node));
        }
    }

    // Special "return" depot
    sp_coords.emplace_back(mp_instance.getCoord(mp_instance.getN() + 1).getX(), mp_instance.getCoord(mp_instance.getN() + 1).getY());
    mp_to_sp_customers[mp_instance.getN() + 1] = sp_n_vertices + 1;
    sp_to_mp_customers[sp_n_vertices + 1] = mp_instance.getN() + 1;
    sp_demand.push_back(mp_instance.getDemand(mp_instance.getN() + 1));
    sp_service_time.push_back(mp_instance.getServiceTime(mp_instance.getN() + 1));

    sp_instance.setProblem(sp_n_vertices, sp_coords, sp_demand, sp_service_time, sp_n_vehicles, false, sp_capacity, sp_max_duration);

    Matrix<double> sp_dist_matrix{sp_n_vertices + 2, 0.0};

    for(auto i = 0; i <= sp_n_vertices + 1; ++i) {
        const auto mp_i = sp_to_mp_customers[i];
        for(auto j = 0; j <= sp_n_vertices + 1; ++j) {
            const auto mp_j = sp_to_mp_customers[j];

            sp_dist_matrix.setElement(i, j, mp_instance.getDist(mp_i, mp_j));
        }
    }

    sp_instance.setDistMatrix(sp_dist_matrix);
}

std::optional<CVRPSolution> SubProblem::solve() {
    if(sp_instance.getN() == 0 || sp_instance.getNVehicles() == 0) {
#ifndef BATCH_MODE
        std::cout << mp_params.outputPrefix << "[Subproblem] Empty subproblem (" << sp_instance.getN() << " customers, " << sp_instance.getNVehicles()
                  << " vehicles). "
                  << "Skipping.\n";
#endif
        return std::nullopt;
    }

#ifndef BATCH_MODE
    std::cout << mp_params.outputPrefix << "[Subproblem] Solving a subproblem (" << sp_instance.getN() << " customers, " << sp_instance.getNVehicles()
              << " vehicles).\n";
#endif

    CVRPSolution sp_sol{sp_instance};
    ParFileParser pfp;
    pfp.loadParFile(cvrp_specific_params);
    pfp.setParameter("MaxIter", "10000");
    pfp.setParameter("DecompositionMethod", std::to_string(mp_params.m_iDecompositionMethod));
    pfp.setParameter("DecompositionIters", std::to_string(mp_params.m_iDecompositionIters));
    pfp.setParameter("DecompositionSize", std::to_string(mp_params.m_iDecompositionSize));
    pfp.setParameter("DecompositionLevel", std::to_string(mp_params.m_iDecompositionLevel));
    pfp.setParameter("OutputPrefix", mp_params.outputPrefix);

    double timeoutSec;
    pfp.getDoubleParameter("TimeoutSec", timeoutSec);
    timeoutSec -= elapsedSec;
    pfp.setParameter("TimeoutSec", std::to_string(timeoutSec));

    PALNS_CVRP solver;
    solver.go(sp_sol, 1 /*thread*/, 1 /*retry*/, false /*print header*/, pfp);

    CVRPSolution new_mp_sol{mp_instance};

    new_mp_sol.clear(sp_sol.getNRoutes());

    for(auto i = 0; i < sp_sol.getNRoutes(); ++i) {
        const auto& sp_route = sp_sol.getRoute(i);
        std::vector<int> mp_route_nodes;

        for(auto node : sp_route.getNodes()) {
            mp_route_nodes.push_back(sp_to_mp_customers[node]);
        }

        new_mp_sol.setRouteNodes(i, mp_route_nodes);
    }

    return new_mp_sol;
}