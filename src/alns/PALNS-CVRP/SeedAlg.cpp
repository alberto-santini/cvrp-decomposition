#include "SeedAlg.h"
#include "CVRPInstance.h"
#include "CVRPSolution.h"

// Assign seed customers to each route in sol. It is assummed that all routes in sol are empty.
// The method first chooses the node furthest from the depot, then the node furthest from the depot
// and the first node, then the node furthest from the depot and the two already selected nodes and so on.

// The method assumes symmetric distances (it wouldn't crash for assymetric distances, but it
// wouldn't make much sense either)
void SeedAlg::assignSeedCustomers(CVRPSolution& sol) {
    // Check that all routes are empty.
    int i, iRoute, iN = m_instance.getN();
    for(i = 0; i < sol.getNRoutes(); i++) {
        if(!sol.getRoute(i).empty())
            error("GreedyParAlg::assignSeedCustomers(...)", "One route was not empty");
    }
    // The set of nodes we want our seed nodes to be far away from.
    vector<int> vecNodeSet;
    vecNodeSet.push_back(0);
    vector<bool> vecIsUnassigned(iN + 1, true);
    // vecDist gives the distance to the nearest node in vecNodeSet for each node
    vector<double> vecDist(iN + 1, DBL_MAX);
    // Now find the seed nodes
    for(iRoute = 0; iRoute < sol.getNRoutes(); iRoute++) {
        int iLastSeedNode = vecNodeSet.back();
        double dLargestDist = -DBL_MAX;
        int iNewSeedNode = -1;
        for(i = 1; i <= iN; ++i) {
            if(vecIsUnassigned[i]) {
                double dDist = minFunc(vecDist[i], m_instance.getDist(iLastSeedNode, i));
                vecDist[i] = dDist;
                if(dDist > dLargestDist) {
                    dLargestDist = dDist;
                    iNewSeedNode = i;
                }
            }
        }
        if(iNewSeedNode == -1)
            error("GreedyParAlg::assignSeedCustomers(...)", "Didn't find a seed node.");
        sol.insertCust(iRoute, iNewSeedNode, 1);
        // iNewSeedNode is now used as a seed node.
        vecIsUnassigned[iNewSeedNode] = false;
        vecNodeSet.push_back(iNewSeedNode);
    }
}
