#include "GreedySeqAlg.h"
#include "CVRPSolution.h"

// From SRP-Utils
#include "VectorUtils.h"

#include <vector>

using namespace std;

struct NodeDist {
    int m_iNode;
    double m_dDist;
    NodeDist(int iNode, double dDist) : m_iNode(iNode), m_dDist(dDist) {}

    // Used for sorting such that objects with large m_dist appears first in the sorted structure.
    bool operator<(const NodeDist& other) const { return m_dDist > other.m_dDist; }
};

// This method should always generate a feasible solution. It can generate solutions with unassigned
// customers though (if the number of vehicles is limited). Such solutions are penalized.
void GreedySeqAlg::go(CVRPSolution& sol, TSRandom& randGen) {
    // We make our own copy of the unassigned customers so we can
    // make changes to the vector. We could keep a reference instead and rely on CVRPSolution
    // to update the vector, but it seems a little dangerous if we change the implementation of
    // CVRPSolution.
    int i, iN = m_instance.getN();
    vector<int> vecUnassigned = sol.getUnassignedCusts();
    vector<bool> vecIsUnassigned;
    makeIsInSetVec(vecUnassigned, vecIsUnassigned, iN + 1);
    vector<NodeDist> vecNodeDepotDist;
    for(i = 0; i < (int)vecUnassigned.size(); i++) {
        int iCustId = vecUnassigned[i];
        double dDepotDist;
        if(m_bRandomSeedCustomers)
            dDepotDist = randGen.getRandomDouble(0, 10);
        else
            dDepotDist = m_instance.getDist(0, iCustId) + m_instance.getDist(iCustId, iN);
        vecNodeDepotDist.push_back(NodeDist(vecUnassigned[i], dDepotDist));
        vecIsUnassigned[iCustId] = true;
    }
    sort(vecNodeDepotDist.begin(), vecNodeDepotDist.end());
    int iCurRoute = 0;
    int iNextSeedNodeIdx = 0;
    if(sol.getNRoutes() == 0) {
        // We are starting from an empty solution. Create an initial route.
        createNewRoute(sol, vecNodeDepotDist, vecIsUnassigned, vecUnassigned, iNextSeedNodeIdx);
    }

    bool bStop = false;
    while(!vecUnassigned.empty() && !bStop) {
        bool bCreateNewRoute = false;
        // Check if we have enough customers to reach the final number of routes.
        if(m_instance.getNVehiclesFixed() && sol.getNRoutes() + (int)vecUnassigned.size() <= m_instance.getNVehicles()) {
            // We should have exactly enough customers left if we get here. Check:
            assert(sol.getNRoutes() + (int)vecUnassigned.size() == m_instance.getNVehicles());
            // Create a new route with just one customer (this is going to be repeated until we
            // have the right number of routes).
            bCreateNewRoute = true;
        }
        // Try to insert a customer into route <iCurRoute>
        if(!bCreateNewRoute && !performBestInsertion(sol, iCurRoute, vecIsUnassigned, vecUnassigned)) {
            // We couldn't insert any customers in the current route. Move to the next route.
            iCurRoute++;
            if(iCurRoute >= sol.getNRoutes()) {
                // The route doesn't exists. We'll have to open a new route if possible
                if(m_instance.getNVehiclesFixed() && sol.getNRoutes() == m_instance.getNVehicles()) {
                    // We cannot create more routes. We'll have to leave some customers unassigned.
                    bStop = true;
                } else {
                    bCreateNewRoute = true;
                    // cout << "Creating new route." << endl;
                    // waitEnter();
                }
            }
        }
        if(bCreateNewRoute)
            createNewRoute(sol, vecNodeDepotDist, vecIsUnassigned, vecUnassigned, iNextSeedNodeIdx);
    }
}

// Evaluates the cost of inserting every node from <vecUnassigned> into every
// possible position in route <iRoute>. Performs the cheapest insertion.
// Returns true if an insertion was found and false otherwise.
bool GreedySeqAlg::performBestInsertion(CVRPSolution& sol, int iRoute, vector<bool>& vecIsUnassigned, vector<int>& vecUnassigned) {
    // @@@ We could do much better in this function by caching costs of insertions. A route only changes at one
    // position every time we perform an insertion, so a lot of computations are repeated.
    // We could perhaps also optimize by recording infeasible moves, but that seems more tricky. Infeasibilities
    // because of capacity could easily be handled, but this is already checked as the first in the
    // <findBestInsertion> method in the Route class so we won't save much by this. The duration check is more
    // complicated as we do not assume that distances satisfy the triangle inequality. This means that insertion
    // of a customer can switch from infeasible to feasible after inserting other customers.
    int i, iBestCust = -1, iBestPos = -1;
    double dBestInc = DBL_MAX;
    int iPos;
    double dInc;
    for(i = 0; i < (int)vecUnassigned.size(); i++) {
        if(sol.findBestInsertion(iRoute, vecUnassigned[i], iPos, dInc)) {
            if(dInc < dBestInc) {
                dBestInc = dInc;
                iBestPos = iPos;
                iBestCust = vecUnassigned[i];
            }
        }
    }
    if(iBestCust != -1) {
        sol.insertCust(iRoute, iBestCust, iBestPos);
        vecIsUnassigned[iBestCust] = false;
        if(!vectorRemoveSingleElement(vecUnassigned, iBestCust))
            error("GreedySeqAlg::performBestInsertion(", "We inserted a customer that wasn't unassigned");
        return true;
    } else
        return false;
}

// Make another route.
void GreedySeqAlg::createNewRoute(CVRPSolution& sol, const vector<NodeDist>& vecNodeDepotDist, vector<bool>& vecIsUnassigned, vector<int>& vecUnassigned,
                                  int iNextSeedNodeIdx) {
    int iCustToInsert = getSeedCustomer(vecIsUnassigned, vecNodeDepotDist, iNextSeedNodeIdx);
    sol.createRoute(iCustToInsert);
    vecIsUnassigned[iCustToInsert] = false;
    if(!vectorRemoveSingleElement(vecUnassigned, iCustToInsert))
        error("GreedySeqAlg::performBestInsertion(", "We inserted a customer that wasn't unassigned");
}

int GreedySeqAlg::getSeedCustomer(const vector<bool>& vecIsUnassigned, const vector<NodeDist>& vecNodeDepotDist, int& iNextSeedNodeIdx) {
    // Find a seed customer
    // Go through the sorted vector of initially unassigned nodes
    // to find the unassigned node that is furthest from the depot.
    while(!vecIsUnassigned[vecNodeDepotDist[iNextSeedNodeIdx].m_iNode]) {
        ++iNextSeedNodeIdx;
        // The following should always be true as <vecUnassigned> is non-empty.
        assert(iNextSeedNodeIdx < (int)vecNodeDepotDist.size());
    }
    return vecNodeDepotDist[iNextSeedNodeIdx].m_iNode;
}
