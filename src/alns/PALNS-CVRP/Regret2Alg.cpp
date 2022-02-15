#include "Regret2Alg.h"
#include "CVRPInstance.h"
#include "CVRPSolution.h"
#include "InsertionCache.h"

// From SRP-Utils
#include "VectorUtils.h"

// <pMatAltCost> If NULL, then use the ordinary cost from <m_instance>. If not NULL, then use this cost matrix instead of the
// standard cost. Currently used by "HistoryRegret2Repair".
void Regret2Alg::go(CVRPSolution& sol, const Matrix<double>* pMatAltCost) {
    // We make our own copy of the unassigned customers so we can
    // make changes to the vector. We could keep a reference instead and rely on CVRPSolution
    // to update the vector, but it seems a little dangerous if we change the implementation of
    // CVRPSolution.
    int i, iN = m_instance.getN();
    vector<int> vecUnassigned = sol.getUnassignedCusts();
    vector<bool> vecIsUnassigned;
    makeIsInSetVec(vecUnassigned, vecIsUnassigned, iN + 1);
    double dRnd = 1;

    InsertionCache insCache(m_instance, sol, vecUnassigned, pMatAltCost);

    bool bStop = false;
    // <iNUnusedRoutes> counts the number of unused routes. This number is important if
    // the number of vehicles is fixed - we should make sure to reach the right number of vehicles
    int iNUnusedRoutes = sol.countUnusedRoutes();

    if(m_instance.getNVehiclesFixed()) {
        if(sol.getNRoutes() != m_instance.getNVehicles()) {
            error("GreedyParAlg::go(...)", string("number of vehicles is fixed, but the solution does not have the right number of ") + string("vehicles"));
        }
        if((int)vecUnassigned.size() < iNUnusedRoutes) {
            error("GreedyParAlg::go(...)",
                  string("Number of vehicles is fixed and number of free") + string(" customers is smaller than the number of empty vehicles"));
        }
    }

    while(!vecUnassigned.empty() && !bStop) {
        // Check if we have enough customers to reach the final number of routes.
        if(m_instance.getNVehiclesFixed() && (int)vecUnassigned.size() <= iNUnusedRoutes) {
            // We should have exactly enough customers left if we get here. Check:
            if((int)vecUnassigned.size() != iNUnusedRoutes)
                error("GreedyParAlg::go(...)", "#vehicles fixed, too few unassigned customers");
            // Insert the remaining customers on the empty routes.
            // First find the empty routes.
            vector<int> vecEmptyRoutes;
            for(i = 0; i < sol.getNRoutes(); i++) {
                if(sol.getRoute(i).size() <= 2)
                    vecEmptyRoutes.push_back(i);
            }
            assert(vecEmptyRoutes.size() == vecUnassigned.size());
            for(i = 0; i < (int)vecUnassigned.size(); ++i) {
                sol.insertCust(vecEmptyRoutes[i], vecUnassigned[i], 1);
            }
            vecUnassigned.clear();
        } else {
            // Go through all unassigned customers and try to insert them on every route.
            // We stop if no insertion is possible (that only happens if the number of vehicles is limited).
            // cout << "dRnd = " << dRnd << "  ";
            double dBest, dSecondBest, dBestRegret = -DBL_MAX;
            bool bOneInsertion = false;
            int iBestIdx = -1;
            int iRoute, iBestRoute = -1;
            for(i = 0; i < (int)vecUnassigned.size(); i++) {
                if(m_dRandFactor > 1)
                    dRnd = m_randGen.getRandomDouble(1.0, m_dRandFactor);

                int iCustId = vecUnassigned[i];
                dBest = insCache.getBestInsertion(iCustId, iRoute);
                dSecondBest = insCache.getNthBestInsertionDelta(iCustId, 1);
                double dRegret;
                if(dBest != DBL_MAX) {
                    if(dSecondBest == DBL_MAX && (!bOneInsertion || dBest * dRnd < dBestRegret)) {
                        // Only one insertion left for <iCustId>
                        bOneInsertion = true;
                        dBestRegret = dBest * dRnd;
                        iBestRoute = iRoute;
                        iBestIdx = i;
                    } else {
                        dRegret = (dSecondBest - dBest) * dRnd;
                        if(dRegret > dBestRegret && !bOneInsertion) {
                            dBestRegret = dRegret;
                            iBestRoute = iRoute;
                            iBestIdx = i;
                        }
                    }
                }
            }

            if(iBestIdx == -1)
                bStop = true;
            else {
                // Insert the customer.
                double dDelta;
                int iPos, iCustId = vecUnassigned[iBestIdx];
                unorderedVectorRemoveAt(vecUnassigned, iBestIdx);
                sol.findBestInsertion(iBestRoute, iCustId, iPos, dDelta);
                sol.insertCust(iBestRoute, iCustId, iPos);
                if(sol.getRoute(iBestRoute).size() == 3)
                    --iNUnusedRoutes;
                if(iNUnusedRoutes == 0 && !m_instance.getNVehiclesFixed()) {
                    // We have no more empty routes in the solution. Add a new one.
                    ++iNUnusedRoutes;
                    sol.createRoute();
                    insCache.updateCache(sol, sol.getNRoutes() - 1, vecUnassigned);
                }
                insCache.updateCache(sol, iBestRoute, vecUnassigned);
            }
        }
    }
    sol.consCalc();
}
