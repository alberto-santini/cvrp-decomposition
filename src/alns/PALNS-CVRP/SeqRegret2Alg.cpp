#include "SeqRegret2Alg.h"
#include "CVRPInstance.h"
#include "CVRPSolution.h"
#include "SeqInsertionCache.h"

// From SRP-Utils
#include "VectorUtils.h"

void SeqRegret2Alg::go(CVRPSolution& sol) {
    // We make our own copy of the unassigned sequences so we can
    // make changes to the vector. We could keep a reference instead and rely on CVRPSolution
    // to update the vector, but it seems a little dangerous if we change the implementation of
    // CVRPSolution.
    int i, iN = m_instance.getN();
    const vector<vector<int>> vecUnassignedSeqs = sol.getUnassignedSeqs();
    // We take over control of the unassinged sequences.
    sol.clearUnassignedSeqs();
    vector<int> vecUnassignedSeqIds;
    makeRange(vecUnassignedSeqIds, 0, (int)vecUnassignedSeqs.size() - 1);
    double dRnd = 1;

    SeqInsertionCache insCache(m_instance, sol, vecUnassignedSeqs);

    bool bStop = false;
    // <iNUnusedRoutes> counts the number of unused routes. This number is important if
    // the number of vehicles is fixed - we should make sure to reach the right number of vehicles
    int iNUnusedRoutes = sol.countUnusedRoutes();

    if(m_instance.getNVehiclesFixed()) {
        if(sol.getNRoutes() != m_instance.getNVehicles()) {
            error("SeqRegret2Alg::go(...)", string("number of vehicles is fixed, but the solution does not have the right number of ") + string("vehicles"));
        }
        if((int)vecUnassignedSeqIds.size() < iNUnusedRoutes) {
            error("SeqRegret2Alg::go(...)",
                  string("Number of vehicles is fixed and number of free") + string(" sequences is smaller than the number of empty vehicles"));
        }
    }

    while(!vecUnassignedSeqIds.empty() && !bStop) {
        // Check if we have enough sequences to reach the final number of routes.
        if(m_instance.getNVehiclesFixed() && (int)vecUnassignedSeqIds.size() <= iNUnusedRoutes) {
            // We should have exactly enough sequences left if we get here. Check:
            if((int)vecUnassignedSeqIds.size() != iNUnusedRoutes)
                error("SeqRegret2Alg::go(...)", "#vehicles fixed, too few unassigned sequences");
            // Insert the remaining sequences on the empty routes.
            // First find the empty routes.
            vector<int> vecEmptyRoutes;
            for(i = 0; i < sol.getNRoutes(); i++) {
                if(sol.getRoute(i).size() <= 2)
                    vecEmptyRoutes.push_back(i);
            }
            assert(vecEmptyRoutes.size() == vecUnassignedSeqIds.size());
            for(i = 0; i < (int)vecUnassignedSeqIds.size(); ++i) {
                sol.insertSeq(vecEmptyRoutes[i], vecUnassignedSeqs[vecUnassignedSeqIds[i]], false, 1);
            }
            vecUnassignedSeqIds.clear();
        } else {
            // Go through all unassigned sequences and try to insert them on every route.
            // We stop if no insertion is possible (that only happens if the number of vehicles is limited).
            // cout << "dRnd = " << dRnd << "  ";
            double dBest, dSecondBest, dBestRegret = -DBL_MAX;
            bool bOneInsertion = false;
            int iBestIdx = -1;
            int iRoute, iBestRoute = -1;
            for(i = 0; i < (int)vecUnassignedSeqIds.size(); i++) {
                if(m_dRandFactor > 1)
                    dRnd = m_randGen.getRandomDouble(1.0, m_dRandFactor);

                int iSeqId = vecUnassignedSeqIds[i];
                dBest = insCache.getBestInsertion(iSeqId, iRoute);
                dSecondBest = insCache.getNthBestInsertionDelta(iSeqId, 1);
                double dRegret;
                if(dBest != DBL_MAX) {
                    if(dSecondBest == DBL_MAX && (!bOneInsertion || dBest * dRnd < dBestRegret)) {
                        // Only one insertion left for <iSeqId>
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
                bool bReversed;
                int iPos, iSeqId = vecUnassignedSeqIds[iBestIdx];
                const vector<int>& vecSeq = vecUnassignedSeqs[iSeqId];
                // cout << "Inserting: ";
                // outputVector(cout, vecSeq, " ", "\n");
                unorderedVectorRemoveAt(vecUnassignedSeqIds, iBestIdx);
                // Notice, the delta calculated by findBestSeqInsertion does not take the length of the sequence into account, just the change in
                // cost related to the detour over the first and last customer of the sequence.
                sol.findBestSeqInsertion(iBestRoute, vecSeq.front(), vecSeq.back(), 0 /* Seq duration, there is no need to calculate this, it won't change
                                                                                      the best insertion position */
                                         ,
                                         iPos, dDelta, bReversed);

                // --- debug code ---:
                double dSeqCost = 0;
                int idx;
                for(idx = 0; idx < vecSeq.size() - 1; ++idx) {
                    dSeqCost += m_instance.getDist(vecSeq[idx], vecSeq[idx + 1]);
                }
                // cout << "Delta: " << dDelta << ", seq cost: " << dSeqCost << " cost of sol: " << sol.getCost();
                double dNewExpectedCost = sol.getCost() + dDelta + dSeqCost - sol.getUnassignedCost() * vecSeq.size();
                // cout << ", expected new solution cost: " << dNewExpectedCost << endl;
                // ------------------
                // Are we inserted into an empty route?
                if(sol.getRoute(iBestRoute).size() == 2)
                    // Yes, decrease the empty route counter.
                    --iNUnusedRoutes;
                sol.insertSeq(iBestRoute, vecSeq, bReversed, iPos);
                // cout << "cost of solution after insertion: " << sol.getCost() << endl;
                if(!epsilonEqual(sol.getCost(), dNewExpectedCost, 1e-5))
                    error("SeqRegret2Alg::go()", "Unexpected solution cost after sequence insertion");

                if(iNUnusedRoutes == 0 && !m_instance.getNVehiclesFixed()) {
                    // We have no more empty routes in the solution. Add a new one.
                    ++iNUnusedRoutes;
                    sol.createRoute();
                    insCache.updateCache(sol, sol.getNRoutes() - 1, iSeqId);
                }
                insCache.updateCache(sol, iBestRoute, iSeqId);
            }
        }
    }
    // If some sequences were left unassigned then insert them back into the solution as unassigned sequences.
    // @@@ We could also try to break the unassigned sequences up into individual customers and insert those, that might produce
    // @@@ feasible solutions in some cases. We could easily do this by calling the appropriate customer oriented insertion algorithm.
    while(!vecUnassignedSeqIds.empty()) {
        sol.addUnassignedSeq(vecUnassignedSeqs[vecUnassignedSeqIds.back()]);
        vecUnassignedSeqIds.pop_back();
    }

    sol.consCalc();
}
