#include "SeqInsertionCache.h"
#include "CVRPInstance.h"
#include "CVRPSolution.h"

SeqInsertionCache::SeqInsertionCache(const CVRPInstance& instance, CVRPSolution& sol, const vector<vector<int>>& vecUnassignedSeqs) : m_instance(instance) {
    initCache(sol, vecUnassignedSeqs);
}

void SeqInsertionCache::initCache(CVRPSolution& sol, const vector<vector<int>>& vecUnassignedSeqs) {
    int i, j, iRoute, iPos, iN = m_instance.getN();
    double dDelta;
    m_vecUnassignedSeqs.assign(vecUnassignedSeqs.begin(), vecUnassignedSeqs.end());
    m_seqIsUnassigned.assign(m_vecUnassignedSeqs.size(), true);
    for(i = 0; i < (int)vecUnassignedSeqs.size(); ++i) {
        const vector<int>& vecSeq = vecUnassignedSeqs[i];
        m_vecSequenceFirstCust.push_back(vecSeq.front());
        m_vecSequenceLastCust.push_back(vecSeq.back());
        int iDemand = 0, iServTime = 0;
        double dDist = 0;
        for(j = 0; j < (int)vecSeq.size(); ++j) {
            int iCustId = vecSeq[j];
            iDemand += m_instance.getDemand(iCustId);
            iServTime += m_instance.getServiceTime(iCustId);
            if(j > 0)
                dDist += m_instance.getDist(vecSeq[j - 1], iCustId);
        }
        m_vecSequenceDemand.push_back(iDemand);
        m_vecSequenceDist.push_back(dDist);
        m_vecSequenceServTime.push_back(iServTime);
    }

    // We reserve room for 5 more routes in case we are adding routes in the course of the insertion.
    // If add more than 5 routes we'll have to extend the size of the matrix.
    m_iNReservedRoutes = sol.getNRoutes() + 5;
    m_vecRouteDelta.resize(vecUnassignedSeqs.size());
    m_matPossible.resizeAndClear((int)vecUnassignedSeqs.size(), m_iNReservedRoutes, true);

    int iMaxLoad = m_instance.getMaxLoad();
    int iMaxDuration = m_instance.getMaxDuration();
    int iSeqId;
    for(iSeqId = 0; iSeqId < (int)vecUnassignedSeqs.size(); ++iSeqId) {
        int iDemand = m_vecSequenceDemand[iSeqId];
        double dDist = m_vecSequenceDist[iSeqId];
        int iServTime = m_vecSequenceServTime[iSeqId];
        int iFirstCust = m_vecSequenceFirstCust[iSeqId];
        int iLastCust = m_vecSequenceLastCust[iSeqId];
        for(iRoute = 0; iRoute < sol.getNRoutes(); iRoute++) {
            // @@@ We only assume that distances are non negative, we do not assume that they satisfy the triangle inequality
            // If we knew that we could mark even more route/sequence combinations as infeasible. Now we cannot do anything (much)
            // because the duration of the route may decrease if more is inserted into it (we only know that it cannot
            // decrease below the total service time on the route, but that is probably a very weak bound).
            if(sol.getRoute(iRoute).getLoad() + iDemand > iMaxLoad)
                m_matPossible.setElement(iSeqId, iRoute, false);
            else {
                bool bFoundSol, bReversed;
                bFoundSol = sol.findBestSeqInsertion(iRoute, iFirstCust, iLastCust, dDist + iServTime, iPos, dDelta, bReversed);
                if(bFoundSol)
                    m_vecRouteDelta[iSeqId].push_back(RouteDelta(iRoute, dDelta + dDist));
            }
        }
        sort(m_vecRouteDelta[iSeqId].begin(), m_vecRouteDelta[iSeqId].end());
    }
}

inline void eraseRouteDelta(vector<RouteDelta>& vecRouteDelta, int iRoute) {
    int i, iSize = (int)vecRouteDelta.size();
    for(i = 0; i < iSize; ++i) {
        if(vecRouteDelta[i].m_iRoute == iRoute) {
            vecRouteDelta.erase(vecRouteDelta.begin() + i);
            break;
        }
    }
}

// updateCache(CVRPSolution &sol, int iChangedRoute, int iSeqIdInserted)
// is called every time a sequence (iSeqIdInserted) has been inserted into a
// Route (iChangedRoute)

void SeqInsertionCache::updateCache(CVRPSolution& sol, int iChangedRoute, int iSeqIdInserted) {
    m_seqIsUnassigned[iSeqIdInserted] = false;
    double dDelta;
    int i, iPos, iMaxLoad = m_instance.getMaxLoad();
    // Did we create a new route? If so, check if we have room for this in our data structures
    // or if we have to extend them.
    if(sol.getNRoutes() >= m_iNReservedRoutes) {
        // Ooops we are out of space.
        // cout << "Out of space! extending cache" << endl;
        m_iNReservedRoutes = sol.getNRoutes() + 5;
        m_matPossible.expand(m_matPossible.getNRows(), m_iNReservedRoutes, true);
    }

    // A sequence has been added to route <iChangedRoute>. Update the cache.
    // @@@ We could be more clever than we are here. We actually only need to
    // look at the two positions around the inserted node, not at all positions
    // in the entire route. There is one exception though because we do not
    // assume that distances obbeys the travel distance: If the duration of the
    // route has decreased after inserting the new node, some insertion positions
    // that were infeasible earlier could have become feasible. In that case we
    // should scan the entire route (I think we can be smarter if the cost and
    // distance matrix is the same (as it is at the moment)).
    for(i = 0; i < (int)m_vecUnassignedSeqs.size(); i++) {
        if(m_seqIsUnassigned[i]) {
            if(m_matPossible.getElement(i, iChangedRoute)) {
                // Find the entry corresponding to <iChangedRoute> in <m_vecRouteDelta[i]>
                // and remove it (if it exists). If insertion into the route is feasible,
                // then we'll insert a new element.
                eraseRouteDelta(m_vecRouteDelta[i], iChangedRoute);
                if(sol.getRoute(iChangedRoute).getLoad() + m_vecSequenceDemand[i] > iMaxLoad)
                    m_matPossible.setElement(i, iChangedRoute, false);
                else {
                    bool bFoundSol, bReversed;
                    bFoundSol = sol.findBestSeqInsertion(iChangedRoute, m_vecSequenceFirstCust[i], m_vecSequenceLastCust[i],
                                                         m_vecSequenceDist[i] + m_vecSequenceServTime[i], iPos, dDelta, bReversed);
                    if(bFoundSol)
                        insertOrderedVec(m_vecRouteDelta[i], RouteDelta(iChangedRoute, dDelta + m_vecSequenceDist[i]));
                }
            }
        }
    }
}
