#include "InsertionCache.h"
#include "CVRPInstance.h"
#include "CVRPSolution.h"

InsertionCache::InsertionCache(const CVRPInstance& instance, CVRPSolution& sol, const vector<int>& vecUnassigned, const Matrix<double>* pMatAltCost)
    : m_instance(instance), m_pMatAltCost(pMatAltCost) {
    initCache(sol, vecUnassigned);
}

void InsertionCache::initCache(CVRPSolution& sol, const vector<int>& vecUnassigned) {
    int iRoute, iPos, iN = m_instance.getN();
    double dDelta;
    // Build the map from the original to compressed id.
    m_vecOrigToCompressedId.assign(iN + 1, -1);
    int iComprId = 0;

    // We reserve room for 5 more routes in case we are adding routes in the course of the insertion.
    // If add more than 5 routes we'll have to extend the size of the matrix.
    m_iNReservedRoutes = sol.getNRoutes() + 5;
    m_vecRouteDelta.resize(vecUnassigned.size());
    m_matPossible.resizeAndClear((int)vecUnassigned.size(), m_iNReservedRoutes, true);

    int iMaxLoad = m_instance.getMaxLoad();
    for(iComprId = 0; iComprId < (int)vecUnassigned.size(); ++iComprId) {
        int iCustId = vecUnassigned[iComprId];
        int iDemand = m_instance.getDemand(iCustId);
        m_vecOrigToCompressedId[iCustId] = iComprId;
        for(iRoute = 0; iRoute < sol.getNRoutes(); iRoute++) {
            // @@@ If we knew that distances satisfy the triangle inequality then we could
            // also use the route duration constraint to mark customer/route combinations as
            // infeasible.
            if(sol.getRoute(iRoute).getLoad() + iDemand > iMaxLoad)
                m_matPossible.setElement(iComprId, iRoute, false);
            else {
                bool bFoundSol;
                if(m_pMatAltCost)
                    bFoundSol = sol.altCostFindBestIns(iRoute, iCustId, iPos, dDelta, *m_pMatAltCost);
                else
                    bFoundSol = sol.findBestInsertion(iRoute, iCustId, iPos, dDelta);
                if(bFoundSol) {
                    m_vecRouteDelta[iComprId].push_back(RouteDelta(iRoute, dDelta));
                }
            }
        }
        sort(m_vecRouteDelta[iComprId].begin(), m_vecRouteDelta[iComprId].end());
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

void InsertionCache::updateCache(CVRPSolution& sol, int iChangedRoute, const vector<int>& vecUnassigned) {
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

    // A node has been added to route <iChangedRoute>. Update the cache.
    // @@@ We could be more clever than we are here. We actually only need to
    // look at the two positions around the inserted node, not at all positions
    // in the entire route. There is one exception though because we do not
    // assume that distances obbeys the travel distance: If the duration of the
    // route has decreased after inserting the new node, some insertion positions
    // that were infeasible earlier could have become feasible. In that case we
    // should scan the entire route (I think we can be smarter if the cost and
    // distance matrix is the same (as it is at the moment)).
    for(i = 0; i < (int)vecUnassigned.size(); i++) {
        int iCustId = vecUnassigned[i];
        int iComprId = m_vecOrigToCompressedId[iCustId];
        if(m_matPossible.getElement(iComprId, iChangedRoute)) {
            // Find the entry corresponding to <iChangedRoute> in <m_vecRouteDelta[iComprId]>
            // and remove it (if it exists). If insertion into the route is feasible,
            // then we'll insert a new element.
            eraseRouteDelta(m_vecRouteDelta[iComprId], iChangedRoute);
            if(sol.getRoute(iChangedRoute).getLoad() + m_instance.getDemand(iCustId) > iMaxLoad)
                m_matPossible.setElement(iComprId, iChangedRoute, false);
            else {
                bool bFoundSol;
                if(m_pMatAltCost)
                    bFoundSol = sol.altCostFindBestIns(iChangedRoute, iCustId, iPos, dDelta, *m_pMatAltCost);
                else
                    bFoundSol = sol.findBestInsertion(iChangedRoute, iCustId, iPos, dDelta);
                if(bFoundSol) {
                    // insertOrderedVec is defined in VectorUtils.h
                    insertOrderedVec(m_vecRouteDelta[iComprId], RouteDelta(iChangedRoute, dDelta));
                }
            }
        }
    }
}
