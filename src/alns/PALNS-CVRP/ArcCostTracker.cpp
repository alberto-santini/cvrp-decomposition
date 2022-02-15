#include "ArcCostTracker.h"
#include "CVRPSolution.h"

void ArcCostTracker::callback(const CVRPSolution& sol, bool bAccepted, bool bImproved, bool bNewBest) {
    double dCost = sol.getCost();
    int iRoute, iPos;
    for(iRoute = 0; iRoute < sol.getNRoutes(); ++iRoute) {
        const Route& route = sol.getRoute(iRoute);
        for(iPos = 0; iPos < route.size() - 1; ++iPos) {
            int iFrom = route[iPos], iTo = route[iPos + 1];
            if(m_matBestSolCost.getElement(iFrom, iTo) > dCost)
                m_matBestSolCost.setElement(iFrom, iTo, dCost);
        }
    }
}