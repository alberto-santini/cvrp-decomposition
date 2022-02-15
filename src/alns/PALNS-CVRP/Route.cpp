// @@@ Notice! the savings algorithm assumes symmetric costs!

#include "Route.h"

// From SRP-Utils
#include "Utils.h"
#include "VectorUtils.h"

Route::Route(const CVRPInstance& instance) : m_pInstance(&instance), m_pMatDists(&(instance.getDists())) {
    int iN = instance.getN();
    m_vecRoute.push_back(0);
    m_vecRoute.push_back(iN + 1);

    m_iLoad = 0;
    m_dCost = m_pMatDists->getElement(0, iN + 1);
    m_dDuration = m_pMatDists->getElement(0, iN + 1);
    m_dToleranceEps = 1E-7;
}

/***********************************************************************
 * bool Route::removeCustomer(int iPos)
 *
 * remove the customer at position <iPos>
 *
 * --- input parameters ---
 * iPos			: The position of the customer to remove
 *                iPos = 0 is invalid (position 0 is occupied by the depot).
 * --- output parameters ---
 * --- return value ---
 ***********************************************************************/

void Route::removeCustomer(int iPos) {
#ifdef _DEBUG
    if(iPos <= 0 || iPos >= (int)m_vecRoute.size() - 1)
        error("Route::removeCustomer(...)", "remove position out of range ");
#endif

    int iCustId = m_vecRoute[iPos];
    double inc = m_pMatDists->getElement(m_vecRoute[iPos - 1], m_vecRoute[iPos + 1]) - m_pMatDists->getElement(m_vecRoute[iPos - 1], iCustId) -
                 m_pMatDists->getElement(iCustId, m_vecRoute[iPos + 1]);

    m_dCost += inc;
    m_iLoad -= m_pInstance->getDemand(iCustId);
    m_dDuration += inc;
    m_dDuration -= m_pInstance->getServiceTime(iCustId);
    m_vecRoute.erase(m_vecRoute.begin() + iPos);
}

// Go through the route and removes customers if vecIsCustRemoved[id] == true
// where <id> is the id of the customer. The method returns the number of customers removed.
void Route::removeCustomers(const vector<bool>& vecIsCustRemoved, vector<int>& vecRemovedCusts) {
    vecRemovedCusts.clear();
    assert((int)vecIsCustRemoved.size() >= (m_pInstance->getN() + 1));
    vector<int> vecNewRoute(m_vecRoute.size());
    vecNewRoute[0] = 0;
    int iOrigIdx, iNewIdx = 1;
    // Don't look at the last element in the route (it's the end depot).
    int iOrigLastIdx = (int)m_vecRoute.size() - 1;
    // dDelta : change in distance.
    double dDelta = 0;
    int iServTimeDelta = 0;
    for(iOrigIdx = 1; iOrigIdx < iOrigLastIdx; ++iOrigIdx) {
        // If the customer at position <iOrigIdx> hasn't been removed, then
        // copy it to the new route.
        if(!vecIsCustRemoved[m_vecRoute[iOrigIdx]]) {
            vecNewRoute[iNewIdx] = m_vecRoute[iOrigIdx];
            ++iNewIdx;
        } else {
            // Customer is removed. Update delta and load.
            int iCustId = m_vecRoute[iOrigIdx];
            dDelta += m_pMatDists->getElement(vecNewRoute[iNewIdx - 1], m_vecRoute[iOrigIdx + 1]) - m_pMatDists->getElement(vecNewRoute[iNewIdx - 1], iCustId) -
                      m_pMatDists->getElement(iCustId, m_vecRoute[iOrigIdx + 1]);
            m_iLoad -= m_pInstance->getDemand(iCustId);
            iServTimeDelta -= m_pInstance->getServiceTime(iCustId);
            vecRemovedCusts.push_back(iCustId);
        }
    }
    m_dCost += dDelta;
    m_dDuration += dDelta + iServTimeDelta;
    vecNewRoute[iNewIdx] = m_pInstance->getN() + 1;
    ++iNewIdx;
    // Shrink the new vector.
    vecNewRoute.resize(iNewIdx);
    m_vecRoute.swap(vecNewRoute);
}

void Route::removeSequence(int iStartPos, int iNToRemove) {
    double dDeltaDist = 0;
    int iServTimeDelta = 0;
    int iDeltaLoad = 0;
    int i;
    for(i = 0; i < iNToRemove; ++i) {
        int iPos = iStartPos + i;
        dDeltaDist -= m_pMatDists->getElement(m_vecRoute[iPos - 1], m_vecRoute[iPos]);
        iServTimeDelta -= m_pInstance->getServiceTime(m_vecRoute[iPos]);
        iDeltaLoad -= m_pInstance->getDemand(m_vecRoute[iPos]);
    }
    dDeltaDist -= m_pMatDists->getElement(m_vecRoute[iStartPos + iNToRemove - 1], m_vecRoute[iStartPos + iNToRemove]);
    dDeltaDist += m_pMatDists->getElement(m_vecRoute[iStartPos - 1], m_vecRoute[iStartPos + iNToRemove]);

    m_dCost += dDeltaDist;
    m_iLoad += iDeltaLoad;
    m_dDuration += dDeltaDist + iServTimeDelta;
    m_vecRoute.erase(m_vecRoute.begin() + iStartPos, m_vecRoute.begin() + (iStartPos + iNToRemove));
}

void Route::concatRoute(const Route& otherRoute) {
    // cout << "concatRoute" << endl;
    int iN = m_pInstance->getN();

    double saving = m_pMatDists->getElement(getLastCustomer(), otherRoute.getFirstCustomer()) - m_pMatDists->getElement(getLastCustomer(), iN + 1) -
                    m_pMatDists->getElement(0, otherRoute.getFirstCustomer());

    m_dCost += otherRoute.getCost() + saving;
    m_iLoad += otherRoute.getLoad();
    m_dDuration += otherRoute.getDuration() + saving;

    m_vecRoute.insert(m_vecRoute.begin() + size() - 1, otherRoute.m_vecRoute.begin() + 1, otherRoute.m_vecRoute.begin() + otherRoute.size() - 1);
}

// concatenate <otherRoute> to this route, inversing the order of <otherRoute>
void Route::concatRouteInv(const Route& otherRoute) {
    // cout << "concatRouteInv" << endl;
    int iN = m_pInstance->getN();

    double saving = m_pMatDists->getElement(getLastCustomer(), otherRoute.getLastCustomer()) - m_pMatDists->getElement(getLastCustomer(), iN + 1) -
                    m_pMatDists->getElement(otherRoute.getLastCustomer(), iN + 1);

    m_dCost += otherRoute.getCost() + saving;
    m_iLoad += otherRoute.getLoad();
    m_dDuration += otherRoute.getDuration() + saving;

    int i;
    m_vecRoute.pop_back();
    for(i = otherRoute.size() - 2; i > 0; i--)
        m_vecRoute.push_back(otherRoute.getNodeId(i));
    m_vecRoute.push_back(iN + 1);
}

// concatenate <otherRoute> to this route, inversing the order of <this> route
void Route::concatRouteInvThis(const Route& otherRoute) {
    // cout << "concatRouteInvThis" << endl;
    int iN = m_pInstance->getN();

    double saving = m_pMatDists->getElement(getFirstCustomer(), otherRoute.getFirstCustomer()) - m_pMatDists->getElement(0, getFirstCustomer()) -
                    m_pMatDists->getElement(0, otherRoute.getFirstCustomer());

    m_dCost += otherRoute.getCost() + saving;
    m_iLoad += otherRoute.getLoad();
    m_dDuration += otherRoute.getDuration() + saving;

    reverse(m_vecRoute.begin(), m_vecRoute.end());
    // Now the first node on m_vecRoute is n+1 and the last node is 0. Correct this:
    m_vecRoute.front() = 0;
    m_vecRoute.back() = iN + 1;
    m_vecRoute.insert(m_vecRoute.begin() + size() - 1, otherRoute.m_vecRoute.begin() + 1, otherRoute.m_vecRoute.begin() + otherRoute.size() - 1);
}

void Route::consCalc() {
    int iPos, iNodeId;
    int iMaxPos = (int)m_vecRoute.size();
    assert(iMaxPos >= 2);
    double dDist = 0;
    int iServTime = 0;
    int iLoad = 0;
    for(iPos = 1; iPos < iMaxPos; iPos++) {
        iNodeId = m_vecRoute[iPos];
        dDist += m_pMatDists->getElement(m_vecRoute[iPos - 1], iNodeId);
        iServTime += m_pInstance->getServiceTime(iNodeId);
        iLoad += m_pInstance->getDemand(iNodeId);
    }
    if(iLoad > m_pInstance->getMaxLoad()) {
        assert(false);
        error("Route::consCalc()", "Load is violated");
    }
    if(dDist + iServTime > m_pInstance->getMaxDuration())
        error("Route::consCalc()", "Duration is violated");
    if(iLoad != m_iLoad)
        error("Route::consCalc()", "Calculated (" + int2String(iLoad) + ") and stored load (" + int2String(m_iLoad) + ") mismatch!");
    if(fabs(dDist - m_dCost) > m_dToleranceEps) {
        assert(false);
        error("Route::consCalc()", "Calculated and stored cost mismatch!");
    }
    if(fabs(dDist + iServTime - m_dDuration) > m_dToleranceEps)
        error("Route::consCalc()", "Calculated and stored duration mismatch!");
    // Use the calculated distance. It is most precise (it's the result of the fewest number
    // of calculations).
    m_dCost = dDist;
    m_dDuration = dDist + iServTime;
}

// Perform two opt local search on the route. Return true iff the route was improved.
// The method only works for symmetric distances.
bool Route::twoOpt(double dImproveEps) {
#ifndef BATCH_MODE
    consCalc();
#endif

    int i, j, iBestI, iBestJ;
    int iRouteLength = (int)m_vecRoute.size();
    bool bImproved = false;
    bool bImprovedThisIter;
    do {
        bImprovedThisIter = false;
        double dBestDelta = DBL_MAX;
        // i indexes the first edge to remove. We are removing the edge (m_vecRoute[i-1], m_vecRoute[i]), that is
        // the edge between the (i-1)'th and i'th node in the tour. As the last node in the tour is node
        // iRouteLength-1 and there has to be one edge between the two edges removed the last edge in the tour that
        // can be the "first edge" is the joining the (iRouteLength-4)'th and (iRouteLength-3)'th nodes as the
        // second edge will joing ((iRouteLength-2)'th and (iRouteLength-1)'th
        for(i = 1; i <= iRouteLength - 3; i++) {
            // i indexes the second edge to remove. We are removing the edge between the (j-1)'th and j'th node in the tour.
            for(j = i + 2; j <= iRouteLength - 1; j++) {
                double dDelta = -m_pMatDists->getElement(m_vecRoute[i - 1], m_vecRoute[i]) - m_pMatDists->getElement(m_vecRoute[j - 1], m_vecRoute[j]) +
                                m_pMatDists->getElement(m_vecRoute[i - 1], m_vecRoute[j - 1]) + m_pMatDists->getElement(m_vecRoute[i], m_vecRoute[j]);

                // DISTANCE MATRIX NOT SYMMETRIC => NEED THIS:
                for(auto k = i; k < j - 1; ++k) {
                    dDelta -= m_pMatDists->getElement(m_vecRoute[k], m_vecRoute[k + 1]);
                    dDelta += m_pMatDists->getElement(m_vecRoute[k + 1], m_vecRoute[k]);
                }

                if(dDelta < dBestDelta) {
                    dBestDelta = dDelta;
                    iBestI = i;
                    iBestJ = j;
                }
            }
        }
        if(dBestDelta < -dImproveEps) {
#ifndef BATCH_MODE
            consCalc();
#endif

            bImproved = true;
            bImprovedThisIter = true;
            do2optMove(iBestI, iBestJ);

#ifndef BATCH_MODE
            consCalc();
#endif
        }
    } while(bImprovedThisIter);
    return bImproved;
}

void Route::do2optMove(int i, int j) {
    assert(i < j);

    double dDelta = -m_pMatDists->getElement(m_vecRoute[i - 1], m_vecRoute[i]) - m_pMatDists->getElement(m_vecRoute[j - 1], m_vecRoute[j]) +
                    m_pMatDists->getElement(m_vecRoute[i - 1], m_vecRoute[j - 1]) + m_pMatDists->getElement(m_vecRoute[i], m_vecRoute[j]);

    // DISTANCE MATRIX NOT SYMMETRIC => NEED THIS:
    for(auto k = i; k < j - 1; ++k) {
        dDelta -= m_pMatDists->getElement(m_vecRoute[k], m_vecRoute[k + 1]);
        dDelta += m_pMatDists->getElement(m_vecRoute[k + 1], m_vecRoute[k]);
    }

    m_dCost += dDelta;
    m_dDuration += dDelta;

    // Update the route. All that the 2 opt move does is reversing the segment
    // (r[i], r[i+1], ... r[j-1]) where r[i] is the node at position i in the route.
    reverse(m_vecRoute.begin() + i, m_vecRoute.begin() + j);
}

// Set the nodes of this route to be the nodes of <vecNodes>
void Route::setNodes(const vector<int>& vecNodes) {
    assert(vecNodes.size() >= 2);
    assert(vecNodes.front() == 0);
    assert(vecNodes.back() == m_pInstance->getN() + 1);

    m_vecRoute.assign(vecNodes.begin(), vecNodes.end());

    int iPos, iNodeId;
    int iMaxPos = (int)m_vecRoute.size();
    double dDist = 0;
    int iServTime = 0;
    m_iLoad = 0;
    for(iPos = 1; iPos < iMaxPos; iPos++) {
        iNodeId = m_vecRoute[iPos];
        dDist += m_pMatDists->getElement(m_vecRoute[iPos - 1], iNodeId);
        iServTime += m_pInstance->getServiceTime(iNodeId);
        m_iLoad += m_pInstance->getDemand(iNodeId);
    }
    m_dCost = dDist;
    m_dDuration = dDist + iServTime;
}

// Set the customers of this route to be those of <vecCustomers>
// (the method is similar to <setNodes(...)>, but in the method below <vecCustomers> should not contain the start and end depot).
void Route::setCustomers(const vector<int>& vecCustomers) {
    assert(vecCustomers.front() != 0);
    assert(vecCustomers.back() != m_pInstance->getN() + 1);

    m_vecRoute.clear();
    m_vecRoute.push_back(0);
    m_vecRoute.insert(m_vecRoute.end(), vecCustomers.begin(), vecCustomers.end());
    m_vecRoute.push_back(m_pInstance->getN() + 1);

    consCalc();
}

ostream& operator<<(ostream& os, const Route& route) {
    os << "{ ";
    outputVector(os, route.m_vecRoute);
    os << " }, cost: " << route.getCost() << ", load: " << route.getLoad() << ", duration: " << route.getDuration();
    return os;
}
