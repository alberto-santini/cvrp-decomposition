#include "CVRPSolution.h"

// From SRP-Utils
#include "GraphObjects.h"
#include "HashNumberGenerator.h"
#include "Utils.h"
#include "VectorUtils.h"

#include <algorithm>
#include <array>
#include <cfloat>
#include <fstream>
#include <limits>
#include <vector>

using namespace std;

CVRPSolution::CVRPSolution(const CVRPInstance& instance) : m_pInstance(&instance) {
    m_dCost = 0;
    m_dUnassignedCost = 10000.0;
    m_dImproveEps = 1E-7;
    m_dToleranceEps = 1E-7;
    m_vecRoutes.reserve(instance.getN());
    m_vecCust2Route = vector<int>(instance.getN() + 1, -1);
    // All customers are initially not assigned to a route.
    makeRange(m_vecUnassigned, 1, instance.getN());
    m_dCost = instance.getN() * m_dUnassignedCost;
}

/***********************************************************************
 * unsigned int CVRPSolution::getHashCode();
 *
 * Returns a hash code that represents the solution. The hash code is calculated
 * as sum(i=0 to |seq|) (hashNumbers[i]*seq[i]). Notice that the hash code is
 * dependent of the hashNumbers. If hashNumbers is changed then so is the hash
 * code.
 *
 * --- input parameters ---
 * --- return values ---
 * The hash code of the solution represented by this object.
 ***********************************************************************/

unsigned int CVRPSolution::getHashCode() const {
    vector<int> vecSeq;
    getHashSequence(vecSeq);
    unsigned int hashKey = 0;
    const vector<unsigned int>& hashNumbers = FixedHashNumberGenerator::getHashNumbers((int)vecSeq.size());
    for(unsigned int i = 0; i < vecSeq.size(); i++)
        hashKey += hashNumbers[i] * vecSeq[i];
    return hashKey;
}

/***********************************************************************
 * void CVRPSolution::getHashSequence(vector<int> &vecOrderId);
 *
 * Returns the order ids of the solution in an unique order. That
 * is if the sequence of two solutions is identical then the two solutions
 * are identical. Used for generating hash codes.
 *
 * --- input parameters ---
 * vecOrderId		: return parameter. A vector of order ids.
 *					  Unique to this solution
 *
 * --- return values ---
 ***********************************************************************/

// Method for comparing vectors by just comparing the first element. Assumes that both vectors are non-empty
bool vectorCompLess(const vector<int>& vecA, const vector<int>& vecB) { return vecA.front() < vecB.front(); }

// @@@ to do: Handle situation when CVRPSolution has unassigned customers.
void CVRPSolution::getHashSequence(vector<int>& vecOrderId) const {
    int iN = m_pInstance->getN();
    vecOrderId.reserve(iN + m_vecRoutes.size());

    int i, iNRoutes = (int)m_vecRoutes.size();
    vector<vector<int>> vecRoutes;
    vecRoutes.reserve(iNRoutes);
    bool bSymmetricCVRP = m_pInstance->isSymmetric();
    for(i = 0; i < iNRoutes; ++i) {
        const Route& route = m_vecRoutes[i];
        vector<int> vecRoute(route.getNodes());
        vecRoute.pop_back();
        vecRoute.erase(vecRoute.begin());
        // If we are solving a symmetric CVRP then we want the smallest node id at the front of the vector.
        if(bSymmetricCVRP && vecRoute.size() >= 2 && vecRoute.back() < vecRoute.front())
            reverseVector(vecRoute);
        if(!vecRoute.empty())
            vecRoutes.push_back(vecRoute);
    }
    // Sort vector by first element.
    sort(vecRoutes.begin(), vecRoutes.end(), vectorCompLess);

    vecOrderId.clear();
    for(i = 0; i < (int)vecRoutes.size(); ++i) {
        vecOrderId.insert(vecOrderId.end(), vecRoutes[i].begin(), vecRoutes[i].end());
        if(i < (int)vecRoutes.size() - 1)
            vecOrderId.push_back(0);
    }
}

void CVRPSolution::createRoute(int iCustId) {
    m_vecRoutes.push_back(Route(*m_pInstance));
    m_vecRoutes.back().insertCust(iCustId, 1);
    m_dCost += m_vecRoutes.back().getCost();
    m_vecCust2Route[iCustId] = (int)m_vecRoutes.size() - 1;
    // <iCustId> are no longer unassigned. Remove him from <m_vecUnassigned>
    vector<int>::iterator it = find(m_vecUnassigned.begin(), m_vecUnassigned.end(), iCustId);
    if(it == m_vecUnassigned.end())
        error("CVRPSolution::createRoute(...)", "Customer was not unassigned!");
    unorderedVectorRemoveAt(m_vecUnassigned, (int)(it - m_vecUnassigned.begin()));
    m_dCost -= m_dUnassignedCost;
}

// Add an empty route.
void CVRPSolution::createRoute() {
    m_vecRoutes.push_back(Route(*m_pInstance));
    m_dCost += m_vecRoutes.back().getCost();
}

void CVRPSolution::createRoutes(int iNRoutes) {
    int i;
    for(i = 0; i < iNRoutes; ++i)
        createRoute();
}

void CVRPSolution::mergeRoutes(int iRoute1, int iRoute2) {
    Route &route1 = m_vecRoutes[iRoute1], &route2 = m_vecRoutes[iRoute2];
    double costBefore = route1.getCost();
    route1.concatRoute(route2);
    m_dCost += route1.getCost() - costBefore - route2.getCost();
    m_vecRoutes.erase(m_vecRoutes.begin() + iRoute2);
    int i;
    for(i = 1; i <= m_pInstance->getN(); i++) {
        if(m_vecCust2Route[i] == iRoute2)
            m_vecCust2Route[i] = iRoute1;
        if(m_vecCust2Route[i] > iRoute2)
            m_vecCust2Route[i]--;
    }
}

void CVRPSolution::mergeRoutesInv(int iRoute1, int iRoute2) {
    Route &route1 = m_vecRoutes[iRoute1], &route2 = m_vecRoutes[iRoute2];
    double costBefore = route1.getCost();
    route1.concatRouteInv(route2);
    m_dCost += route1.getCost() - costBefore - route2.getCost();
    m_vecRoutes.erase(m_vecRoutes.begin() + iRoute2);
    int i;
    for(i = 1; i <= m_pInstance->getN(); i++) {
        if(m_vecCust2Route[i] == iRoute2)
            m_vecCust2Route[i] = iRoute1;
        if(m_vecCust2Route[i] > iRoute2)
            m_vecCust2Route[i]--;
    }
}

void CVRPSolution::mergeRoutesInvThis(int iRoute1, int iRoute2) {
    Route &route1 = m_vecRoutes[iRoute1], &route2 = m_vecRoutes[iRoute2];
    double costBefore = route1.getCost();
    route1.concatRouteInvThis(route2);
    m_dCost += route1.getCost() - costBefore - route2.getCost();
    m_vecRoutes.erase(m_vecRoutes.begin() + iRoute2);
    int i;
    for(i = 1; i <= m_pInstance->getN(); i++) {
        if(m_vecCust2Route[i] == iRoute2)
            m_vecCust2Route[i] = iRoute1;
        if(m_vecCust2Route[i] > iRoute2)
            m_vecCust2Route[i]--;
    }
}

void CVRPSolution::clear() {
    m_vecRoutes.clear();
    m_vecCust2Route.assign(m_pInstance->getN() + 1, -1);
    // no customers is assigned to a route.
    makeRange(m_vecUnassigned, 1, m_pInstance->getN());
    m_dCost = m_pInstance->getN() * m_dUnassignedCost;
}

void CVRPSolution::clear(int iNRoutes) {
    // If the number of vehicles is fixed, then <iNRoutes> should be equal to the number of vehicles.
    assert(!m_pInstance->getNVehiclesFixed() || (iNRoutes == m_pInstance->getNVehicles()));
    m_vecRoutes.assign(iNRoutes, Route(*m_pInstance));
    int iN = m_pInstance->getN();
    m_vecCust2Route.assign(iN + 1, -1);
    // no customers is assigned to any route.
    makeRange(m_vecUnassigned, 1, iN);
    m_dCost = m_pInstance->getN() * m_dUnassignedCost;
    m_dCost += iNRoutes * m_pInstance->getDist(0, iN + 1);
}

// So far, this method only relocates to other routes, not within a route.
// @@@ The method could be speeded up dramatically by using caching.
// Returns true if the solution was improved.
bool CVRPSolution::relocateOpt() {
    int iRoute, iOtherRoute, iPos, iCustId;
    int iMaxDur = m_pInstance->getMaxDuration();
    int iNRoutes = (int)m_vecRoutes.size();
    // Variables storing the best move.
    int iBestFromRoute = -1, iBestToRoute = -1, iBestFromPos = -1, iBestToPos = -1;
    double dBestDelta = 0;
    bool bOk = true;
    bool bImproved = false;
    while(bOk) {
        // This variable is used to indicate if a move was found.
        iBestFromRoute = -1;
        dBestDelta = 0;
        for(iRoute = 0; iRoute < iNRoutes; ++iRoute) {
            Route& route = m_vecRoutes[iRoute];
            // If this route only contains one customer and the number of vehicles is fixed,
            // then don't remove the only customer.
            if(m_pInstance->getNVehiclesFixed() && route.size() <= 3)
                continue;
            for(iPos = 1; iPos < route.size() - 1; ++iPos) {
                iCustId = route.getNodeId(iPos);
                double dRemoveDelta = route.getRemoveDelta(iPos);
                // Is the remove ok, or would it violate the max duration constraint (in case of
                // distances that does not satisfy the triangle inequality)?
                if(dRemoveDelta + route.getDuration() <= iMaxDur) {
                    // It's ok to remove the customer. Try to insert it on all other routes.
                    for(iOtherRoute = 0; iOtherRoute < iNRoutes; ++iOtherRoute) {
                        // Skip the route where the customer originated from.
                        if(iOtherRoute == iRoute)
                            continue;

                        int iLocalBestPos;
                        double dLocalBestInc;
                        if(m_vecRoutes[iOtherRoute].findBestInsertion(iCustId, iLocalBestPos, dLocalBestInc)) {
                            if(dRemoveDelta + dLocalBestInc < dBestDelta - m_dImproveEps) {
                                dBestDelta = dRemoveDelta + dLocalBestInc;
                                iBestFromRoute = iRoute;
                                iBestToRoute = iOtherRoute;
                                iBestFromPos = iPos;
                                iBestToPos = iLocalBestPos;
                            }
                        }
                    }
                }
            }
        }
        // Did we find an improving move?
        if(iBestFromRoute >= 0) {
            // cout << "Cost improved by: " << dBestDelta << endl;
            // Yes! Perform the move.
            bImproved = true;
            iCustId = m_vecRoutes[iBestFromRoute].getNodeId(iBestFromPos);
            // Remove the customer from it's original route.
            m_vecRoutes[iBestFromRoute].removeCustomer(iBestFromPos);
            // Insert it in the new route.
            m_vecRoutes[iBestToRoute].insertCust(iCustId, iBestToPos);
            // Update the cost of the solution:
            m_dCost += dBestDelta;
            // The customer switched route.
            m_vecCust2Route[iCustId] = iBestToRoute;
        } else
            bOk = false;
    }
    // Check if solution is ok after the relocate opt (catching bugs and adjusting the cost
    // to a more precise value).
    consCalc();
    return bImproved;
}

// So far, this method only swaps between different routes, not within a route.
// @@@ The method could be speeded up dramatically by using caching.
// Returns true if the solution was improved.
bool CVRPSolution::swapOpt() {
    const Matrix<double>& matDists = m_pInstance->getDists();
    int iRoute1, iRoute2, iPos1, iPos2, iCustId1, iCustId2;
    int iMaxDur = m_pInstance->getMaxDuration();
    int iMaxLoad = m_pInstance->getMaxLoad();
    int iN = m_pInstance->getN();
    int iNRoutes = (int)m_vecRoutes.size();
    double dRemoveDist1, dRemoveDist2;
    // Variables storing the best move.
    int iBestRoute1 = -1, iBestRoute2 = -1, iBestPos1 = -1, iBestPos2 = -1;
    double dBestDelta;
    bool bOk = true;
    bool bImproved = false;
    // vecRemoveDists[i] stores d(pred(i), i) + d(i, succ(i))
    // where <i> is a customer id, pred(i) is the id of the node preceding i in the current solution,
    // succ(i) is the id of the node succeding i in the current solution and d( , ) is the distance
    // function.
    vector<double> vecRemoveDists(iN + 1, -DBL_MAX);
    while(bOk) {
        // calculate <vecRemoveDists>
        for(iRoute1 = 0; iRoute1 < iNRoutes; ++iRoute1) {
            Route& route1 = m_vecRoutes[iRoute1];
            for(iPos1 = 1; iPos1 < route1.size() - 1; ++iPos1) {
                iCustId1 = route1.getNodeId(iPos1);
                // vecRemoveDists[iCustId1] = d(v[iPos-1], v[iPos-1]
                vecRemoveDists[iCustId1] = -distIJK(matDists, route1[iPos1 - 1], iCustId1, route1[iPos1 + 1]);
            }
        }
        // This variable is used to indicate if a move was found.
        iBestPos1 = -1;
        dBestDelta = 0;
        for(iRoute1 = 0; iRoute1 < iNRoutes; ++iRoute1) {
            Route& route1 = m_vecRoutes[iRoute1];
            // If this route only contains one customer and the number of vehicles is fixed,
            // then don't remove the only customer.
            for(iPos1 = 1; iPos1 < route1.size() - 1; ++iPos1) {
                iCustId1 = route1.getNodeId(iPos1);
                dRemoveDist1 = vecRemoveDists[iCustId1];
                int iDemand1 = m_pInstance->getDemand(iCustId1);
                int iServiceTime1 = m_pInstance->getServiceTime(iCustId1);
                int iLoad1Without = route1.getLoad() - iDemand1;
                // Notice that we avoid testing the same swap twice by letting
                // iRoute2 start from iRoute1+1
                for(iRoute2 = iRoute1 + 1; iRoute2 < iNRoutes; ++iRoute2) {
                    Route& route2 = m_vecRoutes[iRoute2];
                    int iRoute2Load = route2.getLoad();
                    for(iPos2 = 1; iPos2 < route2.size() - 1; ++iPos2) {
                        iCustId2 = route2.getNodeId(iPos2);
                        int iDemand2 = m_pInstance->getDemand(iCustId2);
                        // Check capacities
                        if(iLoad1Without + iDemand2 > iMaxLoad || iRoute2Load - iDemand2 + iDemand1 > iMaxLoad) {
                            // Capacity violation.
                            continue;
                        }
                        double dInsertDist1 = distIJK(matDists, route1[iPos1 - 1], iCustId2, route1[iPos1 + 1]);
                        double dInsertDist2 = distIJK(matDists, route2[iPos2 - 1], iCustId1, route2[iPos2 + 1]);
                        dRemoveDist2 = vecRemoveDists[iCustId2];
                        double dDelta = dInsertDist1 + dInsertDist2 + dRemoveDist1 + dRemoveDist2;
                        if(dDelta >= dBestDelta - m_dImproveEps) {
                            // Delta cost is not interesting
                            continue;
                        }
                        // Check violation of durations
                        int iServiceTime2 = m_pInstance->getServiceTime(iCustId2);
                        if(route1.getDuration() + dInsertDist1 + dRemoveDist1 + iServiceTime2 - iServiceTime1 > iMaxDur ||
                           route2.getDuration() + dInsertDist2 + dRemoveDist2 + iServiceTime1 - iServiceTime2 > iMaxDur) {
                            // Duration is violated.
                            continue;
                        }
                        // If we reach this spot, then the move satisfy the constraints
                        // and the delta cost is better than the best delta observed so far.
                        // Remember this move:
                        dBestDelta = dDelta;
                        iBestRoute1 = iRoute1;
                        iBestRoute2 = iRoute2;
                        iBestPos1 = iPos1;
                        iBestPos2 = iPos2;
                    }
                }
            }
        }
        // Did we find an improving move?
        if(iBestPos1 >= 0) {
            // Yes! Perform it:
            bImproved = true;
            iCustId1 = m_vecRoutes[iBestRoute1].getNodeId(iBestPos1);
            iCustId2 = m_vecRoutes[iBestRoute2].getNodeId(iBestPos2);
            m_vecRoutes[iBestRoute1].removeCustomer(iBestPos1);
            m_vecRoutes[iBestRoute2].removeCustomer(iBestPos2);
            m_vecRoutes[iBestRoute1].insertCust(iCustId2, iBestPos1);
            m_vecRoutes[iBestRoute2].insertCust(iCustId1, iBestPos2);
            // Update cost:
            m_dCost += dBestDelta;
            // update <m_vecCust2Route>:
            m_vecCust2Route[iCustId1] = iBestRoute2;
            m_vecCust2Route[iCustId2] = iBestRoute1;
        } else
            bOk = false;
    }
    consCalc();
    return bImproved;
}

// @@@ Notice: In the sequential search we assume that distances to/from 0 are the same as distance to/from node n+1. This
// allows us to only store one depot in between each tour and makes the code simpler.
void CVRPSolution::getGiantTour(vector<int>& vecGiantTour, vector<int>& vecNodeToPos, vector<int>& vecRouteId, vector<int>& vecDemandBefore,
                                vector<int>& vecDemandAfter, vector<double>& vecDurationBefore, vector<double>& vecDurationAfter) {
    int i, j, iNRoutes = (int)m_vecRoutes.size();
    int iN = m_pInstance->getN();
    vecGiantTour.clear();

    vecGiantTour.resize(iN + iNRoutes);
    vecNodeToPos.resize(iN + 1);
    vecNodeToPos[0] = -1; // Depot occurs multiple times so no position is assigned to the depot

    int idx = 0;
    for(i = 0; i < iNRoutes; ++i) {
        const Route& route = m_vecRoutes[i];
        int iSumDemand = 0;
        double dSumDuration = 0;
        if(!route.empty()) {
            // "route.size()-1": Because we do not add end depot
            for(j = 0; j < route.size() - 1; ++j) {
                int iNodeId = route.getNodeId(j);
                iSumDemand += m_pInstance->getDemand(iNodeId);
                if(j > 0)
                    dSumDuration += m_pInstance->getDist(route.getNodeId(j - 1), iNodeId) + m_pInstance->getServiceTime(iNodeId);
                vecGiantTour[idx] = iNodeId;
                vecRouteId[idx] = i;
                vecDemandBefore[idx] = iSumDemand;
                vecDurationBefore[idx] = dSumDuration;
                vecNodeToPos[iNodeId] = idx;
                ++idx;
            }
            // Process the tour backwards to set <vecDemandAfter> and <vecDurationAfter>
            iSumDemand = 0;
            dSumDuration = 0;
            for(j = route.size() - 2; j >= 0; --j) {
                int iNodeId = route.getNodeId(j);
                dSumDuration += m_pInstance->getDist(iNodeId, route.getNodeId(j + 1)) + m_pInstance->getServiceTime(iNodeId);
                iSumDemand += m_pInstance->getDemand(iNodeId);
                vecDemandAfter[idx] = iSumDemand;
                vecDurationAfter[idx] = dSumDuration;
            }
        }
    }
    // shrink the tour if there were some empty routes:
    vecGiantTour.resize(idx);
}

void CVRPSolution::seqSearchTailSwapOpt() {
    const Matrix<double>& matDists = m_pInstance->getDists();
    int iVehCap = m_pInstance->getMaxLoad();
    double dMaxDur = m_pInstance->getMaxDuration();

    vector<int> vecGiantTour;
    vector<int> vecNodeToPos;
    vector<int> vecRouteId;
    vector<int> vecDemandBefore;
    vector<int> vecDemandAfter;
    vector<double> vecDurationBefore;
    vector<double> vecDurationAfter;

    getGiantTour(vecGiantTour, vecNodeToPos, vecRouteId, vecDemandBefore, vecDemandAfter, vecDurationBefore, vecDurationAfter);
    int iSizeGTour = (int)vecGiantTour.size();
    int iN = m_pInstance->getN();
    int i1, i2, i1Best = -1, i2Best = -1;
    double dGStar = m_dImproveEps;
    for(i1 = 0; i1 < iSizeGTour - 1; ++i1) {
        int t1 = vecGiantTour[i1];
        int t2 = vecGiantTour[i1 + 1];
        double dB1 = matDists.getElement(t1, t2) - m_dImproveEps / 2;
        const vector<int>& vecNeighbors = m_pInstance->getNearestNodes(t2);
        for(int idx = 0; idx < (int)vecNeighbors.size(); ++idx) {
            int t3 = vecNeighbors[idx];
            if(matDists.getElement(t2, t3) >= dB1)
                break;
            else {
                i2 = vecNodeToPos[t3];
                if(t3 == 0) {
                    // Special case when t3 is the depot.
                }
                int t4 = vecGiantTour[i2 + 1];
                double dG = matDists.getElement(t1, t2) - matDists.getElement(t3, t2) + matDists.getElement(t3, t4) - matDists.getElement(t1, t4);
                if(dG > dGStar) {
                    // Check feasibility.
                    // Step 1: Test that there is a depot in between node i1+1 and i2 and between node i2+1 and i1
                    if(vecRouteId[i1 + 1] != vecRouteId[i2] && vecRouteId[i2 + 1] != vecRouteId[i1]) {
                        // Ok. passed that check. Check demands:
                        if(vecDemandBefore[i1] + vecDemandAfter[i2 + 1] <= iVehCap && vecDemandBefore[i2] + vecDemandAfter[i1 + 1] <= iVehCap) {
                            // Passed the demand check as well. Check duration
                            if(vecDurationBefore[i1] + vecDurationAfter[i2 + 1] + matDists.getElement(t1, t4) <= dMaxDur &&
                               vecDurationBefore[i2] + vecDurationAfter[i1 + 1] + matDists.getElement(t3, t2) <= dMaxDur) {
                                // Ok. The move is feasible and better than what we have seen so far. Record it.
                                i1Best = i1;
                                i2Best = i2;
                            }
                        }
                    }
                }
            }
        }
    }
}

enum TailSwapBreakCase { CASE1 = 1, CASE2 = 2, CASE3 = 3, CASE4 = 4 };

struct TailSwapCacheEntry {
    int m_iPos1;
    int m_iPos2;
    //         m_bBestIsFrontBack
    //      True:       |        False:
    //   ---O---O---    |     ---O---O---
    //  /           \   |    /           \ 
// X             X  |   X             X
    //  \           /   |    \           /
    //   ---O---O---    |     ---O---O---
    //                  |
    //     =>           |       =>
    //                  |
    //   ---O   O---    |     ---O   O---
    //  /    \ /    \   |    /   |   |   \ 
// X      |      X  |   X    |   |    X
    //  \    / \    /   |    \   |   |   /
    //   ---O   O---    |     ---O   O---
    bool m_bIsFrontBack;
    TailSwapBreakCase frontBackFalseCase; // How to reverse in case m_bIsFrontBack is false.
};

std::ostream& operator<<(std::ostream& o, const TailSwapCacheEntry& c) {
    o << "{ pos1: " << c.m_iPos1 << ", pos2: " << c.m_iPos2 << ", ";
    o << "frontback: " << std::boolalpha << c.m_bIsFrontBack;

    if(c.m_bIsFrontBack) {
        o << " }";
    } else {
        o << ", case: ";
        switch(c.frontBackFalseCase) {
        case CASE1: o << "1 }"; break;
        case CASE2: o << "2 }"; break;
        case CASE3: o << "3 }"; break;
        case CASE4: o << "4 }"; break;
        default: o << "? }"; break;
        }
    }
    return o;
}

// std::ostream& operator<<(std::ostream& o, const Route& r) {
// 	for(auto c : r.getNodes()) {
// 		o << c << " ";
// 	}
// 	return o;
// }

struct LoadDuration {
    int m_iLoadBefore;
    int m_iLoadAfter;

    int travelTimeBefore;
    int travelTimeAfter;

    int serviceTimeBefore;
    int serviceTimeAfter;

    int durBefore() const { return travelTimeBefore + serviceTimeBefore; }
    int durAfter() const { return travelTimeAfter + serviceTimeAfter; }
};

std::ostream& operator<<(std::ostream& o, const LoadDuration& l) {
    o << "{ tt before: " << l.travelTimeBefore << ", tt after: " << l.travelTimeAfter << ", dur before: " << l.durBefore() << ", dur after: " << l.durAfter()
      << " }";
    return o;
}

// Precalculates the load and duration of individual route pieces (before and after each "breakpoint").
// Provides a nice speedup compared to recalculating this in <evalTailSwap>
// (p05, 1 threads, 75000 iterations, tailSwap in each iter: 70.8 -> 50.2 sec
//  kelly04, 4 cores, 8 threads, 75000 iterations, tailSwap in each iter: 52.9 -> 39.4 sec)
void CVRPSolution::tailSwapPrecompute(const Route& route, vector<LoadDuration>& vecLoadDuration) const {
    int iRouteSize = (int)route.size();
    vecLoadDuration.resize(iRouteSize - 1);
    const Matrix<double>& matDists = m_pInstance->getDists();

    int iLoadBef = route.getLoad();
    int iLoadAft = 0;
    int iTTBef = 0;
    int iTTAft = 0;
    int iSTBef = 0;
    int iSTAft = 0;

    for(auto i = 0; i < route.getNodes().size() - 1; ++i) {
        iTTBef += matDists.getElement(route[i], route[i + 1]);
        iSTBef += m_pInstance->getServiceTime(route[i]);
    }
    iSTBef += m_pInstance->getServiceTime(route.getNodes().back());

    for(int iPos = iRouteSize - 2; iPos >= 0; --iPos) {
        iLoadBef -= m_pInstance->getDemand(route[iPos + 1]);
        iLoadAft += m_pInstance->getDemand(route[iPos + 1]);

        iTTBef -= matDists.getElement(route[iPos], route[iPos + 1]);
        iSTBef -= m_pInstance->getServiceTime(route[iPos + 1]);

        if(iPos < iRouteSize - 2) {
            iTTAft += matDists.getElement(route[iPos + 1], route[iPos + 2]);
            iSTAft += m_pInstance->getServiceTime(route[iPos + 1]);
        }

        LoadDuration& loadDuration = vecLoadDuration[iPos];
        loadDuration.m_iLoadBefore = iLoadBef;
        loadDuration.m_iLoadAfter = iLoadAft;
        loadDuration.travelTimeAfter = iTTAft;
        loadDuration.travelTimeBefore = iTTBef;
        loadDuration.serviceTimeAfter = iSTAft;
        loadDuration.serviceTimeBefore = iSTBef;
    }
}

// @@@ There are several things that can be done here to speed up this local search.
// -We might be able to stop searching for moves between to routes in an early state if we can determine
// that capacity would be violated for remaining moves (or perhaps even duration if we know that distances
// satisfy the triangle inequality).
// We can make a data structure that stores "dDurBef", "dDurAft" for each position in all routes, this could
// speed up the method somewhat compared to now where those values are calculated (in constant time though, but
// still, we access the distance matrix and performs comparisons in order to compute these values.
bool CVRPSolution::tailSwapOpt() {
    int iRoute1, iRoute2, iNRoutes = (int)m_vecRoutes.size();

    // Caching of results.
    // The cache store what the best tail-swap between any pair of routes is. We only
    // need to update this value if one of the two routes involved in the pair has been changed).
    Matrix<double> matDelta(iNRoutes, -DBL_MAX);
    Matrix<TailSwapCacheEntry> matCache(iNRoutes);

    vector<vector<LoadDuration>> vecLoadDuration(iNRoutes);
    for(iRoute1 = 0; iRoute1 < iNRoutes; ++iRoute1)
        tailSwapPrecompute(m_vecRoutes[iRoute1], vecLoadDuration[iRoute1]);

    bool bImproved = false;
    bool bIterImproved = false;
    do {
        bIterImproved = false;
        double dBestDelta = DBL_MAX;
        TailSwapCacheEntry bestMove;
        int iBestRoute1, iBestRoute2;
        bool checkedWithEmpty1 = false;

        for(iRoute1 = 0; iRoute1 < iNRoutes; ++iRoute1) {
            const Route& route1 = m_vecRoutes[iRoute1];
            bool checkedWithEmpty2 = false;

            if(route1.empty()) {
                if(checkedWithEmpty1) {
                    continue;
                } else {
                    checkedWithEmpty1 = true;
                }
            }

            for(iRoute2 = iRoute1 + 1; iRoute2 < iNRoutes; ++iRoute2) {
                double dDelta = matDelta.getElement(iRoute1, iRoute2);
                const Route& route2 = m_vecRoutes[iRoute2];

                if(route2.empty()) {
                    if(checkedWithEmpty2) {
                        continue;
                    } else {
                        checkedWithEmpty2 = true;
                    }
                }

                // Is cache value valid?
                if(dDelta == -DBL_MAX) {
                    // No. Recalculate:
                    dDelta = evalTailSwap(route1, route2, bestMove, vecLoadDuration[iRoute1], vecLoadDuration[iRoute2]);
                    matDelta.setElement(iRoute1, iRoute2, dDelta);
                    matCache.setElement(iRoute1, iRoute2, bestMove);
                }

                if(dDelta < dBestDelta) {
                    iBestRoute1 = iRoute1;
                    iBestRoute2 = iRoute2;
                    dBestDelta = dDelta;
                }
            }
        }

        if(dBestDelta < -m_dImproveEps) {
            bestMove = matCache.getElement(iBestRoute1, iBestRoute2);

            double dCostBefore = m_dCost;
            auto route1before = m_vecRoutes[iBestRoute1], route2before = m_vecRoutes[iBestRoute2];

            doTailSwap(iBestRoute1, iBestRoute2, bestMove);

            assert(epsilonEqual(dCostBefore + dBestDelta, m_dCost, 1e-6));

            bIterImproved = true;
            bImproved = true;

            // All moves involving iBestRoute1 and iBestRoute2 need to be recalculated. Invalidate cached results:
            int iRoute;
            for(iRoute = 0; iRoute < iNRoutes; ++iRoute) {
                if(iBestRoute1 < iRoute) {
                    matDelta.setElement(iBestRoute1, iRoute, -DBL_MAX);
                } else {
                    matDelta.setElement(iRoute, iBestRoute1, -DBL_MAX);
                }

                if(iBestRoute2 < iRoute) {
                    matDelta.setElement(iBestRoute2, iRoute, -DBL_MAX);
                } else {
                    matDelta.setElement(iRoute, iBestRoute2, -DBL_MAX);
                }
            }

            // Update load and duration cache
            tailSwapPrecompute(m_vecRoutes[iBestRoute1], vecLoadDuration[iBestRoute1]);
            tailSwapPrecompute(m_vecRoutes[iBestRoute2], vecLoadDuration[iBestRoute2]);
        }
    } while(bIterImproved);

    consCalc();
    return bImproved;
}

// Returns true if an improved tail-swap was found
double CVRPSolution::evalTailSwap(const Route& route1, const Route& route2, TailSwapCacheEntry& result, const vector<LoadDuration>& vecLoadDurRoute1,
                                  const vector<LoadDuration>& vecLoadDurRoute2) const {
    double dBestDelta = DBL_MAX;
    const Matrix<double>& matDists = m_pInstance->getDists();
    int iMaxLoad = m_pInstance->getMaxLoad();
    int iMaxDur = m_pInstance->getMaxDuration();

    int iR1Pos, iR2Pos;
    for(iR1Pos = route1.size() - 2; iR1Pos >= 0; --iR1Pos) {
        const LoadDuration& loadDur1 = vecLoadDurRoute1[iR1Pos];
        int iR1Prev = route1[iR1Pos];
        int iR1Next = route1[iR1Pos + 1];
        double dR1DistPrevNext = matDists.getElement(iR1Prev, iR1Next);

        for(iR2Pos = route2.size() - 2; iR2Pos >= 0; --iR2Pos) {
            const LoadDuration& loadDur2 = vecLoadDurRoute2[iR2Pos];

            //////// === First Type: it is the same with both symmetric and asymmetric distances ===

            //   ---O---O---
            //  /           \ 
			// X             X
            //  \           /
            //   ---O---O---
            //
            //     =>
            //
            //   ---O   O---
            //  /    \ /    \ 
			// X      |      X
            //  \    / \    /
            //   ---O   O---

            // check if the move is feasible
            if(loadDur1.m_iLoadBefore + loadDur2.m_iLoadAfter <= iMaxLoad && loadDur2.m_iLoadBefore + loadDur1.m_iLoadAfter <= iMaxLoad &&
               loadDur1.durBefore() + loadDur2.durAfter() + matDists.getElement(iR1Prev, route2[iR2Pos + 1]) < iMaxDur &&
               loadDur2.durBefore() + loadDur1.durAfter() + matDists.getElement(route2[iR2Pos], iR1Next) < iMaxDur) {
                // Check distance of the two new routes obtained by removing edges
                // (route1[iR1Pos], route1[iR1Pos+1]) and (route2[iR2Pos], route2[iR2Pos+1])
                // and exchanging them with edges
                // (route1[iR1Pos], route2[iR2Pos+1]) and (route2[iR2Pos], route1[iR1Pos+1])
                double dDelta = -dR1DistPrevNext - matDists.getElement(route2[iR2Pos], route2[iR2Pos + 1]) + matDists.getElement(iR1Prev, route2[iR2Pos + 1]) +
                                matDists.getElement(route2[iR2Pos], iR1Next);
                // Is move better than what we have seen so far (for these two routes)?
                if(dDelta < dBestDelta) {
                    // Duration and load of the new move is okay. If the number of vehicles is fixed,
                    // then we have to make sure that the move doesn't eliminate a route.
                    if(!m_pInstance->getNVehiclesFixed() || ((iR1Pos > 0 || iR2Pos < route2.size() - 2) && (iR2Pos > 0 || iR1Pos < route1.size() - 2))) {
                        // That was okay as well. We should remember this move.
                        result.m_iPos1 = iR1Pos;
                        result.m_iPos2 = iR2Pos;
                        result.m_bIsFrontBack = true;
                        dBestDelta = dDelta;
                    }
                }
            }

            //////// === Second Type: more complex with asymmetric distances! ===

            //   ---O---O---
            //  /           \ 
			// X             X
            //  \           /
            //   ---O---O---
            //
            //     =>
            //
            //   ---O   O---
            //  /   |   |   \ 
			// X    |   |    X
            //  \   |   |   /
            //   ---O   O---
            //

            std::array<int, 4> delta_case{
                {std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}};

            // HERE WE HAVE TO DISTINGUISH DIFFERENT CASES WHEN THE DISTANCE MATRIX IS NOT SYMMETRIC:

            /*	1. Invert 0..R1POS and R1POS+1..0
             *
             *    0 <--- R1POS  R1POS+1 <--- 0
             *    |        A       |         A
             *    |        |       V         |
             *    L----> R2POS  R2POS+1 ----/
             *
             * 	2. Invert 0..R1POS and R2POS+1..0
             *
             *    0 <--- R1POS  R1POS+1 ---> 0
             *    |        A       A         |
             *    |        |       |         |
             *    L----> R2POS  R2POS+1 <---/
             *
             * 	3. Invert 0..R2POS and R1POS+1..0
             *
             *    0 ---> R1POS  R1POS+1 <--- 0
             *    A        |       |         A
             *    |        V       V         |
             *    L----- R2POS  R2POS+1 ----/
             *
             *  4. Invert 0..R2POS and R2POS+1..0
             *
             *    0 ---> R1POS  R1POS+1 ---> 0
             *    A        A       A         |
             *    |        |       |         |
             *    L----- R2POS  R2POS+1 <---/
             */

            // COMMON TO ALL CASES:
            if(loadDur1.m_iLoadBefore + loadDur2.m_iLoadBefore > iMaxLoad || // Left route
               loadDur1.m_iLoadAfter + loadDur2.m_iLoadAfter > iMaxLoad)     // Right route
            {
                continue;
            }

            // CASE 1:

            //	Left route:
            int dur_after1 = 0;
            for(auto k = iR1Pos; k > 0; --k) {
                dur_after1 += matDists.getElement(route1[k], route1[k - 1]);
            }
            const int dur_before2 = loadDur2.travelTimeBefore;
            const int dur_bridge_2_1 = matDists.getElement(route2[iR2Pos], route1[iR1Pos]);

            // 	Right route:
            int dur_before1p1 = 0;
            for(auto k = route1.size() - 1; k > iR1Pos + 1; --k) {
                dur_before1p1 += matDists.getElement(route1[k], route1[k - 1]);
            }
            const int dur_after2p1 = loadDur2.travelTimeAfter;
            const int dur_bridge_1p1_2p1 = matDists.getElement(route1[iR1Pos + 1], route2[iR2Pos + 1]);

            if(dur_before2 + dur_bridge_2_1 + dur_after1 < iMaxDur &&       // Left route
               dur_before1p1 + dur_bridge_1p1_2p1 + dur_after2p1 < iMaxDur) // Right route
            {
                delta_case[0] = -loadDur1.travelTimeBefore + dur_after1 - matDists.getElement(route1[iR1Pos], route1[iR1Pos + 1]) + dur_bridge_2_1 -
                                loadDur1.travelTimeAfter + dur_before1p1 - matDists.getElement(route2[iR2Pos], route2[iR2Pos + 1]) + dur_bridge_1p1_2p1;
            }

            // CASE 2:

            // Left route: as in CASE 1

            // Right route:
            int dur_before2p1 = 0;
            for(auto k = route2.size() - 1; k > iR2Pos + 1; --k) {
                dur_before2p1 += matDists.getElement(route2[k], route2[k - 1]);
            }
            const int dur_after1p1 = loadDur1.travelTimeAfter;
            const int dur_bridge_2p1_1p1 = matDists.getElement(route2[iR2Pos + 1], route1[iR1Pos + 1]);

            if(dur_before2 + dur_bridge_2_1 + dur_after1 < iMaxDur && // Left route
               dur_before2p1 + dur_bridge_2p1_1p1 + dur_after1p1)     // Right route
            {
                delta_case[1] = -loadDur1.travelTimeBefore + dur_after1 - matDists.getElement(route1[iR1Pos], route1[iR1Pos + 1]) + dur_bridge_2_1 -
                                loadDur2.travelTimeAfter + dur_before2p1 - matDists.getElement(route2[iR2Pos], route2[iR2Pos + 1]) + dur_bridge_2p1_1p1;
            }

            // CASE 3:

            // Left route:
            int dur_after2 = 0;
            for(auto k = iR2Pos; k > 0; --k) {
                dur_after2 += matDists.getElement(route2[k], route2[k - 1]);
            }
            const int dur_before1 = loadDur1.travelTimeBefore;
            const int dur_bridge_1_2 = matDists.getElement(route1[iR1Pos], route2[iR2Pos]);

            // Right route: as in CASE 1

            if(dur_before1 + dur_bridge_1_2 + dur_after2 < iMaxDur &&       // Left route
               dur_before1p1 + dur_bridge_1p1_2p1 + dur_after2p1 < iMaxDur) // Right route
            {
                delta_case[2] = -loadDur2.travelTimeBefore + dur_after2 - matDists.getElement(route1[iR1Pos], route1[iR1Pos + 1]) + dur_bridge_1_2 -
                                loadDur1.travelTimeAfter + dur_before1p1 - matDists.getElement(route2[iR2Pos], route2[iR2Pos + 1]) + dur_bridge_1p1_2p1;
            }

            // CASE 4:

            // Left route: as in CASE 3

            // Right route: as in CASE 2

            if(dur_before1 + dur_bridge_1_2 + dur_after2 < iMaxDur && // Left route
               dur_before2p1 + dur_bridge_2p1_1p1 + dur_after1p1)     // Right route
            {
                delta_case[3] = -loadDur2.travelTimeBefore + dur_after2 - matDists.getElement(route1[iR1Pos], route1[iR1Pos + 1]) + dur_bridge_1_2 -
                                loadDur2.travelTimeAfter + dur_before2p1 - matDists.getElement(route2[iR2Pos], route2[iR2Pos + 1]) + dur_bridge_2p1_1p1;
            }

            const auto min = std::min_element(delta_case.begin(), delta_case.end());

            // Is move better than what we have seen so far (for these two routes)?
            if(*min < dBestDelta) {
                // Duration and load of the new move is okay. If the number of vehicles is fixed,
                // then we have to make sure that the move doesn't eliminate a route.
                if(!m_pInstance->getNVehiclesFixed() || ((iR1Pos > 0 || iR2Pos > 0) && (iR1Pos < route1.size() - 2 || iR2Pos < route2.size() - 2))) {
                    result.m_iPos1 = iR1Pos;
                    result.m_iPos2 = iR2Pos;
                    result.m_bIsFrontBack = false;
                    dBestDelta = *min;

                    switch(min - delta_case.begin()) {
                    case 0: result.frontBackFalseCase = CASE1; break;
                    case 1: result.frontBackFalseCase = CASE2; break;
                    case 2: result.frontBackFalseCase = CASE3; break;
                    case 3: result.frontBackFalseCase = CASE4; break;
                    default: assert(false); break;
                    }
                }
            }
        }
    }

    return dBestDelta;
}

void CVRPSolution::doTailSwap(int iRoute1, int iRoute2, const TailSwapCacheEntry& bestMove) {
    int iPos1 = bestMove.m_iPos1, iPos2 = bestMove.m_iPos2;
    bool bBestIsFrontBack = bestMove.m_bIsFrontBack;
    const auto fbfCase = bestMove.frontBackFalseCase;
    int iN = m_pInstance->getN();
    Route& route1 = m_vecRoutes[iRoute1];
    Route& route2 = m_vecRoutes[iRoute2];
    const vector<int>& vecR1 = route1.getNodes();
    const vector<int>& vecR2 = route2.getNodes();
    vector<int> vecNewRoute1, vecNewRoute2;

    if(bBestIsFrontBack) {
        // Make the new routes.
        vecNewRoute1.insert(vecNewRoute1.end(), vecR1.begin(), vecR1.begin() + iPos1 + 1);
        vecNewRoute1.insert(vecNewRoute1.end(), vecR2.begin() + iPos2 + 1, vecR2.end());
        vecNewRoute2.insert(vecNewRoute2.end(), vecR2.begin(), vecR2.begin() + iPos2 + 1);
        vecNewRoute2.insert(vecNewRoute2.end(), vecR1.begin() + iPos1 + 1, vecR1.end());
    } else {
        switch(fbfCase) {
        case CASE1: {
            vecNewRoute1 = std::vector<int>(vecR2.begin(), vecR2.begin() + iPos2 + 1);
            vecNewRoute1.insert(vecNewRoute1.end(), vecR1.begin() + 1, vecR1.begin() + iPos1 + 1);
            std::reverse(vecNewRoute1.end() - iPos1, vecNewRoute1.end());
            vecNewRoute1.push_back(iN + 1);

            vecNewRoute2 = std::vector<int>(vecR1.begin() + iPos1 + 1, vecR1.end());
            const size_t vecNewRoute2_head_sz = vecNewRoute2.size();
            vecNewRoute2.at(vecNewRoute2_head_sz - 1) = 0;
            vecNewRoute2.insert(vecNewRoute2.end(), vecR2.begin() + iPos2 + 1, vecR2.end());
            std::reverse(vecNewRoute2.begin(), vecNewRoute2.begin() + vecNewRoute2_head_sz);
        }

        break;

        case CASE2: {
            vecNewRoute1 = std::vector<int>(vecR2.begin(), vecR2.begin() + iPos2 + 1);
            vecNewRoute1.insert(vecNewRoute1.end(), vecR1.begin() + 1, vecR1.begin() + iPos1 + 1);
            std::reverse(vecNewRoute1.end() - iPos1, vecNewRoute1.end());
            vecNewRoute1.push_back(iN + 1);

            vecNewRoute2 = std::vector<int>(vecR2.begin() + iPos2 + 1, vecR2.end());
            const size_t vecNewRoute2_head_sz = vecNewRoute2.size();
            vecNewRoute2.at(vecNewRoute2_head_sz - 1) = 0;
            vecNewRoute2.insert(vecNewRoute2.end(), vecR1.begin() + iPos1 + 1, vecR1.end());
            std::reverse(vecNewRoute2.begin(), vecNewRoute2.begin() + vecNewRoute2_head_sz);
        }

        break;

        case CASE3: {
            vecNewRoute1 = std::vector<int>(vecR1.begin(), vecR1.begin() + iPos1 + 1);
            vecNewRoute1.insert(vecNewRoute1.end(), vecR2.begin() + 1, vecR2.begin() + iPos2 + 1);
            std::reverse(vecNewRoute1.end() - iPos2, vecNewRoute1.end());
            vecNewRoute1.push_back(iN + 1);

            vecNewRoute2 = std::vector<int>(vecR1.begin() + iPos1 + 1, vecR1.end());
            const size_t vecNewRoute2_head_sz = vecNewRoute2.size();
            vecNewRoute2.at(vecNewRoute2_head_sz - 1) = 0;
            vecNewRoute2.insert(vecNewRoute2.end(), vecR2.begin() + iPos2 + 1, vecR2.end());
            std::reverse(vecNewRoute2.begin(), vecNewRoute2.begin() + vecNewRoute2_head_sz);
        }

        break;

        case CASE4: {
            vecNewRoute1 = std::vector<int>(vecR1.begin(), vecR1.begin() + iPos1 + 1);
            vecNewRoute1.insert(vecNewRoute1.end(), vecR2.begin() + 1, vecR2.begin() + iPos2 + 1);
            std::reverse(vecNewRoute1.end() - iPos2, vecNewRoute1.end());
            vecNewRoute1.push_back(iN + 1);

            vecNewRoute2 = std::vector<int>(vecR2.begin() + iPos2 + 1, vecR2.end());
            const size_t vecNewRoute2_head_sz = vecNewRoute2.size();
            vecNewRoute2.at(vecNewRoute2_head_sz - 1) = 0;
            vecNewRoute2.insert(vecNewRoute2.end(), vecR1.begin() + iPos1 + 1, vecR1.end());
            std::reverse(vecNewRoute2.begin(), vecNewRoute2.begin() + vecNewRoute2_head_sz);
        }

        break;

        default: assert(false); break;
        }

        // Construct route 1:
        // vecNewRoute1.insert(vecNewRoute1.end(), vecR1.begin(), vecR1.begin()+iPos1+1);
        // vecNewRoute1.insert(vecNewRoute1.end(), vecR2.begin(), vecR2.begin()+iPos2+1);
        // Reverse the segment from route two:
        // reverse(vecNewRoute1.begin()+iPos1+1, vecNewRoute1.end());
        // Exchange the node at the end of <vecNewRoute1> with node n+1 (it is node zero at this point).
        // vecNewRoute1.back() = iN+1;

        // Construct route 2:
        // vecNewRoute2.insert(vecNewRoute2.end(), vecR1.begin()+iPos1+1, vecR1.end());
        // Reverse the segment from route one:
        // reverse(vecNewRoute2.begin(), vecNewRoute2.end());
        // Exchange the node at the start of <vecNewRoute1> with node 0 (it is node n+1at this point).
        // vecNewRoute2.front() = 0;
        // Add the segment from route two.
        // vecNewRoute2.insert(vecNewRoute2.end(), vecR2.begin()+iPos2+1, vecR2.end());
    }

    double dCostBefore = route1.getCost() + route2.getCost();
    route1.setNodes(vecNewRoute1);
    route2.setNodes(vecNewRoute2);
    double dCostAfter = route1.getCost() + route2.getCost();
    m_dCost += dCostAfter - dCostBefore;
    int i;
    for(i = 1; i < route1.size() - 1; ++i)
        m_vecCust2Route[route1[i]] = iRoute1;
    for(i = 1; i < route2.size() - 1; ++i)
        m_vecCust2Route[route2[i]] = iRoute2;
}

bool CVRPSolution::route2opt() {
    int iRoute, iNRoutes = (int)m_vecRoutes.size();
    bool bImproved = false;
    for(iRoute = 0; iRoute < iNRoutes; ++iRoute) {
        Route& route = m_vecRoutes[iRoute];

#ifndef BATCH_MODE
        route.consCalc();
#endif

        double dCost = route.getCost();

        if(route.twoOpt(m_dImproveEps)) {
            double dNewCost = route.getCost();
            m_dCost += dNewCost - dCost;
            bImproved = true;
        }

#ifndef BATCH_MODE
        route.consCalc();
#endif
    }
    return bImproved;
}

// Run both the relocate and the swap local search until we run out of improvements.
void CVRPSolution::combinedOpt(bool bRelocate, bool bSwap, bool b2opt, bool bTailSwap) {
    int iLastImproveIter = 0;
    int iIter = 0;
    int iNMethods = 0;
    int iRelocateId, iSwapId, i2optId, iTailSwapId;
    if(b2opt) {
        i2optId = iNMethods;
        ++iNMethods;
    }
    if(bRelocate) {
        iRelocateId = iNMethods;
        ++iNMethods;
    }
    if(bSwap) {
        iSwapId = iNMethods;
        ++iNMethods;
    }
    if(bTailSwap) {
        iTailSwapId = iNMethods;
        ++iNMethods;
    }

    while(iIter - iLastImproveIter < iNMethods) {
        bool bImproved;
        if(b2opt && (iIter % iNMethods == i2optId)) {
            bImproved = route2opt();
            // cout << "Cost after route 2 opt: " << getCost() << endl;
        }
        if(bRelocate && (iIter % iNMethods == iRelocateId)) {
            bImproved = relocateOpt();
            // cout << "Cost after relocateOpt: " << getCost() << endl;
        }
        if(bSwap && (iIter % iNMethods == iSwapId)) {
            bImproved = swapOpt();
            // cout << "Cost after swapOpt: " << getCost() << endl;
        }
        if(bTailSwap && (iIter % iNMethods == iTailSwapId)) {
            // double dCostBefore = getCost();
            bImproved = tailSwapOpt();
            // if (bImproved)
            //	cout << "Cost before tailSwapOpt " << dCostBefore << ", after tailSwapOpt: " << getCost() << endl;
        }
        if(bImproved)
            iLastImproveIter = iIter;
        iIter++;
    }
}

// Method for updating the object (recalculates route costs and so on)
// The method terminates the program with an error message if an inconsistency or a violated
// solution is detected.
void CVRPSolution::consCalc() {
    int iRoute, i;
    double dCost = 0;
    int iN = m_pInstance->getN();
    vector<int> vecCustToRoute(iN + 1, -2);
    for(iRoute = 0; iRoute < (int)m_vecRoutes.size(); ++iRoute) {
        Route& route = m_vecRoutes[iRoute];
        route.consCalc();
        dCost += route.getCost();
        // Check the customers on the route.
        for(i = 1; i < route.size() - 1; ++i) {
            int iCustId = route[i];
            if(iCustId <= 0 || iCustId > iN)
                error("CVRPSolution::consCalc()", "Customer out of range!");
            if(vecCustToRoute[iCustId] != -2)
                error("CVRPSolution::consCalc()", "Customer assigned twice!");
            vecCustToRoute[iCustId] = iRoute;
        }
    }
    // Check that all unassigned customers are in m_vecUnassigned and that they
    // don't appear twice in the vector.
    for(i = 0; i < (int)m_vecUnassigned.size(); ++i) {
        if(vecCustToRoute[m_vecUnassigned[i]] >= 0)
            error("CVRPSolution::consCalc()", "Customer is both assigned and unassigned?!");
        if(vecCustToRoute[m_vecUnassigned[i]] == -1)
            error("CVRPSolution::consCalc()", "Customer appears twice in the vector of unassigned customers?!");
        // Mark that we have seen this customer in <m_vecUnassigned>
        vecCustToRoute[m_vecUnassigned[i]] = -1;
    }
    for(i = 1; i <= iN; ++i) {
        if(vecCustToRoute[i] == -2) {
            error("CVRPSolution::consCalc()", "Customer is neither assigned nor unassigned?!");
        }
    }

    dCost += (int)m_vecUnassigned.size() * m_dUnassignedCost;
    if(fabs(dCost - m_dCost) > m_dToleranceEps) {
        assert(false);
        error("CVRPSolution::consCalc()", "Calculated and stored cost mismatch!");
    }
    // Use the calculated cost. It is most precise (it's the result of the fewest number
    // of calculations).
    m_dCost = dCost;
    // We ought to check if <m_vecUnassigned> is consistent the information in the rest of the
    // solution object
#ifdef _DEBUG
    // Check if m_vecCust2Route is consistent with m_vecRoutes and m_vecUnassigned
    for(iRoute = 0; iRoute < (int)m_vecRoutes.size(); ++iRoute) {
        Route& route = m_vecRoutes[iRoute];
        for(i = 0; i < route.size(); ++i) {
            int iCustId = route[i];
            if(iCustId >= 1 && iCustId <= iN) {
                if(m_vecCust2Route[iCustId] != iRoute) {
                    error("CVRPSolution::consCalc()", "m_vecCust2Route is inconsistent with route information");
                }
            }
        }
    }
    for(i = 0; i < (int)m_vecUnassigned.size(); ++i) {
        if(m_vecCust2Route[m_vecUnassigned[i]] != -1) {
            error("CVRPSolution::consCalc()", "m_vecCust2Route is inconsistent with m_vecUnassigned");
        }
    }
#endif
}

void CVRPSolution::consCalcLight(const vector<bool>& vecRouteChanged) {
    assert(vecRouteChanged.size() == m_vecRoutes.size());
    int iRoute;
    double m_dCost = 0;
    int iN = m_pInstance->getN();
    for(iRoute = 0; iRoute < (int)m_vecRoutes.size(); ++iRoute) {
        Route& route = m_vecRoutes[iRoute];
        if(vecRouteChanged[iRoute]) {
            route.consCalc();
        }
        m_dCost += route.getCost();
    }
}

// Finds the customers assigned to a route and the route they are assigned to.
void CVRPSolution::getAssignedCustomers(vector<int>& vecCustIds, vector<int>& vecRouteIds) {
    int i, idx = 0, iN = m_pInstance->getN();
    vecCustIds.resize(iN - m_vecUnassigned.size());
    vecRouteIds.resize(iN - m_vecUnassigned.size());
    for(i = 1; i <= iN; i++) {
        if(m_vecCust2Route[i] >= 0) {
            vecCustIds[idx] = i;
            vecRouteIds[idx] = m_vecCust2Route[i];
            ++idx;
        }
    }
    assert(idx == (int)vecCustIds.size());
}

/***********************************************************************
 * void CVRPSolution::splitRoute(int iRoute, const vector<int> &vecSplitPos)
 *
 * Splits a route in the solution into two or more routes. The vector
 * <vecSplitPos> indicates where the route should be split. The result of applying the
 * method are routes
 * route[0],..,route[vecSplitPos[1]], n+1
 * 0, route[vecSplitPos[1]+1],..,route[vecSplitPos[2]], n+1
 * 0, route[vecSplitPos[2]+1],..,route[vecSplitPos[3]], n+1
 * ...
 * 0, route[vecSplitPos[k]+1],..,route[last element]]
 *
 * where <route> is the route pointed to by <iRoute>
 * The first of these routes is placed at the position of the old route
 * while the rest are added at the end of the vector of routes.
 *
 * The method is used in the arc removal methods of the LNS.
 *
 * --- input parameters ---
 * iRoute		: The route to split
 * vecSplitPos	: Indicates where the route should be split. The vector should
 *				  be sorted in an increasing way.
 *
 * --- return values ---
 ***********************************************************************/
void CVRPSolution::splitRoute(int iRoute, const vector<int>& vecSplitPos) {
    if(vecSplitPos.empty()) {
        cout << "Warning: CVRPSolution::splitRoute(...): No split points given" << endl;
        return;
    }
    Route& theRoute = m_vecRoutes[iRoute];
    // Reduce the cost of the solution by the cost of the original route. Later
    // on we add the cost of the new routes created.
    m_dCost -= theRoute.getCost();
    const vector<int>& vecNodes = theRoute.getNodes();

    int i, j, iN = m_pInstance->getN();

    vector<vector<int>> vecRoutes(vecSplitPos.size() + 1);
    int iStart = 0;
    for(i = 0; i < (int)vecSplitPos.size(); ++i) {
        vector<int>& vecRoute = vecRoutes[i];

        int iSplitPos = vecSplitPos[i] + 1;
        if(iStart >= iSplitPos || iSplitPos < 2 || iSplitPos > theRoute.size() - 2)
            error("CVRPSolution::splitRoute(...)", "invalid splitposition (bug in program?)");

        if(i != 0)
            vecRoute.push_back(0);
        vecRoute.insert(vecRoute.end(), vecNodes.begin() + iStart, vecNodes.begin() + iSplitPos);
        vecRoute.push_back(iN + 1);
        iStart = iSplitPos;
    }
    // The last route is formed by taking all the nodes after the last split point.
    vector<int>& vecRoute = vecRoutes.back();
    vecRoute.push_back(0);
    vecRoute.insert(vecRoute.end(), vecNodes.begin() + iStart, vecNodes.end());

    // The first route (vecRoutes[0]) we simply insert at the  position of the route we
    // are splitting:
    // Notice: We do not need to update m_vecCust2Route for the customers on this route. They have
    // not been moved around.
    theRoute.setNodes(vecRoutes[0]);
    m_dCost += theRoute.getCost();

    // The other routes are inserted at the back of the solution:
    for(i = 1; i < (int)vecRoutes.size(); ++i) {
        // Create the new route:
        createRoute();
        m_vecRoutes.back().setNodes(vecRoutes[i]);
        m_dCost += m_vecRoutes.back().getCost();
        const vector<int>& vecCurRoute = vecRoutes[i];

        // Update m_vecCust2Route array in this object (CVRPSolution) to reflect that customers
        // have moved.
        for(j = 1; j < (int)vecCurRoute.size() - 1; ++j)
            m_vecCust2Route[vecCurRoute[j]] = (int)m_vecRoutes.size() - 1;
    }
}

// dist computes the L_0 norm between two {0,1,2}-vectors that represents the solutions <this> and <other>
// At the moment the method only works for symmetric instances (dist(i,j) = dist(j,i))
// Think of the conversion from solition to {0,1,2}-vector in the following way:

// The solution is represented using x(i,j) variables where a variable x(i,j) = 1 if the
// solution includes an arc from i to j or from j to i. Since the instance is assumed to be
// symmetric we do not care about the direction of the arc. For the same reason we only include variables
// where i < j.
// A variable x(i,j) can have the value 2 if the edge (i,j) is used twice. This only happens
// for edges adjacent to the depot and corresponds to routes where only a single customer is
// visited.

double CVRPSolution::dist(const CVRPSolution& other) const {
    vector<SimpleArc> vecArcsRepThis, vecArcsRepOther;
    getArcRep(vecArcsRepThis);
    other.getArcRep(vecArcsRepOther);
    sort(vecArcsRepThis.begin(), vecArcsRepThis.end());
    sort(vecArcsRepOther.begin(), vecArcsRepOther.end());
    int idxThis = 0, idxOther = 0, iDist = 0;
    bool bDone = false;
    while(!bDone) {
        SimpleArc& arcThis = vecArcsRepThis[idxThis];
        SimpleArc& arcOther = vecArcsRepOther[idxOther];
        if(arcThis == arcOther) {
            ++idxThis;
            ++idxOther;
        } else {
            iDist++;
            if(arcThis < arcOther)
                ++idxThis;
            else
                ++idxOther;
        }
        if(idxThis >= (int)vecArcsRepThis.size()) {
            bDone = true;
            iDist += (int)vecArcsRepOther.size() - idxOther;
        }
        if(idxOther >= (int)vecArcsRepOther.size()) {
            bDone = true;
            iDist += (int)vecArcsRepThis.size() - idxThis;
        }
    }
    return iDist;
}

void CVRPSolution::getArcRep(vector<SimpleArc>& vecArcs) const {
    vecArcs.clear();
    if(m_pInstance->isSymmetric()) {
        int i, j, iN = m_pInstance->getN();
        for(i = 0; i < (int)m_vecRoutes.size(); ++i) {
            const Route& route = m_vecRoutes[i];
            if(route.size() > 2) {
                int iPrevVisitId = route.getNodeId(0);
                for(j = 1; j < route.size(); ++j) {
                    int iVisitId = route.getNodeId(j);
                    if(iVisitId == iN + 1)
                        iVisitId = 0;
                    assert(iVisitId != iPrevVisitId);
                    if(iPrevVisitId < iVisitId)
                        vecArcs.push_back(SimpleArc(iPrevVisitId, iVisitId));
                    else
                        vecArcs.push_back(SimpleArc(iVisitId, iPrevVisitId));
                    iPrevVisitId = iVisitId;
                }
            }
        }
    } else {
        error("CVRPSolution::getArcRep(...)", "Need to implement code for ACVRP instances");
    }
}

bool CVRPSolution::operator==(const CVRPSolution& other) const {
    vector<int> vecSeqThis, vecSeqOther;
    // Use hash sequence to get a "canonical" representation of the two solutions.
    getHashSequence(vecSeqThis);
    other.getHashSequence(vecSeqOther);
    return (vecSeqThis == vecSeqOther);
}

ostream& operator<<(ostream& os, const CVRPSolution& sol) {
    os << "cost: " << sol.getCost() << ", Routes: " << endl;
    outputVector(os, sol.m_vecRoutes, '\n');
    return os;
}

void CVRPSolution::writeSimpleSol(const string& strFileName) const {
    ofstream ofs(strFileName.c_str());

    ofs << getCost() << "\n";
    int i, j;
    for(i = 0; i < (int)m_vecRoutes.size(); ++i) {
        ofs << "{";
        const Route& route = m_vecRoutes[i];
        for(j = 1; j < route.size() - 1; ++j) {
            ofs << route[j];
            if(j < route.size() - 2)
                ofs << " ";
        }
        ofs << "}\n";
    }
    ofs.flush();
    ofs.close();
}

void CVRPSolution::polish() {
    bool bImproved = true;
    double dCostBefore = m_dCost;
    while(bImproved) {
        tailSwapOpt();
        bImproved = route2opt();
    }
    if(m_dCost < dCostBefore - m_dToleranceEps)
        cout << "Polishing improved solution from " << dCostBefore << " to " << m_dCost << endl;
}

void CVRPSolution::removeSequence(int iRouteId, int iStartPos, int iNToRemove, bool bConsCalc) {
    const vector<int>& vecRoute = m_vecRoutes[iRouteId].getNodes();
    vector<int> vecSeq(vecRoute.begin() + iStartPos, vecRoute.begin() + (iStartPos + iNToRemove));
    m_vecUnassignedSeqs.push_back(vecSeq);
    double dRouteCostBefore = m_vecRoutes[iRouteId].getCost();
    m_vecRoutes[iRouteId].removeSequence(iStartPos, iNToRemove);
    m_dCost += m_vecRoutes[iRouteId].getCost() - dRouteCostBefore;
    m_dCost += iNToRemove * m_dUnassignedCost;
    // update the vector of unassigned customers.
    m_vecUnassigned.insert(m_vecUnassigned.end(), vecSeq.begin(), vecSeq.end());
    // update <m_vecCust2Route>.
    int i;
    for(i = 0; i < (int)vecSeq.size(); ++i)
        m_vecCust2Route[vecSeq[i]] = -1;
    if(bConsCalc)
        m_vecRoutes[iRouteId].consCalc();
}
