#include "DestroyMethods.h"
#include "CVRPSolution.h"
#include "CVRPUtils.h"
#include "Utils.h"

#include <limits>

int calcNToRemove(TSRandom& randGen, int iN, double dMinPercentage = 0.1, double dMaxPercentage = 0.4, int iAbsMin = 10, int iAbsMax = 50) {
    int iMin = maxFunc(iAbsMin, (int)(dMinPercentage * iN));
    iMin = minFunc(iMin, iN);
    int iMax = minFunc(iAbsMax, (int)(dMaxPercentage * iN));
    iMax = minFunc(iMax, iN);

    // In case we have large instances such that n*dMinPercentage > iAbsMax. In that case we would end up with iMin > iMax without the following test.
    if(iMin > iMax)
        iMin = iMax;

    // cout << "calcNToRemove, min: " << iMin << ", max: " << iMax << endl;
    return randGen.getRandom(iMin, iMax);
}

void removeCustomers(CVRPSolution& sol, const vector<bool>& vecIsCustRemoved, const vector<bool>& vecIsRouteChanged) {
#ifndef BATCH_MODE
    sol.consCalc();
#endif

    const CVRPInstance& instance = sol.getInstance();
    int iN = instance.getN();

    int i;
    for(i = 0; i < (int)sol.getNRoutes(); ++i) {
        if(vecIsRouteChanged[i]) {
            sol.removeCustomers(i, vecIsCustRemoved);
            if(sol.getRoute(i).getDuration() > instance.getMaxDuration()) {
                // Removing customers has caused the route duration to increase (distance
                // doesn't satisfy triangle inequality). This should be rare. Remove all customers
                // from the route.
                cout << "ALNS::removeCustomers(...) : duration violated" << endl;
                sol.clearRoute(i);
            }
        }
    }
    // Check that the solution is ok.
    sol.consCalc();
}

void RandomDestroy::destroySolution(CVRPSolution& sol, TSRandom& randGen) {
    int iNToRemove = calcNToRemove(randGen, sol.getInstance().getN());
    // cout << "iNToRemove: " << iNToRemove << endl;
    // Get the customers that are assigned to a route.
    vector<int> vecAssignedCustIds, vecRouteIds;
    sol.getAssignedCustomers(vecAssignedCustIds, vecRouteIds);

    int iNRemoved = 0;
    vector<bool> vecIsCustRemoved(sol.getInstance().getN() + 1, false);
    vector<bool> vecIsRouteChanged(sol.getNRoutes(), false);

    // Select the customers to remove:
    while(iNRemoved < iNToRemove && !vecAssignedCustIds.empty()) {
        int idx = randGen.getRandom(0, (int)vecAssignedCustIds.size() - 1);
        vecIsCustRemoved[vecAssignedCustIds[idx]] = true;
        vecIsRouteChanged[vecRouteIds[idx]] = true;
        unorderedVectorRemoveAt(vecAssignedCustIds, idx);
        unorderedVectorRemoveAt(vecRouteIds, idx);
        ++iNRemoved;
    }
    removeCustomers(sol, vecIsCustRemoved, vecIsRouteChanged);
}

struct RemoveDelta {
    int m_iCust;
    double m_dDelta;
    int m_iRoute;
    int m_iPos;
    RemoveDelta(int iCust, double dDelta, int iRoute, int iPos) : m_iCust(iCust), m_dDelta(dDelta), m_iRoute(iRoute), m_iPos(iPos) {}
    bool operator<(const RemoveDelta& other) const {
        // We want the largest deltas to appear at the front of a list sorted using this operator.
        return m_dDelta > other.m_dDelta;
    }
};

void ExpensiveNodeDestroy::destroySolution(CVRPSolution& sol, TSRandom& randGen) {
    int iNToRemove = calcNToRemove(randGen, sol.getInstance().getN());

    const CVRPInstance& instance = sol.getInstance();
    const Matrix<double>& matDist = instance.getDists();
    int iRoute, i, iN = instance.getN();
    vector<RemoveDelta> vecRemoveDelta;
    vecRemoveDelta.reserve(iN);

    for(iRoute = 0; iRoute < sol.getNRoutes(); ++iRoute) {
        const Route& route = sol.getRoute(iRoute);
        for(i = 1; i < route.size() - 1; ++i) {
            int iCust = route[i];
            double dDelta = matDist.getElement(route[i - 1], iCust) + matDist.getElement(iCust, route[i + 1]) - matDist.getElement(route[i - 1], route[i + 1]);
            vecRemoveDelta.push_back(RemoveDelta(iCust, dDelta, iRoute, i));
        }
    }

    int iNRemoved = 0;
    while(!vecRemoveDelta.empty() && iNRemoved < iNToRemove) {
        int idx = (int)(pow(randGen.getRandomDouble(0, 1), m_dRandFactor) * (double)vecRemoveDelta.size());
        idx = min(idx, (int)vecRemoveDelta.size() - 1); // just a safeguard
        // Find the idx'th element in the sorted vector (without sorting)
        nth_element(vecRemoveDelta.begin(), vecRemoveDelta.begin() + idx, vecRemoveDelta.end());
        RemoveDelta removeDelta = vecRemoveDelta[idx];
        unorderedVectorRemoveAt(vecRemoveDelta, idx);
        sol.removeCustomer(removeDelta.m_iRoute, removeDelta.m_iPos);
        iNRemoved++;
        if(sol.getRoute(removeDelta.m_iRoute).getDuration() > instance.getMaxDuration()) {
            // Removing customers has caused the route duration to increase (distance
            // doesn't satisfy triangle inequality). This should be rare. Remove all customers
            // from the route.
            cout << "ALNS::removeCustomers(...) : duration violated" << endl;
            iNRemoved += sol.getRoute(removeDelta.m_iRoute).size() - 2;
            sol.clearRoute(removeDelta.m_iRoute);
            // Update vecRemoveDelta;
            while(i < (int)vecRemoveDelta.size()) {
                if(vecRemoveDelta[i].m_iRoute == removeDelta.m_iRoute)
                    unorderedVectorRemoveAt(vecRemoveDelta, i);
                else
                    ++i;
            }
        }
        updateVecRemoveDelta(sol, vecRemoveDelta, removeDelta);
    }
    sol.consCalc();
}

void ExpensiveNodeDestroy::updateVecRemoveDelta(const CVRPSolution& sol, vector<RemoveDelta>& vecRemoveDelta, RemoveDelta& removeDelta) const {
    const Matrix<double>& matDist = sol.getInstance().getDists();
    int iRoute = removeDelta.m_iRoute;
    int iRemovePos = removeDelta.m_iPos;
    int i;
    const Route& route = sol.getRoute(iRoute);
    for(i = 0; i < (int)vecRemoveDelta.size(); ++i) {
        if(vecRemoveDelta[i].m_iRoute == iRoute) {
            int iOtherPos = vecRemoveDelta[i].m_iPos;
            // Has the position of this customer changed because of the removal?
            if(iOtherPos > iRemovePos) {
                iOtherPos--;
                vecRemoveDelta[i].m_iPos = iOtherPos;
            }
            // Is this RemoveDelta for one of the nodes that were adjacent with the removed node? If so, we have to recalculate delta values.
            if(iOtherPos == iRemovePos || iOtherPos == iRemovePos - 1) {
                int iCust = route[iOtherPos];
                double dDelta = matDist.getElement(route[iOtherPos - 1], iCust) + matDist.getElement(iCust, route[iOtherPos + 1]) -
                                matDist.getElement(route[iOtherPos - 1], route[iOtherPos + 1]);
                vecRemoveDelta[i].m_dDelta = dDelta;
            }
        }
    }
}

void GeoDestroy::destroySolution(CVRPSolution& sol, TSRandom& randGen) {
#ifndef BATCH_MODE
    sol.consCalc();
#endif

    const CVRPInstance& instance = sol.getInstance();
    int iNToRemove = calcNToRemove(randGen, instance.getN());

    vector<int> vecAssignedCustIds, vecRouteIds;
    sol.getAssignedCustomers(vecAssignedCustIds, vecRouteIds);
    const Matrix<double>& matDist = instance.getDists();

    int iNRemoved = 0;
    vector<bool> vecIsCustRemoved(instance.getN() + 1, false);
    vector<bool> vecIsRouteChanged(sol.getNRoutes(), false);
    vector<int> vecRemoved;
    if(vecAssignedCustIds.empty())
        return;

    // Chose a random customer as the "seed"
    int i, iIdx = randGen.getRandom(0, (int)vecAssignedCustIds.size() - 1);
    int iCustId = vecAssignedCustIds[iIdx];
    vecRemoved.push_back(iCustId);
    vecIsCustRemoved[iCustId] = true;
    vecIsRouteChanged[vecRouteIds[iIdx]] = true;
    unorderedVectorRemoveAt(vecAssignedCustIds, iIdx);
    unorderedVectorRemoveAt(vecRouteIds, iIdx);

    ++iNRemoved;
    // Select the customers to remove:
    while(iNRemoved < iNToRemove && !vecAssignedCustIds.empty()) {
        vector<NodeData> vecNodeDist;
        // Select a random node from the set of removed customers:
        iCustId = vecRemoved[randGen.getRandom(0, (int)vecRemoved.size() - 1)];
        // Get the distance from this node to all other assigned nodes
        for(i = 0; i < (int)vecAssignedCustIds.size(); ++i) {
            vecNodeDist.push_back(NodeData(i, matDist.getElement(iCustId, vecAssignedCustIds[i])));
        }

        int iOffset = (int)(pow(randGen.getRandomDouble(0, 1), m_dRandFactor) * (vecNodeDist.size() - 1));
        if(iOffset > (int)vecNodeDist.size() - 1 || iOffset < 0)
            cout << "Offset out of range" << endl;
        nth_element(vecNodeDist.begin(), vecNodeDist.begin() + iOffset, vecNodeDist.end(), sortIncreasing);
        int iChosenIdx = vecNodeDist[iOffset].m_iNode;
        vecRemoved.push_back(vecAssignedCustIds[iChosenIdx]);

        vecIsCustRemoved[vecAssignedCustIds[iChosenIdx]] = true;
        vecIsRouteChanged[vecRouteIds[iChosenIdx]] = true;
        unorderedVectorRemoveAt(vecAssignedCustIds, iChosenIdx);
        unorderedVectorRemoveAt(vecRouteIds, iChosenIdx);
        ++iNRemoved;
    }
    // Physically remove the customers.
    removeCustomers(sol, vecIsCustRemoved, vecIsRouteChanged);
}

void SolHistoryDestroy::destroySolution(CVRPSolution& sol, TSRandom& randGen) {
#ifndef BATCH_MODE
    sol.consCalc();
#endif

    int iNToRemove = calcNToRemove(randGen, sol.getInstance().getN());

    const CVRPInstance& instance = sol.getInstance();
    int iRoute, i, iN = instance.getN();
    vector<RemoveDelta> vecRemoveDelta;
    vecRemoveDelta.reserve(iN);

    for(iRoute = 0; iRoute < sol.getNRoutes(); ++iRoute) {
        const Route& route = sol.getRoute(iRoute);
        for(i = 1; i < route.size() - 1; ++i) {
            int iCust = route[i];
            double dCost = (m_pArcCostTracker->getSolCost(route[i - 1], iCust) + m_pArcCostTracker->getSolCost(iCust, route[i + 1])) / 2;
            vecRemoveDelta.push_back(RemoveDelta(iCust, dCost, iRoute, i));
        }
    }

    int iNRemoved = 0;
    while(!vecRemoveDelta.empty() && iNRemoved < iNToRemove) {
        int idx = (int)(pow(randGen.getRandomDouble(0, 1), m_dRandFactor) * (double)vecRemoveDelta.size());
        idx = min(idx, (int)vecRemoveDelta.size() - 1); // just a safeguard
        // Find the idx'th element in the sorted vector (without sorting)
        nth_element(vecRemoveDelta.begin(), vecRemoveDelta.begin() + idx, vecRemoveDelta.end());
        RemoveDelta removeDelta = vecRemoveDelta[idx];
        unorderedVectorRemoveAt(vecRemoveDelta, idx);
        sol.removeCustomer(removeDelta.m_iRoute, removeDelta.m_iPos);
        iNRemoved++;
        if(sol.getRoute(removeDelta.m_iRoute).getDuration() > instance.getMaxDuration()) {
            // Removing customers has caused the route duration to increase (distance
            // doesn't satisfy triangle inequality). This should be rare. Remove all customers
            // from the route.
            cout << "ALNS::removeCustomers(...) : duration violated" << endl;
            iNRemoved += sol.getRoute(removeDelta.m_iRoute).size() - 2;
            sol.clearRoute(removeDelta.m_iRoute);
            // Update vecRemoveDelta;
            while(i < (int)vecRemoveDelta.size()) {
                if(vecRemoveDelta[i].m_iRoute == removeDelta.m_iRoute)
                    unorderedVectorRemoveAt(vecRemoveDelta, i);
                else
                    ++i;
            }
        }
        updateVecRemoveDelta(sol, vecRemoveDelta, removeDelta);
    }
    sol.consCalc();
}

void SolHistoryDestroy::updateVecRemoveDelta(const CVRPSolution& sol, vector<RemoveDelta>& vecRemoveDelta, RemoveDelta& removeDelta) const {
    int iRoute = removeDelta.m_iRoute;
    int iRemovePos = removeDelta.m_iPos;
    int i;
    const Route& route = sol.getRoute(iRoute);
    for(i = 0; i < (int)vecRemoveDelta.size(); ++i) {
        if(vecRemoveDelta[i].m_iRoute == iRoute) {
            int iOtherPos = vecRemoveDelta[i].m_iPos;
            // Has the position of this customer changed because of the removal?
            if(iOtherPos > iRemovePos) {
                iOtherPos--;
                vecRemoveDelta[i].m_iPos = iOtherPos;
            }
            // @@@ As an experiment, try not to update the delta value. Doing so had a positive impact
            // in the PALNS-VRPTW code (see that).

            // Is this RemoveDelta for one of the nodes that were adjacent with the removed node? If so, we have to recalculate delta values.

            if(iOtherPos == iRemovePos || iOtherPos == iRemovePos - 1) {
                int iCust = route[iOtherPos];
                double dDelta = m_pArcCostTracker->getSolCost(route[iOtherPos - 1], iCust) + m_pArcCostTracker->getSolCost(iCust, route[iOtherPos + 1]);
                vecRemoveDelta[i].m_dDelta = dDelta;
            }
        }
    }
}

LocalRandomDestroy::LocalRandomDestroy(const LocalRandomDestroy& other)
    : m_iNeighborhoodSize(other.m_iNeighborhoodSize), m_vecNodeCount(other.m_vecNodeCount) {}

LocalRandomDestroy::LocalRandomDestroy(int iNeighborhoodSize, int iN) : m_iNeighborhoodSize(iNeighborhoodSize), m_vecNodeCount(iN + 1, 0) {}

void LocalRandomDestroy::destroySolution(CVRPSolution& sol, TSRandom& randGen) {
    // Find the customer that should be the "center of the neighborhood":
    const CVRPInstance& instance = sol.getInstance();
    int iN = instance.getN();
    int i, iNMinCount = std::numeric_limits<int>::max();
    vector<int> vecMinimizers;
    for(i = 1; i < (int)m_vecNodeCount.size(); ++i) {
        if(m_vecNodeCount[i] < iNMinCount) {
            vecMinimizers.clear();
            iNMinCount = m_vecNodeCount[i];
        }
        if(m_vecNodeCount[i] == iNMinCount)
            vecMinimizers.push_back(i);
    }
    int iCenterNode = vecMinimizers[randGen.getRandom(0, (int)vecMinimizers.size() - 1)];
    const vector<int>& vecNearestCust = instance.getNearestNodes(iCenterNode);
    // Get the <m_iNeighborhoodSize> customers nearest to <iCenterNode> (or all if m_iNeighborhoodSize > iN).
    vector<int> vecAvailNodes(vecNearestCust.begin(), vecNearestCust.begin() + min(m_iNeighborhoodSize, (int)vecNearestCust.size()));
    vector<bool> vecIsAvail;
    makeIsInSetVec(vecAvailNodes, vecIsAvail, iN + 1);

    int iNToRemove = calcNToRemove(randGen, iN);
    // cout << "iNToRemove: " << iNToRemove << endl;
    // Get the customers that are assigned to a route.
    vector<int> vecAssignedCustIds, vecRouteIds;
    sol.getAssignedCustomers(vecAssignedCustIds, vecRouteIds);
    // Filter out all the customers that are not in the neighborhood of the center customer:
    i = 0;
    while(i < (int)vecAssignedCustIds.size()) {
        if(vecIsAvail[vecAssignedCustIds[i]])
            ++i;
        else {
            unorderedVectorRemoveAt(vecAssignedCustIds, i);
            unorderedVectorRemoveAt(vecRouteIds, i);
        }
    }

    int iNRemoved = 0;
    vector<bool> vecIsCustRemoved(iN + 1, false);
    vector<bool> vecIsRouteChanged(sol.getNRoutes(), false);

    // Select the customers to remove:
    while(iNRemoved < iNToRemove && !vecAssignedCustIds.empty()) {
        int idx = randGen.getRandom(0, (int)vecAssignedCustIds.size() - 1);
        vecIsCustRemoved[vecAssignedCustIds[idx]] = true;
        vecIsRouteChanged[vecRouteIds[idx]] = true;
        m_vecNodeCount[vecAssignedCustIds[idx]]++;
        unorderedVectorRemoveAt(vecAssignedCustIds, idx);
        unorderedVectorRemoveAt(vecRouteIds, idx);
        ++iNRemoved;
    }
    removeCustomers(sol, vecIsCustRemoved, vecIsRouteChanged);
}

struct SequenceInfo {
    SequenceInfo(int iRouteId, int iStartPos, int iLength) : m_iRouteId(iRouteId), m_iStartPos(iStartPos), m_iLength(iLength) {}

    int m_iRouteId;
    int m_iStartPos;
    int m_iLength;
    // Operator that ensures that sorting will be in increasing order of first route id and then in decreasing order of start position.
    // This makes it easy to remove the sequences later on.
    bool operator<(const SequenceInfo& other) const {
        return m_iRouteId < other.m_iRouteId || (m_iRouteId == other.m_iRouteId && m_iStartPos >= other.m_iStartPos);
    }
};

RandomSequenceDestroy::RandomSequenceDestroy() {
    // dMinPercentage = 0.1, double dMaxPercentage = 0.4, int iAbsMin = 10, int iAbsMax = 50
    m_iMinSeqSize = 2;
    m_iMaxSeqSize = 2;
    m_iMinAbsNSeq = 5;
    m_iMaxAbsNSeq = 50;
    m_dMinRelNSeq = 0.001;
    m_dMaxRelNSeq = 0.3;
}

void RandomSequenceDestroy::destroySolution(CVRPSolution& sol, TSRandom& randGen) {
    int iN = sol.getInstance().getN();
    int iNSeqsToRemove = calcNToRemove(randGen, iN, m_dMinRelNSeq, m_dMaxRelNSeq, m_iMinAbsNSeq, m_iMaxAbsNSeq);

    // Find max route length:
    int i, j, iMaxRouteLength = -1;
    vector<vector<int>> vecSequences;
    vector<int> vecAssignedCusts;
    vector<int> vecCustToSeq(iN + 1, -1);
    // vecSeqToRouteId[i] is the route that sequence i originated from
    vector<int> vecSeqToRouteId;
    // vecSeqStartRoutePos[i] is the position in the original route that the first element in sequence i occupied.
    // As we remove sequences we start to split the original sequences and this number "becomes important".
    vector<int> vecSeqStartRoutePos;

    for(i = 0; i < sol.getNRoutes(); ++i) {
        int iSize = sol.getRoute(i).size();
        if(iSize > iMaxRouteLength)
            iMaxRouteLength = iSize;
        if(iSize > 2) {
            const vector<int>& vecVisits = sol.getVisits(i);
            vecSequences.push_back(vector<int>(vecVisits.begin() + 1, vecVisits.begin() + (vecVisits.size() - 1)));
            vecSeqToRouteId.push_back(i);
            vecSeqStartRoutePos.push_back(1);
            int iSeqId = (int)vecSequences.size() - 1;
            for(j = 1; j < (int)vecVisits.size() - 1; ++j) {
                vecCustToSeq[vecVisits[j]] = iSeqId;
                vecAssignedCusts.push_back(vecVisits[j]);
            }
        }
    }
    if(((int)vecAssignedCusts.size()) != (iN - sol.getUnassignedCusts().size()))
        error("RandomSequenceDestroy::destroySolution(...)", "#unassigned customers mismatch");

    // Subtract 2 because of start and end depot.
    iMaxRouteLength -= 2;

    int iNRemoved = 0;
    vector<bool> vecRouteChanged(sol.getNRoutes(), false);
    vector<SequenceInfo> vecSequencesToRemove;
    while(iNRemoved < iNSeqsToRemove && !vecAssignedCusts.empty()) {
        // Select an assigned customer (compared to selecting a random route, this should make it more probable that we take sequences from routes with many
        // customers)
        int iIdx = randGen.getRandom(0, (int)vecAssignedCusts.size() - 1);
        // Note, we might not even remove <iCustId>, we just use it to select the sequence to remove from.
        int iCustId = vecAssignedCusts[iIdx];
        // Find the sequence that this customer belongs to
        int iSeqId = vecCustToSeq[iCustId];
        // It is on purpose that we take a copy here and do not just use a reference. Later on we add to
        // vecSequences and that could make the reference invalid if STL needs to reallocate the array. The code
        // could be rearranged so it would not cause problems, but it seems to be a source for future errors, so
        // it seems worthwhile to waste a few clock cycles here.
        vector<int> vecSeq(vecSequences[iSeqId]);
        int iSize = (int)vecSequences[iSeqId].size();
        int iMaxToRemove = min(iSize, m_iMaxSeqSize);
        int iMinToRemove = min(iSize, m_iMinSeqSize);
        int iNToRemove = randGen.getRandom(iMinToRemove, iMaxToRemove);
        int iStartIdx = randGen.getRandom(0, iSize - iNToRemove);
        vecSequencesToRemove.push_back(SequenceInfo(vecSeqToRouteId[iSeqId], iStartIdx + vecSeqStartRoutePos[iSeqId], iNToRemove));

        // Update <vecCustToSeq> and <vecAssignedCusts>
        vector<bool> vecIsRemoved(iN + 1, false);
        for(i = iStartIdx; i < iStartIdx + iNToRemove; ++i) {
            vecCustToSeq[vecSeq[i]] = -1;
            vecIsRemoved[vecSeq[i]] = true;
        }
        i = 0;
        while(i < vecAssignedCusts.size()) {
            if(vecIsRemoved[vecAssignedCusts[i]])
                unorderedVectorRemoveAt(vecAssignedCusts, i);
            else
                ++i;
        }

        // Update vecSequences. Did we cut the sequence in two pieces?
        if(iStartIdx > 0 && (iStartIdx + iNToRemove < vecSeq.size())) {
            // Yes, Copy the last piece to a new sequence:
            vecSequences.push_back(vector<int>(vecSeq.begin() + (iStartIdx + iNToRemove), vecSeq.end()));
            const vector<int>& vecNewSeq = vecSequences.back();
            int iNewSeqId = (int)vecSequences.size() - 1;
            for(i = 0; i < (int)vecNewSeq.size(); ++i)
                vecCustToSeq[vecNewSeq[i]] = iNewSeqId;
            vecSeqToRouteId.push_back(vecSeqToRouteId[iSeqId]);
            vecSeqStartRoutePos.push_back(vecSeqStartRoutePos[iSeqId] + iStartIdx + iNToRemove);
            // erase the removed sequence as well as the last piece from the original sequence:
            vecSequences[iSeqId].erase(vecSequences[iSeqId].begin() + iStartIdx, vecSequences[iSeqId].end());
        } else {
            // No, either there is only one piece left or the whole sequence was removed:
            if(iStartIdx > 0)
                // Remove end:
                vecSequences[iSeqId].erase(vecSequences[iSeqId].begin() + iStartIdx, vecSequences[iSeqId].end());
            else {
                // Remove start?
                if(iStartIdx == 0 && iNToRemove < vecSeq.size()) {
                    vecSequences[iSeqId].erase(vecSequences[iSeqId].begin(), vecSequences[iSeqId].begin() + iStartIdx + iNToRemove);
                    vecSeqStartRoutePos[iSeqId] += iNToRemove;
                } else {
                    // Whole sequence was removed
                    assert(iStartIdx == 0 && iNToRemove == vecSeq.size());
                    vecSequences[iSeqId].clear();
                    vecSeqStartRoutePos[iSeqId] = -1;
                }
            }
        }

        ++iNRemoved;
    }

    // Now we actually remove the sequences:
    // After sorting <vecSequencesToRemove> will be in increasing order of first route id and then in decreasing order of start position.
    // This makes it easy to remove the sequences: If we process vecSequencesToRemove in order then we will do not have to
    // worry about positions changing because of removals.
    sort(vecSequencesToRemove.begin(), vecSequencesToRemove.end());
    for(i = 0; i < (int)vecSequencesToRemove.size(); ++i) {
        const SequenceInfo& SeqInfo = vecSequencesToRemove[i];
        vecRouteChanged[SeqInfo.m_iRouteId] = true;
        sol.removeSequence(SeqInfo.m_iRouteId, SeqInfo.m_iStartPos, SeqInfo.m_iLength, false);
    }

    sol.consCalcLight(vecRouteChanged);
#ifdef _DEBUG
    sol.consCalc();
#endif
}
