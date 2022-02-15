#include "PALNS-CVRP.h"
#include "CVRPSolution.h"
#include "DestroyMethods.h"
#include "GreedySeqAlg.h"
#include "PALNS.h"
#include "RepairMethods.h"

class CVRPInitialSolutionCreator : public InitialSolutionCreator {
public:
    CVRPInitialSolutionCreator(const CVRPInstance& instance) : m_greedySeqAlg(instance, true) {}
    virtual void createInitialSolution(CVRPSolution& sol, TSRandom& randGen) {
        sol.clear();
        m_greedySeqAlg.go(sol, randGen);
    }

protected:
    GreedySeqAlg m_greedySeqAlg;
};

void PALNS_CVRP::go(CVRPSolution& sol, int iNThreads, int iNRetries, bool bHeader, const ParFileParser& pfp) const {
    PALNS palns(sol.getInstance());

    ArcCostTracker arcCostTracker(sol.getInstance());
    palns.setSolTrackerCallBack(arcCostTracker);

    vector<int> vecNodeDestroy, vecSeqDestroy;
    RandomDestroy randomDestroy;
    int iRandomDestroy = palns.addDestroyMethod(randomDestroy, "Random destroy");
    vecNodeDestroy.push_back(iRandomDestroy);

    // ExpensiveNodeDestroy expensiveNodeDestroy(4.0);
    // palns.addDestroyMethod(expensiveNodeDestroy, "Expensive arc destroy");

    GeoDestroy geoDestroy(5.0);
    int iGeoDestroy = palns.addDestroyMethod(geoDestroy, "Geo destroy");
    vecNodeDestroy.push_back(iGeoDestroy);

    SolHistoryDestroy solHistoryDestroy(5.0);
    int iHistoryDestroy = palns.addDestroyMethod(solHistoryDestroy, "Solution history destroy");
    vecNodeDestroy.push_back(iHistoryDestroy);

    // LocalRandomDestroy localRandomDestroy(150, sol.getInstance().getN());
    // palns.addDestroyMethod(localRandomDestroy, "Local random destroy");

    int iSeqDestroyRepair = 0, iRandSeqDestroyId = -1;
    pfp.getIntParameter("SeqDestroyRepair", iSeqDestroyRepair);

    RandomSequenceDestroy randSeqDestroy;
    if(iSeqDestroyRepair > 0) {
        int iMinSeqSize = 2, iMaxSeqSize = 2, iMinAbsNSeq = 5, iMaxAbsNSeq = 50;
        double dMinRelNSeq = 0.001, dMaxRelNSeq = 0.4;
        pfp.getIntParameter("MinSeqSize", iMinSeqSize);
        pfp.getIntParameter("MaxSeqSize", iMaxSeqSize);
        pfp.getIntParameter("MinAbsNSeq", iMinAbsNSeq);
        pfp.getIntParameter("MaxAbsNSeq", iMaxAbsNSeq);
        pfp.getDoubleParameter("MinRelNSeq", dMinRelNSeq);
        pfp.getDoubleParameter("MaxRelNSeq", dMaxRelNSeq);

        randSeqDestroy.setParamAbsNSeq(iMinAbsNSeq, iMaxAbsNSeq);
        randSeqDestroy.setParamRelNSeq(dMinRelNSeq, dMaxRelNSeq);
        randSeqDestroy.setParamSeqSize(iMinSeqSize, iMaxSeqSize);

        iRandSeqDestroyId = palns.addDestroyMethod(randSeqDestroy, "Random sequence destroy");
        vecSeqDestroy.push_back(iRandSeqDestroyId);
    }

    // ------------------------------------------------------------------------------------------------
    double dInterRouteImprProb;
    pfp.getDoubleParameter("RepairInterRouteImprProb", dInterRouteImprProb);
    Regret2Repair regret2RepairNormal(1.0, dInterRouteImprProb);
    Regret2Repair regret2RepairRandomized(1.5, dInterRouteImprProb);
    palns.addRepairMethod(regret2RepairNormal, "regret-2 normal", false, vecNodeDestroy);
    palns.addRepairMethod(regret2RepairRandomized, "regret-2 randomized", false, vecNodeDestroy);
    // HistoryRegret2Repair historyRegret2Repair(1.5, dInterRouteImprProb);
    // palns.addRepairMethod(historyRegret2Repair, "History regret-2");

    SequenceRegretRepair seqRegret2RepairNormal(1.0, dInterRouteImprProb);
    SequenceRegretRepair seqRet2RepairRandomized(1.5, dInterRouteImprProb);
    if(iSeqDestroyRepair > 0) {
        palns.addRepairMethod(seqRegret2RepairNormal, "sequence regret-2 normal", false, vecSeqDestroy);
        palns.addRepairMethod(seqRet2RepairRandomized, "sequence regret-2 randomized", false, vecSeqDestroy);
    }

    CVRPInitialSolutionCreator cvrpSC(sol.getInstance());
    PALNSParameters params(pfp);

    string strName = sol.getInstance().getName();
    palns.go(sol, params, &cvrpSC);
}