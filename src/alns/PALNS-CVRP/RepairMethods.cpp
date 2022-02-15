#include "RepairMethods.h"
#include "ArcCostTracker.h"
#include "Regret2Alg.h"
#include "SeqRegret2Alg.h"

Regret2Repair::Regret2Repair(double dRandFactor, double dInterRouteImprProb) : m_dRandFactor(dRandFactor), m_dInterRouteImprProb(dInterRouteImprProb) {}

void Regret2Repair::repairSolution(CVRPSolution& sol, TSRandom& randGen) {
#ifndef BATCH_MODE
    sol.consCalc();
#endif

    Regret2Alg regret2Alg(sol.getInstance(), randGen);
    regret2Alg.setRandFactor(m_dRandFactor);
    regret2Alg.go(sol);

#ifndef BATCH_MODE
    sol.consCalc();
#endif

    sol.route2opt();

#ifndef BATCH_MODE
    sol.consCalc();
#endif

    if((m_dInterRouteImprProb > 0) && randGen.getRandomDouble(0, 1) <= m_dInterRouteImprProb) {
        bool bImproved = true;
        while(bImproved) {
            bImproved = sol.tailSwapOpt();

#ifndef BATCH_MODE
            sol.consCalc();
#endif

            if(bImproved) {
                bImproved = sol.route2opt();

#ifndef BATCH_MODE
                sol.consCalc();
#endif
            }
        }
    }

#ifndef BATCH_MODE
    sol.consCalc();
#endif
    // sol.combinedOpt();
}

void HistoryRegret2Repair::repairSolution(CVRPSolution& sol, TSRandom& randGen) {
#ifndef BATCH_MODE
    sol.consCalc();
#endif

    Regret2Alg regret2Alg(sol.getInstance(), randGen);
    regret2Alg.setRandFactor(m_dRandFactor);
    regret2Alg.go(sol, &(m_pArcCostTracker->getCostMatrix()));
    sol.route2opt();
    if((m_dInterRouteImprProb > 0) && randGen.getRandomDouble(0, 1) <= m_dInterRouteImprProb) {
        bool bImproved = true;
        while(bImproved) {
            bImproved = sol.tailSwapOpt();
            if(bImproved)
                bImproved = sol.route2opt();
        }
    }

#ifndef BATCH_MODE
    sol.consCalc();
#endif
}

SequenceRegretRepair::SequenceRegretRepair(double dRandFactor, double dInterRouteImprProb)
    : m_dRandFactor(dRandFactor), m_dInterRouteImprProb(dInterRouteImprProb) {}

void SequenceRegretRepair::repairSolution(CVRPSolution& sol, TSRandom& randGen) {
#ifndef BATCH_MODE
    sol.consCalc();
#endif

    SeqRegret2Alg seqRegret2Alg(sol.getInstance(), randGen);
    seqRegret2Alg.setRandFactor(m_dRandFactor);
    seqRegret2Alg.go(sol);
    sol.route2opt();
    // cout << "Cost after sequence repair: " << sol.getCost() << endl;
    if((m_dInterRouteImprProb > 0) && randGen.getRandomDouble(0, 1) <= m_dInterRouteImprProb) {
        bool bImproved = true;
        while(bImproved) {
            bImproved = sol.tailSwapOpt();
            if(bImproved)
                bImproved = sol.route2opt();
        }
    }

#ifndef BATCH_MODE
    sol.consCalc();
#endif
}