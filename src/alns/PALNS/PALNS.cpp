//
// Created by Alberto Santini on 04/06/2021.
//

#include "PALNS.h"
#include "AcceptanceCriteria/RecordDeviation.h"

void PALNS::resetLocalParameters(CVRPSolution& sol, PALNSParameters& par) {
    m_iIter = 0;
    m_initialSol = sol;
    m_bestSol = sol;
    m_curSol = sol;
    m_dEpsilon = 1e-7;
    m_par = par;

    m_vecDestroyWeights.clear();
    m_vecDestroyWeights.resize(m_vecDestroyMethods.size(), 1);
    m_vecRepairWeights.clear();
    m_vecRepairWeights.resize(m_vecRepairMethods.size(), 1);

    m_dStartTime = wallClock();

    m_acceptanceCriterion = std::unique_ptr<AcceptanceCriterion>(new RecordDeviation(m_par));

    if(m_par.m_bUseDecomposition) {
        m_decomposition = get_decomposition_method(*this);
        m_iNextDecompositionIter = m_par.m_iDecompositionIters;
    } else {
        m_decomposition = nullptr;
        m_iNextDecompositionIter = 0;
    }

    vfreq = std::vector<std::vector<double>>(m_instance.getN(), std::vector<double>(m_instance.getN(), 0.0));
    afreq = std::vector<std::vector<double>>(m_instance.getN(), std::vector<double>(m_instance.getN(), 0.0));
    pfreq = PathFrequencyMatrix{};
}

void PALNS::go(CVRPSolution& sol, PALNSParameters par, InitialSolutionCreator* pISC) {
    TSRandom rng;
    pISC->createInitialSolution(sol, rng);

    resetLocalParameters(sol, par);

    m_acceptanceCriterion->initialise();

    start();

    m_bestSol.postProcess();

    sol = m_bestSol;
}

int PALNS::addDestroyMethod(AbstractDestroyMethod& destroyMethod, const string& strDescr) {
    m_vecDestroyMethods.push_back(&destroyMethod);
    m_vecDestroyDescr.push_back(strDescr);
    return (int)m_vecDestroyMethods.size() - 1;
}

void PALNS::addRepairMethod(AbstractRepairMethod& repairMethod, const string& strDescr, bool bAllCompatible, const vector<int>& vecCompatibleDestroyIds) {
    m_vecRepairMethods.push_back(&repairMethod);
    m_vecRepairDescr.push_back(strDescr);
    m_vecRepairDestroyAllCompatible.push_back(bAllCompatible);
    m_vecRepairCompatibleDestroyIds.push_back(vecCompatibleDestroyIds);
}

void PALNS::start() {
    TSRandom randGen{};
    randGen.robustInitialise(m_iMasterSeed);

    for(auto& dm : m_vecDestroyMethods) {
        dm->setSolTracker(m_pSolTrackerCallBack);
    }

    for(auto& rm : m_vecRepairMethods) {
        rm->setSolTracker(m_pSolTrackerCallBack);
    }

    std::vector<std::vector<int>> vecCompatibleRepair(m_vecDestroyMethods.size());
    for(auto i = 0u; i < m_vecRepairMethods.size(); ++i) {
        if(!m_vecRepairDestroyAllCompatible[i]) {
            const std::vector<int>& vecCompatibleDestroyIds = m_vecRepairCompatibleDestroyIds[i];
            for(const auto& compatibleDestroyId : vecCompatibleDestroyIds) {
                vecCompatibleRepair[compatibleDestroyId].push_back(i);
            }
        } else {
            for(auto iDestroyIdx = 0u; iDestroyIdx < m_vecDestroyMethods.size(); ++iDestroyIdx) {
                vecCompatibleRepair[iDestroyIdx].push_back(i);
            }
        }
    }

#ifndef BATCH_MODE
    m_curSol.consCalc();
#endif

    if(m_pSolTrackerCallBack) {
        m_pSolTrackerCallBack->callback(m_curSol, true, true, true);
    }

    // Main loop
    while(++m_iIter < m_par.m_iMaxIter) {
        double elapsedSec = wallClock() - m_dStartTime;

        if(elapsedSec > m_par.m_dTimeouSec) {
            break;
        }

        double bestObj = m_bestSol.getCost();
        double curObj = m_curSol.getCost();

#ifndef BATCH_MODE
        if((m_iIter % m_iPrintOutPutEveryNIterations) == 0) {
            std::cout << m_par.outputPrefix << m_iIter << " (remains: " << m_par.m_dTimeouSec - elapsedSec << "), best solution: " << bestObj
                      << ", current solution: " << m_curSol.getCost() << "\n";
        }
#endif

        int iChosenDestroy = rouletteWheelSelection(m_vecDestroyWeights, randGen);
        int iChosenRepair = rouletteWheelSelection(vecCompatibleRepair[iChosenDestroy], m_vecRepairWeights, randGen);

#ifndef BATCH_MODE
        m_curSol.consCalc();
#endif

        CVRPSolution tempSol(m_curSol);

        m_vecDestroyMethods[iChosenDestroy]->destroySolution(tempSol, randGen);

#ifndef BATCH_MODE
        tempSol.consCalc();
#endif

        m_vecRepairMethods[iChosenRepair]->repairSolution(tempSol, randGen);

#ifndef BATCH_MODE
        tempSol.consCalc();
#endif

        double newObj = tempSol.getCost();

        m_acceptanceCriterion->updateParameters(elapsedSec);

        if(m_par.m_bUseDecomposition && m_iIter >= m_iNextDecompositionIter) {
            if(!m_decomposition) {
                throw std::runtime_error("Using decomposition, but m_decomposition is nullptr!\n");
            }

            // Decomposition phase if on the main thread
            m_iNextDecompositionIter += m_par.m_iDecompositionIters;

#ifndef BATCH_MODE
            cout << m_par.outputPrefix << "Decomposition at iter " << m_iIter << " - current solution cost before: " << tempSol.getCost() << "\n";
#endif

            m_par.outputPrefix += "\t";
            m_par.m_iDecompositionLevel += 1;
            tempSol = m_decomposition->decompose(tempSol, elapsedSec);
            newObj = tempSol.getCost();
            m_par.m_iDecompositionLevel -= 1;
            m_par.outputPrefix.pop_back();

#ifndef BATCH_MODE
            cout << m_par.outputPrefix << "Decomposition at iter " << m_iIter << " - current solution cost after: " << tempSol.getCost() << "\n";
#endif
        }

        bool bGlobalBest = false;
        bool bImproved = false;
        bool bAccepted = m_acceptanceCriterion->shouldAccept(bestObj, newObj, m_dEpsilon);

        if(bAccepted) {
            m_curSol = tempSol;

            if(newObj < curObj - m_dEpsilon) {
                bImproved = true;
            }

            if(newObj < bestObj - m_dEpsilon) {
                bestObj = newObj;
                m_bestSol = tempSol;
                bGlobalBest = true;
            }

            // Accepted solution: a good moment to record "good" solutions' characteristics:
            for(auto r = 0; r < tempSol.getNRoutes(); ++r) {
                const auto& route = tempSol.getRoute(r);
                if(route.empty()) {
                    continue;
                }

                // --- --- --- UPDATE VFREQ
                for(auto i : route.getNodes()) {
                    if(i == 0 || i == m_instance.getN() + 1) {
                        continue;
                    }

                    for(auto j : route.getNodes()) {
                        if(j == 0 || j == m_instance.getN() + 1) {
                            continue;
                        }

                        vfreq[i - 1][j - 1] += 1.0;
                    }
                }

                // --- --- --- UPDATE AFREQ
                for(auto idx = 1; idx < route.getNodes().size() - 2; ++idx) {
                    const auto i = route.getNodes()[idx];
                    const auto j = route.getNodes()[idx + 1];

                    assert(i > 0);
                    assert(i <= m_instance.getN());
                    assert(j > 0);
                    assert(j <= m_instance.getN());

                    assert(i - 1 < afreq.size());
                    assert(j - 1 < afreq[i - 1].size());

                    afreq[i - 1][j - 1] += 1.0;
                }

                // --- --- --- UPDATE PFREQ
                for(auto it = route.getNodes().begin(); it <= route.getNodes().end() - m_par.m_iPathLength; ++it) {
                    const auto path = std::vector<int>(it, it + m_par.m_iPathLength);
                    const auto trieIt = pfreq.find(path);

                    if(trieIt != pfreq.end()) {
                        ++(trieIt->second);
                    } else {
                        pfreq[path] = 1;
                    }
                }
            }
        }

#ifndef BATCH_MODE
        if(bGlobalBest) {
            cout << m_par.outputPrefix << m_iIter << " (remains: " << m_par.m_dTimeouSec - elapsedSec << "), New best solution: " << tempSol.getCost() << ". "
                 << "Destroy method: " << m_vecDestroyDescr[iChosenDestroy] << ". Repair method: " << m_vecRepairDescr[iChosenRepair] << "\n";
        }
#endif

        if(m_pSolTrackerCallBack) {
            m_pSolTrackerCallBack->callback(tempSol, bAccepted, bImproved, bGlobalBest);
        }

        double dScore = calcScore(bAccepted, bImproved, bGlobalBest);

        m_vecDestroyWeights[iChosenDestroy] = m_vecDestroyWeights[iChosenDestroy] * m_par.m_dALNSDecay + dScore * (1 - m_par.m_dALNSDecay);
        m_vecRepairWeights[iChosenRepair] = m_vecRepairWeights[iChosenRepair] * m_par.m_dALNSDecay + dScore * (1 - m_par.m_dALNSDecay);
    }
}

double PALNS::calcScore(bool bAccepted, bool bImproved, bool bNewGlobalBest) const {
    double dScore = 1;

    if(bAccepted) {
        dScore = std::max(dScore, m_par.m_dALNSScoreAccepted);
    }

    if(bImproved) {
        dScore = std::max(dScore, m_par.m_dALNSScoreImproved);
    }

    if(bNewGlobalBest) {
        dScore = std::max(dScore, m_par.m_dALNSScoreGlobalBest);
    }

    return dScore;
}

int PALNS::rouletteWheelSelection(const vector<double>& vecWeights, TSRandom& randGen) const {
    double dSumWeights = accumulate(vecWeights.begin(), vecWeights.end(), 0.0);
    const double dRand = randGen.getRandomDouble(0.0, dSumWeights);

    double dSummedWeights = 0;
    int iMethod = -1; // This will store the chosen index

    for(int i = 0; i < (int)vecWeights.size() - 1; ++i) {
        dSummedWeights += vecWeights[i];

        if(dRand <= dSummedWeights) {
            iMethod = i;
            break;
        }
    }
    if(iMethod < 0) {
        iMethod = (int)vecWeights.size() - 1;
    }

    return iMethod;
}

int PALNS::rouletteWheelSelection(const vector<int>& vecCompatibleMethods, const vector<double>& vecWeights, TSRandom& randGen) const {
    // Like the previous method, but here only certain elements can be chosen, because they
    // correspond to repair methods that are compatible with the chosen destroy method.

    double dSumWeights = 0.0;

    for(const auto& compMethodId : vecCompatibleMethods) {
        dSumWeights += vecWeights[compMethodId];
    }

    const double dRand = randGen.getRandomDouble(0.0, dSumWeights);

    double dSummedWeights = 0;
    int iMethod = -1;

    for(int i = 0; i < (int)vecCompatibleMethods.size() - 1; ++i) {
        dSummedWeights += vecWeights[vecCompatibleMethods[i]];

        if(dRand <= dSummedWeights) {
            iMethod = vecCompatibleMethods[i];
            break;
        }
    }

    if(iMethod < 0) {
        iMethod = vecCompatibleMethods.back();
    }

    return iMethod;
}