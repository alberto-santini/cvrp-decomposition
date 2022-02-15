#ifndef PALNS_H
#define PALNS_H

#include "PALNSParameters.h"
#include "PALNSHelpers.h"

#include "PALNSAcceptanceCriterion.h"
#include "PALNSDecompositionMethod.h"

#include "../PALNS-CVRP/CVRPInstance.h"
#include "../PALNS-CVRP/CVRPSolution.h"

#include "Utils.h"
#include "VectorUtils.h"
#include "ParFileParser.h"
#include "StatUtils.h"

#include <vector>
#include <cfloat>
#include <memory>
#include <stdexcept>
#include <numeric>
#include <ext/pb_ds/assoc_container.hpp>

// In the following, we assume that the template class CVRPSolution implements:
// - The following methods: 
//			double getCost() const
//			void postProcess() [eventually empty, if no post-processing on the solution is needed]
// - A constructor that takes a <CVRPInstance> object as input
// - A copy constructor CVRPSolution(const CVRPSolution &sol) [eventually the compiler-generated one]
// - An assignment operator [eventually the compiler-generated one]

struct PALNS {
  using FrequencyMatrix = std::vector<std::vector<double>>;
  using RecordedPath = std::vector<int>;

  struct AccTraits {
    using size_type = RecordedPath::size_type;
    using key_type = RecordedPath;
    using key_const_reference = const RecordedPath&;
    using e_type = std::size_t;
    using const_iterator = RecordedPath::const_iterator;
    enum {
      max_size = 1100
    };

    inline static const_iterator begin(key_const_reference k) { return k.begin(); }
    inline static const_iterator end(key_const_reference k) { return k.end(); }
    inline static std::size_t e_pos(e_type e) { return e; }
  };

  using PathFrequencyMatrix = __gnu_pbds::trie<RecordedPath, std::uint32_t, AccTraits>;

  explicit PALNS(const CVRPInstance& instance) :
      m_instance(instance),
      m_bestSol(instance),
      m_curSol(instance),
      m_initialSol(instance),
      m_pSolTrackerCallBack(nullptr),
      m_iMasterSeed(std::time(nullptr)) {}

  int addDestroyMethod(AbstractDestroyMethod& destroyMethod, const string& strDescr);
  void addRepairMethod(AbstractRepairMethod& repairMethod, const string& strDescr, bool bAllCompatible = true, const vector<int>& vecCompatibleDestroyIds = vector<int>());

  void setSolTrackerCallBack(AbstractSolTrackerCallBack& solTrackerCallBack) { m_pSolTrackerCallBack = &solTrackerCallBack; }

  void go(CVRPSolution& sol, PALNSParameters par, InitialSolutionCreator* pISC);
  void start();

  double calcScore(bool bAccepted, bool bImproved, bool bNewGlobalBest) const;
  int rouletteWheelSelection(const vector<double>& vecWeights, TSRandom& randGen) const;
  int rouletteWheelSelection(const vector<int>& vecCompatibleMethods, const vector<double>& vecWeights, TSRandom& randGen) const;

  void resetLocalParameters(CVRPSolution& solm, PALNSParameters& par);

  static const int m_iPrintOutPutEveryNIterations = 1000;

  AbstractSolTrackerCallBack* m_pSolTrackerCallBack; // Eventual callback on the solution
  vector<AbstractDestroyMethod*> m_vecDestroyMethods; // List of destroy methods
  vector<AbstractRepairMethod*> m_vecRepairMethods; // List of repair methods
  CVRPSolution m_initialSol; // Initial solution for the current run

  FrequencyMatrix vfreq; // vfreq[i-1][j-1] = How many times i and j were in the same route in good solutions.
  FrequencyMatrix afreq; // afreq[i-1][j-1] = How many times arc (i,j) was in a good solution.
  PathFrequencyMatrix pfreq; // pfreq[path] = How many times path was in a good solution.

  CVRPSolution m_curSol; // Current solution
  CVRPSolution m_bestSol; // Best solution so far
  int m_iIter; // Number of iterations
  int m_iNextDecompositionIter; // Next iter when decomposition will happen
  std::unique_ptr<AcceptanceCriterion> m_acceptanceCriterion; // The acceptance criterion to be used
  std::unique_ptr<DecompositionMethod> m_decomposition; // Decomposition method
  std::vector<double> m_vecDestroyWeights; // Destroy heuristics weights
  std::vector<double> m_vecRepairWeights; // Repair heuristics weights

  const CVRPInstance m_instance; // Instance of the problem
  PALNSParameters m_par; // Program parameters
  
  double m_dEpsilon; // EPS to be used in floating point comparisons
  int m_iMasterSeed; // General random seed
  double m_dStartTime; // Start time of the algorithm.

  std::vector<string> m_vecDestroyDescr; // Names of destroy heuristics
  std::vector<string> m_vecRepairDescr; // Names of repair heuristics
  std::vector<bool> m_vecRepairDestroyAllCompatible; // Are repair heuristics able to repair from any destroy method?
  std::vector<std::vector<int>> m_vecRepairCompatibleDestroyIds; // Compatibility table for repair and destroy methods
};

#endif
