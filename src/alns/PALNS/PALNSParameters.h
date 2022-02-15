#ifndef PALNSPARAMETERS
#define PALNSPARAMETERS

#include <string>

class ParFileParser;

struct PALNSParameters {
  struct RecordDeviationParameters {
    double m_dStartDeviation;

    RecordDeviationParameters() = default;

    RecordDeviationParameters(double sd) : m_dStartDeviation(sd) {}
  };

  RecordDeviationParameters m_deviationPars;

  // Maximum number of iterations
  int m_iMaxIter;

  // Timeout in seconds
  double m_dTimeouSec;

  // Parameters for adaptive weight adjustment
  double m_dALNSScoreAccepted;
  double m_dALNSScoreImproved;
  double m_dALNSScoreGlobalBest;
  double m_dALNSDecay;

  // Decomposition parameters
  int m_iDecompositionMethod;
  bool m_bUseDecomposition = false;
  int  m_iDecompositionIters = 10000;
  int  m_iDecompositionSize = 200;
  int  m_iDecompositionLevel = 0;
  int  m_iPathLength = 3;

  std::string outputPrefix = "";

  PALNSParameters() { setDefault(); }

  explicit PALNSParameters(const ParFileParser& pfp);

  void setDefault();
};

#endif