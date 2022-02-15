#include "PALNSParameters.h"

// From SRP-Utils
#include "ParFileParser.h"

#include <iostream>
#include <limits>

void PALNSParameters::setDefault() {
    // Default parameters:

    m_deviationPars = RecordDeviationParameters(0.1);

    m_iMaxIter = std::numeric_limits<int>::max();
    m_dTimeouSec = 1200;
    m_dALNSScoreAccepted = 2;
    m_dALNSScoreImproved = 4;
    m_dALNSScoreGlobalBest = 10;
    m_dALNSDecay = 0.99;

    m_iDecompositionMethod = 1;
    m_bUseDecomposition = false;
    m_iDecompositionIters = 10000;
    m_iDecompositionSize = 200;
    m_iDecompositionLevel = 0;
    m_iPathLength = 3;
    outputPrefix = "";
}

PALNSParameters::PALNSParameters(const ParFileParser& pfp) {
    // Start out with the default parameters. Modify values using the entries in ParFileParser object
    setDefault();

    // Generic ALNS options
    pfp.getIntParameter("MaxIter", m_iMaxIter);
    pfp.getDoubleParameter("TimeoutSec", m_dTimeouSec);

    // Record Deviation
    pfp.getDoubleParameter("StartDeviation", m_deviationPars.m_dStartDeviation);

    // ALNS specific (heuristic selection)
    pfp.getDoubleParameter("ALNSScoreAccepted", m_dALNSScoreAccepted);
    pfp.getDoubleParameter("ALNSScoreImproved", m_dALNSScoreImproved);
    pfp.getDoubleParameter("ALNSScoreGlobalBest", m_dALNSScoreGlobalBest);
    pfp.getDoubleParameter("ALNSDecay", m_dALNSDecay);

    // Decomposition options
    pfp.getBoolParameter("UseDecomposition", m_bUseDecomposition);
    pfp.getIntParameter("DecompositionMethod", m_iDecompositionMethod);
    pfp.getIntParameter("DecompositionSize", m_iDecompositionSize);
    pfp.getIntParameter("DecompositionIters", m_iDecompositionIters);
    pfp.getIntParameter("DecompositionLevel", m_iDecompositionLevel);
    pfp.getParameter("OutputPrefix", outputPrefix);
}