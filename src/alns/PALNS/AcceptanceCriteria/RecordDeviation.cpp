#include "RecordDeviation.h"

void RecordDeviation::initialise() {
    m_dDeviation = this->m_par.m_deviationPars.m_dStartDeviation;
    m_dDeviationStep = this->m_par.m_deviationPars.m_dStartDeviation / this->m_par.m_dTimeouSec;
}

void RecordDeviation::updateParameters(double elapsedSec) { m_dDeviation = this->m_par.m_deviationPars.m_dStartDeviation - elapsedSec * m_dDeviationStep; }

bool RecordDeviation::shouldAccept(double bestObj, double newObj, double eps) {
    const double gap = (newObj - bestObj) / newObj;
    return (gap < m_dDeviation - eps);
}