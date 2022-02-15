#ifndef ARCCOSTTRACKER_H
#define ARCCOSTTRACKER_H

#include "CVRPInstance.h"

// From PALNS
#include "PALNS.h"

// From SRP-Utils
#include "Matrix.h"

#include <cfloat>

class CVRPSolution;

class ArcCostTracker : public AbstractSolTrackerCallBack
{
public:
	ArcCostTracker(const CVRPInstance &instance) : m_instance(instance), m_matBestSolCost(instance.getN()+2, DBL_MAX/(10*instance.getN()))
	{}
	ArcCostTracker(const ArcCostTracker &other) : m_instance(other.m_instance), m_matBestSolCost(other.m_matBestSolCost)
	{}
	AbstractSolTrackerCallBack* clone() const			{ return new ArcCostTracker(*this); }
	
	void callback(const CVRPSolution &sol, bool bAccepted, bool bImproved, bool bNewBest);
	
	inline double getSolCost(int i, int j) const					{ return m_matBestSolCost.getElement(i,j); }
	const Matrix<double> &getCostMatrix() const						{ return m_matBestSolCost; }

protected:
	const CVRPInstance &m_instance;
	// For each arc (i,j) m_matBestSolCost stores the cost of the best solution observed that contained (i,j).
	Matrix<double> m_matBestSolCost;
};

#endif