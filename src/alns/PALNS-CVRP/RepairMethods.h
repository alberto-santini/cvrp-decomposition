#ifndef REPAIRMETHODS_H
#define REPAIRMETHODS_H

#include "PALNS.h"
#include "CVRPSolution.h"

class ArcCostTracker;

class Regret2Repair : public AbstractRepairMethod
{
public:
	// Copy constructor:
	Regret2Repair(const Regret2Repair &other) 
		: m_dRandFactor(other.m_dRandFactor), m_dInterRouteImprProb(other.m_dInterRouteImprProb) 
	{}
	Regret2Repair(double dRandFactor, double dInterRouteImprProb);
	
	void repairSolution(CVRPSolution &sol, TSRandom &randGen);
	virtual AbstractRepairMethod* clone() const			{ return new Regret2Repair(*this); }

protected:
	double m_dRandFactor; 
	double m_dInterRouteImprProb;
};

class HistoryRegret2Repair : public AbstractRepairMethod
{
public:
	// Copy constructor:
	HistoryRegret2Repair(const HistoryRegret2Repair &other) 
		: m_dRandFactor(other.m_dRandFactor), m_pArcCostTracker(other.m_pArcCostTracker), m_dInterRouteImprProb(other.m_dInterRouteImprProb) 
	{}
	HistoryRegret2Repair(double dRandFactor, double dInterRouteImprProb) 
		: m_dRandFactor(dRandFactor), m_dInterRouteImprProb(dInterRouteImprProb)
	{ m_pArcCostTracker = NULL; }
	
	void repairSolution(CVRPSolution &sol, TSRandom &randGen);
	virtual AbstractRepairMethod* clone() const				{ return new HistoryRegret2Repair(*this); }
	void setSolTracker(const AbstractSolTrackerCallBack *pST)	{ m_pArcCostTracker = (ArcCostTracker*) pST; }

protected:
	double m_dRandFactor; 
	const ArcCostTracker *m_pArcCostTracker;
	double m_dInterRouteImprProb;
};

class SequenceRegretRepair : public AbstractRepairMethod
{
public:
	// Copy constructor:
	SequenceRegretRepair(const SequenceRegretRepair &other) 
		: m_dRandFactor(other.m_dRandFactor), m_dInterRouteImprProb(other.m_dInterRouteImprProb)
	{}
	SequenceRegretRepair(double dRandFactor, double dInterRouteImprProb);
	
	void repairSolution(CVRPSolution &sol, TSRandom &randGen);
	virtual AbstractRepairMethod* clone() const			{ return new SequenceRegretRepair(*this); }

protected:
	double m_dRandFactor; 
	double m_dInterRouteImprProb;
};

#endif
