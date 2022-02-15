#ifndef DESTROYMETHODS_H
#define DESTROYMETHODS_H

#include "PALNS.h"

#include "CVRPSolution.h"
#include "ArcCostTracker.h"

class RandomDestroy : public AbstractDestroyMethod
{
public:
	// Copy constructor:
	RandomDestroy(const RandomDestroy &other)							{};
	RandomDestroy()														{};

	void destroySolution(CVRPSolution &sol, TSRandom &randGen);
	virtual AbstractDestroyMethod* clone() const			{ return new RandomDestroy(*this); }
};

class LocalRandomDestroy : public AbstractDestroyMethod
{
public:
	// Copy constructor:
	LocalRandomDestroy(const LocalRandomDestroy &other);	
	LocalRandomDestroy(int iNeighborhoodSize, int iN);
	

	void destroySolution(CVRPSolution &sol, TSRandom &randGen);
	virtual AbstractDestroyMethod* clone() const			{ return new LocalRandomDestroy(*this); }
protected:
	int m_iNeighborhoodSize;
	vector<int> m_vecNodeCount;
};

struct RemoveDelta;

class ExpensiveNodeDestroy : public AbstractDestroyMethod
{
public:
	// Copy constructor:
	ExpensiveNodeDestroy(const ExpensiveNodeDestroy &other) : m_dRandFactor(other.m_dRandFactor) {}
	ExpensiveNodeDestroy(double dRandFactor) : m_dRandFactor(dRandFactor) {}
	
	void destroySolution(CVRPSolution &sol, TSRandom &randGen);
	virtual AbstractDestroyMethod* clone() const			{ return new ExpensiveNodeDestroy(*this); }
protected:
	void updateVecRemoveDelta(const CVRPSolution &sol, vector<RemoveDelta> &vecRemoveDelta, RemoveDelta &removeDelta) const;

	const double m_dRandFactor;
};

class GeoDestroy : public AbstractDestroyMethod
{
public:
	// Copy constructor:
	GeoDestroy(const GeoDestroy &other) : m_dRandFactor(other.m_dRandFactor) {};
	GeoDestroy(double dRandFactor) : m_dRandFactor(dRandFactor) {}
	
	void destroySolution(CVRPSolution &sol, TSRandom &randGen);
	virtual AbstractDestroyMethod* clone() const			{ return new GeoDestroy(*this); }
; 
protected:
	const double m_dRandFactor;
};

class SolHistoryDestroy : public AbstractDestroyMethod
{
public:
	// Copy constructor:
	SolHistoryDestroy(const SolHistoryDestroy &other) : m_dRandFactor(other.m_dRandFactor), m_pArcCostTracker(other.m_pArcCostTracker) {}
	SolHistoryDestroy(double dRandFactor) : m_dRandFactor(dRandFactor)		{ m_pArcCostTracker = NULL; }
	
	void destroySolution(CVRPSolution &sol, TSRandom &randGen);
	virtual AbstractDestroyMethod* clone() const				{ return new SolHistoryDestroy(*this); }
	void setSolTracker(const AbstractSolTrackerCallBack *pST)	{ m_pArcCostTracker = (ArcCostTracker*) pST; }

protected:
	void updateVecRemoveDelta(const CVRPSolution &sol, vector<RemoveDelta> &vecRemoveDelta, RemoveDelta &removeDelta) const;

	const double m_dRandFactor;
	const ArcCostTracker *m_pArcCostTracker;
};

class RandomSequenceDestroy : public AbstractDestroyMethod
{
public:
	// Copy constructor:
	RandomSequenceDestroy(const RandomDestroy &other)					{}
	RandomSequenceDestroy();

	void destroySolution(CVRPSolution &sol, TSRandom &randGen);
	virtual AbstractDestroyMethod* clone() const			{ return new RandomSequenceDestroy(*this); }

	void setParamSeqSize(int iMinSeqSize, int iMaxSeqSize)				{ m_iMaxSeqSize = iMaxSeqSize; m_iMinSeqSize = iMinSeqSize; }
	void setParamAbsNSeq(int iMinAbsNSeq, int iMaxAbsNSeq)				{ m_iMinAbsNSeq = iMinAbsNSeq; m_iMaxAbsNSeq = iMaxAbsNSeq; }
	void setParamRelNSeq(double dMinRelNSeq, double dMaxRelNSeq)		{ m_dMinRelNSeq = dMinRelNSeq; m_dMaxRelNSeq = dMaxRelNSeq; }
private:
	int m_iMinSeqSize;
	int m_iMaxSeqSize;
	int m_iMinAbsNSeq;
	int m_iMaxAbsNSeq;
	double m_dMinRelNSeq;
	double m_dMaxRelNSeq;
};


#endif