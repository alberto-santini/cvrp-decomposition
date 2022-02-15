#ifndef ROUTE_H
#define ROUTE_H

#include "CVRPInstance.h"

// From SRP-Utils
#include "Utils.h"
#include "Matrix.h"

#include <vector>
#include <float.h>

using namespace std;


// First element in a route should always be node 0 and last element should always be
// node n+1.

class Route
{
public:
	Route(const CVRPInstance &instance);

	int size() const					{ return (int) m_vecRoute.size(); }
	int getLoad() const					{ return m_iLoad; }
	double getDuration() const			{ return m_dDuration; }
	double getCost() const				{ return m_dCost; }

	int getFirstCustomer() const		{ return m_vecRoute[1]; }
	int getLastCustomer() const			{ return m_vecRoute[(int) m_vecRoute.size()-2]; }
	int getNodeId(int iPos) const		{ return m_vecRoute[iPos]; }
	int operator[](int iPos) const		{ return m_vecRoute[iPos]; }
	bool empty() const					{ return m_vecRoute.size() == 2; }
	const vector<int> &getNodes() const	{ return m_vecRoute; }
	void setNodes(const vector<int> &vecNodes);
	void setCustomers(const vector<int> &vecCustomers);
	inline void clear();

	void concatRoute(const Route &otherRoute);
	void concatRouteInv(const Route &otherRoute);
	void concatRouteInvThis(const Route &otherRoute);

	// Local search
	bool twoOpt(double dImproveEps);

	inline bool getInsertCost(int iCustId, int iPos, double &cost) const;
	inline bool findBestInsertion(int iCustId, int &iPos, double &dCost);
	inline bool altCostFindBestIns(int iCustId, int &iBestPos, double &dBestInc, const Matrix<double> &matAltCost);
	inline bool findBestSeqInsertion(int iFirstCust, int iLastCust, double dSeqDuration, int &iBestPos, double &dBestInc, bool &bReversed);
	inline double insertCust(int iCustId, int iPos);
	inline double insertSeq(const vector<int> &vecSeq, int iPos);
	inline double getRemoveDelta(int iPos) const;
	void removeCustomer(int iPos);
	void removeCustomers(const vector<bool> &vecIsCustRemoved, vector<int> &vecRemovedCusts);
	void removeSequence(int iStartPos, int iNToRemove);

	// Updates the datastructures on the route. 
	// The method terminates with an error message if an inconsistency is found.
	void consCalc();
	
	
	friend ostream& operator << (ostream& os, const Route &route);

private:
	void do2optMove(int i, int j);

	vector<int> m_vecRoute;
	
	int m_iLoad;
	double m_dCost;
	double m_dDuration;
	
	// How much do we accept the calculated route cost and stored cost to vary when we call the 
	// consCalc method? (they can vary because of rounding errors).
	double m_dToleranceEps;

	const CVRPInstance *m_pInstance;
	const Matrix<double> *m_pMatDists;
};

// inline methods (to speed up program):

/***********************************************************************
 * bool Route::getInsertCost(int iCustId, int iPos, double &deltaCost)
 *
 * Calculates the cost of inserting a customer into the route.
 *
 * --- input parameters ---
 * iCustId		: The customer to insert.
 * iPos			: The position to insert at. The customer is inserted before the visit at position <iPos>. 
 *                iPos = 0 is invalid (position 0 is occupied by the depot).
 * --- output parameters ---
 * deltaCost	: How much more expensive does the route get after inserting the customer
 * --- return value ---
 * true if and only if the insertion is possible.
 ***********************************************************************/

inline bool Route::getInsertCost(int iCustId, int iPos, double &deltaCost) const
{
#ifdef _DEBUG
	if (iPos <= 0 || iPos >= (int) m_vecRoute.size())
		error("Route::getInsertCost(...)", "insertion position out of range ");
#endif
	double inc = 
		m_pMatDists->getElement(m_vecRoute[iPos-1], iCustId) +
		m_pMatDists->getElement(iCustId, m_vecRoute[iPos]) -
		m_pMatDists->getElement(m_vecRoute[iPos-1], m_vecRoute[iPos]);

	if (m_dDuration + inc + m_pInstance->getServiceTime(iCustId) > m_pInstance->getMaxDuration())
		return false;
	if (m_iLoad + m_pInstance->getDemand(iCustId) > m_pInstance->getMaxLoad())
		return false;

	deltaCost = inc;
	return true;
}

/***********************************************************************
 * bool Route::findBestInsertion(int iCustId, int &iBestPos, double &dBestInc)
 *
 * Find the best insertion of customer <iCustId> into this route.
 *
 * --- input parameters ---
 * iCustId		: The customer to insert.
 * --- output parameters ---
 * iBestPos		: The best position for inserting the customer (before the node at position <iBestPos>.
 * dBestInc		: How much the cost increase by inserting the customer at the best position.
 *				  If no insertion is feasible this variable has the value DBL_MAX.
 * --- return value ---
 * true if and only if an insertion is possible.
 ***********************************************************************/

inline bool Route::findBestInsertion(int iCustId, int &iBestPos, double &dBestInc)
{
	dBestInc = DBL_MAX;
	if (m_iLoad + m_pInstance->getDemand(iCustId) > m_pInstance->getMaxLoad())
		return false;
	
	int iPos;
	int iMaxPos = (int) m_vecRoute.size();
	assert(iMaxPos >= 2);
	// Creating variables for previous and next customer instead of looking up in <m_vecRoute> all the time created a nice speedup
	// (~122 sec vs. ~142s on Kelly02 [two processors Core 2 6300, 1.86GHz]).
	int iPrevCust = m_vecRoute[0];
	int iNextCust;
	for (iPos = 1; iPos < iMaxPos; ++iPos)
	{
		iNextCust = m_vecRoute[iPos];
		double dInc = 
			m_pMatDists->getElement(iPrevCust, iCustId) +
			m_pMatDists->getElement(iCustId, iNextCust) -
			m_pMatDists->getElement(iPrevCust, iNextCust);
		iPrevCust = iNextCust;
		if (dInc < dBestInc)
		{
			dBestInc = dInc;
			iBestPos = iPos;
		}
	}
	// Check duration
	if (m_dDuration + dBestInc + m_pInstance->getServiceTime(iCustId) > (double) m_pInstance->getMaxDuration())
	{
		dBestInc = DBL_MAX;
		return false;
	}
	else
		// There should always be one insertion position (see assert in the beginning 
		// of the method)).
		return true;
}

/***********************************************************************
 * bool Route::findBestSeqInsertion(int iFirstCust, int iLastCust, double dSeqDuration, int &iBestPos, double &dBestInc, bool &bReversed)
 *
 * Find the best insertion of a sequence of customers into this route. The sequence
 * is represented by its first and last customer and duration. We "forget" about the intermediate
 * customers. We also assume that the sequence fits within the route "capacity wise" - that 
 * should have been checked up front.
 *
 * @@@ The method assumes that the distance matrix is symmetric
 *
 * --- input parameters ---
 * iFirstCust		: The first customer in the sequence
 * iLastCust		: last customer in the sequence
 * dSeqDuration		: duration of the sequence
 * --- output parameters ---
 * iBestPos		: The best position for inserting the sequence (before the node at position <iBestPos>.
 * dBestInc		: How much the cost increase by inserting the sequence at the best position.
 *				  If no insertion is feasible this variable has the value DBL_MAX.
 * bReversed    : True if the sequence has to be reversed before insertion.
 * --- return value ---
 * true if and only if an insertion is possible.
 ***********************************************************************/

inline bool Route::findBestSeqInsertion(int iFirstCust, int iLastCust, double dSeqDuration, int &iBestPos, double &dBestInc, bool &bReversed)
{
	dBestInc = DBL_MAX;

	int iPos;
	int iMaxPos = (int) m_vecRoute.size();
	assert(iMaxPos >= 2);
	// Creating variables for previous and next customer instead of looking up in <m_vecRoute> all the time created a nice speedup
	int iPrevCust = m_vecRoute[0];
	int iNextCust;
	for (iPos = 1; iPos < iMaxPos; ++iPos)
	{
		iNextCust = m_vecRoute[iPos];
		double dInc1 = 
			m_pMatDists->getElement(iPrevCust, iFirstCust) +
			m_pMatDists->getElement(iLastCust, iNextCust) -
			m_pMatDists->getElement(iPrevCust, iNextCust);
		// cost if sequence is reversed before insertion:
		double dInc2 = 
			m_pMatDists->getElement(iPrevCust, iLastCust) +
			m_pMatDists->getElement(iFirstCust, iNextCust) -
			m_pMatDists->getElement(iPrevCust, iNextCust);
		iPrevCust = iNextCust;
		if (dInc1 < dBestInc)
		{
			dBestInc = dInc1;
			iBestPos = iPos;
			bReversed = false;
		}
		if (dInc2 < dBestInc)
		{
			dBestInc = dInc2;
			iBestPos = iPos;
			bReversed = true;
		}
	}
	// Check duration
	if (m_dDuration + dBestInc + dSeqDuration > (double) m_pInstance->getMaxDuration())
	{
		dBestInc = DBL_MAX;
		return false;
	}
	else
		// There should always be one insertion position (see assert in the beginning 
		// of the method)).
		return true;
}

/***********************************************************************
 * bool Route::altCostFindBestIns(int iCustId, int &iBestPos, double &dBestInc, const Matrix<double> &matAltCost)
 *
 * Find the best insertion of customer <iCustId> into this route. Use an alternative cost matrix (given by matAltCost).
 *
 * --- input parameters ---
 * iCustId		: The customer to insert.
 * matAltCost	: The alternative cost marix.
 * --- output parameters ---
 * iBestPos		: The best position for inserting the customer (before the node at position <iBestPos>.
 * dBestInc		: How much the cost increase by inserting the customer at the best position.
 *				  If no insertion is feasible this variable has the value DBL_MAX.
 * --- return value ---
 * true if and only if an insertion is possible.
 ***********************************************************************/

inline bool Route::altCostFindBestIns(int iCustId, int &iBestPos, double &dBestInc, const Matrix<double> &matAltCost)
{
	dBestInc = DBL_MAX;
	if (m_iLoad + m_pInstance->getDemand(iCustId) > m_pInstance->getMaxLoad())
		return false;
	
	int iPos;
	iBestPos = -1;
	int iMaxPos = (int) m_vecRoute.size();
	assert(iMaxPos >= 2);
	double dBaseDur = m_dDuration + m_pInstance->getServiceTime(iCustId);
	double dMaxDur = (double) m_pInstance->getMaxDuration();
	// Creating variables for previous and next customer instead of looking up in <m_vecRoute> all the time 
	// also created a nice speedup here (see findBestInsertion).
	int iPrevCust = m_vecRoute[0];
	int iNextCust;
	for (iPos = 1; iPos < iMaxPos; ++iPos)
	{
		iNextCust = m_vecRoute[iPos];
		double dInc = 
			matAltCost.getElement(iPrevCust, iCustId) +
			matAltCost.getElement(iCustId, iNextCust) -
			matAltCost.getElement(iPrevCust, iNextCust);
		double dDurInc = 
			m_pMatDists->getElement(iPrevCust, iCustId) +
			m_pMatDists->getElement(iCustId, iNextCust) -
			m_pMatDists->getElement(iPrevCust, iNextCust);
		iPrevCust = iNextCust;
		if (dBaseDur + dDurInc <= dMaxDur && dInc < dBestInc)
		{
			dBestInc = dInc;
			iBestPos = iPos;
		}
	}
	if (iBestPos == -1)
		return false;
	else
		return true;
}

/***********************************************************************
 * double Route::insertCust(int iCustId, int iPos)
 *
 * insert customer <iCustId> on position <iPos>
 *
 * --- input parameters ---
 * iCustId		: The customer to insert.
 * iPos			: The position to insert at. The customer is inserted before the visit at position <iPos>. 
 *                iPos = 0 is invalid (position 0 is occupied by the depot).
 * --- output parameters ---
 * --- return value ---
 * the delta cost of performing the insertion.
 ***********************************************************************/

inline double Route::insertCust(int iCustId, int iPos)
{
	m_iLoad += m_pInstance->getDemand(iCustId);
	// Load should not be violated:
	assert(m_iLoad <= m_pInstance->getMaxLoad());
	assert(0 < iPos && iPos < (int) m_vecRoute.size());

	double dInc = 
		m_pMatDists->getElement(m_vecRoute[iPos-1], iCustId) +
		m_pMatDists->getElement(iCustId, m_vecRoute[iPos]) -
		m_pMatDists->getElement(m_vecRoute[iPos-1], m_vecRoute[iPos]);

	m_dDuration += dInc + m_pInstance->getServiceTime(iCustId);
	// Duration should not be violated.
	assert(m_dDuration <= (double) m_pInstance->getMaxDuration());
	m_dCost += dInc;
	// Perform the actual insertion:
	m_vecRoute.insert(m_vecRoute.begin()+iPos,iCustId);
	return dInc;
}

inline double Route::insertSeq(const vector<int> &vecSeq, int iPos)
{
	assert(0 < iPos && iPos < (int) m_vecRoute.size());

	double dInc = 
		m_pMatDists->getElement(m_vecRoute[iPos-1], vecSeq.front()) +
		m_pMatDists->getElement(vecSeq.back(), m_vecRoute[iPos]) -
		m_pMatDists->getElement(m_vecRoute[iPos-1], m_vecRoute[iPos]);

	int i;
	int iServTime = 0;
	int iLoad = 0;
	for (i=0; i < (int) vecSeq.size(); ++i)
	{
		int iCustId = vecSeq[i];
		if (i > 0)
			dInc += m_pMatDists->getElement(vecSeq[i-1], iCustId);
		iServTime += m_pInstance->getServiceTime(iCustId);
		iLoad += m_pInstance->getDemand(iCustId);
	}

	m_iLoad += iLoad;
	// Load should not be violated:
	assert(m_iLoad <= m_pInstance->getMaxLoad());
		
	m_dDuration += dInc + iServTime;
	// Duration should not be violated.
	assert(m_dDuration <= (double) m_pInstance->getMaxDuration());
	m_dCost += dInc;
	// Perform the actual insertion:
	m_vecRoute.insert(m_vecRoute.begin()+iPos, vecSeq.begin(), vecSeq.end());
#ifdef _DEBUG
	consCalc();
#endif
	return dInc;
}

inline double Route::getRemoveDelta(int iPos) const
{
	// It has to be a customer that we are removing.
	assert(0 < iPos && iPos < (int) m_vecRoute.size()-1);

	return 
		m_pMatDists->getElement(m_vecRoute[iPos-1], m_vecRoute[iPos+1]) 
		- m_pMatDists->getElement(m_vecRoute[iPos-1], m_vecRoute[iPos]) 
		- m_pMatDists->getElement(m_vecRoute[iPos], m_vecRoute[iPos+1]);
}

inline void Route::clear()
{
	int iN = m_pInstance->getN();
	m_vecRoute.resize(2);	
	m_vecRoute[0] = 0;
	m_vecRoute[1] = iN+1;

	m_iLoad = 0;
	m_dCost = m_pMatDists->getElement(0, iN+1);
	m_dDuration = m_pMatDists->getElement(0, iN+1);
}


#endif
