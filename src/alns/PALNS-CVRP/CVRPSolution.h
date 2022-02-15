#ifndef CVRPSOLUTION_H
#define CVRPSOLUTION_H

#include "Route.h"

// From SRP-Utils
#include "VectorUtils.h"

#include <vector>
#include <algorithm>

using namespace std;

class CVRPInstance;
class SimpleArc;
struct TailSwapCacheEntry;
struct LoadDuration;

class CVRPSolution
{
public:
	typedef double TObj;

	CVRPSolution() {}
	CVRPSolution(const CVRPInstance &instance);
	CVRPSolution(const CVRPInstance* instance) : m_pInstance(instance) {}
	void createRoute(int iCustId);
	void createRoute();
	void createRoutes(int iNRoutes);
	void mergeRoutes(int iRoute1, int iRoute2);
	void mergeRoutesInv(int iRoute1, int iRoute2);
	void mergeRoutesInvThis(int iRoute1, int iRoute2);
	double getCost() const							{ return m_dCost; }
	int getNRoutes() const							{ return (int) m_vecRoutes.size(); }
	const Route &getRoute(int iRouteId)	const		{ return m_vecRoutes[iRouteId]; }
	Route& getRoute(int iRouteId) { return m_vecRoutes[iRouteId]; }
	int getRouteId(int iNodeId) const				{ return m_vecCust2Route[iNodeId]; }
	const vector<int> &getUnassignedCusts() const	{ return m_vecUnassigned; }
	const CVRPInstance &getInstance() const			{ return *m_pInstance; }
	const vector<vector<int> > &getUnassignedSeqs() const	{ return m_vecUnassignedSeqs; }
	void clearUnassignedSeqs()						{ m_vecUnassignedSeqs.clear(); }
	void addUnassignedSeq(const vector<int> &vecSeq)	{ m_vecUnassignedSeqs.push_back(vecSeq); }
	const vector<int> &getVisits(int iRoute) const	{ return m_vecRoutes[iRoute].getNodes(); }
	double getUnassignedCost() const				{ return m_dUnassignedCost; }

	void setRouteNodes(int i_nRoute, const std::vector<int> v_routeNodes) {
	  const int iN = m_pInstance->getN();

	  m_vecRoutes[i_nRoute].clear();
	  m_vecRoutes[i_nRoute].setNodes(v_routeNodes);
	  for(auto node : v_routeNodes) {
	    if(node == 0 || node == iN + 1) { continue; }

	    m_vecCust2Route[node] = i_nRoute;
	    m_vecUnassigned.erase(std::remove(m_vecUnassigned.begin(), m_vecUnassigned.end(), node), m_vecUnassigned.end());
	    m_dCost -= m_dUnassignedCost;
	  }

	  m_dCost -= m_pInstance->getDist(0, iN + 1);
	  m_dCost += m_vecRoutes[i_nRoute].getCost();

	  consCalc();
	}
	
	void writeSimpleSol(const string &strFileName) const;
	inline bool findBestInsertion(int iRouteId, int iCustId, int &iBestPos, double &dBestInc);
	inline bool findBestSeqInsertion(int iRoute, int iFirstCust, int iLastCust, double dSeqDuration, int &iBestPos, double &dBestInc, bool &bReversed);
	inline bool altCostFindBestIns(int iRouteId, int iCustId, int &iBestPos, double &dBestInc, const Matrix<double> &matAltCost);
	inline void insertCust(int iRoute, int iBestCust, int iBestPos);
	// Notice: insertSeq does not remove the inserted sequence from <m_vecUnassignedSeqs>
	inline void insertSeq(int iRoute, const vector<int> &vecSeq, bool bReversed, int iPos);
	inline void removeCustomers(int iRoute, const vector<bool> &vecIsCustRemoved);	
	inline void removeCustomer(int iRoute, int iPos);
	inline void clearRoute(int iRoute);
	void clear();
	void clear(int iNRoutes);
	inline int countUnusedRoutes() const;
	void getAssignedCustomers(vector<int> &vecCustIds, vector<int> &vecRouteIds);

  // Hash
  unsigned int getHashCode() const;
  void getHashSequence(vector<int>& vecOrderId) const;
  double dist(const CVRPSolution& other) const;
	
	// For sequence oriented removal:
	void removeSequence(int iRouteId, int iStartPos, int iNToRemove, bool bConsCalc = true);

	// Method used for arc removal in LNS method
	void splitRoute(int iRoute, const vector<int> &vecSplitPos);

	// Demanded by PALNS framework. Does not do anything.
	void postProcess() {}

	// Simple local search methods
	bool relocateOpt();
	bool swapOpt();
	bool tailSwapOpt();
	bool route2opt();
	void combinedOpt(bool bRelocate = true, bool bSwap = true, bool b2opt = true, bool bTailSwap = true);
	void polish();
	void seqSearchTailSwapOpt();

	// Method for updating the object (recalculates route costs and so on)
	// The method terminates the program with an error message if an inconsistency or a violated
	// solution is detected.
	void consCalc();
	// Like consCalc but only validates a subset of the routes and does not check the sanity. This is used by the removal methods such that only the routes from which orders were
	// removed are validated.
	void consCalcLight(const vector<bool> &vecRouteChanged);

	bool operator== (const CVRPSolution &other) const;

	friend ostream& operator << (ostream& os, const CVRPSolution &sol);

private:
	inline double distIJK(const Matrix<double> &matDists, int iNodeI, int iNodeJ, int iNodeK) const;
	void tailSwapPrecompute(const Route &route, vector<LoadDuration> &vecLoadDuration) const;
	double evalTailSwap(const Route &route1, const Route &route2, TailSwapCacheEntry &Result, 
						const vector<LoadDuration> &vecLoadDurRoute1, const vector<LoadDuration> &vecLoadDurRoute2) const;
	void doTailSwap(int iRoute1, int iRoute2, const TailSwapCacheEntry& bestMove);
	void getArcRep(vector<SimpleArc> &vecArcs) const;

	void getGiantTour(vector<int> &vecGiantTour, vector<int> &vecNodeToPos, vector<int> &vecRouteId,
								vector<int> &vecDemandBefore, vector<int> &vecDemandAfter,
								vector<double> &vecDurationBefore, vector<double> &vecDurationAfter);

	vector<Route> m_vecRoutes;
	vector<int> m_vecCust2Route;
	// Customers that are not served:
	vector<int> m_vecUnassigned;
	// Sequences of customers that are not served (should be inserted again as sequences)
	vector<vector<int> > m_vecUnassignedSeqs;

	double m_dUnassignedCost;
	double m_dCost;
	// By how much should a move improve a solution before it's accepted. Used to avoid 
	// cycling because of numerical unstability.
	double m_dImproveEps;
	// How much do we accept the calculated route cost and stored cost to vary when we call the 
	// consCalc method? (they can vary because of rounding errors).
	double m_dToleranceEps;

	const CVRPInstance *m_pInstance;
};

inline double CVRPSolution::distIJK(const Matrix<double> &matDists, 
									int iNodeI, int iNodeJ, int iNodeK) const
{
	return matDists.getElement(iNodeI, iNodeJ) + matDists.getElement(iNodeJ, iNodeK);
}

inline bool CVRPSolution::findBestInsertion(int iRouteId, int iCustId, int &iBestPos, double &dBestInc)
{
	assert(0 <= iRouteId && iRouteId < (int) m_vecRoutes.size());
	return m_vecRoutes[iRouteId].findBestInsertion(iCustId, iBestPos, dBestInc);
}

inline bool CVRPSolution::findBestSeqInsertion(int iRouteId, int iFirstCust, int iLastCust, double dSeqDuration, int &iBestPos, double &dBestInc, bool &bReversed)
{
	assert(0 <= iRouteId && iRouteId < (int) m_vecRoutes.size());
	return m_vecRoutes[iRouteId].findBestSeqInsertion(iFirstCust, iLastCust, dSeqDuration, iBestPos, dBestInc, bReversed);
}

inline bool CVRPSolution::altCostFindBestIns(int iRouteId, int iCustId, int &iBestPos, double &dBestInc, const Matrix<double> &matAltCost)
{
	assert(0 <= iRouteId && iRouteId < (int) m_vecRoutes.size());
	return m_vecRoutes[iRouteId].altCostFindBestIns(iCustId, iBestPos, dBestInc, matAltCost);
}

inline void CVRPSolution::insertCust(int iRouteId, int iCustId, int iPos)
{
	assert(0 <= iRouteId && iRouteId < (int) m_vecRoutes.size());
	m_dCost += m_vecRoutes[iRouteId].insertCust(iCustId, iPos);
	m_dCost -= m_dUnassignedCost;
	m_vecCust2Route[iCustId] = iRouteId;
	if (!vectorRemoveSingleElement(m_vecUnassigned, iCustId))
	{
		error("CVRPSolution::insertCust(...)", 
			"Customer" + int2String(iCustId) + " was already assigned to a route.");
	}
}

// Notice: insertSeq does not remove the inserted sequence from <m_vecUnassignedSeqs>
inline void CVRPSolution::insertSeq(int iRouteId, const vector<int> &vecSeq, bool bReversed, int iPos)
{
	assert(0 <= iRouteId && iRouteId < (int) m_vecRoutes.size());
	if (bReversed)
	{
		vector<int> vecSeqRev(vecSeq);
		reverse(vecSeqRev.begin(), vecSeqRev.end());
		m_dCost += m_vecRoutes[iRouteId].insertSeq(vecSeqRev, iPos);
	}
	else
		m_dCost += m_vecRoutes[iRouteId].insertSeq(vecSeq, iPos);
	m_dCost -= m_dUnassignedCost*vecSeq.size();
	int i;
	vector<bool> vecIsInSeq(m_pInstance->getN()+1, false);
	for (i=0; i < (int) vecSeq.size(); ++i)
	{
		m_vecCust2Route[vecSeq[i]] = iRouteId;
		vecIsInSeq[vecSeq[i]] = true;
	}
	i = 0;
	while (i < (int) m_vecUnassigned.size())
	{
		if (vecIsInSeq[m_vecUnassigned[i]])
			unorderedVectorRemoveAt(m_vecUnassigned, i);
		else
			++i;
	}
}

inline int CVRPSolution::countUnusedRoutes() const
{
	int i, iCount = 0;
	for (i = 0; i < getNRoutes(); ++i)
		if (getRoute(i).empty())
			++iCount;
	return iCount;
}

// <vecIsCustRemoved> should contain (at least) iN+1 elements.
inline void CVRPSolution::removeCustomers(int iRoute, const vector<bool> &vecIsCustRemoved)
{ 
	assert((int) vecIsCustRemoved.size() >= (m_pInstance->getN()+1));

	vector<int> vecRemovedCusts;
	m_dCost -= m_vecRoutes[iRoute].getCost();
	m_vecRoutes[iRoute].removeCustomers(vecIsCustRemoved, vecRemovedCusts); 
	m_dCost += m_vecRoutes[iRoute].getCost() + m_dUnassignedCost * vecRemovedCusts.size();
	// update the vector of unassinged customers.
	m_vecUnassigned.insert(m_vecUnassigned.end(), vecRemovedCusts.begin(), vecRemovedCusts.end());
	// update <m_vecCust2Route>.
	int i;
	for (i = 0; i < (int) vecRemovedCusts.size(); ++i)
		m_vecCust2Route[vecRemovedCusts[i]] = -1;
}

inline void CVRPSolution::removeCustomer(int iRoute, int iPos)
{
	assert(0 < iPos && iPos < m_vecRoutes[iRoute].size()-1);
	int iCust = m_vecRoutes[iRoute][iPos];
	m_dCost -= m_vecRoutes[iRoute].getCost();
	m_vecRoutes[iRoute].removeCustomer(iPos); 
	m_dCost += m_vecRoutes[iRoute].getCost() + m_dUnassignedCost;
	// update the vector of unassinged customers.
	m_vecUnassigned.push_back(iCust);
	// update <m_vecCust2Route>.
	m_vecCust2Route[iCust] = -1;
}

// The route is cleared (all customers are removed).
inline void CVRPSolution::clearRoute(int iRoute)
{
	const vector<int> &vecNodes = m_vecRoutes[iRoute].getNodes();
	// the first and last node in <vecNodes> is the depot. We should skip these
	// when we insert into the unassigned nodes.
	m_vecUnassigned.insert(m_vecUnassigned.end(), vecNodes.begin()+1,
		vecNodes.begin()+vecNodes.size()-1);
	// update <m_vecCust2Route>.
	int i;
	for (i = 1; i < (int) vecNodes.size()-1; ++i)
		m_vecCust2Route[vecNodes[i]] = -1;

	// Now perform the actual update of the route.
	m_dCost -= m_vecRoutes[iRoute].getCost();
	m_vecRoutes[iRoute].clear();
	m_dCost += m_vecRoutes[iRoute].getCost();
}


#endif 
