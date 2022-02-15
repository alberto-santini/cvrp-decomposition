#ifndef INSERTIONCACHE_H
#define INSERTIONCACHE_H

#include <vector>
#include <float.h>

// From SRP-Utils
#include "NSMatrix.h"

using namespace std;

class CVRPInstance;
class CVRPSolution;

template<class T>
class Matrix;

struct RouteDelta
{
	int m_iRoute;
	double m_dDelta;

	RouteDelta(int iRoute, double dDelta)
		: m_iRoute(iRoute), m_dDelta(dDelta)
	{}

	bool operator < (const RouteDelta &other) const
	{
		return m_dDelta < other.m_dDelta;
	}
};

class InsertionCache
{
public:
	InsertionCache(const CVRPInstance &instance, CVRPSolution &sol, const vector<int> &vecUnassigned, const Matrix<double> *pMatAltCost = NULL);
	inline double getNthBestInsertionDelta(int iCustId, int iNth) const;
	inline double getBestInsertion(int iCustId, int &iRoute) const;
	void updateCache(CVRPSolution &sol, int iChangedRoute, const vector<int> &vecUnassigned);
	
private:
	void initCache(CVRPSolution &sol, const vector<int> &vecUnassigned);

	const CVRPInstance &m_instance;
	const Matrix<double> *m_pMatAltCost;
	// Cache members variables:
	// As the cache often will be used where only a subset of the customers are to be inserted
	// we renumber the customers to have more compact datastructures. The vector
	// m_vecOrigToCompressedId maps from the original id's to the "compressed" ids. That is,
	// m_vecOrigToCompressedId[i] gives the compressed id of customer i.
	vector<int> m_vecOrigToCompressedId;
	// m_vecRouteDelta[i] contains a vector of RouteDelta objects for the customer with compressed
	// id <i>. The first RouteDelta object in the vector is the one corresponding to the best possible
	// insertion, the next one corresponds to the second best insertion, etc. That is each vector
	// is sorted according to increasing delta value.
	vector<vector<RouteDelta> > m_vecRouteDelta;

	// m_matPossible[i,j] is false if and only if it is impossible to insert the customer with 
	// compressed id <i> into route j (this would be because of the capacity).
	NSMatrix<bool> m_matPossible;
	// How many routes have we reserved space for in our data structures?
	int m_iNReservedRoutes;
};

// Get the delta cost of the nth best insertion (the best insertion is insertion number 0).
// Returns DBL_MAX if there is no nth best insertion (the customer can not be inserted in that
// many routes).
inline double InsertionCache::getNthBestInsertionDelta(int iCustId, int iNth) const
{
	const vector<RouteDelta> &vec = m_vecRouteDelta[m_vecOrigToCompressedId[iCustId]];
	if (iNth >= (int) vec.size())
		return DBL_MAX;
	else
		return vec[iNth].m_dDelta;
}

// Get the delta cost and the route id associated with the best insertion.
// Returns DBL_MAX if the customer cannot be inserted in any route.
inline double InsertionCache::getBestInsertion(int iCustId, int &iRoute) const
{
	const vector<RouteDelta> &vec = m_vecRouteDelta[m_vecOrigToCompressedId[iCustId]];
	if (vec.empty())
		return DBL_MAX;
	else
	{
		const RouteDelta &rd = vec.front();
		iRoute = rd.m_iRoute;
		return rd.m_dDelta;
	}
}


#endif
