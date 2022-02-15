#ifndef SEQINSERTIONCACHE_H
#define SEQINSERTIONCACHE_H

#include <vector>
#include <float.h>

// To get RouteDelta class
#include "InsertionCache.h"

// From SRP-Utils
#include "NSMatrix.h"

using namespace std;

class CVRPInstance;
class CVRPSolution;

template<class T>
class Matrix;

class SeqInsertionCache
{
public:
	SeqInsertionCache(const CVRPInstance &instance, CVRPSolution &sol, const vector<vector<int> > &vecUnassignedSeqs);
	inline double getNthBestInsertionDelta(int iSeqId, int iNth) const;
	inline double getBestInsertion(int iSeqId, int &iRoute) const;
	void updateCache(CVRPSolution &sol, int iChangedRoute, int iSeqIdInserted);
	
private:
	void initCache(CVRPSolution &sol, const vector<vector<int> > &vecUnassignedSeqs);

	const CVRPInstance &m_instance;

	vector<vector<int> > m_vecUnassignedSeqs;
	vector<int> m_vecSequenceDemand;
	vector<double> m_vecSequenceDist;
	vector<int> m_vecSequenceServTime;
	vector<int> m_vecSequenceFirstCust;
	vector<int> m_vecSequenceLastCust;
	// m_seqIsUnassigned keeps track of the sequences that still are unassigned.
	vector<bool> m_seqIsUnassigned;

	// m_vecRouteDelta[i] contains a vector of RouteDelta objects for the sequence with 
	// id <i>. The first RouteDelta object in the vector is the one corresponding to the best possible
	// insertion, the next one corresponds to the second best insertion, etc. That is, each vector
	// is sorted according to increasing delta value.
	vector<vector<RouteDelta> > m_vecRouteDelta;

	// m_matPossible[i,j] is false if and only if it is impossible to insert the sequence with 
	// id <i> into route j (this would be because of the capacity).
	NSMatrix<bool> m_matPossible;
	// How many routes have we reserved space for in our data structures?
	int m_iNReservedRoutes;
};

// Get the delta cost of the nth best insertion (the best insertion is insertion number 0).
// Returns DBL_MAX if there is no nth best insertion (the customer can not be inserted in that
// many routes).
inline double SeqInsertionCache::getNthBestInsertionDelta(int iSeqId, int iNth) const
{
	const vector<RouteDelta> &vec = m_vecRouteDelta[iSeqId];
	if (iNth >= (int) vec.size())
		return DBL_MAX;
	else
		return vec[iNth].m_dDelta;
}

// Get the delta cost and the route id associated with the best insertion.
// Returns DBL_MAX if the customer cannot be inserted in any route.
inline double SeqInsertionCache::getBestInsertion(int iSeqId, int &iRoute) const
{
	const vector<RouteDelta> &vec = m_vecRouteDelta[iSeqId];
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
