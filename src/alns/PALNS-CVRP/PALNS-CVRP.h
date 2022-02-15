#ifndef PALNSCVRP_H
#define PALNSCVRP_H

#include "CVRPSolution.h"
#include "CVRPInstance.h"

#include <vector>

#include "../PALNS/PALNS.h"

using namespace std;

class PALNS_CVRP
{
public:
	void go(CVRPSolution &sol, int iNThreads, int iNRetries, bool bHeader, const ParFileParser &pfp) const;
	void setResultPrefix(const string &strResultPrefix) { m_strResultPrefix = strResultPrefix; }
	
protected:
	string m_strResultPrefix;
};

#endif