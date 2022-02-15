#ifndef REGRET2ALG_H
#define REGRET2ALG_H

// Get NULL
#include <stdlib.h>

class CVRPInstance;
class CVRPSolution;
class InsertionCache;
class TSRandom;

template<class T>
class Matrix;

class Regret2Alg
{
public:
	Regret2Alg(const CVRPInstance &instance, TSRandom &randGen) 
		: m_instance(instance), m_randGen(randGen)   
	{
		m_dRandFactor = 1;
	}

	void go(CVRPSolution &sol, const Matrix<double> *pMatAltCost = NULL);
	void setRandFactor(double dRandFactor)		{ m_dRandFactor = dRandFactor; }

private:
	const CVRPInstance &m_instance;
	TSRandom &m_randGen;

	// a number greater than or equal to 1. 
	// 1 means no randomness and values above zero means more randomness 
	// To be more precise, the regret value is multiplied by RND(1, m_dRandFactor) everytime
	// it is evalutated (RND(a,b) returns a random number in the interval [a,b]).
	double m_dRandFactor;

};


#endif
