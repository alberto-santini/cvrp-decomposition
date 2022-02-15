// Sequence based regret-2 algorithms 

#ifndef SEQREGRET2ALG_H
#define SEQREGRET2ALG_H

// Get NULL
#include <stdlib.h>

class CVRPInstance;
class CVRPSolution;
class InsertionCache;
class TSRandom;

template<class T>
class Matrix;

class SeqRegret2Alg
{
public:
	SeqRegret2Alg(const CVRPInstance &instance, TSRandom &randGen) 
		: m_instance(instance), m_randGen(randGen)   
	{
		m_dRandFactor = 1;
	}

	void go(CVRPSolution &sol);
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
