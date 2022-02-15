#ifndef _HASH_NUMBER_GENERATOR_H
#define _HASH_NUMBER_GENERATOR_H

#include <vector>
#include "Utils.h"

class HashNumberGenerator
{
private:
	static std::vector<unsigned int> m_hashNumbers;
	
public:
	static const std::vector<unsigned int> &getHashNumbers(int n);
	static void clear();
};

class FixedHashNumberGenerator
{
private:
	static std::vector<unsigned int> m_hashNumbers;
	
public:
	// Returns a vector with at least <n> random numbers (the same random numbers every call)
	static const std::vector<unsigned int> &getHashNumbers(int n) 
	{
		if ((int) m_hashNumbers.size() < n)
			error("FixedHashNumberGenerator::getHashNumbers(...)","Not enough pregenerated numbers");
		return m_hashNumbers;
	}

	static void clear()
	{
		m_hashNumbers.clear();
	}
	
	static void init(int iSize)
	{
		int i;
		m_hashNumbers.resize(iSize);
		for (i=0; i < iSize; ++i)
			m_hashNumbers[i] = Random::getRandom(0, 2000000000);
	}

	static int size()
	{
		return (int) m_hashNumbers.size();
	}
};

#endif
