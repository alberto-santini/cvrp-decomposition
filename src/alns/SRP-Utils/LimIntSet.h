#ifndef LIMINTSET_H
#define LIMINTSET_H

#include "Utils.h"

// This class represents a set of integers that are not to different in an efficient way for testing 
// inclusion in the set. 
//
// The type of the elements stored is a template paramter, but it should be an integer type like:
// int, long, unsigned short, char, etc.
//
// The only method the class offers (currently) is <isInSet> which answers if a certain integer is in
// the set. 
//
// The internal representation is a vector of bool, one per element in the range iMin, ..., iMax where
// iMin and iMax are the smallest and largest elements in the set. Consequently the difference between 
// the largest and smallest element should not be too big.


template<class T>
class LimIntSet
{
public:
	template<class InputIterator>
	LimIntSet(InputIterator first, InputIterator last);
	inline bool isInSet(T id) const;

private:
	T m_min;
	T m_max;
	std::vector<bool> m_vecIsInSet;
	bool m_bEmpty;
};

template<class T>
bool LimIntSet<T>::isInSet(T id) const
{
	if (!m_bEmpty && m_min <= id && id <= m_max)
	{
		return m_vecIsInSet[id - m_min];
	}
	else
		return false;
}

template<class T>
template<class InputIterator>
LimIntSet<T>::LimIntSet(InputIterator first, InputIterator last)
{
	if (first == last)
		m_bEmpty = true;
	else 
	{
		m_bEmpty = false;

		m_min = *first; m_max = *first;
		InputIterator iter;
		for (iter = first ; iter != last; iter++)
		{
			if (*iter < m_min)
				m_min = *iter;
			if (*iter > m_max)
				m_max = *iter;
		}

		m_vecIsInSet.assign(m_max-m_min+1, false);
		for (iter = first ; iter != last; iter++)
		{
			m_vecIsInSet[*iter-m_min] = true;
		}
	}
}

#endif

