#ifndef VECTORUTILS_H
#define VECTORUTILS_H

#include <vector>
#include <assert.h>
#include "Utils.h"
#include "LimIntSet.h"

using namespace std;

int subset(const vector<int> &vecA, const vector<int> &vecB, int iMinElem, int iMaxElem);
bool compareVecVecInt(const vector<vector<int> > &vecVecInt1, const vector<vector<int> > &vecVecInt2);
void makeRandomIntVec(vector<int> &vec, int iNElems, int iMin, int iMax);

/***********************************************************************
 * template <class T>
 * T findMin(const vector<T> &vec)
 *
 * Find the minimum element in <vec>. It is assumed that the vector is non-empty.
 *
 * --- input parameters ---
 * vec			: The vector to search.
 * --- output parameters ---
 * --- return values ---
 * The minimum value.
 ***********************************************************************/

template <class T>
T findMin(const vector<T> &vec) 
{
	if (vec.empty())
		error("findMin(...)", "vector is empty");
	int i;
	T min = vec[0];
	for (i=1; i < (int) vec.size(); ++i)
	{
		if (vec[i] < min)
			min = vec[i];
	}
	return min;
}

template <class T>
void findMinMax(const vector<T> &vec, T &min, T &Max) 
{
	if (vec.empty())
		error("findMinMax(...)", "vector is empty");
	int i;
	min = vec[0];
	max = vec[0];
	for (i=1; i < (int) vec.size(); ++i)
	{
		if (vec[i] < min)
			min = vec[i];
		if (vec[i] > max)
			max = vec[i];
	}
}

/***********************************************************************
 * template <class T>
 * T findMin(const vector<T> &vec, int &iMinIdx)
 *
 * Find the minimum element in <vec>. It is assumed that the vector is non-empty.
 * This version of the method also returns the index of the minimum element
 *
 * --- input parameters ---
 * vec			: The vector to search.
 * --- output parameters ---
 * iMinIdx		: the index of the minimum element.
 * --- return values ---
 * The minimum value.
 ***********************************************************************/

template <class T>
T findMin(const vector<T> &vec, int &iMinIdx) 
{
	if (vec.empty())
		error("findMin(...)", "vector is empty");
	int i;
	// @@@ If T is a very large class, then it's not efficient to declare a new T object and
	// copy the T object every time a new better object has been found. In that case it would
	// be better to keep an index to the current best element.
	T min = vec[0];
	iMinIdx = 0;
	for (i=1; i < (int) vec.size(); ++i)
	{
		if (vec[i] < min)
		{
			min = vec[i];
			iMinIdx = i;
		}
	}
	return min;
}

/***********************************************************************
 * template <class T>
 * T vectorSum(const vector<T> &vec)
 *
 * Find the sum of the elements in <vec>. It is assumed that the vector is non-empty.
 *
 * --- input parameters ---
 * vec			: The vector to calculate sum of.
 * --- output parameters ---
 * --- return values ---
 * The sum of the elements in <vec>.
 ***********************************************************************/

template <class T>
T vectorSum(const vector<T> &vec) 
{
	if (vec.empty())
	{
		error("vectorSum(...)", "vector is empty");
		return vec.front();
	}
	else
	{
		int i;
		T sum = vec[0];
		for (i=1; i < (int) vec.size(); i++)
			sum += vec[i];
		return sum;
	}
}

/***********************************************************************
 * template <class T>
 * void vectMult(vector<T> &vec, T multiplier) 
 *
 * multiplies each element of vec by <multiplier>
 *
 * --- input parameters ---
 * vec			: The vector to operate on.
 * multiplier	: The value to multiply with
 * --- output parameters ---
 * --- return values ---
 ***********************************************************************/

template <class T>
void vecMult(vector<T> &vec, T multiplier) 
{
	int i;
	int iSize = (int) vec.size();
	for (i=0; i < iSize; ++i)
		vec[i] *= multiplier;
}

/***********************************************************************
 * template <class T>
 * void makeRange(vector<T> &vec, T min, T max)
 *
 * Constructs a vector with the elements (min, min+1, min+2, ..., max-1, max)
 * T should be a type where "<=" and "++" have reasonable meanings (int is the prime
 * example).
 *
 * --- input parameters ---
 * min, max		: the range of numbers to generate.
 * --- output parameters ---
 * vec			: the vector with the desired range
 * --- return values ---
 ***********************************************************************/

// Constructs a vector with the elements (min, min+1, min+2, ..., max-1, max)
template <class T>
void makeRange(vector<T> &vec, T min, T max)
{
	vec.clear();
	T i;
	for (i = min; i <= max; i++)
		vec.push_back(i);	
}

template <class T>
void vectorRemoveElement(std::vector<T> &vec, const T &element)
{
	typename std::vector<T>::iterator start = vec.begin();

	typename std::vector<T>::iterator last = remove(start, vec.end(), element);
	vec.resize(last-start);
}

// Returns true if the element was deleted. Returns false if no such element was found.
template <class T>
bool vectorRemoveSingleElement(std::vector<T> &vec, const T &element)
{
	typename std::vector<T>::iterator it = find(vec.begin(), vec.end(), element);
	if (it != vec.end())
	{
		vec.erase(it);
		return true;
	}
	else
		return false;
}

/***********************************************************************
 * template <class T>
 * void unorderedVectorRemoveAt(std::vector<T> &vec, int iPos)
 *
 * Removes an element from a vector. The ordering of the vector 
 * can be changed by the operation.
 *
 * --- input parameters ---
 * vec			: The vector to remove from
 * pos			: The index of the element to remove
 * --- return values ---
 * none
 ***********************************************************************/

template <class T>
void unorderedVectorRemoveAt(std::vector<T> &vec, int iPos)
{
	vec[iPos] = vec.back();
	vec.pop_back();
}

/***********************************************************************
 * template <class T>
 * void unorderedVectorRemoveElem(std::vector<T> &vec, const T &elem)
 *
 * Removes the first occurence of an element from a vector. The ordering of the vector 
 * can be changed by the operation. The method should be faster that 
 * <vectorRemoveElement> on average, but it only removes one element and it does
 * not keep the ordering.
 *
 * --- input parameters ---
 * vec			: The vector to remove from
 * elem			: The element to remove
 * --- return values ---
 * none
 ***********************************************************************/

template <class T>
void unorderedVectorRemoveElem(std::vector<T> &vec, const T &elem)
{
	typename std::vector<T>::iterator start = vec.begin();
	typename std::vector<T>::iterator end = vec.end();
	typename std::vector<T>::iterator iter;
	iter = find(start, end, elem);
	if (iter != end)
	{
		typename vector<T>::iterator::difference_type idx = iter-start;
		vec[idx] = vec.back();
		vec.pop_back();
	}
}


/***********************************************************************
 * template<class IntType>
 * bool vectorsIdentical(vector<IntType> vec1, vector<IntType> vec2, IntType iMinElem, IntType iMaxElem)
 *
 * Examines if two vectors contains the same elements. The type of the elements
 * stored in the vectors should be integer-like (int, shor, long, char, etc.).
 * The run time complexity is O(max{vec1.size(), iMaxElem-iMinElem})
 *
 * --- input parameters ---
 * vec1, vec2           : The two vectors we are comparing.
 * iMinElem             : The smallest element in either vector. Warning: We do not check if
 *                        this really is the smallest element. The method might (will) crash
 *                        if the smallest element really is smaller than indicated.
 * iMaxElem             : The largest element in either vector.  Warning: We do not check if
 *                        this really is the largest element. The method might (will) crash
 *                        if the largest element really is larger than indicated.
 * --- output parameters ---
 * --- return value ---
 * true => The two vectors were identical.
 ***********************************************************************/

template<class IntType>
bool vectorsIdentical(const vector<IntType> &vec1, const vector<IntType> &vec2, IntType iMinElem, IntType iMaxElem)
{
    if (vec1.size() != vec2.size())
        return false;
    
    if (iMinElem > iMaxElem)
        error("vectorsIdentical(...)", "iMinElem > iMaxElem");
    
    vector<unsigned int> vecElemInVec1Count(iMaxElem-iMinElem+1, 0);
    vector<unsigned int> vecElemInVec2Count(iMaxElem-iMinElem+1, 0);
    typename vector<IntType>::const_iterator cIter;
    for (cIter = vec1.begin(); cIter != vec1.end(); cIter++)
	{
		assert(iMinElem <= (*cIter) && (*cIter) <= iMaxElem);
        vecElemInVec1Count[(*cIter)-iMinElem]++;
	}

    for (cIter = vec2.begin(); cIter != vec2.end(); cIter++)
	{
		assert(iMinElem <= (*cIter) && (*cIter) <= iMaxElem);
        vecElemInVec2Count[(*cIter)-iMinElem]++;
	}

    for (unsigned int i=0; i < vecElemInVec1Count.size(); i++)
    {
        if (vecElemInVec1Count[i] != vecElemInVec2Count[i])
            return false;
    }
    return true;
}

//To do: make a version of the above method that do not need max and min elem (can be found in linear time).

/***********************************************************************
 * template<class SortType>
 * bool vectorsIdenticalSortType(vector<SortType> vec1, vector<SortType> vec2)
 *
 * Examines if two vectors contains the same elements. The type of the elements
 * stored in the vectors should be sortable (that is, they should have a meaningful
 * "<" operator. 
 * The method makes a copy of the two vectors, sort these copies and finally
 * it compares the sorted vectors element by element.
 * The time complexity is O(n log n) where n=max{|vec1|,|vec2|}
 * Thus if the two vectors contain integers then <vectorsIdentical> is faster.
 *
 * --- input parameters ---
 * --- output parameters ---
 * --- return value ---
 * true => The two vectors were identical.
 ***********************************************************************/

template<class SortType>
bool vectorsIdentical(const vector<SortType> &vec1, const vector<SortType> &vec2)
{
    if (vec1.size() != vec2.size())
        return false;
    
	vector<SortType> vec1Sorted(vec1), vec2Sorted(vec2);
	sort(vec1Sorted.begin(), vec1Sorted.end());
	sort(vec2Sorted.begin(), vec2Sorted.end());

	int i, iSize;
	iSize = (int) vec1.size();

	for (i=0; i < iSize; ++i)
	{
		if (vec1Sorted[i] != vec2Sorted[i])
			return false;
	}

    return true;
}

/***********************************************************************
 * template<class T>
 * void reverseVector(vector<T> &vec);
 *
 * reverses the vector <vec>. The method makes O(n) calls to T's copy constructor, where n is the length
 * of the vector.
 * 
 * --- input parameters ---
 * vec			: The vector to be reversed. After calling the method the reversed vector is returned through this
 *				  parameter.
 * --- output parameters ---
 * --- return value ---
 ***********************************************************************/

template<class T>
void reverseVector(vector<T> &vec)
{
	T temp;
	int i, iLength = (int) vec.size();
	int iMax = iLength/2;
	for (i = 0; i < iMax; i++)
	{
		temp = vec[i];
		vec[i] = vec[iLength-1-i];
		vec[iLength-1-i] = temp;
	}
}

/***********************************************************************
 * template <class T>
 * void vecToString ()
 *
 *
 * --- input parameters ---
 * --- return values ---
 ***********************************************************************/

template <class T>
void vecToString (	string &str, const vector<T> &vec, string (*f) (T) , 
					const string &strDelim = " ", const string &strEnd = "")
{
    int i;
    for (i=0; i < (int) vec.size(); i++)
	{
        str += f(vec[i]);
		if (i < (int) vec.size() - 1)
			str += strDelim;
	}
	str += strEnd;
}

/***********************************************************************
 * template <class T>
 * ostream &outputVector (ostream &os, const std::vector<T> &vec, const string &strDelim = " ", const string &strEnd = "")
 *
 * Writes the content of a vector to the outstream <os>
 *
 * --- input parameters ---
 * vec			: The vector to output
 * os			: The outstream to write on.
 * --- return values ---
 * the modified outstream.
 ***********************************************************************/

template <class T>
ostream &outputVector (ostream &os, const std::vector<T> &vec, const string &strDelim = " ", const string &strEnd = "")
{
    int i;
    for (i=0; i < (int) vec.size(); i++)
	{
        os << vec[i];
		if (i < (int) vec.size()-1)
			os << strDelim;
	}
	os << strEnd;
    return os;
}

/***********************************************************************
 * template <class T>
 * ostream &outputVector (ostream &os, const std::vector<T> &vec, const string &strDelim, 
 *					      const string &strStart, const string &strEnd)
 *
 * Writes the content of a vector to the outstream <os>
 *
 * --- input parameters ---
 * vec			: The vector to output
 * os			: The outstream to write on.
 * strDelim		: string used to separate elements in the vector.
 * strStart		: string printed before starting to print elements from the vector.
 * strEnd		: string printed after printing all elements from the vector.
 * --- return values ---
 * the modified outstream.
 ***********************************************************************/

template <class T>
ostream &outputVector (ostream &os, const std::vector<T> &vec, const string &strDelim, 
					   const string &strStart, const string &strEnd)
{
    int i;
	os << strStart;
    for (i=0; i < (int) vec.size(); i++)
	{
        os << vec[i];
		if (i < (int) vec.size()-1)
			os << strDelim;
	}
	os << strEnd;
    return os;
}

/***********************************************************************
 * template <class T>
 * ostream &outputVector (ostream &os, const std::vector<T> &vec)
 *
 * Writes the content of a vector to the outstream <os>
 *
 * --- input parameters ---
 * vec			: The vector to output
 * os			: The outstream to write on.
 * --- return values ---
 * the modified outstream.
 ***********************************************************************/

template <class T>
ostream &outputVector (ostream &os, const std::vector<T> &vec, char delim, const string &strEnd = "")
{
	int i;
    for (i=0; i < (int) vec.size(); i++)
	{
        os << vec[i];
		if (i < (int) vec.size()-1)
			os << delim;
	}
	os << strEnd;
    return os;
}

/***********************************************************************
 * template <class T>
 * ostream &outputVecNum (ostream &os, const std::vector<T> &vec, const string &strDelim = ", ", const string &strEnd = "", const string &strIdxDelim = ": ")
 *
 * Writes the content of a vector to the outstream <os>. Each element is prefixed with its index in the vector.
 * that is if we pass the vector (4, 1, 6, 9) the method will output:
 * "0: 4, 1: 1, 2: 6, 3: 9"
 * By changing the parameters the method can also output:
 * 0: 4
 * 1: 1
 * 2: 6
 * 3: 9
 *
 * --- input parameters ---
 * vec			: The vector to output
 * os			: The outstream to write on.
 * strDelim		: The string to output after each element.
 * strEnd		: The string to output after the entire vector.
 * strIdxDelim	: The string to output after each index.
 * --- return values ---
 * the modified outstream.
 ***********************************************************************/

template <class T>
ostream &outputVecNum (ostream &os, const std::vector<T> &vec, const string &strDelim = ", ", const string &strEnd = "", const string &strIdxDelim = ": ")
{
	int i;
    for (i=0; i < (int) vec.size(); i++)
	{
		os << i << strIdxDelim << vec[i];
		if (i < (int) vec.size()-1)
			os << strDelim;
	}
	os << strEnd;
    return os;
}

/***********************************************************************
 * template <class T>
 * int subset(const vector<T> &vecA, const vector<T> &vecB)
 *
 * Determines if vecA is a subset of vecB or if the opposite is true. 
 * The method asumes that sensible "<", "==" and "!=" operators are defined on T.
 * The method takes O(n log n) where n=max(vecA.size, vecB.size()).
 *
 * --- input parameters ---
 * vecA           :
 * vecB           : The two vectors to compare.
 * --- output parameters ---
 * --- return value ---
 * 1        : vecB is a subset of vecA ("vecA-vecB > 0")
 * 0        : vecA and VecB are identical (when seen as sets) 
 * -99      : neither vecA or vecB are subsets of each other
 * -1       : vecA is a subset of vecB ("vecA-vecB < 0")
 ***********************************************************************/
template <class T>
int subset(const vector<T> &vecA, const vector<T> &vecB)
{
	// Some special cases when one or both sets are empty:
	if (vecA.empty() && vecB.empty())
		return 0;
	if (vecA.empty() && !vecB.empty())
		return -1;
	if (vecB.empty() && !vecA.empty())
		return 1;
	// We now know that both sets are non-empty.

	vector<T> vec1(vecA), vec2(vecB);
	sort(vec1.begin(), vec1.end());
	sort(vec2.begin(), vec2.end());
	
	bool bASubsetB = false;
	bool bBSubsetA = false;
	bool bIdentical = true;
	bool bStop = false;
	typename vector<T>::const_iterator cIter1 = vec1.begin(), cIter2 = vec2.begin();
	typename vector<T>::const_iterator cIter1Next, cIter2Next;

	while (!bStop)
	{
		cIter1Next = cIter1; cIter2Next = cIter2;
		cIter1Next++; cIter2Next++;
		// Skip past duplicates in both vectors.
		while (cIter1Next != vec1.end() && *cIter1 == *cIter1Next)
		{
			cIter1++; cIter1Next++;
		}
		while (cIter2Next != vec2.end() && *cIter2 == *cIter2Next)
		{
			cIter2++; cIter2Next++;
		}
		if (*cIter1 != *cIter2)
		{
			bIdentical = false;
			if (*cIter1 < *cIter2)
			{
				// vec1 seems to contain more elements than vec2.
				bBSubsetA = true;
				cIter1++;
			}
			else
			{
				// vec2 seems to contain more elements than vec1.
				bASubsetB = true;
				cIter2++;
			}
		}
		else
		{
			cIter1++; cIter2++;
		}
		if (cIter1Next == vec1.end() || cIter2Next == vec2.end() || (bBSubsetA && bASubsetB))
			bStop = true;
	}
	// Let's figure if A is a subset of B or vice versa
	if (bIdentical && cIter1Next != vec1.end())
		// vec2 is a subset of vec1 (they have been identical untill we met the end of vec2, 
		// but vec1 still has some elements).
		return 1;
	if (bIdentical && cIter2Next != vec2.end())
		// vec1 is a subset of vec2 (they have been identical untill we met the end of vec1, 
		// but vec2 still has some elements).
		return -1;
	if (bASubsetB && cIter1Next != vec1.end())
		// It seemed like A was a subset of B, but vec1 still contains some elements and we reached the end of
		// vec2. Consequently neither vec is a subset of the other.
		return -99;
	if (bBSubsetA && cIter2Next != vec2.end())
		// It seemed like B was a subset of A, but vec2 still contains some elements and we reached the end of
		// vec1. Consequently neither vec is a subset of the other.
		return -99;
	if (bIdentical)	
		// If we get to this point, then we have reached the end of both vectors and they are identical.
		return 0;
	if (bASubsetB && bBSubsetA)
		// Neither is a subset of the other;
		return -99;
	if (bASubsetB)
		return -1;
	if (bBSubsetA)
		return 1;
	// We should never reach this part.
	error("int subset(...)", "Reached part of code that should be unreachable");
	return 10000;
}

/***********************************************************************
 * template<class T>
 * void makeIsInSetVec(const vector<T> &vecSet, vector<bool> &vecIsInSet, int iMaxId)
 *
 * Given a set of integers (vecSet) in the range [0, iMaxId-1] the method creates a vector 
 * of bools (vecIsInSet) with <iMaxId> elements. vecIsInSet[x] is true if and 
 * only if x is in vecSet.
 *
 * The type T must have integer type (eg. int, char, long, unsigned int, etc.)
 * 
 * --- input parameters ---
 * vecSet			: The set of integers we wish to process
 * iMaxId			: No element in the set can have value equal to or higher than <iMaxId>
 * --- output parameters ---
 * vecIsInSet		: vecIsInSet[x] is true if and only if x is in <vecSet>
 * --- return value ---
 ***********************************************************************/
template<class T>
void makeIsInSetVec(const vector<T> &vecSet, vector<bool> &vecIsInSet, int iMaxId)
{
	vecIsInSet.clear();
	vecIsInSet.resize(iMaxId, false);
	int i;
	for (i=0; i < (int) vecSet.size(); i++)
	{
		assert(0 <= vecSet[i] && vecSet[i] < iMaxId);
		vecIsInSet[vecSet[i]] = true;
	}
}

/***********************************************************************
 * template<class T>
 * void addToIsInSetVec(const vector<T> &vecNewSet, vector<bool> &vecIsInSet)
 *
 * Adds to a "IsInSet" bool vector (see <makeIsInSetVec>).
 *
 * The type T must have integer type (eg. int, char, long, unsigned int, etc.)
 * 
 * --- input parameters ---
 * vecSet			: The set of integers we wish to add. 
 * --- output parameters ---
 * vecIsInSet		: vecIsInSet[x] is true if and only if x is in <vecNewSet> 
 *					  or in the original set.
 * --- return value ---
 ***********************************************************************/
template<class T>
void addToIsInSetVec(const vector<T> &vecNewSet, vector<bool> &vecIsInSet)
{
	int i;
	for (i=0; i < (int) vecNewSet.size(); ++i)
	{
		assert(0 <= vecNewSet[i] && vecNewSet[i] < (int) vecIsInSet.size());
		vecIsInSet[vecNewSet[i]] = true;
	}
}

inline void boolVecIntersect(vector<bool> &vecInOut, const vector<bool> &vecOther)
{
	assert(vecInOut.size() == vecOther.size());
	vector<bool>::size_type i;
	for (i=0; i < vecInOut.size(); ++i)
		vecInOut[i] = vecInOut[i] && vecOther[i];
}

inline void boolVecUnion(vector<bool> &vecInOut, const vector<bool> &vecOther)
{
	assert(vecInOut.size() == vecOther.size());
	vector<bool>::size_type i;
	for (i=0; i < vecInOut.size(); ++i)
		vecInOut[i] = vecInOut[i] || vecOther[i];
}


/***********************************************************************
 * template<class T>
 * T *vectorToArray(const vector<T> &vec)
 *
 * Copies the elements in the vector <vec> to a classic array. A pointer to
 * the array is returned. It is up to the caller to delete the array
 * (with delete []).
 * --- input parameters ---
 * vec			: The vector to copy.
 * --- output parameters ---
 * --- return value ---
 * A pointer to the created array.
 ***********************************************************************/
template<class T>
T *vectorToArray(const vector<T> &vec)
{
	T *pArray = new T[vec.size()];
	int i;
	for (i=0; i < (int) vec.size(); i++)
		pArray[i] = vec[i];
	return pArray;
}

/***********************************************************************
 * template<class T>
 * void permuteRandom(vector<T> &vec)
 *
 * Permute <vec> in a random fashion (that is, elements are reordered).
 * The method is not as efficient as it could be, especially if the elements
 * stored in the vector are very large.
 *
 * --- input parameters ---
 * vec			: The vector to permute-
 * --- output parameters ---
 * --- return value ---
 ***********************************************************************/
template<class T>
void permuteRandom(vector<T> &vec)
{
	vector<T> vecTemp;
	vecTemp.reserve(vec.size());
	vector<int> vecAvailIdx;
	makeRange(vecAvailIdx, 0, (int) vec.size()-1);
	while (!vecAvailIdx.empty())
	{
		int idx1 = Random::getRandom(0, (int) vecAvailIdx.size()-1);
		int idx2 = vecAvailIdx[idx1];
		vecTemp.push_back(vec[idx2]);
		unorderedVectorRemoveAt(vecAvailIdx, idx1);
	}
	// Would it be faster to use vec.swap(vecTemp) below???
	vec = vecTemp;
}

/***********************************************************************
 * void unorderedFilterVec(vector<int> &vec, const vector<bool> &vecForbiddenElements)
 *
 * Remove certain elements from the vector <vec>.
 * vecForbiddenElements[i] is true if element <i> should be removed from <vec>.
 * It is assumed that all elements in <vec> should is in the interval
 * [0, vecForbiddenElements.size()]. WARNING: This is only checked in debug mode.
 *
 * --- input parameters ---
 * vec						: The vector to filter (is modified by the method)
 * vecForbiddenElements		: Specifies the forbidden elements (see above)
 * --- output parameters ---
 * --- return values ---
 ***********************************************************************/

inline void unorderedFilterVec(vector<int> &vec, const vector<bool> &vecForbiddenElements)
{
	int i = 0;
	while (i < (int) vec.size())
	{
		assert(0 <= vec[i] && vec[i] < (int) vecForbiddenElements.size());
		if (vecForbiddenElements[vec[i]])
			unorderedVectorRemoveAt(vec, i);
		else 
			i++;
	}
}

/***********************************************************************
 * void unorderedFilterVec(vector<int> &vec, const vector<bool> &vecAllowedElements)
 *
 * Remove certain elements from the vector <vec>.
 * vecAllowedElements[i] is false if element <i> should be removed from <vec>.
 * It is assumed that all elements in <vec> should is in the interval
 * [0, vecAllowedElements.size()]. WARNING: This is only checked in debug mode.
 *
 * --- input parameters ---
 * vec						: The vector to filter (is modified by the method)
 * vecAllowedElements		: Specifies the allowed elements (see above)
 * --- output parameters ---
 * --- return values ---
 ***********************************************************************/

inline void unorderedPosFilterVec(vector<int> &vec, const vector<bool> &vecAllowedElements)
{
	int i = 0;
	while (i < (int) vec.size())
	{
		assert(0 <= vec[i] && vec[i] < (int) vecAllowedElements.size());
		if (vecAllowedElements[vec[i]])
			i++;
		else
			unorderedVectorRemoveAt(vec, i);
	}
}

/***********************************************************************
 * void unorderedRemoveDuplicates(vector<T> &vec)
 *
 * Remove duplicate elements from vec. The sequence of vec will be 
 * rearanged with this call. 
 * The time complexity is O(n log n) where n is the number of elements.
 *
 * --- input parameters ---
 * vec						: The vector to remove duplicates from
 * --- output parameters ---
 * --- return values ---
 ***********************************************************************/

template<class T> 
void unorderedRemoveDuplicates(vector<T> &vec)
{
	sort(vec.begin(), vec.end());
	int i=1;
	while (i < (int) vec.size())
	{
		if (vec[i-1] == vec[i])
			unorderedVectorRemoveAt(vec, i);
		else
			++i;
	}
}

/***********************************************************************
 * template<class T> 
 * void insertOrderedVec(vector<T> &vec, const &T elem)
 *
 * Inserts <elem> into the ordered vector <vec>. Uses binary search
 * to find the right insertion position. The time complexity is O(n)
 * where n is the number of elements in the vector (the insertion itself
 * is what causes the linear time complexity, not the search operation).
 * The class T must have an operator "<" that defines a total ordering.
 * 
 * --- input parameters ---
 * vec						: The sorted vector.
 * --- output parameters ---
 * --- return values ---
 ***********************************************************************/

template<class T> 
void insertOrderedVec(vector<T> &vec, const T &elem)
{
	int iMin = 0, iMax = (int) vec.size();
	while (iMin < iMax)
	{
		int iMiddle = iMin + (iMax-iMin)/2;
		if (vec[iMiddle] < elem)
			iMin = iMiddle+1;
		else
			iMax = iMiddle;
	}
	vec.insert(vec.begin()+iMax, elem);
}

/***********************************************************************
 * template<class T> 
 * void insertUnique(vector<T> &vec, const &T elem)
 *
 * Inserts <elem> into the vector <vec>. The element is only inserted
 * if it does not already exist in the vector. The method uses the 
 * "==" operator to determine if two elements are identical. The T class
 * should therefore provide an implementation of this operator.
 *
 * --- input parameters ---
 * vec						: The vector to insert into.
 * elem						: The element we wish to insert.
 * --- output parameters ---
 * --- return values ---
 ***********************************************************************/
template<class T>
void insertUnique(vector<T> &vec, const T &elem)
{
	int i, iSize = (int) vec.size();
	bool bFound = false;
	for (i=0; i < iSize; ++i)
	{
		if (vec[i] == elem)
		{
			bFound = true;
			break;
		}
	}
	if (!bFound)
		vec.push_back(elem);
}

/***********************************************************************
 * template<class IntType>
 * void complementSet(const vector<IntType> &vecS, vector<IntType> &vecSBar, IntType iMinElem, IntType iMaxElem)
 *
 * Given a set S and V the method computes V \setminus S. The set V is specified by
 * its smallest and largest elements. It is assumed that S \subseteq V.
 *
 * The running time is O(|V| + |vecS|).
 *
 * --- input parameters ---
 * vecS					: The set we want to complement.
 * iMinElem             : The smallest element in the base set V. Remark: this MUST be less than or equal to the smallest
 *						  element in S as well.
 * iMaxElem             : The largest element in the base set V. Remark: this MUST be larger than or equal to the largest
 *						  element in S as well.
 * --- output parameters ---
 * vecSBar				: The complement of vecS in V. 
 * --- return value ---
 ***********************************************************************/

template<class IntType>
void complementSet(const vector<IntType> &vecS, vector<IntType> &vecSBar, IntType iMinElem, IntType iMaxElem)
{
	int i;
	vecSBar.clear();
	vector<bool> vecIsInS;
	makeIsInSetVec(vecS, vecIsInS, iMaxElem-iMinElem+1);
	for (i = 0; i <= iMaxElem-iMinElem; ++i)
		if (!vecIsInS[i])
			vecSBar.push_back(i+iMinElem);
}


// Calculates vecRes = S1 \setminus S2;
// 
// time complexity:
// O(k + |S1|) where k = max{i : i in S2}-min{i : i in S2}.
template<class IntType>
void setMinus(const vector<IntType> &vecS1, const vector<IntType> &vecS2, vector<IntType> &vecRes)
{
	int i;
	vecRes.clear();
	LimIntSet<IntType> lisS2(vecS2.begin(), vecS2.end());

	for (i=0; i < (int) vecS1.size(); ++i)
		if (!lisS2.isInSet(vecS1[i]))
			vecRes.push_back(vecS1[i]);
}

/***********************************************************************
 * template<class T> 
 * double dist(const vector<T> &vec1, const vector<T> &vec2)
 *
 * Calculates the euclidean distance between the two vectors vec1 and vec2
 *
 * --- input parameters ---
 * vec1			
 * vec2
 * 
 * --- output parameters ---
 * --- return value ---
 * The euclidean distance between the two vectors.
 ***********************************************************************/

template<class T> 
double dist(const vector<T> &vec1, const vector<T> &vec2)
{
	int i, iSize = (int) vec1.size();
	if ((int) vec2.size() != iSize)
		error("dist(...)", "vector size mismatch");
	double dSum = 0;
	T temp;
	for (i=0; i < iSize; ++i)
	{
		temp = vec1[i]-vec2[i];
		dSum += temp*temp;
	}
	return sqrt(dSum);
}

/***********************************************************************
 * template<class T> 
 * void removeElements(const vector<T> &vec, const vector<int> &vecIdsToRemoveConst)
 *
 * removes elements from <vec>, while keeping the order of the vector. The elements
 * to remove are indicated by their indeces in <vecIdsToRemove>.
 * The method has an assymptotic running time of 
 * O(max{ |vec|, |vecIdsToRemove|*log(|vecIdsToRemove|) } )
 *
 * which is better than a naive implentation which would require:
 * O(|vec|^2).
 *
 * It would be possible to do the removal in O(|vec|) by making a boolean vector
 * that indicates membership in |vecIdsToRemove| and then through |vec|
 * and check for membership in vecIdsToRemove for every element in <vec>.
 * (or simply use algorithm::remove_if with a predicate using the boolean vector).
 *
 * In practice, when vecIdsToRemove is relatively small, I expect the 
 * current implementation to be more efficient because of a low number of 
 * comparisons.
 *
 * --- input parameters ---
 * vec				: The vector to remove from.
 * vecIdsToRemove	: The indeces of the elements to remove. It is assumed that this vector doesn't contain duplicates.
 * 
 * --- output parameters ---
 * --- return value ---
 * 
 ***********************************************************************/

template<class T> 
void removeElements(vector<T> &vec, const vector<int> &vecIdsToRemoveConst)
{
	vector<int> vecIdsToRemove(vecIdsToRemoveConst);
	sort(vecIdsToRemove.begin(), vecIdsToRemove.end());
	int i = 0;
	int j;
	int iNRemoved = 0;
	while (i < (int) vecIdsToRemove.size())
	{
		for (j=i+1; j < (int) vecIdsToRemove.size() && (vecIdsToRemove[j] == vecIdsToRemove[j-1]+1); ++j);
		// j now points to the first element in <vecIdsToRemove> that isn't following consequtively from vecIdsToRemove[i].
		// That is, if vecIdsToRemove = [ 2, 4, 5, 6, 8, 12] and i=1 then we would have j = 4.
		if (j < (int) vecIdsToRemove.size())
			copy(vec.begin()+vecIdsToRemove[j-1]+1, vec.begin()+vecIdsToRemove[j], vec.begin()+vecIdsToRemove[i]-iNRemoved);
		else
			copy(vec.begin()+vecIdsToRemove[j-1]+1, vec.end(), vec.begin()+vecIdsToRemove[i]-iNRemoved);
		iNRemoved += vecIdsToRemove[j-1]+1-vecIdsToRemove[i];
		i = j;
	}
	vec.resize(vec.size()-iNRemoved);
}

/***********************************************************************
 * template<class T> 
 * void removeElementsByValue(vector<T> &vec, const vector<T> &vecElementsToRemove)
 *
 * removes elements from <vec>, while keeping the order of the vector. The elements
 * to remove are indicated by their values in  <vecElementsToRemove>.
 * The method has an assymptotic running time of 
 * O(max{ |vec| + |vecElementsToRemove| + (max(vecElementsToRemove)-min(vecElementsToRemove)) } 
 *
 * T is assumed to be an "integer-like" type.
 *
 * --- input parameters ---
 * vec					: The vector to remove from.
 * vecElementsToRemove	: The elements to remove. It is okay if the vector contains duplicates.
 * 
 * --- output parameters ---
 * --- return value ---
 * 
 ***********************************************************************/
template<class T> 
void removeElementsByValue(vector<T> &vec, const vector<T> &vecElementsToRemove)
{
	LimIntSet<T> lis(vecElementsToRemove.begin(), vecElementsToRemove.end());
	int i;
	vector<T> vecCopy;
	vecCopy.reserve(vec.size());
	for (i=0; i < (int) vec.size(); ++i)
	{
		if (!lis.isInSet(vec[i]))
			vecCopy.push_back(vec[i]);
	}
	vec = vecCopy;
}

/***********************************************************************
 * template<class TInt> 
 * void removeSortedAndDecrease(vector<TInt> &vec, const vector<int> &vecIdsToRemoveConst)
 *
 * removes elements from <vec>, while keeping the order of the vector.
 * IT IS ASSUMED THAT VEC IS SORTED IN AN INCREASING FASHION. The elements
 * to remove are indicated by their indeces in <vecIdsToRemove>.
 * The method has an assymptotic running time of 
 * O(max{ |vec|, |vecIdsToRemove|*log(|vecIdsToRemove|) } )
 * (similar to <removeElements>)
 * 
 * The method also modifies the elements of <vec> that remains. If k elements have
 * been removed in front of an element then the element is decreased by k. 
 * This can be useful if <vec> stores a set of indexes into a vector and
 * some elements from the original vector is removed. It was used in the
 * BCP-2idx project. It's not entirely clear if this method really will be
 * useful for other purposes, but now it is here.
 *
 * Example:
 * input: vec = (3, 7, 9, 15, 18, 24), vecIdsToRemoveConst = (1, 4)
 * vec after call (the elements 7 and 18 are removed):
 * (3, 8, 14, 22)
 * 
 * --- input parameters ---
 * vec				: The vector to remove from. MUST BE SORTED IN INCREASING FAHSION.
 * vecIdsToRemove	: The indeces of the elements to remove. It is assumed that this vector doesn't contain duplicates.
 * 
 * --- output parameters ---
 * --- return value ---
 * 
 ***********************************************************************/
template<class TInt> 
void removeSortedAndDecrease(vector<TInt> &vec, const vector<int> &vecIdsToRemoveConst)
{
	vector<int> vecIdsToRemove(vecIdsToRemoveConst);
	sort(vecIdsToRemove.begin(), vecIdsToRemove.end());
	int iSrcIdx = 0, iDestIdx = 0;
	int i;
	int iNRemoved = 0;
	for (i=0; i < (int) vecIdsToRemove.size(); ++i)
	{
		int iNextRemIdx = vecIdsToRemove[i];
		while (iSrcIdx < iNextRemIdx)
		{
			vec[iDestIdx] = vec[iSrcIdx]-iNRemoved;
			++iDestIdx;
			++iSrcIdx;
		}
		++iNRemoved;
		// Increment src idx so we skip the element to be deleted.
		++iSrcIdx;
	}
	// Copy the rest of the elements (the ones after the last id to remove):
	for ( ; iSrcIdx < (int) vec.size(); ++iSrcIdx)
	{
		vec[iDestIdx] = vec[iSrcIdx]-iNRemoved;
		iDestIdx++;
	}
	vec.resize(vec.size()-iNRemoved);
}

/***********************************************************************
 * template<class TInt> 
 * void sortedVecDecrease(vector<TInt> &vec, const vector<int> &vecRemovedConst)
 *
 * <vec> is a vector of indices, representing some subset of a larger, ordered set. 
 * When elements are removed from the larger set, the indices in <vec> must be updated
 * because indices of the elements in the larger set are shifted to "fill up holes".
 * The indices of elements that were removed are given by <vecRemovedConst>.
 * - It is assummed that the intersection between vec and vecRemovedConst is empty.
 * - It is assumed that <vec> is sorted in an increasing fashion. 
 *
 * The method has an assymptotic running time of 
 * O(max{ |vec|, |vecRemovedConst|*log(|vecRemovedConst|) } )
 * 
 * Method is used in the BCP-2idx project.
 *
 * Example:
 * input: vec = (3, 7, 9, 15, 18, 24), vecRemovedConst = (5, 16)
 * vec after call 
 * (3, 6, 8, 14, 16, 22)
 * 
 * --- input parameters ---
 * vec				: The vector of indeces. MUST BE SORTED IN INCREASING FAHSION.
 * vecIdsToRemove	: The indeces of the elements that was removed. It is assumed that this vector doesn't contain duplicates.
 *					  It is assumed that the intersection between vec and vecRemovedConst is empty.
 * 
 * --- output parameters ---
 * --- return value ---
 * 
 ***********************************************************************/

template<class TInt> 
void sortedVecDecrease(vector<TInt> &vec, const vector<int> &vecRemovedConst)
{
	vector<int> vecRemoved(vecRemovedConst);
	sort(vecRemoved.begin(), vecRemoved.end());
	int i;
	int iRemovedIdx = 0;
	int iNRemovedBefore = 0;
	for (i=0; i < (int) vec.size() && iRemovedIdx < (int) vecRemoved.size(); ++i)
	{
		while (iRemovedIdx < (int) vecRemoved.size() && vecRemoved[iRemovedIdx] < vec[i])
		{
			++iNRemovedBefore;
			++iRemovedIdx;
		}
		vec[i] = vec[i]-iNRemovedBefore;
	}
	// Decrease the last ids (we get to this point when we are done with <vecRemoved>).
	for (; i < (int) vec.size(); ++i)
		vec[i] = vec[i]-iNRemovedBefore;
}

/***********************************************************************
 * template<class T> 
 * void vecMoveToBack(vector<T> &vec, int iStart, int iEnd)
 *
 * Moves a segment of a vector to the back of the vector. 
 * Parameters are only checked in debug mode.
 * Used in BCP-EVRPTW. 
 *
 * Example:
 * input: vec = (3, 7, 9, 15, 18, 24), iStart = 2, iEnd = 3
 * result: vec = (3, 7, 18, 24, 9, 15)
 * 
 * --- input parameters ---
 * vec				: vector to permute
 * iStart			: start of sequence to move. 0 <= iStart < vec.size()
 * iEnd				: end of sequence to move. 0 <= iStart <= iEnd < vec.size()
 * --- output parameters ---
 * --- return value ---
 * 
 ***********************************************************************/
template<class T> 
void vecMoveToBack(vector<T> &vec, int iStart, int iEnd)
{
	// Could probably be implemented with some STL magic, but one have to be careful.
	// For example, inserting the sequence in question at the end of the vector and then 
	// deleting the sequence  is bound to give problems as the insert operation may cause the
	// vector to be reallocated.
	// We take a simple and safe approach.
	assert(0 <= iStart && iStart <= iEnd && iEnd < vec.size());
	vector<T> vecCopy(vec.begin()+iStart, vec.begin()+iEnd+1);
	vec.erase(vec.begin()+iStart, vec.begin()+iEnd+1);
	vec.insert(vec.end(), vecCopy.begin(), vecCopy.end());
}

#endif

