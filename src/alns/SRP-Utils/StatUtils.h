#ifndef STATUTILS_H
#define STATUTILS_H

#include <vector>

#include "Utils.h"

using namespace std;

/***********************************************************************
 * template <class T>
 * T findAvg(const vector<T> &vec)
 *
 * Find the average of the elements in <vec>. It is assumed that the vector is non-empty.
 *
 * --- input parameters ---
 * vec			: The vector to calculate average of.
 * --- output parameters ---
 * --- return values ---
 * The average value (calculated using type T, thus if T is integer, then the
 * value is going to be rounded).
 ***********************************************************************/

template <class T>
T findAvg(const vector<T> &vec) 
{
	if (vec.empty())
		error("findAvg(...)", "vector is empty");
	int i;
	T sum = vec[0];
	for (i=1; i < (int) vec.size(); i++)
	{
		sum += vec[i];
	}
	return sum / vec.size();
}

// As above, but result is returned in a double. This mean that the average of
// a vector of integers will be computed correctly.

template <class T>
double findAvgDouble(const vector<T> &vec) 
{
	if (vec.empty())
		error("findAvg(...)", "vector is empty");
	int i;
	T sum = vec[0];
	for (i=1; i < (int) vec.size(); i++)
	{
		sum += vec[i];
	}
	return (double) sum / (double) vec.size();
}

// if the vector contains zero or one element the result of sample standard deviation is undefined. 
// In this case we simply return 0.
template <class T>
double computeSampleStandardDeviation(const vector<T> &vec, double dMean) 
{
	if (vec.size() <= 1)
		return 0;
	int i;
	double dSum = 0;
	for (i=0; i < (int) vec.size(); i++)
		dSum += (vec[i]-dMean)*(vec[i]-dMean);
	return sqrt(1.0/(vec.size() - 1) * dSum);
}

template <class T>
double computeSampleStandardDeviation(const vector<T> &vec) 
{
	return computeSampleStandardDeviation(vec, findAvgDouble(vec));
}

#endif