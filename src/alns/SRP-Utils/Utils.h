#ifndef _UTILS_H
#define _UTILS_H

#include <string>
#include <ctime>
#include <algorithm>
#include <vector>
#include <iostream>
#include <assert.h>
#include <math.h>

template<class TCoord>
class Coordinate;

void error(std::string where, std::string what);

clock_t my_clock(void);
double elapsedSeconds();
double wallClock();
bool setAffinity(int iVal);

double trunc(double val);
double trunc(double val, unsigned int iPrecision);
double round(double val);
double round(double val, unsigned int iPrecision);
double roundUp(double val);
double roundUp(double val, unsigned int iPrecision);

int gcd(int a, int b);

double truncToNSignifDigits(double dVal, int iNDigits);

void generateFeasiblePos(const std::vector<int> &vecItemsLengths, int iMaxLength, std::vector<int> &vecFeasiblePos);

void parseNumbers(const std::string &str, std::vector<double> &vecNumbers);
void parseNumbers(const std::string &str, std::vector<int> &vecNumbers);
void splitString(const std::string &str, const std::string &strDelim, std::vector<std::string> &vecTokens, bool bRepeatedDelimsCountAsOne = true);
void stringReplace(const std::string &strIn, std::string &strOut, const std::string &strSearchFor, const std::string &strReplaceBy);

std::string secToHHMMSS(int sec);

#ifndef _WIN32
	#define __max(x,y) ((x) > (y) ? (x) : (y))
	#define __min(x,y) ((x) < (y) ? (x) : (y))
#endif

double calcDist(double dX1, double dY1, double dX2, double dY2);
double calcDist(const Coordinate<double> &c1, const Coordinate<double> &c2);

template<class T>
inline T maxFunc(T x,T y) 
{
	return (x > y ? x : y);
}

template<class T>
inline T minFunc(T x,T y) 
{
	return (x < y ? x : y);
}

template<class T>
inline T clamp(T x, T min, T max) 
{
	if (x < min)
		return min;
	if (x > max)
		return max;
	return x;
}


void waitEnter();

std::string bool2String(bool b, const std::string &strIsTrue = "true", const std::string &strIsFalse = "false");
std::string int2String(int i);
std::string double2String(double d, int iPrecision = 4);

int getNeededPrecision(double dVal, int iNDecimals);

void randomPermutation(std::vector<int> &vecInts, int iLower, int iUpper);
void firstPerm(std::vector<int> &vecPerm, int iN);
bool nextPerm(std::vector<int> &vecPerm, int iN);
std::string repeatString(int iNRepeats, const std::string &str);

void skipPath(const std::string &strFullPath, std::string &strFilename);
std::string getDirectory(const std::string &strFullPath);
std::string stripExtension(const std::string &strFileName);

// e.g. dropTrailingCharacters(str, " \t")  drops trailing whitespace
void dropTrailingCharacters(std::string &str, const std::string &dropChars);
// e.g. dropCharacters("super duper", "pu") would return "ser der"
std::string dropCharacters(const std::string &str, const std::string &dropChars);

template <class T>
inline T dotProduct(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
	if (vec1.size() != vec2.size())
		error("dotProduct(...)","std::vectors have different sizes");
	if (vec1.empty())
		error("dotProduct(...)","std::vectors are empty");
	int i, iSize = (int) vec1.size();
	
	T res(vec1.front() * vec2.front());
	for (i=1; i < iSize; ++i)
		res += vec1[i] * vec2[i];
	
	return res;
}

template <class T>
inline double norm(const std::vector<T> &vec)
{
	if (vec.empty())
		error("norm(...)","std::vector is empty");
	int i, iSize = (int) vec.size();
	
	T temp(vec.front() * vec.front());
	for (i=1; i < iSize; ++i)
		temp += vec[i] * vec[i];

	return sqrt((double) temp);
}

// returns cosinus to the angle between the hyperplanes defined by the normals std::vectors <vec1> and <vec2>
// The "real angle" can be obtained by using <acos> on the result.
template <class T>
inline double angle(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
	if (vec1.empty())
		error("hyperPlaneAngle(...)","std::vector is empty");
	if (vec1.size() != vec2.size())
		error("hyperPlaneAngle(...)","std::vectors have different length");
	
	return fabs((double) dotProduct(vec1, vec2)) / (norm(vec1) * norm(vec2));
}

inline bool epsilonEqual(double dNum1, double dNum2, double dEpsilon)
{
	return (fabs(dNum1-dNum2) <= dEpsilon);
}

/***********************************************************************
 * int factorial(int val)
 *
 * Calculates val!, that is is if val = 5 then the method calculates 
 * 5! = 5*4*3*2*1 = 120. <val> must be non-negative.
 *
 * --- input parameters ---
 * val			: The value we want to know the factorial of.
 * --- return values ---
 * The factorial of val.
 ***********************************************************************/

inline int factorial(int val)
{
	if (val == 0)
		return 1;
	if (val < 0)
	{
		error("factorial(int val)", "val < 0");
		return -1;
	}
	else
	{
		int res(val);
		val--;
		while(val > 1)
		{
			res *= val;
			val--;
		}
		return res;
	}
}

// *********************************************************************
// *** Methods for reading from ifstream
// *********************************************************************

void skip(std::ifstream &ifs, const std::string &strSkipString);

// *********************************************************************
// *********************************************************************
// *********************************************************************


/***********************************************************************
 * template <class T>
 * void resizeArray(int newSize, T **oldArray, int numberToCopy)
 *
 * Resizes the array pointed to by *oldArray. The method expects that
 * newSize >= numberToCopy. The function can be used to both enlarge and
 * shrink an array.
 *
 * --- input parameters ---
 * newSize			: The size of the new array
 * oldArray			: A pointer to a pointer that points to the first 
 *					  element in the array. This enables us to make the
 *					  original pointer point somewhere else (is this poetic or pure nonsense ???) ;-)
 * numberToCopy		: The number of elements to copy from the old array to 
 *					  the new.
 * defaultValue		: if newSize > numberToCopy then the rest of the elements is set to <defaultValue>.
 * 
 * --- return values ---
 * none
 ***********************************************************************/

template <class T>
void resizeArray(int newSize, T **oldArray, int numberToCopy, T defaultValue)
{
	T *newArray = new T[newSize];
	int i;
	for (i=0; i<numberToCopy; i++)
		newArray[i] = (*oldArray)[i];
	for (i = numberToCopy; i<newSize; i++)
		newArray[i] = defaultValue;
	delete [] (*oldArray);
	*oldArray = newArray;	
}

// Returns <end1> if no match was found.
// Returns an iterator in [begin1, end1) if a match was found. The iterator
// points to the beginning of the match.
template<class _FwdIt1, class _FwdIt2>
_FwdIt1 findSubString(_FwdIt1 begin1, _FwdIt1 end1, _FwdIt2 begin2, _FwdIt2 end2)
{	// Is [begin2, end2) a "substd::string" of [begin1, end1)
	bool bMatch = false;
	for (; begin1 != end1 && !bMatch; ++begin1)
	{
		bMatch = true;
		int i;
		for (i = 0; i < end2-begin2 && bMatch; i++)
		{
			if (begin1+i == end1)
			{
				bMatch = false;
			}
			else
				if (*(begin1+i) != *(begin2+i))
					bMatch = false;
		}
	}
	if (bMatch)
		return begin1;
	else
		return end1;
}

inline double clamp(double dVal, double dLB, double dUB) 
{  
	if (dVal < dLB)
		return dLB;
	if (dVal > dUB)
		return dUB;
	return dVal;
}

class Random
{	
private:
	static unsigned int _h48, _l48;
	//static int count;

	// Random generator, code "stolen" from David Pisingers Knapsack codes.
	static void srand48x(int s)
	{
		_h48 = s;
		_l48 = 0x330E;
	}

	static int lrand48x()
	{
		_h48 = (_h48 * 0xDEECE66D) + (_l48 * 0x5DEEC);
		_l48 = _l48 * 0xE66D + 0xB;
		_h48 = _h48 + (_l48 >> 16);
		_l48 = _l48 & 0xFFFF;
		return (_h48 >> 1);
	}

	static double myRand()		{ return (double) lrand48x() / 0x7FFFFFFF; }

public:
	static void initialize(int seed)
	{
		//count = 0;
		if (seed < 0)
			srand48x((int) time(0));
		else
			srand48x(seed);
	}
	
	// Perform a more robust initialisation.
	static void robustInitialise(int iSeed)
	{
		int i, iNReseeds = 5;
		for (i=0; i < iNReseeds; i++)
		{
			initialize(iSeed);
			if (i < iNReseeds-1)
				iSeed = Random::getRandom(0, 0x7FFFFFFF);
		}
	}

	// returns an integer between i and j, both are included. It is assumed that i <= j
	static int getRandom(int i, int j)					{ return i + lrand48x() % (j-i+1); }
	static double getRandomDouble(int i, int j)			{ return (double) i + myRand() * (double) (j-i); }
	static double getRandomDouble(double i, double j)	{ return i + myRand() * (j-i); }	
};

// Thread safe version of the Random class. 
// Does no longer use static members.
class TSRandom
{	
private:
	unsigned int _h48, _l48;
	int count;

	// Random generator, code "stolen" from David Pisingers Knapsack codes.
	void srand48x(int s)
	{
		_h48 = s;
		_l48 = 0x330E;
	}

	int lrand48x()
	{
		_h48 = (_h48 * 0xDEECE66D) + (_l48 * 0x5DEEC);
		_l48 = _l48 * 0xE66D + 0xB;
		_h48 = _h48 + (_l48 >> 16);
		_l48 = _l48 & 0xFFFF;
		return (_h48 >> 1);
	}

	double myRand()		{ return (double) lrand48x() / 0x7FFFFFFF; }

public:
	void initialize(int seed)
	{
		count = 0;
		if (seed < 0)
			srand48x((int) time(0));
		else
			srand48x(seed);
	}
	
	// Perform a more robust initialisation.
	void robustInitialise(int iSeed)
	{
		int i, iNReseeds = 5;
		for (i=0; i < iNReseeds; i++)
		{
			initialize(iSeed);
			if (i < iNReseeds-1)
				iSeed = getRandom(0, 0x7FFFFFFF);
		}
	}

	// returns an integer between i and j, both are included.
	int getRandom(int i, int j)					{ return i + lrand48x() % (j-i+1); }
	double getRandomDouble(int i, int j)		{ return (double) i + myRand() * (double) (j-i); }
	double getRandomDouble(double i, double j)	{ return i + myRand() * (j-i); }	
};

#endif
