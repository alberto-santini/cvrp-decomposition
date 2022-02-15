#include "VectorUtils.h"
#include "Utils.h"

#include <algorithm>

/***********************************************************************
 * int subset(const vector<int> &vecA, const vector<int> &vecB, int iMinElem, int iMaxElem)
 *
 * Determines if vecA is a subset of vecB or if the opposite is true. The method uses
 * memory proportional with iMaxElem - iMinElem. The running time is
 * O(iMaxElem-iMinElem + vecA.size() + vecB.size())
 *
 * The method could be made more efficient if we assume that the two vectors do not contain
 * duplicates (not by improving big O notation, but by improving constants).
 *
 * --- input parameters ---
 * vecA           :
 * vecB           : The two vectors to compare.
 * iMinElem       : The smallest element in both vectors should be larger than iMinElem
 * iMaxElem       : The largest element in both vectors should be smaller than iMaxElem
 * --- output parameters ---
 * --- return value ---
 * 1        : vecB is a subset of vecA ("vecA-vecB > 0")
 * 0        : vecA and VecB are identical (when seen as sets)
 * -99      : neither vecA or vecB are subsets of each other
 * -1       : vecA is a subset of vecB ("vecA-vecB < 0")
 ***********************************************************************/

int subset(const vector<int>& vecA, const vector<int>& vecB, int iMinElem, int iMaxElem) {
    vector<bool> vecElemInA(iMaxElem - iMinElem + 1, false);
    vector<bool> vecElemInB(iMaxElem - iMinElem + 1, false);
    int i;
    // Number of elements from A that also is in B.
    int iFromAInB = 0;
    int iUniqueCountA = 0, iUniqueCountB = 0;
    for(i = 0; i < (int)vecA.size(); i++) {
        assert(iMinElem <= vecA[i] && vecA[i] <= iMaxElem);

        if(!vecElemInA[vecA[i] - iMinElem]) {
            vecElemInA[vecA[i] - iMinElem] = true;
            iUniqueCountA++;
        }
    }
    for(i = 0; i < (int)vecB.size(); i++) {
        assert(iMinElem <= vecB[i] && vecB[i] <= iMaxElem);

        if(!vecElemInB[vecB[i] - iMinElem]) {
            vecElemInB[vecB[i] - iMinElem] = true;
            iUniqueCountB++;
            if(vecElemInA[vecB[i] - iMinElem])
                iFromAInB++;
        }
    }

    int iSubset = -99;
    if(iUniqueCountB == iUniqueCountA && iFromAInB == iUniqueCountA)
        // the two sets are identical
        iSubset = 0;
    else {
        if(iUniqueCountB < iUniqueCountA && iFromAInB == iUniqueCountB)
            // B is a subset of A
            iSubset = 1;
        if(iUniqueCountA < iUniqueCountB && iFromAInB == iUniqueCountA)
            // A is a subset of B
            iSubset = -1;
    }
    return iSubset;
}

// Used by <compareVecVecInt>. Copies each vector in <vecVecInt> to <vecSortedCopy>
// then sorts the each copied vector and finaly removes any duplicates from
// each vector.
// Used to transform each vector to a well defined definition of a set.
void sortAndTrim(const vector<vector<int>>& vecVecInt, vector<vector<int>>& vecSortedCopy) {
    int i, j;
    for(i = 0; i < (int)vecVecInt.size(); i++) {
        vecSortedCopy.push_back(vecVecInt[i]);
        vector<int>& vecInt = vecSortedCopy[i];
        sort(vecInt.begin(), vecInt.end());
        j = 0;
        while(j < (int)vecInt.size() - 1) {
            if(vecInt[j] == vecInt[j + 1])
                vecInt.erase(vecInt.begin() + j);
            else
                j++;
        }
    }
}

// Used by <compareVecVecInt>. Calculates a hash code for an order vector of integers.
int calcHashCode(const vector<int>& vec, const vector<int>& vecHashNumbers) {
    int i, iHashCode = 0;
    for(i = 0; i < (int)vec.size(); i++)
        iHashCode += vec[i] * vecHashNumbers[i];
    return iHashCode;
}

// Used by <compareVecVecInt>.
struct HashIdxPair {
    int m_iHashCode;
    int m_iIdx;
    HashIdxPair(int iHashCode, int iIdx) {
        m_iHashCode = iHashCode;
        m_iIdx = iIdx;
    }
    HashIdxPair() {}

    bool operator<(const HashIdxPair& other) const { return (m_iHashCode < other.m_iHashCode); }

    // To structures are different if their hash codes are different.
    bool operator!=(const HashIdxPair& other) const { return (m_iHashCode != other.m_iHashCode); }
};

/***********************************************************************
 * bool compareVecVecInt(const vector<vector<int> > &vecVecInt1, const vector<vector<int> > &vecVecInt2)
 *
 * Compares two vectors of vectors. The vectors are interpretted as two unordered list of sets.
 * The function compares if the two unordered list of sets are identical, that is,
 * if the each set in the first list can be coupled with an identical set in the other
 * list.
 *
 * The function is not too fast. It sorts every vector in the two input vectors
 * and removes duplicates (removing duplicates takes O(v^2) where <v> is the
 * length of the vector.
 *
 * --- input parameters ---
 * vecVecInt1       : Interpretted as an unordered list of sets.
 * vecVecInt		:
 * --- return value ---
 * true if the two unordered list of sets are identical.
 * --- protection status ---
 ***********************************************************************/

bool compareVecVecInt(const vector<vector<int>>& vecVecInt1, const vector<vector<int>>& vecVecInt2) {
    if(vecVecInt1.size() != vecVecInt2.size())
        return false;
    vector<vector<int>> vecSortedCopy1, vecSortedCopy2;
    sortAndTrim(vecVecInt1, vecSortedCopy1);
    sortAndTrim(vecVecInt2, vecSortedCopy2);
    int i;

    int iMaxSize = 0;
    for(i = 0; i < (int)vecSortedCopy1.size(); i++) {
        // if (vecSortedCopy1.size() != vecSortedCopy2.size())
        //	return false;
        if((int)vecSortedCopy1[i].size() > iMaxSize)
            iMaxSize = (int)vecSortedCopy1[i].size();
        if((int)vecSortedCopy2[i].size() > iMaxSize)
            iMaxSize = (int)vecSortedCopy2[i].size();
    }

    // Generate <iMaxSize> random numbers:
    vector<int> vecRandom(iMaxSize);
    for(i = 0; i < iMaxSize; i++)
        vecRandom[i] = Random::getRandom(0, 2000000000);

    vector<HashIdxPair> vecHashCodes1(vecSortedCopy1.size()), vecHashCodes2(vecSortedCopy2.size());
    for(i = 0; i < (int)vecSortedCopy1.size(); i++) {
        vecHashCodes1[i] = HashIdxPair(calcHashCode(vecSortedCopy1[i], vecRandom), i);
        vecHashCodes2[i] = HashIdxPair(calcHashCode(vecSortedCopy2[i], vecRandom), i);
    }

    sort(vecHashCodes1.begin(), vecHashCodes1.end());
    sort(vecHashCodes2.begin(), vecHashCodes2.end());

    // Quick check:
    for(i = 0; i < (int)vecHashCodes1.size(); i++)
        if(vecHashCodes1[i] != vecHashCodes2[i])
            return false;

    // Now we are pretty sure that the sets are identical, but we have to check carefully, the hash codes
    // cannot be trusted.
    for(i = 0; i < (int)vecHashCodes1.size(); i++) {
        if(vecSortedCopy1[vecHashCodes1[i].m_iIdx] != vecSortedCopy2[vecHashCodes2[i].m_iIdx])
            return false;
    }

    // Okay. The sets were identical.
    return true;
}

/***********************************************************************
 * void makeRandomIntVec(vector<int> &vec, int iNElems, int iMin, int iMax)
 *
 * Make a vector with <iNElems> random integers in the interval [iMin; iMax]
 *
 * --- input parameters ---
 * iNElems		: Number of elements in the vector
 * iMin, iMax	: The interval of the random numbers
 *
 * --- output parameters ---
 * vec				 : The random vector constructed by this method.
 * --- return value ---
 * true if the two unordered list of sets are identical.
 * --- protection status ---
 ***********************************************************************/

void makeRandomIntVec(vector<int>& vec, int iNElems, int iMin, int iMax) {
    vec.resize(iNElems);
    int i;
    for(i = 0; i < iNElems; i++)
        vec[i] = Random::getRandom(iMin, iMax);
}
