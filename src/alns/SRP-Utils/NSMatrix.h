// NSMatrix = non-square matrix

// @@@ It might be an idea to create a specialized version of this class if the type is 
// char, short or even int. In that case we could probably make the "+=" and "<=" operators a lot faster
// by implementing them using SSE intrinsic instructions.


#ifndef NSMATRIX_H
#define NSMATRIX_H

#include "Utils.h"

#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;

template<class T>
class NSMatrix
{
public:
	NSMatrix();
	NSMatrix(int iNRows, int iNCols);
	NSMatrix(int iNRows, int iNCols, T initValue);
	NSMatrix(const NSMatrix<T> &other);

	void setElement(int iRow, int iCol, T value) 
	{ 
		assert( (0 <= iRow && iRow < m_iNRows) && (0 <= iCol && iCol < m_iNCols));
		m_vecNSMatrix[iRow*m_iNCols + iCol] = value; 
	}
	void addToElement(int iRow, int iCol, T value) 
	{ 
		assert( (0 <= iRow && iRow < m_iNRows) && (0 <= iCol && iCol < m_iNCols));
		m_vecNSMatrix[iRow*m_iNCols + iCol] += value; 
	}
	inline T getElement(int iRow, int iCol) const 
	{ 
		assert( (0 <= iRow && iRow < m_iNRows) && (0 <= iCol && iCol < m_iNCols));
		return m_vecNSMatrix[iRow*m_iNCols + iCol]; 
	}
	const T &getConstRefToElem(int iRow, int iCol) const 
	{
		assert( (0 <= iRow && iRow < m_iNRows) && (0 <= iCol && iCol < m_iNCols));
		return m_vecNSMatrix[iRow*m_iNCols + iCol]; 
	}
	T &getRefToElem(int iRow, int iCol) 
	{
		assert( (0 <= iRow && iRow < m_iNRows) && (0 <= iCol && iCol < m_iNCols));
		return m_vecNSMatrix[iRow*m_iNCols + iCol]; 
	}
	int getNCols() const { return m_iNCols; }
	int getNRows() const { return m_iNRows; }
	
	void maxWith(const NSMatrix &other, NSMatrix &matMax) const;
	void maxWith(const NSMatrix &other);
	inline void maxWith(const T &val);

	// Finds the largest element in this matrix.
	T getMaxElement() const;
	// Find the column number of the smallest element in a row
	int findMinElementRow(int iRow) const;
	// Find the row number of the smallest element in a column
	int findMinElementColumn(int iColumn) const;

	void coverCount(const NSMatrix &other, int &iNCovered, int &iNWasted) const;

	template<class TDistType>
	TDistType dist(const NSMatrix &other, const TDistType &zero) const;
	bool allEqualTo(const T &compareVal) const;

	bool subtractLargerEqual(const NSMatrix &subtract, const NSMatrix &compare) const;

	// Add <other> to <this>.
	const NSMatrix &operator += (const NSMatrix &other);
	// Subtract <other> from <this>.
	const NSMatrix &operator -= (const NSMatrix &other);
	// Is the elements in this matrix less than or equal to the elements in <other> 
	// (element by element)?
	bool operator <= (const NSMatrix &other) const;
	// operator == : Is this matrix identical to <other> (element by element)
	bool operator == (const NSMatrix &other) const;

	void resizeAndClear(int iNRows, int iNCols, T initValue);
	void resize(int iNRows, int iNCols);
	// Expands the matrix, while keeping the original content.
	void expand(int iNRows, int iNColumns, T initValue);
	void outputMatrix(ostream &os, const string &strRowStart, const string &strRowEnd, const string &strRowDelim, const string &strElemDelim) const;

	template<class U>
	friend ostream& operator << (ostream& os, const NSMatrix<U> &mat); 
	
	template<class U>
	friend istream& operator>>(istream& is, NSMatrix<U>& mat);
	

private:
	vector<T> m_vecNSMatrix;
	int m_iNCols;
	int m_iNRows;
};

template<class T>
NSMatrix<T>::NSMatrix()
{
	m_iNCols = m_iNRows = -999;
}
template<class T>
NSMatrix<T>::NSMatrix(int iNRows, int iNCols)
{
	m_iNCols = iNCols;
	m_iNRows = iNRows;
	m_vecNSMatrix.resize(iNCols*iNRows);
}

template<class T>
NSMatrix<T>::NSMatrix(int iNRows, int iNCols, T initValue)
{
	m_iNCols = iNCols;
	m_iNRows = iNRows;
	m_vecNSMatrix.resize(iNCols*iNRows, initValue);
}
template<class T>
NSMatrix<T>::NSMatrix(const NSMatrix<T> &other)
	: m_vecNSMatrix(other.m_vecNSMatrix)
{
	m_iNCols = other.m_iNCols;
	m_iNRows = other.m_iNRows;
}

template<class T>
void NSMatrix<T>::resizeAndClear(int iNRows, int iNCols, T initValue)
{
	m_iNCols = iNCols;
	m_iNRows = iNRows;
	m_vecNSMatrix.clear();
	m_vecNSMatrix.resize(iNCols*iNRows, initValue);
}

template<class T>
void NSMatrix<T>::resize(int iNRows, int iNCols)
{
	m_iNCols = iNCols;
	m_iNRows = iNRows;
	m_vecNSMatrix.resize(iNCols*iNRows);
}

// Expands the matrix, while keeping the original content.
template<class T>
void NSMatrix<T>::expand(int iNRows, int iNColumns, T initValue)
{
	if (iNRows < m_iNRows || iNColumns < m_iNCols)
	{
		error("NSMatrix<T>::expand(...)", 
			"new matrix dimensions are smaller than the original dimensions"); 
	}
	// Is there anything to do at all:
	if (iNRows > m_iNRows || iNColumns > m_iNCols)
	{
		// Yes, we have to do some work.
		// If we only add rows, then the expand operation is pretty simple:
		if (m_iNCols == iNColumns)
			m_vecNSMatrix.resize(iNColumns*iNRows, initValue);
		else
		{
			// we are adding columns and perhaps also rows, this is a little more tricky.

			// Create a new vector with room for the expanded matrix:
			vector<T> vecTemp(iNColumns*iNRows, initValue);
			// Copy elements from original matrix.:
			int iCol, iRow;
			int iOrigOffset = 0, iNewOffset;
			for (iRow = 0; iRow < m_iNRows; ++iRow)
			{
				iNewOffset = iRow*iNColumns;
				// We could optimize here. There is actually no need to to have the 
				// iCol variable, we could either loop over iNewOffset or 
				// iOrigOffset, but that would be a little harder to understand.
				for (iCol = 0; iCol < m_iNCols; ++iCol)
				{
					vecTemp[iNewOffset] = m_vecNSMatrix[iOrigOffset];
					++iNewOffset;
					++iOrigOffset;
				}
			}
			// Use the new matrix and update the dimensions:
			m_vecNSMatrix.swap(vecTemp);
		}
		m_iNCols = iNColumns;
		m_iNRows = iNRows;
	}
}

/***********************************************************************
 * template<class T>
 * T NSMatrix<T>::dist(const NSMatrix &other) const 
 *
 * Calculate the distance between this matrix (A) and <other> (B). The distance is defined as
 * \Sum_{i,j} |A_ij - B_ij|. It is assumed that A and B have the same dimensions.
 *
 * --- input parameters ---
 * other			: The other matrix that we using in the calculation
 * zero				: The zero element for type T (in case the class is used with some weird abstract class).
 * --- output parameters ---
 * --- return value ---
 * the distance as defined above.
 * --- Protection status ---
 * public.
 ***********************************************************************/

template<class T>
template<class TDistType>
TDistType NSMatrix<T>::dist(const NSMatrix &other, const TDistType &zero) const
{
	TDistType dist(zero);
	assert(m_iNCols == other.m_iNCols && m_iNRows == other.m_iNRows);
	int i, iMax = m_iNCols*m_iNRows;
	for (i=0; i < iMax; i++)
		dist += abs((T) (m_vecNSMatrix[i]-other.m_vecNSMatrix[i]));
	return dist;
}

/***********************************************************************
 * template<class T>
 * void NSMatrix<T>::maxWith(const NSMatrix &other, NSMatrix &matMax) const
 *
 * computes a new matrix whose entries are the maximum of the corresponding entries in 
 * the <this> and <other> matrices. If <this> and <other> are A and B and <matMax> is C, then
 *
 * C_ij = max { A_ij, B_ij}  for all i,j
 *
 * --- input parameters ---
 * other			: The other matrix that we using in the calculation
 * --- output parameters ---
 * matMax			: The "maximum" matrix computed as stated above.
 * --- return value ---
 * --- Protection status ---
 * public.
 ***********************************************************************/
template<class T>
void NSMatrix<T>::maxWith(const NSMatrix &other, NSMatrix &matMax) const
{
	assert(m_iNCols == other.m_iNCols && m_iNRows == other.m_iNRows);
	matMax.resize(m_iNCols, m_iNRows);
	int i, iMax = m_iNCols*m_iNRows;
	for (i=0; i < iMax; i++)
		matMax.m_vecNSMatrix[i] = maxFunc(m_vecNSMatrix[i], other.m_vecNSMatrix[i]);
}

/***********************************************************************
 * template<class T>
 * void NSMatrix<T>::maxWith(const NSMatrix &other) 
 *
 * Updates the entries in this matrix such that each element is equal to the max
 * of that element and the corresponding element in <other>.
 * If <this> and <other> are A and B, then
 *
 * A_ij = max { A_ij, B_ij}  for all i,j
 *
 * --- input parameters ---
 * other			: The other matrix that we using in the calculation
 * --- output parameters ---
 * --- return value ---
 * --- Protection status ---
 * public.
 ***********************************************************************/
template<class T>
void NSMatrix<T>::maxWith(const NSMatrix &other)
{
	assert(m_iNCols == other.m_iNCols && m_iNRows == other.m_iNRows);
	int i, iMax = m_iNCols*m_iNRows;
	for (i=0; i < iMax; i++)
		m_vecNSMatrix[i] = maxFunc(m_vecNSMatrix[i], other.m_vecNSMatrix[i]);
}

/***********************************************************************
 * template<class T>
 * void NSMatrix<T>::maxWith(const T &elem) 
 *
 * Updates the entries in this matrix such that each element is equal to the max
 * of that element and <elem>
 * If <this> is A, then
 *
 * A_ij = max { A_ij, elem}  for all i,j
 *
 * --- input parameters ---
 * elem				: The element we compare agains. <elem> is the smallest element in the
 *					  resulting matrix.
 * --- output parameters ---
 * --- return value ---
 * --- Protection status ---
 * public.
 ***********************************************************************/
template<class T>
inline void NSMatrix<T>::maxWith(const T &val)
{
	int i, iMax = m_iNCols*m_iNRows;
	for (i=0; i < iMax; i++)
		m_vecNSMatrix[i] = maxFunc(m_vecNSMatrix[i], val);
}

/***********************************************************************
 * template<class T>
 * const NSMatrix &NSMatrix<T>::operator += (const NSMatrix &other) 
 *
 * Add <other> to <this>.
 *
 * --- input parameters ---
 * other			: The other matrix that we using in the calculation
 * --- output parameters ---
 * --- return value ---
 * --- Protection status ---
 * public.
 ***********************************************************************/

template<class T>
const NSMatrix<T> &NSMatrix<T>::operator += (const NSMatrix<T> &other)
{
	assert(m_iNCols == other.m_iNCols && m_iNRows == other.m_iNRows);
	int i, iMax = m_iNCols*m_iNRows;
	for (i=0; i < iMax; i++)
		m_vecNSMatrix[i] += other.m_vecNSMatrix[i];
	return *this;
}

/***********************************************************************
 * template<class T>
 * const NSMatrix &NSMatrix<T>::operator -= (const NSMatrix &other) 
 *
 * Subtract <other> from <this>.
 *
 * --- input parameters ---
 * other			: The other matrix that we using in the calculation
 * --- output parameters ---
 * --- return value ---
 * --- Protection status ---
 * public.
 ***********************************************************************/
template<class T>
const NSMatrix<T> &NSMatrix<T>::operator -= (const NSMatrix &other)
{
	assert(m_iNCols == other.m_iNCols && m_iNRows == other.m_iNRows);
	int i, iMax = m_iNCols*m_iNRows;
	for (i=0; i < iMax; i++)
		m_vecNSMatrix[i] -= other.m_vecNSMatrix[i];
	return *this;
}

/***********************************************************************
 * template<class T>
 * bool NSMatrix<T>::operator <= (const NSMatrix &other) const
 *
 * Is the elements in this matrix less than or equal to the elements in <other> (element by element)?
 *
 * --- input parameters ---
 * other			: The other matrix we are comparing to.
 * --- output parameters ---
 * --- return value ---
 * true if the elements in this matrix is less than or equal to the elements in the <other>
 * matrix.
 * --- Protection status ---
 * public.
 ***********************************************************************/
template<class T>
bool NSMatrix<T>::operator <= (const NSMatrix &other) const
{
	assert(m_iNCols == other.m_iNCols && m_iNRows == other.m_iNRows);
	int i, iMax = m_iNCols*m_iNRows;
	bool bLEQ = true;
	for (i=0; i < iMax && bLEQ; i++)
		bLEQ = m_vecNSMatrix[i] <= other.m_vecNSMatrix[i];
	return bLEQ;
}
	
/***********************************************************************
 * template<class T>
 * bool NSMatrix<T>::operator == (const NSMatrix &other) const
 *
 * Is the elements in this matrix equal to the elements in <other> (element by element)?
 * If the two matrices have different dimensions then they are different and
 * the method returns false.
 * 
 * --- input parameters ---
 * other			: The other matrix we are comparing to.
 * --- output parameters ---
 * --- return value ---
 * true if the elements in this matrix is equal to the elements in the <other>
 * matrix.
 * --- Protection status ---
 * public.
 ***********************************************************************/
template<class T>
bool NSMatrix<T>::operator == (const NSMatrix &other) const
{
	
	int i, iMax = m_iNCols*m_iNRows;
	// If the two matrices are to be equal, then they must have the same number of 
	// elements.
	bool bEQ = (m_iNCols == other.m_iNCols && m_iNRows == other.m_iNRows);
	for (i=0; i < iMax && bEQ; i++)
		bEQ = m_vecNSMatrix[i] == other.m_vecNSMatrix[i];
	return bEQ;
}

/***********************************************************************
 * template<class T>
 * bool NSMatrix<T>::allEqualTo(const T &compareVal) const
 *
 * Checks if all elements in the matrix are equal to <compareVal>.
 *
 * --- input parameters ---
 * compareVal		: The value to compare against
 * --- output parameters ---
 * --- return value ---
 * true if all elements in the matrix are equal to <compareVal>.
 * matrix.
 * --- Protection status ---
 * public.
 ***********************************************************************/
template<class T>
bool NSMatrix<T>::allEqualTo(const T &compareVal) const
{
	int i, iMax = m_iNCols*m_iNRows;
	bool bAllEqualTo = true;
	for (i=0; i < iMax && bAllEqualTo; i++)
		bAllEqualTo = (m_vecNSMatrix[i] == compareVal);
	return bAllEqualTo;
}

/***********************************************************************
 * template<class T>
 * bool NSMatrix<T>::subtractLargerEqual(const NSMatrix &subtract, const NSMatrix &compare) const
 *
 * Tests if [this] - [subtract] >= [compare]. The >= returns true if the matrix
 * on the LHS. is greater or equal to the matrix on the RHS, element by element.
 *
 * --- input parameters ---
 * subtract			: The matrix to subtract from this matrix.
 * compare			: The matrix to compare with.
 * --- output parameters ---
 * --- return value ---
 * true if [this] - [subtract] >= [compare].
 * --- Protection status ---
 * public.
 ***********************************************************************/

template<class T>
bool NSMatrix<T>::subtractLargerEqual(const NSMatrix &subtract, const NSMatrix &compare) const
{
	int i, iMax = m_iNCols*m_iNRows;
	bool bSubtractLargerEqual = true;
	for (i=0; i < iMax && bSubtractLargerEqual; i++)
	{
		bSubtractLargerEqual = 
			(m_vecNSMatrix[i] - subtract.m_vecNSMatrix[i] >= compare.m_vecNSMatrix[i]);
	}
	return bSubtractLargerEqual;
}

/***********************************************************************
 * template<class T>
 * void NSMatrix<T>::coverCount(const NSMatrix &other, int &iCoverCount, int &iWasteCount) const
 *
 * Compares two matrices (<this> vs. <other>). Let A_ij = this and B_ij = other.
 * Assume A_ij >= 0 and B_ij >= 0 for i,j. The method computes
 * iCoverCount = sum(i,j) { min {A_ij, B_ij} }
 * iWasteCount = sum(i,j) { max {0, B_ij-A_ij } }
 *
 * Example 
 * This:              other
 * 0 1 1              1 2 1 
 * 2 3 1              1 0 0 
 *
 * iCoverCount = 3
 * iWasteCount = 2
 *
 * The method is used in the Roadef2007 project to see how well a technician matches a 
 * certain task.
 *
 * --- input parameters ---
 * other			: The other matrix to use in the comparison
 * --- output parameters ---
 * iCoverCount		: Formal definition, see above. Informal: Counts how well <other> covers <this>
 * iWasteCount		: Formal definition, see above. Informal: Counts how many entries in <other>
 *					  we waste if we "cover" this by "other" (this makes sense in the Roadef2007
 *					  environment).
 * --- return value ---
 * --- Protection status ---
 * public.
 ***********************************************************************/

template<class T>
void NSMatrix<T>::coverCount(const NSMatrix &other, int &iCoverCount, int &iWasteCount) const
{
	iCoverCount = 0; iWasteCount = 0;
	assert(m_iNCols == other.m_iNCols && m_iNRows == other.m_iNRows);
	typename vector<T>::const_iterator cIter = m_vecNSMatrix.begin();
	typename vector<T>::const_iterator cIterOther = other.m_vecNSMatrix.begin();
	typename vector<T>::const_iterator cIterEnd = m_vecNSMatrix.end();
	for (; cIter != cIterEnd; cIter++)
	{
		int iDiff = *cIter - *cIterOther;
		if (iDiff >= 0)
			iCoverCount += *cIterOther;
		else
		{
			iCoverCount += *cIter;
			iWasteCount += -iDiff;
		}
		cIterOther++;
	}
}

/***********************************************************************
 * template<class T>
 * T NSMatrix<T>::getMaxElement() const
 *
 * Finds the largest element in this matrix. The matrix
 * must be non-empty.
 *
 * --- input parameters ---
 * --- output parameters ---
 * --- return value ---
 * The largest value in the matrix
 * --- Protection status ---
 * public.
 ***********************************************************************/

template<class T>
T NSMatrix<T>::getMaxElement() const
{
	assert(!m_vecNSMatrix.empty());

	// @@@ if the T object is very large (and slow to copy) then this method is not
	// so clever, then we should rather keep a pointer to the largest element so far
	// and only do a copy of the T object when we return the largest value.
	
	// If T is simple, then the current code should be fine as our comparisons are going
	// to be faster.
	T maxElem = m_vecNSMatrix[0];
	int i;
	for (i=1; i < (int) m_vecNSMatrix.size(); i++)
		if (m_vecNSMatrix[i] > maxElem)
			maxElem = m_vecNSMatrix[i];

	return maxElem;
}

// Find the column number of the smallest element in a row
template<class T>
int NSMatrix<T>::findMinElementRow(int iRow) const
{
	if (m_iNCols == 0)
		error("NSMatrix<T>::findMinElementRow(...)", "No columns");
	
	int i;
	int iStart = iRow*m_iNCols;
	int iEnd = (iRow+1)*m_iNCols;
	// @@@ If T is a very large class, then it's not efficient to declare a new T object and
	// copy the T object every time a new better object has been found. In that case it would
	// be better to keep an index to the current best element.
	T min = m_vecNSMatrix[iStart];
	int iMinIdx = iStart;
	
	for (i=iStart+1; i < iEnd; ++i)
	{
		if (m_vecNSMatrix[i] < min)
		{
			min = m_vecNSMatrix[i];
			iMinIdx = i;
		}
	}
	return iMinIdx-iStart;
}

// Find the row number of the smallest element in a column
template<class T>
int NSMatrix<T>::findMinElementColumn(int iColumn) const
{
	if (m_iNRows == 0)
		error("NSMatrix<T>::findMinElementColumn(...)", "No rows");

	int i;
	int iStart = iColumn;
	int iEnd = m_iNRows*m_iNCols;
	// @@@ If T is a very large class, then it's not efficient to declare a new T object and
	// copy the T object every time a new better object has been found. In that case it would
	// be better to keep an index to the current best element.
	T min = m_vecNSMatrix[iStart];
	int iMinIdx = iStart;
	
	for (i=iStart+m_iNCols; i < iEnd; i += m_iNCols)
	{
		if (m_vecNSMatrix[i] < min)
		{
			min = m_vecNSMatrix[i];
			iMinIdx = i;
		}
	}
	return iMinIdx / m_iNCols;
}

/***********************************************************************
 * template<class T>
 * void NSMatrix<T>::outputMatrix(ostream &os, const string &strRowStart, const string &strRowEnd, const string &strRowDelim, const string &strElemDelim) const
 *
 * "pretty prints" the matrix. 
 *
 * --- input parameters ---
 * os			: stream to write on.
 * strRowStart	: Written at the start of each row
 * strRowEnd	: Written at the end of each row
 * strRowDelim	: Written between each row. Not written after the last row.
 * strElemDelim	: Written in between each column. Not written after the last element in each row.
 * --- output parameters ---
 * --- return value ---
 * --- Protection status ---
 * public.
 ***********************************************************************/

template<class T>
void NSMatrix<T>::outputMatrix(ostream &os, const string &strRowStart, const string &strRowEnd, const string &strRowDelim, const string &strElemDelim) const
{
	int iRow, iCol;
	for (iRow = 0; iRow < m_iNRows; iRow++)
	{
		os << strRowStart;
		for (iCol = 0; iCol < m_iNCols; iCol++)
		{
			os << getElement(iRow, iCol);
			if (iCol < m_iNCols-1) 
				os << strElemDelim; 
		}
		os << strRowEnd;
		if (iRow < m_iNRows-1)
			os << strRowDelim;
	}	
}

template<class U>
ostream& operator<<(ostream& os, const NSMatrix<U>& mat) 
{
	int iRow, iCol;
	for (iRow = 0; iRow < mat.m_iNRows; iRow++)
	{
		for (iCol = 0; iCol < mat.m_iNCols; iCol++)
		{
			os << mat.getElement(iRow, iCol) << " "; 
		}
		if (iRow < mat.m_iNRows-1)
			os << endl;
	}	
	return os;
}

// Reads the content of a matrix from a stream. The mehtod simply reads m_iNRows * m_iNCols values
// from <is> (that is, the dimension of the matrix should be known before using this method.
template<class U>
istream& operator>>(istream& is, NSMatrix<U>& mat) 
{
	int iRow, iColumn;
	for (iRow = 0; iRow < mat.m_iNRows; iRow++)
	{
		for (iColumn = 0; iColumn < mat.m_iNCols; iColumn++)
		{
			U val;
			is >> val;
			mat.setElement(iRow, iColumn, val);
		}
	}
	return is;
}

#endif
