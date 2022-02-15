#ifndef MATRIX_H
#define MATRIX_H

#include "Utils.h"

#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;

template<class T>
class Matrix
{
public:
	Matrix();
	Matrix(int iDim);
	Matrix(int iDim, T initValue);
	Matrix(const Matrix<T> &other);
	Matrix<T>& operator=(const Matrix<T>& other) { m_vectorMatrix = other.m_vectorMatrix; m_iDim = other.m_iDim; return *this; }

#ifdef _WIN32
	__forceinline 
#endif 
	void setElement(int i, int j, T value)			
	{ 
		assert(0 <= i && i < m_iDim && 0 <= j && j < m_iDim);
		m_vectorMatrix[i*m_iDim + j] = value; 
	}

	void addToElement(int i, int j, T value)		
	{ 
		assert(0 <= i && i < m_iDim && 0 <= j && j < m_iDim);
		m_vectorMatrix[i*m_iDim + j] += value; 
	}

#ifdef _WIN32
	__forceinline 
#endif 
	T getElement(int i, int j) const				
	{ 
		assert(0 <= i && i < m_iDim && 0 <= j && j < m_iDim);
		return m_vectorMatrix[i*m_iDim + j]; 
	}

#ifdef _WIN32
	__forceinline 
#endif 
	T &getRefToElem(int i, int j)					
	{ 
		assert(0 <= i && i < m_iDim && 0 <= j && j < m_iDim);
		return m_vectorMatrix[i*m_iDim + j]; 
	}

#ifdef _WIN32
	__forceinline 
#endif 
	const T &getConstRefToElem(int i, int j) const	
	{ 
		assert(0 <= i && i < m_iDim && 0 <= j && j < m_iDim);
		return m_vectorMatrix[i*m_iDim + j]; 
	}

#ifdef _WIN32
	__forceinline 
#endif 
	int getDim() const								{ return m_iDim; }
	
	inline bool operator <= (const Matrix<T> &other) const;
	void clear(T initValue)							{ fill(m_vectorMatrix.begin(), m_vectorMatrix.end(), initValue); }
	void resizeAndClear(int iDim, T initValue);
	// Expands the matrix, while keeping the original content.
	void expand(int iDim, T initValue);
	// Shrink the matrix, keep the content that is not touched by the shrinking operation
	void shrink(int iDim);
	
	template<class T2> friend ostream& operator<<(ostream& os, const Matrix<T2>& mat);
	template<class T2> friend istream& operator>>(istream& is, Matrix<T2>& mat);

private:
	vector<T> m_vectorMatrix;
	int m_iDim;
};

template<class T>
inline bool Matrix<T>::operator <= (const Matrix<T> &other) const
{
	assert(m_iDim == other.m_iDim);
	int i;
	bool bLEQ = true;
	int iNElements = m_iDim*m_iDim;
	for (i=0; i < iNElements && bLEQ; i++)
		bLEQ = m_vectorMatrix[i] <= other.m_vectorMatrix[i];
	return bLEQ;
}

template<class T>
Matrix<T>::Matrix()
{
	m_iDim = -999;
}
template<class T>
Matrix<T>::Matrix(int iDim)
{
	m_iDim = iDim;
	m_vectorMatrix.resize(iDim*iDim);
}

template<class T>
Matrix<T>::Matrix(int iDim, T initValue)
{
	m_iDim = iDim;
	m_vectorMatrix.resize(iDim*iDim, initValue);
}
template<class T>
Matrix<T>::Matrix(const Matrix<T> &other)
	: m_vectorMatrix(other.m_vectorMatrix)
{
	m_iDim = other.m_iDim;
}

template<class T>
void Matrix<T>::resizeAndClear(int iDim, T initValue)
{
	m_iDim = iDim;
	m_vectorMatrix = vector<T>(iDim*iDim, initValue);
}

// Expands the matrix, while keeping the original content.
template<class T>
void Matrix<T>::expand(int iDim, T initValue)
{
	if (iDim < m_iDim)
	{
		error("NSMatrix<T>::expand(...)", 
			"new matrix dimension is smaller than the original dimension"); 
	}
	// Is there anything to do at all:
	if (iDim > m_iDim)
	{
		// Yes, we have to do some work.
		
		// Create a new vector with room for the expanded matrix:
		vector<T> vecTemp(iDim*iDim, initValue);
		// Copy elements from original matrix.:
		int iCol, iRow;
		int iOrigOffset = 0, iNewOffset;
		for (iRow = 0; iRow < m_iDim; ++iRow)
		{
			iNewOffset = iRow*iDim;
			// We could optimize here. There is actually no need to have the 
			// iCol variable, we could either loop over iNewOffset or 
			// iOrigOffset, but that would be a little harder to understand.
			for (iCol = 0; iCol < m_iDim; ++iCol)
			{
				vecTemp[iNewOffset] = m_vectorMatrix[iOrigOffset];
				++iNewOffset;
				++iOrigOffset;
			}
		}
		// Use the new matrix and update the dimension:
		m_vectorMatrix.swap(vecTemp);
		m_iDim = iDim;
	}
}

// Shrink the matrix, keep the content that is not touched by the shrinking operation
template<class T>
void Matrix<T>::shrink(int iDim)
{
	if (iDim > m_iDim)
	{
		error("NSMatrix<T>::shrink(...)", 
			"new matrix dimension is larger than the original dimension"); 
	}
	// Is there anything to do at all:
	if (iDim < m_iDim)
	{
		// Yes, we have to do some work.	

		// Copy elements from original matrix.
		// Notice, no need to copy first row.
		int iRow, iOrigOffset, iNewOffset = iDim;
		for (iRow = 1; iRow < iDim; ++iRow)
		{
			int iNewOffsetEnd = iNewOffset + iDim;
			iOrigOffset = iRow*m_iDim;
			
			for ( ; iNewOffset < iNewOffsetEnd; ++iNewOffset)
			{
				m_vectorMatrix[iNewOffset] = m_vectorMatrix[iOrigOffset];
				++iOrigOffset;
			}
		}
		// Update the dimension:
		m_iDim = iDim;
		// Shrink the vector that stores the matrix:
		m_vectorMatrix.resize(m_iDim*m_iDim);

	}
}

template<class T2>
ostream& operator<<(ostream& os, const Matrix<T2>& mat) 
{
	int iRow, iColumn;
	for (iRow = 0; iRow < mat.getDim(); iRow++)
	{
		for (iColumn = 0; iColumn < mat.getDim(); iColumn++)
			os << mat.getElement(iRow, iColumn) << " ";
		os << endl;
	}
	return os;
}

// Reads the content of a matrix from a stream. The mehtod simply reads iDim * iDim values
// from <is> (that is, the dimension of the matrix should be known before using this method.
template<class T2>
istream& operator>>(istream& is, Matrix<T2>& mat) 
{
	int iRow, iColumn;
	for (iRow = 0; iRow < mat.getDim(); iRow++)
	{
		for (iColumn = 0; iColumn < mat.getDim(); iColumn++)
		{
			T2 val;
			is >> val;
			mat.setElement(iRow, iColumn, val);
		}
	}
	return is;
}

#endif
