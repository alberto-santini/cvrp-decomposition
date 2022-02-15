#ifndef GRAPHOBJECTS_H
#define GRAPHOBJECTS_H

#include <iostream>

using namespace std;

class SimpleArc;

class SimpleEdge
{
public:
	SimpleEdge()
	{
		m_iNode1 = -1;
		m_iNode2 = -1;
	}

	SimpleEdge(int iNode1, int iNode2)
	{
		m_iNode1 = iNode1;
		m_iNode2 = iNode2;
	}

	int getNode1() const { return m_iNode1; }
	int getNode2() const { return m_iNode2; }

	bool adjacent(const SimpleEdge &other) const
	{
		return (m_iNode1 == other.m_iNode1 || m_iNode1 == other.m_iNode2 ||
			m_iNode2 == other.m_iNode1 || m_iNode2 == other.m_iNode2); 
	}

	bool operator < (const SimpleEdge &other) const
	{
		return m_iNode1 < other.m_iNode1 || 
			(m_iNode1 == other.m_iNode1 && m_iNode2 < other.m_iNode2);
	}
	bool operator == (const SimpleEdge &other) const
	{
		return (m_iNode1 == other.m_iNode1 && m_iNode2 == other.m_iNode2);
	}
	bool operator != (const SimpleEdge &other) const
	{
		return (m_iNode1 != other.m_iNode1 || m_iNode2 != other.m_iNode2);
	}
	
	const SimpleEdge &operator = (const SimpleEdge &other) 
	{
		if (this != &other)
		{
			m_iNode1 = other.m_iNode1;
			m_iNode2 = other.m_iNode2;
		}
		return *this;
	}
	
	friend ostream& operator << (ostream& os, const SimpleEdge &edge);

private:
	int m_iNode1;
	int m_iNode2;
};

// ******************************************************************************
// ******************************************************************************
// *************************   class SimpleArc   ********************************
// ******************************************************************************
// ******************************************************************************

// For now <SimpleArc> is the same as <SimpleEdge>, but we ought to change 
// <SimpleEdge>, so it represents an edge in an undirected graph.

class SimpleArc
{
public:
	SimpleArc()
	{
		m_iNode1 = -1;
		m_iNode2 = -1;
	}

	SimpleArc(int iNode1, int iNode2)
	{
		m_iNode1 = iNode1;
		m_iNode2 = iNode2;
	}

	int getNode1() const { return m_iNode1; }
	int getNode2() const { return m_iNode2; }
	
	bool operator < (const SimpleArc &other) const
	{
		return m_iNode1 < other.m_iNode1 || 
			(m_iNode1 == other.m_iNode1 && m_iNode2 < other.m_iNode2);
	}
	bool operator == (const SimpleArc &other) const
	{
		return (m_iNode1 == other.m_iNode1 && m_iNode2 == other.m_iNode2);
	}
	bool operator != (const SimpleArc &other) const
	{
		return (m_iNode1 != other.m_iNode1 || m_iNode2 != other.m_iNode2);
	}
	
	friend ostream& operator << (ostream& os, const SimpleArc &edge);

private:
	int m_iNode1;
	int m_iNode2;
};

// ******************************************************************************
// ******************************************************************************
// ***************************   class Edge   ***********************************
// ******************************************************************************
// ******************************************************************************

template <class T>
class Edge
{
public:
	Edge(int iNode1, int iNode2, T weight)
	{
		m_iNode1 = iNode1;
		m_iNode2 = iNode2;
		m_weight = weight;
	}

	int getNode1() const { return m_iNode1; }
	int getNode2() const { return m_iNode2; }
	const T &getWeight() const { return m_weight; }
	bool operator < (const Edge<T> &other) const;
	SimpleEdge getSimpleEdge() const;
	SimpleArc getSimpleArc() const;

	template<class U> friend ostream& operator << (ostream& os, const Edge<U> &edge);

private:
	int m_iNode1;
	int m_iNode2;
	T m_weight;
};

template<class U>
ostream& operator<<(ostream& os, const Edge<U>& edge) {
	os << "(" << edge.getNode1() << ", " << edge.getNode2() << "), weight = " << edge.m_weight;
	return os;
}

// An edge1 < edge2 if the weight of edge1 is less than the weight of edge2 or if the weights are equal and 
// the lexigographic ordering says that edge1 < edge2.
template <class T>
bool Edge<T>::operator < (const Edge<T> &other) const
{
	return 
		m_weight < other.m_weight || 
		(
			m_weight == other.m_weight && 
			(m_iNode1 < other.m_iNode1 || (m_iNode1 == other.m_iNode1 && m_iNode2 < other.m_iNode2))
		);
}

template <class T>
SimpleEdge Edge<T>::getSimpleEdge() const
{
	SimpleEdge simpleEdge(m_iNode1, m_iNode2);
	return simpleEdge;
}
template <class T>
SimpleArc Edge<T>::getSimpleArc() const
{
	SimpleArc simpleArc(m_iNode1, m_iNode2);
	return simpleArc;
}
// For sorting in decreasing order. UDGrater(e1,e2) = (e2 < e1)
template<class T>
bool EdgeGreater ( const Edge<T> &edge1, const Edge<T> &edge2 )
{
	return edge1.getWeight() > edge2.getWeight() || 
		(
			edge1.getWeight() == edge2.getWeight() && 
			(edge1.getNode1() > edge2.getNode1() || (edge1.getNode1() == edge2.getNode1() && edge1.getNode2() > edge2.getNode2()))
		);
}

// ******************************************************************************
// ******************************************************************************
// ***************************  class MyArc   ***********************************
// ******************************************************************************
// ******************************************************************************

// We do not use the naem "Arc" as it could cause conflict with some Microsoft functions.

template <class T>
class MyArc
{
public:
	MyArc()
	{
		m_iNode1 = -1;
		m_iNode2 = -1;
	};
	
	MyArc(int iNode1, int iNode2, T weight)
	{
		m_iNode1 = iNode1;
		m_iNode2 = iNode2;
		m_weight = weight;
	}

	int getNode1() const { return m_iNode1; }
	int getNode2() const { return m_iNode2; }
	const T &getWeight() const { return m_weight; }
	bool operator < (const MyArc<T> &other) const;
	SimpleEdge getSimpleEdge() const;
	SimpleArc getSimpleArc() const;

	template<class U> friend ostream& operator << (ostream& os, const MyArc<U> &edge);

private:
	int m_iNode1;
	int m_iNode2;
	T m_weight;
};

template<class U>
ostream& operator<<(ostream& os, const MyArc<U>& edge) {
	os << "(" << edge.getNode1() << ", " << edge.getNode2() << "), weight = " << edge.m_weight;
	return os;
}

// An edge1 < edge2 if the weight of edge1 is less than the weight of edge2 or if the weights are equal and 
// the lexigographic ordering says that edge1 < edge2.
template <class T>
bool MyArc<T>::operator < (const MyArc<T> &other) const
{
	return 
		m_weight < other.m_weight || 
		(
			m_weight == other.m_weight && 
			(m_iNode1 < other.m_iNode1 || (m_iNode1 == other.m_iNode1 && m_iNode2 < other.m_iNode2))
		);
}

template <class T>
SimpleEdge MyArc<T>::getSimpleEdge() const
{
	SimpleEdge simpleEdge(m_iNode1, m_iNode2);
	return simpleEdge;
}

template <class T>
SimpleArc MyArc<T>::getSimpleArc() const
{
	SimpleArc simpleArc(m_iNode1, m_iNode2);
	return simpleArc;
}

// ******************************************************************************
// ******************************************************************************
// *************************** class NodeData ***********************************
// ******************************************************************************
// ******************************************************************************

struct NodeData
{
	int m_iNode;
	double m_dData;
	NodeData(int iNode, double dData) : m_iNode(iNode), m_dData(dData)	{}
};

// Used for sorting a container in decreasing order
inline bool sortDecreasing ( const NodeData &elem1, const NodeData &elem2 )
{
   return elem1.m_dData > elem2.m_dData;
}

// Used for sorting a container in increasing order
inline bool sortIncreasing  ( const NodeData &elem1, const NodeData &elem2 )
{
   return elem1.m_dData < elem2.m_dData;
}

#endif

