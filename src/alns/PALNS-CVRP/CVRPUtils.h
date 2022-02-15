#ifndef CVRPUTILS
#define CVRPUTILS

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