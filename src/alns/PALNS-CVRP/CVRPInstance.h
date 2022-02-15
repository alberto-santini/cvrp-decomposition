#ifndef CVRPPROB_H
#define CVRPPROB_H

// From SRP-Utils
#include "Coordinate.h"
#include "Matrix.h"

#include <string>
#include <vector>

using namespace std;

class DistCache;
class Customer;

class CVRPInstance
{
public:
	CVRPInstance() { m_bNVehiclesFixed = false; }
	~CVRPInstance();

	int getN() const									{ return m_iN; }
	const Matrix<double> &getDists() const				{ return m_matDists; }
	inline int getMaxLoad() const						{ return m_iMaxLoad; }
	int getMaxDuration() const							{ return m_iMaxDuration; }
	inline int getDemand(int iCustId) const				{ return m_vecDemand[iCustId]; }
	int getServiceTime(int iCustId) const				{ return m_vecServiceTime[iCustId]; }
	Coordinate<double> getCoord(int iCustId) const		{ return m_vecCoords[iCustId]; }
	bool getNVehiclesFixed() const						{ return m_bNVehiclesFixed; }
	int getNVehicles() const							{ return m_iK; }
	double getDist(int i, int j) const					{ return m_matDists.getElement(i,j); }
	const string &getName() const						{ return m_strName; }
	const vector<int> &getNearestNodes(int iNode) const	{ return m_vecNearestNodes[iNode]; }
	// For now we only handle symmetric CVRP instances:
	bool isSymmetric() const							{ return true; }
	string getInstanceInfo() const						{ return int2String(m_iN)+ ","; }
	string getInstanceInfoHeader() const				{ return "n,"; }
	int getInstanceSize() const							{ return m_iN; }

	void loadCordeauVRP(const string &filename, int nVehicles);
	void loadCVRP_GWKC(const string &filename, int nVehicles);
	void loadCVRP_TSPLIB(const string &filename, int iNVehicles = -1, bool bNVehicleFixed = true);
	void loadTailard(const string &filename);

	void writeCordeauInstance();

	// methods for "loading" the problem as an alternative to the "load from file" methods
	void setProblem(int iN, const vector<Coordinate<double> > &vecCoords, const vector<int> &vecDemand, 
					const vector<int> &vecServiceTime, int iNVeh, bool bNVehFixed, int iMaxLoad, int iMaxDur);
	void setDistMatrix(const Matrix<double> &matDists);
			
private:
	void initDistMatrix();
	bool findNextOccurence(istream &is, const string &strSearchString, string &strFound);
	void findNearestNodes();

	vector<Coordinate<double> > m_vecCoords;
	vector<int> m_vecDemand;
	vector<int> m_vecServiceTime;

	// m_vecNearestNodes[i] : All customers, sorted in order of increasing distance to node i.
	vector<vector<int> > m_vecNearestNodes;

	// m_iN : Number of customers.
	int m_iN;
	// m_iK : Number of vehicles.
	int m_iK;
	bool m_bNVehiclesFixed;

	int m_iMaxLoad;
	int m_iMaxDuration;

	Matrix<double> m_matDists;
	string m_strName;
};

#endif
