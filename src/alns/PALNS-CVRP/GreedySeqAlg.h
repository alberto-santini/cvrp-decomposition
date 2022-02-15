// Greedy sequential algorithm
#ifndef	GREEDYSEQALG_H	
#define	GREEDYSEQALG_H	

#include <vector>

using namespace std;


class CVRPInstance;
class CVRPSolution;
class TSRandom;
struct NodeDist;

class GreedySeqAlg
{
public:
	GreedySeqAlg(const CVRPInstance &instance, bool bRandomSeedCustomers = false) 
		: m_instance(instance), m_bRandomSeedCustomers(bRandomSeedCustomers)
	{}

	void go(CVRPSolution &sol, TSRandom &randGen);

private:
	int getSeedCustomer(const vector<bool> &vecIsUnassigned, 
						const vector<NodeDist> &vecNodeDepotDist,
						 int &iNextSeedNodeIdx);
	bool performBestInsertion(CVRPSolution &sol, int iRoute, vector<bool> &vecIsUnassigned, vector<int> &vecUnassigned);
	void createNewRoute(CVRPSolution &sol, const vector<NodeDist> &vecNodeDepotDist, 
						vector<bool> &vecIsUnassigned, vector<int> &vecUnassigned, int iNextSeedNodeIdx);

	const CVRPInstance &m_instance;
	bool m_bRandomSeedCustomers;
};


#endif

