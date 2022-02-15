#ifndef SEEDALG_H
#define SEEDALG_H

class CVRPSolution;
class CVRPInstance;

class SeedAlg
{
public:
	SeedAlg(const CVRPInstance &instance) 
		: m_instance(instance)   
	{}

	void go(CVRPSolution &sol);
	void assignSeedCustomers(CVRPSolution &sol);

private:
	const CVRPInstance &m_instance;

};

#endif
