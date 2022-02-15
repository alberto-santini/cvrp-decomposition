#ifndef PALNS_HELPERS_H
#define PALNS_HELPERS_H

#include "../PALNS-CVRP/CVRPSolution.h"

// Solution tracker that can be used to collect data about the solution process
struct AbstractSolTrackerCallBack {
    // Gives a pointer to the solution, and says if:
    //  - it's been accepted by the acceptance criterion
    //  - it improves on the current solution
    //  - it improves on the global best solution
	virtual void callback(const CVRPSolution& sol, bool bAccepted, bool bImproved, bool bNewBest) = 0;
    
    // Clones the callback object
	virtual AbstractSolTrackerCallBack* clone() const = 0;
	virtual ~AbstractSolTrackerCallBack() = default;
};

// Destroy method
struct AbstractDestroyMethod {
	virtual void destroySolution(CVRPSolution &sol, TSRandom &randGen) = 0;
	virtual void setSolTracker(const AbstractSolTrackerCallBack*) {};
	virtual AbstractDestroyMethod* clone() const = 0;
	virtual ~AbstractDestroyMethod() = default;
};

// Repair method
struct AbstractRepairMethod {
	virtual void repairSolution(CVRPSolution &sol, TSRandom &randGen) = 0;
	virtual void setSolTracker(const AbstractSolTrackerCallBack*) {};
	virtual AbstractRepairMethod* clone() const = 0;
	virtual ~AbstractRepairMethod() = default;
};

// Constructive heuristic to get an initial solution
struct InitialSolutionCreator {
	virtual void createInitialSolution(CVRPSolution &sol, TSRandom &randGen) = 0;
	virtual ~InitialSolutionCreator() = default;
};

#endif