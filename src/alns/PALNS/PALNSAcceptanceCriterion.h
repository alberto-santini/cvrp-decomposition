#ifndef PALNS_ACCEPTANCE_CRITERIA_H
#define PALNS_ACCEPTANCE_CRITERIA_H

#include "Utils.h"
#include "PALNSParameters.h"
#include "../PALNS-CVRP/CVRPSolution.h"

/*! \brief This class models a base, generic acceptance criterion */
struct AcceptanceCriterion {
    /*! Reference to parameters object */
    PALNSParameters& m_par;
    
    /*! Basic constructor */
    AcceptanceCriterion(PALNSParameters& par) : m_par(par) {}
    
    // Methods each subclass should implement:
    
    /*! This method initialise method-specific parameters */
    virtual void initialise() = 0;
    
    /*! This method updates the acceptance criterion's parameters
     */
    virtual void updateParameters(double elapsedSec) = 0;
    
    /*! This method returns true iff the solution should be accepted according to the acceptance criterion
     *
     *  @param bestObj is the objective value of the best solution (globally) found so far
     *  @param newObj is the objective value of the new solution (which we must determine wether to accept or not)
     *  @param eps is an epsilon to be used in numerical comaprisons
     */
    virtual bool shouldAccept(double bestObj, double newObj, double eps) = 0;
};

#endif