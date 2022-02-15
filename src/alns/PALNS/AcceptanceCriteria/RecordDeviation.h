#ifndef RECORD_DEVIATION_H
#define RECORD_DEVIATION_H

#include "../PALNSAcceptanceCriterion.h"
#include "../../PALNS-CVRP/CVRPSolution.h"

/*! \brief This class models the Record Deviation acceptance criterion for ALNS */
struct RecordDeviation : public AcceptanceCriterion {
    /*! Deviation */
    double m_dDeviation;
    
    /*! Deviation step, when decrease is linear */
    double m_dDeviationStep;
    
    /*! Basic constructor */
    RecordDeviation(PALNSParameters& par) : AcceptanceCriterion(par) {}
    
    /*! This method initialise method-specific parameters */
    void initialise();
    
    /*! This method updates the acceptance criterion's parameters, based on running info
     */
    void updateParameters(double elapsedSec);

   /*! This method returns true iff the solution should be accepted according to the acceptance criterion
    *
    *  @param bestObj is the objective value of the best solution (globally) found so far
    *  @param newObj is the objective value of the new solution (which we must determine wether to accept or not)
    *  @param eps is an epsilon to be used in numerical comaprisons
    */
   bool shouldAccept(double bestObj, double newObj, double eps);
};

#endif