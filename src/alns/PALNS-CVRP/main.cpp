#include "CVRPInstance.h"
#include "CVRPSolution.h"
#include "GreedySeqAlg.h"
#include "PALNS-CVRP.h"
#include "Regret2Alg.h"
#include "SeedAlg.h"

// From SRP-Utils
#include "ArgParser.h"
#include "HashNumberGenerator.h"
#include "ParFileParser.h"

#include <fstream>
#include <iostream>

void printHelp() {
    cout << "-f <arg>    : Load problem given by <arg>." << endl
         << "-t <arg>    : problem type to load: " << endl
         << "              1 = JFC CVRP (default)" << endl
         << "              2 = GWKC" << endl
         << "              3 = TSPLIB format (used by exact algs.)" << endl
         << "              4 = Tailard" << endl
         << "              5 = TSPLIB format, but #vehicles is free (for Uchoa et al. instances)" << endl
         << "-p <arg>    : Use parameter file given by <arg>." << endl
         << "-P <arg>    : override parameters from par file with <arg>." << endl
         << "-n <arg>    : number of threads to use." << endl
         << "-v <arg>	: Number of vehicles (for type 3 instances)" << endl
         << "-r <arg>	: Number of retries" << endl
         << "-h          : This help" << endl
         << "-H          : Add header line to output file" << endl
         << "-o <arg>    : prefix of summary files" << endl
         << endl
         << "Example:" << endl
         << "..." << endl;
}

void writeHeaderLine(fstream& fsInfoOut) { cout << "name,method,sol. cost,time (s)" << endl; }

int main(int argc, char* argv[]) {
    ArgParser argParser;
    argParser.parseCmdLine("f:hHn:o:p:P:r:t:v:", argc, argv);

    if(argParser.hasArg('h')) {
        printHelp();
        exit(0);
    }

    int iNThreads = 2;
    argParser.getIntArg('n', iNThreads);

    // *** Write header line if desired ***
    bool bHeader = argParser.hasArg('H');

    string strResultPrefix = "";
    argParser.getArg('o', strResultPrefix);

    string strProblemFileName;
    // *** get problem file name ***
    if(!argParser.getArg('f', strProblemFileName)) {
        printHelp();
        exit(0);
    }
    string strParFile = "params.txt";
    argParser.getArg('p', strParFile);
    ParFileParser pfp;
    pfp.loadParFile(strParFile);

    if(argParser.hasArg('P')) {
        string strNewParFileParams;
        argParser.getArg('P', strNewParFileParams);
        pfp.addParamsFromString(strNewParFileParams);
    }

    int iNVehicles = -1;

    int iNRetries = 1;
    argParser.getIntArg('r', iNRetries);
    // *** get problem type ***
    int iProblemType = 1;
    argParser.getIntArg('t', iProblemType);

    // *** load problem ***
    CVRPInstance instance;
    switch(iProblemType) {
    case 1: instance.loadCordeauVRP(strProblemFileName, -1); break;
    case 2: instance.loadCVRP_GWKC(strProblemFileName, -1); break;
    case 3:
        argParser.getIntArg('v', iNVehicles);
        instance.loadCVRP_TSPLIB(strProblemFileName, iNVehicles);
        break;
    case 4: instance.loadTailard(strProblemFileName); break;
    case 5: instance.loadCVRP_TSPLIB(strProblemFileName, -1, false); break;
    default:
        cout << "Unknown problem type" << endl;
        printHelp();
        exit(0);
    }
    //	cout << "Instance contains " << instance.getN() << " customers." << endl;
    FixedHashNumberGenerator::init(10 * instance.getN());

    // instance.writeCordeauInstance();
    // exit(0);

    CVRPSolution sol(instance);

    TSRandom randGen;
    GreedySeqAlg greedySeqAlg(instance);
    double dTimeStart = elapsedSeconds();
    greedySeqAlg.go(sol, randGen);
    //	cout << "Cost of solution found by greedy sequential algorithm: " << sol.getCost() << endl;
    //	cout << "Time spent by greedy sequential algorithm: " << elapsedSeconds() - dTimeStart << endl;

    /*
    int i;
    for (i=0; i < 10; ++i)
    {
        sol.clear();
        GreedySeqAlg greedySeqAlg2(instance, true);
        greedySeqAlg2.go(sol, randGen);
        cout << "Cost of solution found by randomized greedy sequential algorithm: " << sol.getCost() << endl;
    }
    */

    int iNRoutes;
    if(instance.getNVehiclesFixed())
        iNRoutes = instance.getNVehicles();
    else
        // Use the number of routes found by the greedy algorithm if the number of routes is not fixed.
        iNRoutes = sol.getNRoutes();

    /*TSRandom randGen;
    // Regret-2 algortihm
    Regret2Alg regret2Alg(instance, randGen);
    sol.clear(iNRoutes);
    dTimeStart = elapsedSeconds();
    SeedAlg seedAlg(instance);
    seedAlg.assignSeedCustomers(sol);
    regret2Alg.go(sol);
    cout << "Cost of solution found by regret-2 algorithm: " << sol.getCost() << endl;
    cout << "Time spent by regret-2 algorithm: " << elapsedSeconds() - dTimeStart << endl;
    */
    dTimeStart = elapsedSeconds();
    double dWallClock = wallClock();
    PALNS_CVRP alns;
    alns.setResultPrefix(strResultPrefix);
    alns.go(sol, iNThreads, iNRetries, bHeader, pfp);
    // cout << "Time spent by ALNS (CPU): " << elapsedSeconds() - dTimeStart << endl;
    // cout << "Time spent by ALNS (wall clock): " << wallClock() - dWallClock << endl;
    cout << sol.getCost() << endl;
}
