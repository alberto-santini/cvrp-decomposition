#include "CVRPInstance.h"
#include "CVRPUtils.h"

// From SRP-Utils
#include "Coordinate.h"
#include "Utils.h"

#include <fstream>
#include <limits.h>

CVRPInstance::~CVRPInstance() {}

void CVRPInstance::loadCordeauVRP(const string& filename, int nVehicles) {
    ifstream input(filename.c_str());
    if(!input.is_open())
        error("CVRPInstance::loadCordeauVRP(...)", "Specified file could not be opened!");
    string strNoPath;
    skipPath(filename, strNoPath);
    m_strName = stripExtension(strNoPath);

    int iType, iTemp;
    input >> iType;
    if(iType != 1)
        error("void CVRPInstance::LoadCordeauVRP(...)", "Unexpected type");

    input >> m_iK;
    input >> m_iN;
    input >> iTemp; // Read number of days, should be equal to one.
    if(iTemp != 1)
        error("void CVRPInstance::LoadCordeauVRP(...)", "Number of days != 1");
    input >> m_iMaxDuration;
    input >> m_iMaxLoad;

    // It is possible to override the number of vehicles by specifying it as a paramter to this method:
    if(nVehicles > 0)
        m_iK = nVehicles;

    if(m_iMaxDuration == 0)
        m_iMaxDuration = INT_MAX; // this ought to be +infinity

    if(m_iMaxLoad == 0)
        m_iMaxLoad = INT_MAX; // This ought to be +infinity

    // Read the data from the file:
    double x, y;
    Coordinate<double> coord;
    int demand;
    int serviceTime;
    int i;
    m_vecCoords.clear();
    m_vecDemand.clear();
    m_vecServiceTime.clear();
    for(i = 0; i <= m_iN; i++) {
        // cout << "Reading row " << i << endl;
        input >> iTemp; // Read customer no.
        input >> x;
        input >> y;
        coord = Coordinate<double>(x, y);
        m_vecCoords.push_back(coord);
        input >> serviceTime;
        m_vecServiceTime.push_back(serviceTime);
        input >> demand;
        m_vecDemand.push_back(demand);
        input >> iTemp; // frequency of visit (ignored).
        input >> iTemp; // number of possible visit combinations (ignored).
        int iDummy;
        for(int j = 0; j < iTemp; j++)
            input >> iDummy; // read and discard visit combination.
    }
    m_vecCoords.push_back(m_vecCoords.front());
    m_vecServiceTime.push_back(m_vecServiceTime.front());
    m_vecDemand.push_back(m_vecDemand.front());
    input.close();
    initDistMatrix();
    findNearestNodes();
}

void CVRPInstance::loadTailard(const string& filename) {
    ifstream input(filename.c_str());
    if(!input.is_open())
        error("CVRPInstance::loadTailard(...)", "Specified file could not be opened!");
    string strNoPath;
    skipPath(filename, strNoPath);
    m_strName = stripExtension(strNoPath);

    input >> m_iN;
    m_iK = -1;
    m_bNVehiclesFixed = false;
    double dTemp;
    input >> dTemp; // Read Best known objective
    input >> m_iMaxLoad;
    m_iMaxDuration = INT_MAX; // this ought to be +infinity

    if(m_iMaxLoad == 0)
        m_iMaxLoad = INT_MAX; // This ought to be +infinity

    double x, y;
    int demand;
    int i, iTemp;
    m_vecCoords.clear();
    m_vecDemand.clear();
    m_vecServiceTime.clear();
    // Read depot coordinates
    input >> x >> y;
    m_vecCoords.push_back(Coordinate<double>(x, y));
    m_vecDemand.push_back(0);
    m_vecServiceTime.assign(m_iN + 2, 0);
    for(i = 1; i <= m_iN; i++) {
        input >> iTemp; // Read customer no.
        input >> x;
        input >> y;
        m_vecCoords.push_back(Coordinate<double>(x, y));
        input >> demand;
        m_vecDemand.push_back(demand);
    }
    m_vecCoords.push_back(m_vecCoords.front());
    m_vecDemand.push_back(m_vecDemand.front());
    input.close();
    initDistMatrix();
    findNearestNodes();
}

/***********************************************************************
 * void CVRPInstance::loadCVRP_GWKC(const string &filename, int nVehicles)
 *
 * Reads a CVRP in the Golden, Wasil, Kelly and Chao format.
 ***********************************************************************/

void CVRPInstance::loadCVRP_GWKC(const string& filename, int nVehicles) {
    ifstream input(filename.c_str());
    if(!input.is_open())
        error("CVRPInstance::loadCVRP_GWKC(...)", "Specified file could not be opened!");
    string strNoPath;
    skipPath(filename, strNoPath);
    m_strName = stripExtension(strNoPath);

    int iTemp;
    double dTemp;

    input >> m_iN;
    input >> m_iMaxLoad;
    input >> m_iMaxDuration;

    // This next entry should always be zero:
    input >> iTemp;
    if(iTemp != 0)
        error("void VRPPDProblem::loadCVRP_GWKC(...)", "Fourth entry in file should be zero!");

    // Read upper bound.
    input >> dTemp;

    if(nVehicles > 0)
        m_iK = nVehicles;
    else
        m_iK = m_iN;

    if(m_iMaxDuration == 0)
        m_iMaxDuration = INT_MAX; // this ought to be +infinity

    if(m_iMaxLoad == 0)
        m_iMaxLoad = INT_MAX; // This ought to be +infinity

    // Read the data from the file:
    double x, y;
    Coordinate<double> coord;
    int demand;
    int serviceTime = 0;
    int i;
    m_vecCoords.clear();
    m_vecDemand.clear();
    m_vecServiceTime.clear();
    for(i = 0; i <= m_iN; i++) {
        // cout << "Reading row " << i << endl;
        input >> x;
        input >> y;
        coord = Coordinate<double>(x, y);
        m_vecCoords.push_back(coord);
        // Notice that for some of the original problems, the demand was not listed for the depot.
        // I have changed these files such that the demand (0) is always listed for the depot.
        // This makes the files much easier to read.
        input >> demand;
        m_vecDemand.push_back(demand);
        m_vecServiceTime.push_back(serviceTime);
    }
    assert(m_iN > 0);
    m_vecCoords.push_back(m_vecCoords.front());
    m_vecServiceTime.push_back(m_vecServiceTime.front());
    m_vecDemand.push_back(m_vecDemand.front());

    input.close();
    initDistMatrix();
    findNearestNodes();
}

bool CVRPInstance::findNextOccurence(istream& is, const string& strSearchString, string& strFound) {
    // cout << "findNextOccurence: " << strSearchString << endl;
    string str;
    bool bFound = false;

    while(!bFound && !is.eof()) {
        getline(is, str);
        // skip comments:
        if(str[0] != '#') {
            int iPos = (int)str.find(strSearchString, 0);
            if(iPos != string::npos) {
                iPos = (int)str.find_first_not_of("\t :", iPos + strSearchString.length());
                // If the last part of the string doesn't contain any characters different from tab, space and colon then
                // find_first_not_of returns npos. This can happen if the string only contains the search string.
                if(iPos == string::npos)
                    strFound = "";
                else
                    strFound = str.substr(iPos);
                // cout << "*** Match found! ***" << endl;
                bFound = true;
                // cout << "Return string: " << strFound << endl;
            }
        }
    }
    return bFound;
}

/***********************************************************************
 * void CVRPInstance::loadCVRP_TSPLIB(const string &filename, int iNVehicles)
 *
 * Reads CVRP instance in the "TSPLIB" format (the format used for the instances
 * solved by the exact community).
 * see www.branchandcut.org
 ***********************************************************************/

void CVRPInstance::loadCVRP_TSPLIB(const string& filename, int iNVehicles, bool bNVehicleFixed) {
    m_bNVehiclesFixed = bNVehicleFixed;
    //	cout << "Loading CVRP problem..." << endl;

    ifstream input(filename.c_str());
    if(!input.is_open()) {
        string strError = "Specified file could not be opened! filename=" + filename;
        error("void CVRPInstance::loadCVRP_TSP(...)", strError);
    }
    string strNoPath;
    skipPath(filename, strNoPath);
    m_strName = stripExtension(strNoPath);

    m_iK = 0;
    m_iN = 0;
    // No duration in these instances
    m_iMaxDuration = INT_MAX;

    string strTemp;
    if(bNVehicleFixed) {
        if(!findNextOccurence(input, "VEHICLES", strTemp)) {
            if(iNVehicles <= 0)
                error("CVRPInstance::loadCVRP_TSP(...)", "#Vehicles not specified by caller or in instance file");
            m_iK = iNVehicles; // number of vehicles not specified.
            input.clear();
            input.seekg(0); // reset the stream
        } else
            m_iK = atoi(strTemp.c_str());
    }

    if(!findNextOccurence(input, "DIMENSION", strTemp))
        error("CVRPInstance::loadCVRP_TSP(...)", "DIMENSION key not found in file.");
    m_iN = atoi(strTemp.c_str()) - 1;

    if(!findNextOccurence(input, "CAPACITY", strTemp))
        error("CVRPInstance::loadCVRP_TSP(...)", "CAPACITY key not found in file.");
    m_iMaxLoad = atoi(strTemp.c_str());

    if(!findNextOccurence(input, "MAXDIST", strTemp)) {
        // cout << "MAXDIST not found" << endl;
        input.clear();
        input.seekg(0);
    } else {
        error("CVRPInstance::loadCVRP_TSP(...)", "MAXDIST specified, this is not supported.");
    }

    // cout << "iN = " << m_iN << endl;
    if(!findNextOccurence(input, "NODE_COORD_SECTION", strTemp))
        error("CVRPInstance::loadCVRP_TSP(...)", "NODE_COORD_SECTION key not found in file.");

    // Read the data from the file:
    m_vecDemand.clear();
    m_vecCoords.clear();

    double x, y;
    int demand;
    int i, j, iId;
    for(i = 0; i < m_iN + 1; i++) {
        // cout << "Reading row " << i << endl;
        input >> iId;
        input >> x;
        input >> y;
        m_vecCoords.push_back(Coordinate<double>(x, y));
    }
    if(!findNextOccurence(input, "DEMAND_SECTION", strTemp))
        error("CVRPInstance::loadCVRP_TSP(...)", "DEMAND_SECTION key not found in file.");
    for(i = 0; i < m_iN + 1; i++) {
        // cout << "Reading row " << i << endl;
        input >> iId; // Read customer no.
        input >> demand;
        m_vecDemand.push_back(demand);
        // We are reading a CVRP problem, so some VRP2DBin fields should be set to "dummy values".
    }
    m_vecServiceTime.assign(m_iN + 2, 0);
    // Node n+1 is a copy of the depot.
    m_vecDemand.push_back(m_vecDemand.front());
    m_vecCoords.push_back(m_vecCoords.front());
    // we ignore the depot section, we assume that node 1 is the depot.
    input.close();

    m_matDists.resizeAndClear(m_iN + 2, 0);
    for(i = 0; i <= m_iN + 1; i++) {
        for(j = 0; j <= m_iN + 1; j++) {
            double dDist = calcDist(m_vecCoords[i], m_vecCoords[j]);
            m_matDists.setElement(i, j, round(dDist));
        }
    }
    findNearestNodes();
}

void CVRPInstance::initDistMatrix() {
    m_matDists.resizeAndClear(m_iN + 2, 0);

    int i, j;
    for(i = 0; i < m_iN + 2; i++) {
        for(j = 0; j < m_iN + 2; j++) {
            m_matDists.setElement(i, j, calcDist(m_vecCoords[i], m_vecCoords[j]));
        }
    }
}

void CVRPInstance::writeCordeauInstance() {
    ofstream ofs((m_strName + "_C.txt").c_str());
    ofs << 1 << " " << m_iN << " " << m_iN << " " << 1 << "\n";
    if(m_iMaxDuration == INT_MAX)
        ofs << 0;
    else
        ofs << m_iMaxDuration;
    ofs << " " << m_iMaxLoad << "\n";
    int i;
    for(i = 0; i <= m_iN; i++) {
        ofs << i << " " << m_vecCoords[i].getX() << " " << m_vecCoords[i].getY() << " " << m_vecServiceTime[i] << " " << m_vecDemand[i];
        if(i == 0)
            ofs << " 0 0 \n";
        else
            ofs << " 1 1 1\n";
    }
    ofs.flush();
    ofs.close();
}

void CVRPInstance::setProblem(int iN, const vector<Coordinate<double>>& vecCoords, const vector<int>& vecDemand, const vector<int>& vecServiceTime, int iNVeh,
                              bool bNVehFixed, int iMaxLoad, int iMaxDur) {
    m_iN = iN;
    assert(iN == (int)vecCoords.size() - 2);
    assert(vecCoords.size() == vecDemand.size());
    assert(vecCoords.size() == vecServiceTime.size());
    m_vecCoords = vecCoords;
    m_vecDemand = vecDemand;
    m_vecServiceTime = vecServiceTime;
    m_iK = iNVeh;
    m_bNVehiclesFixed = bNVehFixed;
    m_iMaxLoad = iMaxLoad;
    m_iMaxDuration = iMaxDur;
}

void CVRPInstance::setDistMatrix(const Matrix<double>& matDists) {
    if(matDists.getDim() != m_iN + 2)
        error("CVRPInstance::setDistMatrix(...)", "Distance matrix dimension mis-match");
    m_matDists = matDists;
    findNearestNodes();
}

void CVRPInstance::findNearestNodes() {
    m_vecNearestNodes.resize(m_iN + 1, vector<int>(m_iN + 1));
    int i, j;
    for(i = 0; i <= m_iN; ++i) {
        vector<NodeData> vecNodeData;
        for(j = 0; j <= m_iN; ++j) {
            vecNodeData.push_back(NodeData(j, m_matDists.getElement(i, j)));
        }
        sort(vecNodeData.begin(), vecNodeData.end(), sortIncreasing);
        for(j = 0; j <= m_iN; ++j)
            m_vecNearestNodes[i][j] = vecNodeData[j].m_iNode;
    }
}