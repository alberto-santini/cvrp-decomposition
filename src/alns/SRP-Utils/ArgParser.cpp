// ***************************************************************************************
// * ArgParser.cpp
// *
// * Author: Stefan Ropke
// ***************************************************************************************

#include "ArgParser.h"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

// ***************************************************************************************
// * void ArgParser::parseCmdLine(const string &strExpectedArgs, int argc, char *argv[ ])
// *
// * Parses the "command line" given as parameters <argc> and <argv> (it is expected that
// * <argc> and <argv> are the std. "main(...)" parameters), and stores the result internally
// * in the object.
// *
// * The method is inspired by the "getopts" method in Perl. The function parses
// * command lines like:
// * myprogram -f 10 -p5 -F -r help.txt
// *
// * --- input parameters ---
// * strExpectedArgs		: A string that specifies the flags or switches to look for. In the example
// *						  <strExpectedArgs> could be "f:p:r:F" which is interpretted as:
// *						  look for "-f" followed by a parameter, look for "-p" followed by a parameter
// *						  look for "-r" followed by a parameter, look for "-F" _not_ followed by a parameter
// *						  (notice that there are no ":" after F in "f:p:r:F").
// *						  The method is quite stupid so you should probably be carefull with "-" characters
// *						  in the parameters following a switch.
// *						  If the method is called with <strExpectedArgs> = "f:p:r:F" and <argc> and <argv>
// *						  is defined by the example command line shown above then an internal map is built
// *						  that would contain the following bindings:
// *
// *						  'f' -> "10"
// *						  'p' -> "5"
// *						  'F' -> ""
// *						  'r' -> "help.txt"
// *
// *						  Notice that the switches must be single-characters.
// * --- return values ---
// ***************************************************************************************

void ArgParser::parseCmdLine(const string& strExpectedArgs, int argc, char* argv[]) {
    mapArguments.clear();
    int i = 0;
    while(i < (int)strExpectedArgs.size()) {
        bool bParameterExpected = false;
        char cSwitch = strExpectedArgs[i];
        i++;
        if(i < (int)strExpectedArgs.size() && strExpectedArgs[i] == ':') {
            bParameterExpected = true;
            i++;
        }
        // search after the switch:
        for(int j = 1; j < argc; j++) {
            // Notice that we do not run into problems when argv[j] only contain one char (we might fear array-out-of-bound errors)
            // The reason is that the argv[j] strings are zero terminated.
            if((argv[j])[0] == '-' && (argv[j])[1] == cSwitch) {
                string strPar = "";
                // switch found:
                if(bParameterExpected) {
                    if((argv[j])[2] != 0)
                        strPar = string(argv[j] + 2);
                    else {
                        if(j + 1 < argc)
                            strPar = string(argv[j + 1]);
                    }
                }
                mapArguments.insert(pair<char, string>(cSwitch, strPar));
            }
        }
    }
}

void ArgParser::writeParams(ostream& os, int argc, char* argv[]) const {
    for(int j = 1; j < argc; j++)
        os << argv[j] << " ";
}

// ***************************************************************************************
// * bool ArgParser::getArg(char cArg, string &strValue) const
// *
// * gets a parameter from the parsed command line.
// *
// * --- input parameters ---
// * cArg		: The switch we are looking for
// * --- output parameters ---
// * strValue	: The parameter is returned through this parameter (nice sentence eh?).
// * --- return value: ---
// * true => the switch was found.
// * false => switch not found.
// ***************************************************************************************

bool ArgParser::getArg(char cArg, string& strValue) const {
    map<char, string>::const_iterator itMap = mapArguments.find(cArg);
    if(itMap != mapArguments.end()) {
        strValue = itMap->second;
        return true;
    } else
        return false;
}

// ***************************************************************************************
// * bool ArgParser::getArg(char cArg, int &iValue) const
// *
// * gets an integer parameter from the parsed command line.
// *
// * --- input parameters ---
// * cArg		: The switch we are looking for
// * --- output parameters ---
// * iValue		: The parameter is returned through this parameter (nice sentence eh?).
// * --- return value: ---
// * true => the parameter was found.
// * false => parameter not found.
// ***************************************************************************************

bool ArgParser::getIntArg(char cArg, int& iValue) const {
    string str;
    if(getArg(cArg, str)) {
        iValue = atoi(str.c_str());
        return true;
    } else
        return false;
}

// ***************************************************************************************
// * bool ArgParser::getArg(char cArg, double &dValue) const
// *
// * gets a double parameter from the parsed command line.
// *
// * --- input parameters ---
// * cArg		: The switch we are looking for
// * --- output parameters ---
// * dValue		: The parameter is returned through this parameter.
// * --- return value: ---
// * true => the parameter was found.
// * false => parameter not found.
// ***************************************************************************************

bool ArgParser::getDoubleArg(char cArg, double& dValue) const {
    string str;
    if(getArg(cArg, str)) {
        dValue = atof(str.c_str());
        return true;
    } else
        return false;
}

// ***************************************************************************************
// * bool ArgParser::hasArg(char cArg) const
// *
// * Checks if a specific argument was used on the command line.
// *
// * --- input parameters ---
// * cArg		: The switch we are looking for
// * --- output parameters ---
// * --- return value: ---
// * true => the parameter was found.
// * false => parameter not found.
// ***************************************************************************************

bool ArgParser::hasArg(char cArg) const {
    string str;
    if(getArg(cArg, str))
        return true;
    else
        return false;
}
