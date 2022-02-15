// ***************************************************************************************
// * ArgParser.h 
// * 
// * Author: Stefan Ropke
// ***************************************************************************************

#ifndef _ARGPARSER_H
#define _ARGPARSER_H

#include <map>
#include <string>

using namespace std;

class ArgParser
{
public:
	void parseCmdLine(const string &strExpectedArgs, int argc, char *argv[ ]);
	void writeParams(ostream &os, int argc, char *argv[ ]) const;

	bool getArg(char cArg, string &strValue) const;
	bool getIntArg(char cArg, int &iValue) const;
	bool getDoubleArg(char cArg, double &dValue) const;
	bool hasArg(char cArg) const;

	int getNArgsRead() const { return (int) mapArguments.size(); }

private:
	// Key : switch character.
	// Value : argument or "" if no argument followed the switch.
	map<char, string> mapArguments;
};

#endif
