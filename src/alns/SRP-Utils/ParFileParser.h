#ifndef PARFILEPARSER_H
#define PARFILEPARSER_H

#include <string>
#include <map>

using namespace std;

class ParFileParser
{
public:
	ParFileParser() {};
	bool loadParFile(const string &strFileName);
	void addParamsFromString(const string &strPars);
	void outputParams(ostream &os) const;

	bool hasParameter(const string &strParName) const;
	bool getParameter(const string &strParName, string &strRes) const;
	bool getBoolParameter(const string &strParName, bool &bValue) const;
	bool getIntParameter(const string &strParName, int &iValue) const;
	bool getDoubleParameter(const string &strParName, double &dValue) const;
	void setParameter(const string& strParName, string value);

private:
	bool scanLine(const string &str, string &strKey, string &strValue);
	void parseError(const string &str);

	map<string, string> m_mapParameters;

};

#endif
