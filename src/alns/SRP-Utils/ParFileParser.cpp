// #define USE_GRETA

#include "ParFileParser.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

#ifdef USE_GRETA
#include "regexpr2.h"
using namespace regex;

bool ParFileParser::loadParFile(const string& strFileName) {
    fstream fsParFile;
    fsParFile.open(strFileName.c_str());
    if(!fsParFile.good())
        return false;

    char buffer[10000];
    int count = 9900;
    // fsParFile.read(buffer, count);
    while(!fsParFile.eof()) {
        fsParFile.getline(buffer, count);
        string strParFileContent(buffer);

        // Remove comments from the line (comments start with "//").
        rpattern removeComments("//.*", "", GLOBAL);
        subst_results results;
        removeComments.substitute(strParFileContent, results);

        // cout << strParFileContent << endl;

        // See if the line matches a definition
        // the following pattern matches something like
        // * = {*}
        // eg.: abc = { def }
        rpattern pat("(\\w*)\\s*=\\s*\\{\\s*(\\w*)\\s*\\}");
        match_results matchResults;
        match_results::backref_type br = pat.match(strParFileContent, matchResults);
        if(br.matched) {
            const vector<match_results::backref_type>& vec = matchResults.all_backrefs();
            /*cout << "match success: " << br << endl;
            for (int i=0; i <  (int) vec.size(); i++)
            {
                cout << i << " : " << vec[i].str() << endl;
            }
            */
            string strKey = vec[1].str();
            string strValue = vec[2].str();

            m_mapParameters.insert(make_pair(strKey, strValue));
        } else {
            // Do nothing, try next line
            // cout << "match failed!" << endl;
        }
    }
    return true;
}

#else

bool ParFileParser::loadParFile(const string& strFileName) {
    fstream fsParFile;
    fsParFile.open(strFileName.c_str());
    if(!fsParFile.good())
        return false;

    char buffer[10000];
    int count = 9900;
    // fsParFile.read(buffer, count);
    while(!fsParFile.eof()) {
        fsParFile.getline(buffer, count);
        string strParFileContent(buffer);
        string strKey, strValue;
        if(scanLine(strParFileContent, strKey, strValue))
            m_mapParameters.insert(make_pair(strKey, strValue));
    }
    return true;
}

#endif

/***********************************************************************
 * void ParFileParser::addParamsFromString(const string &strPars)
 *
 * Add one or more parameters. The extra parameters are passed in strPars.
 * The new parameters overrides existing ones.
 *
 * --- input parameters ---
 * strPars		:
 * --- output parameters ---
 * --- return values ---
 * --- protection status ---
 * public
 ***********************************************************************/

typedef basic_string<char>::size_type size_type;

void ParFileParser::addParamsFromString(const string& strPars) {
    bool bDone = false;
    size_type pos = 0;
    static const basic_string<char>::size_type npos = -1;
    while(true) {
        bDone = true;

        pos = strPars.find_first_not_of(" \t", pos);
        if(pos == npos) {
            break;
        }

        size_type pos2 = strPars.find_first_of("= \t", pos + 1);
        if(pos2 == npos) {
            break;
        }

        string strKey = strPars.substr(pos, pos2 - pos);
        size_type pos3 = strPars.find_first_of("{", pos2 + 1);
        if(pos3 == npos) {
            break;
        }

        pos3 = strPars.find_first_not_of(" \t", pos3 + 1);
        if(pos3 == npos) {
            break;
        }

        size_type pos4 = strPars.find_first_of("} \t", pos3 + 1);
        if(pos4 == npos) {
            break;
        }

        size_type pos5 = strPars.find_first_of("}", pos4);
        if(pos5 == npos) {
            break;
        }

        string strVal = strPars.substr(pos3, pos4 - pos3);

        // Make sure the map doesn't already contain this parameter. This is necessary because "insert" doesn't overwrite.
        m_mapParameters.erase(strKey);
        m_mapParameters.insert(make_pair(strKey, strVal));

        pos = pos5 + 1;
    }
}

/***********************************************************************
 * bool ParFileParser::scanLine(const string &str, string &strKey, string &strValue)
 *
 * Scans a line, looking for a key and value pair of the form
 * key = { value }
 *
 * The method terminates the program if a bad line is read.
 *
 * --- input parameters ---
 * str			: the line to scan
 * --- output parameters ---
 * strKey		: the key found
 * strValue		: the value found
 * --- return values ---
 * true if a (key, value) pair was found, false if an empty line or a comment was read.
 * --- protection status ---
 * private
 ***********************************************************************/

bool ParFileParser::scanLine(const string& str, string& strKey, string& strValue) {
    static const basic_string<char>::size_type npos = (basic_string<char>::size_type) - 1;
    int i = 0;
    int strSize = (int)str.size();

    i = (int)str.find_first_not_of(" \t\n\r", i);
    if(i == npos)
        return false;
    if(i < strSize && str[i] == '/') {
        if(i + 1 < strSize && str[i + 1] == '/')
            // The rest of the string is a comment.
            return false;
        else
            parseError(str);
    } else {
        // Read the key.
        int iEndOfKey = (int)str.find_first_of(" \t=", i);
        if(iEndOfKey == npos)
            parseError(str);
        strKey = str.substr(i, iEndOfKey - i);
        // Find "=" symbol.
        i = (int)str.find_first_of('=', iEndOfKey);
        if(iEndOfKey == npos)
            parseError(str);
        // Find "{" symbol.
        i = (int)str.find_first_of('{', i);
        if(iEndOfKey == npos)
            parseError(str);
        // skip the '{' symbol.
        i++;
        // Skip whitespace.
        i = (int)str.find_first_not_of(" \t", i);
        // Read the value
        // Find closing "}".
        int iEndOfValue = (int)str.find_first_of("}", i);
        if(iEndOfValue == npos)
            parseError(str);
        // Skip white space before closing "}".
        iEndOfValue = (int)str.find_last_not_of(" \t}", iEndOfValue) + 1;
        strValue = str.substr(i, iEndOfValue - i);
        // Skip whitespace.
        i = (int)str.find_first_not_of(" \t", iEndOfValue);
        // make sure that we got '}';
        if(str[i] != '}')
            parseError(str);
    }
    // key and value read.
    return true;
}

void ParFileParser::parseError(const string& str) {
    cout << "Error while reading parameter file" << endl;
    cout << "malformed line \"" << str << "\"" << endl;
    exit(0);
}

bool ParFileParser::getParameter(const string& strParName, string& strRes) const {
    map<string, string>::const_iterator iter = m_mapParameters.find(strParName);
    if(iter != m_mapParameters.end()) {
        strRes = iter->second;
        return true;
    } else {
        return false;
    }
}

bool ParFileParser::hasParameter(const string& strParName) const {
    map<string, string>::const_iterator iter = m_mapParameters.find(strParName);
    if(iter != m_mapParameters.end())
        return true;
    else
        return false;
}

bool ParFileParser::getBoolParameter(const string& strParName, bool& bValue) const {
    string strDummy;
    bool bOk = getParameter(strParName, strDummy);
    if(bOk)
        bValue = atoi(strDummy.c_str()) != 0;
    return bOk;
}

bool ParFileParser::getIntParameter(const string& strParName, int& iValue) const {
    string strDummy;
    bool bOk = getParameter(strParName, strDummy);
    if(bOk)
        iValue = atoi(strDummy.c_str());
    return bOk;
}

bool ParFileParser::getDoubleParameter(const string& strParName, double& dValue) const {
    string strDummy;
    bool bOk = getParameter(strParName, strDummy);
    if(bOk)
        dValue = atof(strDummy.c_str());
    return bOk;
}

void ParFileParser::outputParams(ostream& os) const {
    map<string, string>::const_iterator iter = m_mapParameters.begin();
    while(iter != m_mapParameters.end()) {
        os << iter->first << " = { " << iter->second << " }" << endl;
        iter++;
    }
}

void ParFileParser::setParameter(const string& strParName, string value) { m_mapParameters[strParName] = std::move(value); }