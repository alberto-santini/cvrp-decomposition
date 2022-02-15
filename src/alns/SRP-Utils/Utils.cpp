#include "Utils.h"
#include "Coordinate.h"
#include "NSMatrix.h"
#include "VectorUtils.h"

#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>

#if(defined _MSC_VER || defined __CYGWIN__) && !defined COMP_LANG_C_PEDANTIC
// For wall clock function
#include <sys/timeb.h>
#include <time.h>
#else
#include <sys/time.h>
#endif

using namespace std;

unsigned int Random::_h48 = 0;
unsigned int Random::_l48 = 0x330E;
// int Random::count;

void error(string where, string what) {
    cout << "<" << where << "> : " << what << "." << endl;
    /*
    cout << "Press 1 and enter" << endl;
    int x;
    cin >> x;
    */
    exit(1);
}

/***********************************************************************
 * double round(double val, unsigned int iPrecision)
 *
 * Rounds a double precision number to <iPrecision> decimals. This
 * function does what it is supposed to do,  but it might loose some
 * precision in the conversion process. Could it be written in a smarter way?
 ***********************************************************************/

double round(double val, unsigned int iPrecision) { return floor(val * pow(10.0, (int)iPrecision) + 0.5) / pow(10.0, (int)iPrecision); }

double round(double val) { return floor(val + 0.5); }

/***********************************************************************
 * double trunc(double val, unsigned int iPrecision)
 *
 * Truncates a double precision number to <iPrecision> decimals. This
 * function does what it is supposed to do,  but it might loose some
 * precision in the conversion process. Could it be written in a smarter way?
 ***********************************************************************/

double trunc(double val, unsigned int iPrecision) { return floor(val * pow(10.0, (int)iPrecision)) / pow(10.0, (int)iPrecision); }

double trunc(double val) { return floor(val); }

/***********************************************************************
 * double roundUp(double val, unsigned int iPrecision)
 *
 * rounds <val> to <iPrecision> decimals, the method always rounds up. This
 * function does what it is supposed to do,  but it might loose some
 * precision in the conversion process. Could it be written in a smarter way?
 ***********************************************************************/

double roundUp(double val, unsigned int iPrecision) { return ceil(val * pow(10.0, (int)iPrecision)) / pow(10.0, (int)iPrecision); }

double roundUp(double val) { return ceil(val); }

// It is assumed that a and b are positive. If we want to lift that assumption we can make a
// "gateway" function that passes absolute values on to the gcd function. Maybe the method
// below actually works for negative numbers - I have not looked into this. However, one has
// to be careful with the % operator. In older versions of C++ its behaviour is implementation
// dependent if I remember correctly.
int gcd(int a, int b) {
    if(b == 0)
        return a;
    else
        return gcd(b, a % b);
}

double calcDist(double dX1, double dY1, double dX2, double dY2) { return sqrt((dX1 - dX2) * (dX1 - dX2) + (dY1 - dY2) * (dY1 - dY2)); }

double calcDist(const Coordinate<double>& c1, const Coordinate<double>& c2) { return calcDist(c1.getX(), c1.getY(), c2.getX(), c2.getY()); }

/* for more accurate timing under Windows NT */
/* TODO: better error checking */
/* TODO: fix potential overflow converting FILETIME to __int64 under MSVC */
#if(defined _MSC_VER || defined __CYGWIN__) && !defined COMP_LANG_C_PEDANTIC
/* wishful thinking */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#ifdef _MSC_VER
/* MSVC doesn't support conversion from unsigned __int64
   to double, so __int64 it is */
typedef __int64 my_clock_64;
#else
typedef unsigned long long my_clock_64;
#endif
clock_t my_clock(void) {
    static HANDLE hProcess = 0;
    static int use_GetProcessTimes = 0;
    FILETIME Ignored1, Ignored2, KernelTime, UserTime;
    clock_t current_time;

    if(hProcess == 0) {
        hProcess = GetCurrentProcess();
        if(GetProcessTimes(hProcess, &Ignored1, &Ignored2, &KernelTime, &UserTime)) {
            use_GetProcessTimes = 1;
            // fprintf(stderr, "my_clock: using GetProcessTimes() for timing\n");
        } else { // fprintf(stderr, "my_clock: using clock() for timing\n");
        }
    }

    if(use_GetProcessTimes) {
        if(GetProcessTimes(hProcess, &Ignored1, &Ignored2, &KernelTime, &UserTime)) {
            my_clock_64 Kernel = ((my_clock_64)KernelTime.dwHighDateTime << 32) | KernelTime.dwLowDateTime;
            my_clock_64 User = ((my_clock_64)UserTime.dwHighDateTime << 32) | UserTime.dwLowDateTime;

            current_time = (clock_t)((double)(Kernel + User) * CLOCKS_PER_SEC / 10000000);
        } else {
            current_time = -1;
        }
    } else {
        current_time = clock();
    }

    return current_time;
}
#else
clock_t my_clock(void) { return clock(); }
#endif

#if(defined _MSC_VER || defined __CYGWIN__) && !defined COMP_LANG_C_PEDANTIC
#else
#include <sys/times.h>
#include <unistd.h>
#endif

#ifndef _WIN32
namespace ClockVars {
    long lSC_CLK_TCK = -1;
}
#endif

double elapsedSeconds() {
#if(defined _MSC_VER || defined __CYGWIN__) && !defined COMP_LANG_C_PEDANTIC
    // windows version: just call my_clock and divide by CLOCKS_PER_SEC
    return my_clock() / (double)CLOCKS_PER_SEC;
#else
    // linux version: Use the "times" function. This allows us to measure time intervals
    // that are longer than about 70 minutes which are the maximum time span that can
    // be measured using "clock".
    tms timeInfo;
    times(&timeInfo);
    // avoid calling sysconf(...) more than once. Calling sysconf(...) might be
    // costly (honestly I don't know if this is the case... It seems like a simple
    // function that just return the right constant.)
    if(ClockVars::lSC_CLK_TCK == -1)
        ClockVars::lSC_CLK_TCK = sysconf(_SC_CLK_TCK);

    return timeInfo.tms_utime / (double)ClockVars::lSC_CLK_TCK;
#endif
}

// Set process affinity (fixes a process to a subset of cores)
// iVal should be seen as a binary vector so 1 means core A, 2 means core B and 3 means core A and B.
// Apparantly it is not well-defined how each bit is assigned to physical/logical cores.
bool setAffinity(int iVal) {
#if(defined _MSC_VER || defined __CYGWIN__) && !defined COMP_LANG_C_PEDANTIC
    // Gets the current process handle
    HANDLE hProc = GetCurrentProcess();
    DWORD newMask = iVal;
    BOOL res = SetProcessAffinityMask(hProc, (DWORD_PTR)newMask);
    return (bool)res;
#else
    return false;
#endif
}

double wallClock() {
#if(defined _MSC_VER || defined __CYGWIN__) && !defined COMP_LANG_C_PEDANTIC
    // windows version:
    __timeb64 timeStruct;
    _ftime64_s(&timeStruct);
    return ((double)timeStruct.time + timeStruct.millitm / 1000.0);

#else
    // Linux version
    timeval tim;
    gettimeofday(&tim, NULL);
    return ((double)tim.tv_sec + (tim.tv_usec / 1000000.0));
#endif
}

string int2String(int i) {
    // a positive number i requires (int) log10(i)+1 characters to be represented in decimal notation.
    // (eg (int) (log(99) + 1) = (int) (1.9956 + 1) = (int) 2.9956 = 2)
    // We need two extra characters as we must reserve space for the string termination character
    // and a minus sign.
    int iStringSize;
    if(i == 0)
        iStringSize = 2;
    else
        iStringSize = (int)log10((float)abs(i)) + 3;
    char* buffer = new char[iStringSize];
    sprintf(buffer, "%d", i);
    string strRes(buffer);
    delete[] buffer;
    return strRes;
}

string double2String(double d, int iPrecision) {
    char buffer[100];
    string stringFormat("%." + int2String(iPrecision) + "f");
    sprintf(buffer, stringFormat.c_str(), d);
    string strRes(buffer);
    return strRes;
}

string bool2String(bool b, const string& strIsTrue, const string& strIsFalse) {
    if(b)
        return strIsTrue;
    else
        return strIsFalse;
}

void waitEnter() {
    cout << "Press enter!" << endl;
    getchar();
}

/***********************************************************************
 * void firstPerm(vector<int> &vecPerm, int iN)
 *
 * Generate an initial permutation of the set {0,...,iN-1}
 *
 * --- input parameters ---
 * iN       : The number of elements in the permutation
 * --- output parameters ---
 * vecPerm  : The initial permuation
 * --- return value ---
 ***********************************************************************/

void firstPerm(vector<int>& vecPerm, int iN) {
    vecPerm.clear();
    vecPerm.resize(iN, -1);
    for(int i = 0; i < iN; i++)
        vecPerm[i] = i;
}

/* Test code:
    vector<int> vecPerm;
    int iN = 4;
    firstPerm(vecPerm,iN);
    bool bDone = false;
    while (!bDone)
    {
        for (int i=0; i < iN; i++)
            cout << vecPerm[i] << " ";
        cout << endl;
        bDone = !nextPerm(vecPerm,iN);
    }
*/

/***********************************************************************
 * bool nextPerm(vector<int> &vecPerm, int iN)
 *
 * Generate the next permutation of the set {0,...,iN-1} given
 * a previous permutation
 *
 * --- input parameters ---
 * iN       : The number of elements in the permutation
 * vecPerm  : The old permutation
 * --- output parameters ---
 * vecPerm  : The new permuation
 * --- return value ---
 * true => a new permutation was found.
 * false => no more permutations exists.
 ***********************************************************************/

bool nextPerm(vector<int>& vecPerm, int iN) {
    if(iN <= 0)
        return false;
    vector<bool> vecUsed(iN, false);
    int i, j;
    for(i = 0; i < iN - 1; i++)
        vecUsed[vecPerm[i]] = true;

    int iPos = iN - 2;
    bool bChanged = false;
    while(!bChanged && iPos >= 0) {
        vecUsed[vecPerm[iPos]] = false;
        for(i = vecPerm[iPos] + 1; i < iN && !bChanged; i++) {
            if(!vecUsed[i]) {
                vecPerm[iPos] = i;
                vecUsed[i] = true;
                bChanged = true;
            }
        }
        if(!bChanged)
            iPos--;
    }
    /*    cout << "iPos = " << iPos << ", vecUsed = ";
        for (i=0; i < (int) vecUsed.size(); i++)
            cout << vecUsed[i] << " ";
            cout << endl;*/
    if(bChanged) {
        // Find values for the positions after iPos
        for(i = iPos + 1; i < iN; i++) {
            bool bValAssigned = false;
            for(j = 0; j < iN && !bValAssigned; j++) {
                if(!vecUsed[j]) {
                    vecPerm[i] = j;
                    vecUsed[j] = true;
                    bValAssigned = true;
                }
            }
        }
    }
    return bChanged;
}

string bool2String(bool bBool) {
    if(bBool)
        return "true";
    else
        return "false";
}

/***********************************************************************
 * void parseNumbers(const string &str, vector<double> &vecNumbers)
 *
 * Converts a string into a series of numbers. If the string looks like this
 * 234 55   234  1 15
 * then the vector <vecNumbers> will contain (234, 55, 234, 1, 15). Numbers are
 * assumed to be separated by spaces or tabs. If the string contains letters,
 * then the behaviour depends on the function <atof>. That is, the behaviour is
 * undefined!
 *  The method can hopefully be useful when parsing data files.
 *
 * --- input parameters ---
 * str				: The string to parse.
 * vecNumbers		: The numbers found in the string.
 * --- output parameters ---
 * --- return value ---
 ***********************************************************************/
void parseNumbers(const string& str, vector<double>& vecNumbers) {
    vecNumbers.clear();
    string::size_type idxCurrent = 0;
    bool bDone = false;
    string::size_type npos = (string::size_type)-1;
    while(!bDone) {
        string::size_type idxStart = str.find_first_not_of(" \t", idxCurrent);
        if(idxStart == npos)
            bDone = true;
        else {
            string::size_type idxEnd = str.find_first_of(" \t", idxStart);
            string strNumber;
            if(idxEnd == npos) {
                strNumber = str.substr(idxStart);
                bDone = true;
            } else
                strNumber = str.substr(idxStart, idxEnd - idxStart);
            vecNumbers.push_back(atof(strNumber.c_str()));
            idxCurrent = idxEnd;
        }
    }
}

void parseNumbers(const string& str, vector<int>& vecNumbers) {
    vecNumbers.clear();
    string::size_type idxCurrent = 0;
    bool bDone = false;
    string::size_type npos = (string::size_type)-1;
    string::size_type idxStart, idxEnd;
    int iVal;
    while(!bDone) {
        idxStart = str.find_first_not_of(" \t", idxCurrent);
        if(idxStart == npos)
            bDone = true;
        else {
            idxEnd = str.find_first_of(" \t", idxStart);
            if(idxEnd == npos) {
                iVal = atoi(str.substr(idxStart).c_str());
                bDone = true;
            } else
                iVal = atoi(str.substr(idxStart, idxEnd - idxStart).c_str());
            vecNumbers.push_back(iVal);
            idxCurrent = idxEnd;
        }
    }
}

/***********************************************************************
 * void splitString(const string &str, const string &strDelim, vector<string> &vecTokens, bool bRepeatedDelimsCountAsOne)
 *
 * Splits the string <str> into a number of substrings. The delimiters in <strDelim>
 * decides where the string is split. Examples
 *
 * str = "Splits the string <str> into a number of substrings", strDelim = " ".
 * result: vecTokens = "Splits","the","string","<str>","into","a","number","of","substrings"
 *
 * The method has two modes of operation controlled by <bRepeatedDelimsCountAsOne>. if bRepeatedDelimsCountAsOne = true (default)
 * then repeated delimiters are counted as one so for the following example
 * str = "1,2,3, 4,,5,6", strDelim = " ,"
 * we would get
 * result: vecTokens = "1","2","3","4","5","6"
 *
 * If bRepeatedDelimsCountAsOne = false the same example would give (notice the empty strings)
 * vecTokens = "1","2","3","","4","","5","6"
 *
 * Notice that if bRepeatedDelimsCountAsOne = false and strDelim = "," then we get the following results
 * str = ","       =>  vecTokens = "",""
 * str = ",,"      =>  vecTokens = "","",""
 * str = "test,"   =>  vecTokens = "test",""
 * str = ",test,"  =>  vecTokens = "","test",""
 * str = ",,test," =>  vecTokens = "","","test",""
 *
 * --- input parameters ---
 * str				: The string to parse.
 * strDelim			: Delimiters
 * --- output parameters ---
 * vecTokens        : The tokens found
 * --- return value ---
 ***********************************************************************/
void splitString(const string& str, const string& strDelim, vector<string>& vecTokens, bool bRepeatedDelimsCountAsOne) {
    assert(!strDelim.empty());
    vecTokens.clear();
    string::size_type idxCurrent = 0;
    bool bDone = false;
    string::size_type npos = (string::size_type)-1;
    while(!bDone) {
        string::size_type idxStart = str.find_first_not_of(strDelim, idxCurrent);
        if(idxStart == npos) {
            if(!bRepeatedDelimsCountAsOne) {
                int iCountDelims = (int)(str.size() - idxCurrent);
                if(iCountDelims > 0) {
                    // Special case. If the string contains no tokens, but only delimiters, then we want an extra token corresponding
                    // to the entry to the left of the first delimiter. E.g. if "," is the delimiter then the string
                    // "," should give two empty strings as tokens and
                    // ",,," should give four empty strings as tokens and so on.
                    if(vecTokens.empty())
                        iCountDelims++;
                    vector<string> vecRepeatEmptyStrings(iCountDelims, "");
                    vecTokens.insert(vecTokens.end(), vecRepeatEmptyStrings.begin(), vecRepeatEmptyStrings.end());
                }
            }
            bDone = true;
        } else {
            if(!bRepeatedDelimsCountAsOne) {
                int iCountDelims = (int)(idxStart - idxCurrent);
                // Special case of a starting delimiter. In that case we want a empty string as the first token
                if(vecTokens.empty() && iCountDelims >= 1)
                    iCountDelims++;
                if(iCountDelims > 1) {
                    vector<string> vecRepeatEmptyStrings(iCountDelims - 1, "");
                    vecTokens.insert(vecTokens.end(), vecRepeatEmptyStrings.begin(), vecRepeatEmptyStrings.end());
                }
            }
            string::size_type idxEnd = str.find_first_of(strDelim, idxStart);
            string strToken;
            if(idxEnd == npos) {
                strToken = str.substr(idxStart);
                bDone = true;
            } else
                strToken = str.substr(idxStart, idxEnd - idxStart);
            vecTokens.push_back(strToken);
            idxCurrent = idxEnd;
        }
    }
}

/***********************************************************************
 * string repeatString(int iNRepeats, const string &str)
 *
 * Repeats the string <str> <iNRepeats> times.
 *
 * --- input parameters ---
 * iNRepeats		: The number of time to repeat the string.
 * str				: The string to repeat.
 * --- output parameters ---
 * --- return value ---
 * A string containing <str> repeated <iNRepeats> times.
 ***********************************************************************/

string repeatString(int iNRepeats, const string& str) {
    string str2;
    int i;
    for(i = 0; i < iNRepeats; i++)
        str2 += str;
    return str2;
}

/***********************************************************************
 * void skipPath(const string &strFullPath, string &strFilename)
 *
 * Strips the leading path from a "full path".
 * (the function is quite inefficient but who cares).
 *
 * --- input parameters ---
 * strFullPath		: a full path like "/home/sropke/c101.txt" or
 *					  "c:\data\r207.txt".
 * --- output parameters ---
 * strFilename		: The filename from the input path. In the two examples
 *					  above filename would be "c101.txt" and "r207.txt"
 *					  respectively.
 * --- return value ---
 ***********************************************************************/

void skipPath(const string& strFullPath, string& strFilename) {
    // Get the problem filename (skip the path)
    string::const_reverse_iterator rIter;
    strFilename = "";
    for(rIter = strFullPath.rbegin(); rIter != strFullPath.rend(); rIter++) {
        if(*rIter == '/' || *rIter == '\\')
            break;
        else
            strFilename = *rIter + strFilename;
    }
}

/***********************************************************************
 * string stripExtension(const string &strFileName)
 *
 * Remove the file extension from strFileName. If the function is called
 * with "testfile.txt" it will return "testfile".
 * If the function is called with "testfile" (that is, no extension)
 * it will return "testfile".
 * If the function is called with "testfile.tar.gz"  it will return
 * "testfile.tar".
 *
 * --- input parameters ---
 * strFileName		: The file name we should process.
 * --- output parameters ---
 * --- return value ---
 * the file name stripped of the extension
 ***********************************************************************/
string stripExtension(const string& strFileName) {
    // Remove the file extension
    string::const_reverse_iterator rIter;
    string strStrippedFilename = "";
    int iCount = 0;
    for(rIter = strFileName.rbegin(); rIter != strFileName.rend(); rIter++) {
        iCount++;
        if(*rIter == '.')
            break;
    }
    // We copy the filename until the last dot ("."). This could probably be done
    // in a more elegant way, but as long as it works I'm happy.
    if(rIter != strFileName.rend())
        strStrippedFilename.append(strFileName.begin(), strFileName.begin() + strFileName.size() - iCount);
    else
        // No extension
        strStrippedFilename = strFileName;
    return strStrippedFilename;
}

/***********************************************************************
 * string getDirectory(const string &strFullPath)
 *
 * Strips the directory from a full path. That is if given
 * "data/test/file1.txt" as paramter it returns "data/test".
 * If given "data\test\file1.txt" as paramter it returns "data\test".
 *
 * --- input parameters ---
 * strFullPath		: a full path like "/home/sropke/c101.txt" or
 *					  "c:\data\r207.txt".
 * --- output parameters ---
 * --- return value ---
 * the directory part of the full path.
 ***********************************************************************/
string getDirectory(const string& strFullPath) {
    string::const_iterator cIter;
    string strDirectory = "";
    for(cIter = strFullPath.begin(); cIter != strFullPath.end(); cIter++) {
        if(*cIter == '/' || *cIter == '\\')
            break;
        else
            strDirectory += *cIter;
    }
    return strDirectory;
}

void dropTrailingCharacters(string& str, const string& dropChars) {
    string::size_type endPos = str.find_last_not_of(dropChars);

    str = str.substr(0, endPos + 1);
}

string dropCharacters(const string& str, const string& strDropChars) {
    string::size_type i;
    string strRes;
    for(i = 0; i < str.size(); ++i) {
        if(strDropChars.find(str[i]) == string::npos)
            strRes.push_back(str[i]);
    }
    return strRes;
}

// Note: If you the strings your are searching for and replacing by are both single chars then use the STL algorithm "replace" instead.
void stringReplace(const string& strIn, string& strOut, const string& strSearchFor, const string& strReplaceBy) {
    string::size_type iPos = 0;
    strOut.clear();
    while(iPos <= strIn.size() - strSearchFor.size()) {
        // Does strIn[iPos..(iPos+L-1)] match strSearchFor (where L is the length of strSearchFor).
        bool bMatch = true;
        string::size_type iPos2 = 0;
        while(bMatch && iPos2 < strSearchFor.size()) {
            if(strIn[iPos + iPos2] != strSearchFor[iPos2])
                bMatch = false;
            iPos2++;
        }
        if(bMatch) {
            strOut.append(strReplaceBy);
            iPos += strSearchFor.size();
        } else {
            strOut.push_back(strIn[iPos]);
            iPos++;
        }
    }
}

/***********************************************************************
 * void randomPermutation(vector<int> &vecInts, int iLower, int iUpper)
 *
 * Makes a random permutation of the numbers {iLower,..., iUpper}
 *
 * --- input parameters ---
 * iLower		: Lowest number in the sequence.
 * iUpper		: Highest number in the sequence-
 * --- output parameters ---
 * vecInts		: The random permutation
 * --- return value ---
 ***********************************************************************/

void randomPermutation(vector<int>& vecInts, int iLower, int iUpper) {
    assert(iLower <= iUpper);
    int i, iSize = iUpper - iLower + 1;
    vector<int> vecTemp(iSize);
    for(i = 0; i < iSize; i++)
        vecTemp[i] = iLower + i;

    vecInts.clear();
    vecInts.resize(iSize);
    for(i = 0; i < iSize; i++) {
        int iRnd = Random::getRandom(0, (int)vecTemp.size() - 1);
        int iVal = vecTemp[iRnd];
        unorderedVectorRemoveAt(vecTemp, iRnd);
        vecInts[i] = iVal;
    }
}

/***********************************************************************
 * int getNeededPrecision(double dVal, int iNDecimals)
 *
 * We want to output the value <dVal> using the ios_base class (e.g. cout).
 * We want a fixed number of decimals to be shown (<iNDecimals>). We can
 * use the method <ios_base::precision>. But what parameter should we pass
 * for precision?? This method provides the answer.
 *
 * --- input parameters ---
 * dVal			: The number we wish to output.
 * iNDecimals	: The number of decimals we want to output.
 * --- output parameters ---
 * --- return value ---
 * The necessary parameter for <ios_base::precision>.
 ***********************************************************************/

int getNeededPrecision(double dVal, int iNDecimals) {
    double dLog = log10(dVal);
    // Number of digits to the left of the decimal separator
    // 95.1234	should give iNDigits = 2	[max(1, floor(1.978+1)) = min(1,2) = 2]
    // 100		should give iNDigits = 3	[max(1, floor(2+1)) = min(1,3) = 3]
    // 0.54		should give iNDigits = 1	[max(1, floor(-0.268+1)) = min(1,0.732) = 1]
    int iNDigits = (int)maxFunc(1.0, floor(dLog + 1));
    return iNDigits + iNDecimals;
}

// *********************************************************************
// *** Methods for reading from ifstream
// *********************************************************************

/***********************************************************************
 * void skip(ifstream &ifs, string skipString)
 *
 * keeps reading from <ifs> as long as we encounters characters from
 * <skipString> (useful for skipping whitespace).
 *
 * The method is weak on error handling!!!
 *
 * --- input parameters ---
 * ifs			: The ifstream to read from.
 * skipString	: the characters to skip.
 * --- output parameters ---
 * --- return value ---
 ***********************************************************************/

void skip(ifstream& ifs, const string& skipString) {
    char nextChar;
    bool bOk = true;
    int i;
    while(bOk) {
        nextChar = ifs.peek();
        bool bSkip = false;
        for(i = 0; i < (int)skipString.size() && !bSkip; i++) {
            if(nextChar == skipString[i])
                bSkip = true;
        }
        if(bSkip)
            ifs.get();
        else
            bOk = false;
    }
}

// *********************************************************************
// *********************************************************************
// *********************************************************************

/***********************************************************************
 * void generateFeasiblePos(const vector<int> &vecItemsLengths, int iMaxLength, vector<int> &vecFeasiblePos)
 *
 * Given the item lengths in <vecItemsLengths> find the feasible positions that can be reached by placing
 * the items after each other, with no gaps. The maximum position we consider is <MaxLength>.
 * One can also think of the "lengths" as "demands" and the maximum length as a "capacity". With this terminology
 * the method finds the total demands that can be realized when packing the items in a bin
 * with the specified capacity.
 *
 * --- input parameters ---
 * vecItemsLengths		: Item lengths (demands)
 * iMaxLength			: Maximum length to consider (capacity of bin)
 *
 * --- output parameters ---
 * vecFeasiblePos		: The positions that can be reached with the items (the total demands that can be realized).
 *						  The vector is sorted according to increasing position. The vector always start with 0 (realizable by not
 *						  using any items).
 * --- return value ---
 ***********************************************************************/

void generateFeasiblePos(const vector<int>& vecItemsLengths, int iMaxLength, vector<int>& vecFeasiblePos) {
    int iNItems = (int)vecItemsLengths.size();
    cout << endl << "Item lengths: ";
    outputVector(cout, vecItemsLengths, " ", "\n");
    vecFeasiblePos.clear();
    vecFeasiblePos.push_back(0);
    NSMatrix<bool> table(iNItems, iMaxLength + 1, false);
    int iPos, iItem, i;
    // This code could be implemented much more efficiently using bit-parallelity. Especially (*)
    for(iPos = 1; iPos <= iMaxLength; ++iPos) {
        bool bPosOk = false;
        for(iItem = 0; iItem < iNItems; ++iItem) {
            int iItemLength = vecItemsLengths[iItem];
            bool bOk = false;
            if(iItemLength == iPos)
                bOk = true;
            else {
                int iOtherPos = iPos - iItemLength;
                if(iOtherPos > 0) {
                    bOk = false;
                    // Notice we use "|" operator on purpose in order to avoid overhead of shortcircuiting operator.
                    // I wonder if it's worthwhile?

                    // (*) this part could be speed up by using bit-parallelity.
                    for(i = 0; i < iItem; ++i)
                        bOk = bOk | table.getElement(i, iOtherPos);
                }
            }
            if(bOk) {
                table.setElement(iItem, iPos, true);
                bPosOk = true;
            }
        }
        if(bPosOk)
            vecFeasiblePos.push_back(iPos);
    }
    cout << "Feasible positions" << endl;
    outputVector(cout, vecFeasiblePos, " ", "\n");
    // cout << "Dyn. Prog. table:" << endl;
    // cout << table << endl;
    // waitEnter();
}

// This is a quick hack... There must be a more educated way of doing this.
// Assumes that dVal > 0
double truncToNSignifDigits(double dVal, int iNDigits) {
    if(dVal <= 0)
        error("truncToNSignifDigits", "input is expected to be > 0");
    if(iNDigits <= 0)
        error("truncToNSignifDigits", "#digits is expected to be > 0");
    double dLog = log10(dVal);
    double dFloorLog = floor(dLog);
    return floor(dVal * pow(10.0, iNDigits - 1 - dFloorLog)) / pow(10.0, iNDigits - 1 - dFloorLog);
}

string secToHHMMSS(int sec) { return int2String(sec / 3600) + ":" + int2String((sec % 3600) / 60) + ":" + int2String(sec % 60); }
