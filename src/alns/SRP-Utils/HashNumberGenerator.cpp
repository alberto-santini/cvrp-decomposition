#include "HashNumberGenerator.h"
#include "Utils.h"

using namespace std;

vector<unsigned int> HashNumberGenerator::m_hashNumbers = vector<unsigned int>();
vector<unsigned int> FixedHashNumberGenerator::m_hashNumbers = vector<unsigned int>();

const vector<unsigned int>& HashNumberGenerator::getHashNumbers(int n) {
    if((int)m_hashNumbers.size() < n) {
        //		cout << "Extending hash number vector!" << endl;
        m_hashNumbers.reserve(n);
        while((int)m_hashNumbers.size() < n)
            m_hashNumbers.push_back(Random::getRandom(0, 2000000000));
    }
    return m_hashNumbers;
}

void HashNumberGenerator::clear() { m_hashNumbers.clear(); }
