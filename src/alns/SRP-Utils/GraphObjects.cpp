#include "GraphObjects.h"
#include <iostream>

using namespace std;

ostream& operator<<(ostream& os, const SimpleEdge& edge) {
    os << "(" << edge.getNode1() << ", " << edge.getNode2() << ")";
    return os;
}

ostream& operator<<(ostream& os, const SimpleArc& arc) {
    os << "(" << arc.getNode1() << ", " << arc.getNode2() << ")";
    return os;
}
