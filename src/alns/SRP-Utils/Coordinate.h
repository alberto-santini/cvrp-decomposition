#ifndef COORDINATE_H
#define COORDINATE_H

#include "Utils.h"

#include <ostream>
#include <vector>

using namespace std;

/***********************************************************************
 * class Coordinate
 *
 * A class representing a coordinate, not much to this for now.
 ***********************************************************************/
template<class TCoord>
class Coordinate
{
public:
	Coordinate()					{ m_x = 0; m_y = 0; }
	Coordinate(TCoord x, TCoord y)	{ m_x = x; m_y = y; }

	TCoord getX() const				{ return m_x; }
	TCoord getY() const				{ return m_y; }
	void setX(const TCoord &x)		{ m_x = x; }
	void setY(const TCoord &y)		{ m_y = y; }
	
	bool operator == (const Coordinate &other) const { return m_x == other.m_x && m_y == other.m_y; }
	bool operator != (const Coordinate &other) const { return m_x != other.m_x || m_y != other.m_y; }

	void print(ostream & os) const
	{
		os << "(";
		os.width (4);
		os << getX();
		os << ",";
		os.width (4);
		os << getY();
		os << ")";
	}

	template<class T> friend ostream& operator<<(ostream& os, const Coordinate<T>& coord);

private:
	TCoord m_x;
	TCoord m_y;
};

template<class TCoord>
void makeRandom(Coordinate<TCoord> &coord, 
				const TCoord &minX, const TCoord &minY, const TCoord &maxX, const TCoord &maxY)
{
	coord.setX((TCoord) Random::getRandomDouble(minX, maxX));
	coord.setY((TCoord) Random::getRandomDouble(minY, maxY));
}

/***********************************************************************
 * template<class TCoord>
 * bool isWithin(const Coordinate<TCoord> &coord, const Coordinate<TCoord> &ll, const Coordinate<TCoord> &ur) 
 *
 * Determines whether <coord> is within the rectangle defined by <ll>
 * and <ur>.
 *
 * --- input parameters ---
 * ll				: Lower-Left corner of the rectangle.
 * ur				: Upper-Right corner of the rectangle.
 *
 * --- return values ---
 * true => The coordinate is within the rectangle.
 * false => The coordinate is outside the rectangle.
 ***********************************************************************/

template<class TCoord>
bool isWithin(const Coordinate<TCoord> &coord, const Coordinate<TCoord> &ll, const Coordinate<TCoord> &ur) 
{
	return ((ll.getX() <= coord.getX()) && (ll.getY() <= coord.getY()) && 
			(coord.getX() <= ur.getX()) && (coord.getY() <= ur.getY()));
}

template<class U>
ostream& operator<<(ostream& os, const Coordinate<U>& coord) {
	coord.print(os);
	return os;
}

template<class T>
double calcDist(const Coordinate<T> &c1, const Coordinate<T> &c2)
{
	return calcDist(c1.getX(), c1.getY(), c2.getX(), c2.getY());
}

#endif
