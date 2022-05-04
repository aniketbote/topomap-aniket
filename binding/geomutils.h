#ifndef GEOMUTILS_H
#define GEOMUTILS_H

#include <vector>

struct Point {
    double x,y;

    Point(){}
    Point(double x, double y):x(x),y(y){}
};

typedef std::vector<Point> Polygon;

void computeConvexHull(Polygon &pts, Polygon &chull);

#endif // GEOMUTILS_H
