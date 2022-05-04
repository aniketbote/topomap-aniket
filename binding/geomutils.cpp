#include "geomutils.h"

#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)


void computeConvexHull(Polygon &pts, Polygon &chull) {
    chull.clear();
    if(pts.size() == 1) {
        chull.push_back(pts[0]);
        chull.push_back(pts[0]);
        return;
    } else if(pts.size() == 2) {
        chull.push_back(pts[0]);
        chull.push_back(pts[1]);
        chull.push_back(pts[0]);
        return;
    }

    typedef boost::tuple<double, double> point;
    typedef boost::geometry::model::multi_point<point> mpoints;
    typedef boost::geometry::model::polygon<point> polygon;

    mpoints mpts;

    for(int i = 0;i < pts.size();i ++) {
        boost::geometry::append(mpts,point(pts[i].x,pts[i].y));
    }
    polygon hull;

    // Polygon is closed
    boost::geometry::convex_hull(mpts, hull);
    for(auto pt : hull.outer()) {
        chull.push_back(Point(pt.get<0>(), pt.get<1>()));
    }
}
