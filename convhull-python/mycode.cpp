#include <vector>
//#include "geomutils.h"
#include "mycode.h"
#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

using namespace std;

struct Point {
    double x, y;

    Point() {}
    Point(double x, double y) :x(x), y(y) {}
};

typedef std::vector<Point> Polygon;

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian);


void computeConvexHull(Polygon& pts, Polygon& chull) {
    chull.clear();
    if (pts.size() == 1) {
        chull.push_back(pts[0]);
        chull.push_back(pts[0]);
        return;
    }
    else if (pts.size() == 2) {
        chull.push_back(pts[0]);
        chull.push_back(pts[1]);
        chull.push_back(pts[0]);
        return;
    }

    typedef boost::tuple<double, double> point;
    typedef boost::geometry::model::multi_point<point> mpoints;
    typedef boost::geometry::model::polygon<point> polygon;

    mpoints mpts;

    for (int i = 0; i < pts.size(); i++) {
        boost::geometry::append(mpts, point(pts[i].x, pts[i].y));
    }
    polygon hull;

    // Polygon is closed
    boost::geometry::convex_hull(mpts, hull);
    for (auto pt : hull.outer()) {
        chull.push_back(Point(pt.get<0>(), pt.get<1>()));
    }
}

vector< vector<double> > customComputeConvexHull(vector< vector<double> > i_matrix) {
    // cout << "\nDone1.1";
    Polygon custompts, customhull;
    for (int r = 0; r < i_matrix.size(); r++) {
        custompts.push_back(Point(i_matrix[r][0], i_matrix[r][1]));
    }
    computeConvexHull(custompts, customhull);
    vector< vector<double> > res;
    for (int i = 0; i < customhull.size(); i++) 
    {
        vector<double> v1{ customhull[i].x,  customhull[i].y };
        res.push_back(v1);
    }
    return res;
}

// int main()
// {
//     // Create an empty vector
//     vector< vector<double> > mat, mat2;

//     vector<double> myRow1{0, 0};
//     mat.push_back(myRow1);

//     vector<double> myRow2{ 7.61, 9.48 };
//     mat.push_back(myRow2);

//     vector<double> myRow3{ 0, 9.48 };
//     mat.push_back(myRow3);
//     cout << "Done1";

//     mat2 = customComputeConvexHull(mat);
//     cout << "Done2";
//     for(int i = 0; i < mat2.size(); i++)
//     {
//         for(int j = 0; j < 2; j++)
//         {
//             cout << mat2[i][j] << " ";
//         }
//         cout << endl;
//     }

//     return 0;
// }