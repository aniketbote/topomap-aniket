#include <vector>
#include "geomutils.h"
#include "mycode.h"

using namespace std;

vector< vector<double> > computeConvexHull(vector< vector<double> > i_matrix){
  Polygon custompts, customhull;
  for (int r = 0; r < i_matrix.size(); r++){
    custompts.push_back(Point(i_matrix[r][0], i_matrix[r][1]));
  }
  computeConvexHull(custompts, customhull);
  vector< vector<double> > res;
  for(int i = 0;i < customhull.size();i ++) {
        res[i][0] = customhull[i].x;
        res[i][1] = customhull[i].y;
  }
  return res;
}


// #include <vector>
// #include "code.h"

// using namespace std;

// vector<double> average (vector< vector<double> > i_matrix) {

//   // Compute average of each row..
//   vector <double> averages;
//   for (int r = 0; r < i_matrix.size(); r++){
//     double rsum = 0.0;
//     double ncols= i_matrix[r].size();
//     for (int c = 0; c< i_matrix[r].size(); c++){
//       rsum += i_matrix[r][c];
//     }
//     averages.push_back(rsum/ncols);
//   }
//   return averages;
// }