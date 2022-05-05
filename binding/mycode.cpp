#include <vector>
#include "geomutils.h"
#include "mycode.h"
#include <iostream>

using namespace std;

vector< vector<double> > customComputeConvexHull(vector< vector<double> > i_matrix){
  cout <<"\nDone1.1";
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

int main()
{
    // Create an empty vector
    vector< vector<double> > mat,mat2;
    
    vector<double> myRow1(0,0);
    mat.push_back(myRow1);

    vector<double> myRow2(7.61,9.48);
    mat.push_back(myRow2);
    
    vector<double> myRow3(0,9.48);
    mat.push_back(myRow3);
    cout <<"Done1";

    mat2 = customComputeConvexHull(mat);
    cout <<"Done2";

    return 0;
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