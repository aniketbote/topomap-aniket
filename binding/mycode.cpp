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
  // vector< vector<double> > res;
  vector<vector<double>> res( customhull.size() , vector<double> (2));
  for(int i = 0;i < customhull.size();i ++) {
        res[i][0] = customhull[i].x;
        res[i][1] = customhull[i].y;
  }
  return res;
}

void print_polygon(Polygon &h, int name){
    std::cout << "\nHull in "<< name << ": \n"<<"[";
    for(int i = 0;i < h.size();i ++) {
        std::cout << "("<< h[i].x<< ", "<< h[i].y<<"), ";
    }
    std::cout <<"]\n";
}

void get_convex_hull_custom(){
    Polygon custompts;
    Polygon customhull;
    custompts.push_back(Point(0,0));
    custompts.push_back(Point(4.58,7.14));
    custompts.push_back(Point(0,7.14));
    computeConvexHull(custompts, customhull);
    print_polygon(customhull, -99999);
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
    cout <<"Done1\n";

    get_convex_hull_custom();

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