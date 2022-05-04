#include "geomutils.h"

#include <iostream>

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

int main(){
    get_convex_hull_custom();
    return 0;
}


