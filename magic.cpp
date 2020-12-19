#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include "vec3.h"

using namespace std;

int main ()
{
    int nx = 1500;
    int ny = 800;
    string filePath = "C:\\abhishek\\RayTracing\\images\\ppmFiles\\ppmtest.ppm";
    string fileName = "ppmtest";
    ofstream myfile;
    myfile.open (filePath);
    myfile<<"P3\n"<<nx<<" "<<ny<<"\n255\n";

    for(int j = ny-1; j >= 0; j--)
    {
        for(int i = 0; i < nx; i++)
        {
            vec3 col (float(i) / float(nx), float(j) / float(ny), 0);
            int ir = int(255.99*col[0]);
            int ig = int(255.99*col[1]);
            int ib = int(255.99*col[2]);
            myfile<<ir<<" "<<ig<<" "<<ib<<"\n";
        }
    }
    return 1;
}