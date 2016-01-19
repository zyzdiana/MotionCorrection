#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <complex>
#include <vector>
#include <math.h>
#include <string>

#define ZZ 32
#define YY 32
#define XX 32

using namespace std;
typedef vector<vector<vector<vector<double> > > > vector4D;
typedef vector<vector<vector<double> > > vector3D;
typedef vector<vector<double> > vector2D;
typedef vector<double> vector1D;

#define PI 3.1415926535897932

vector3D read_file(string filename) {
    //std::ifstream file_input("/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0_slice_0.dat",std::ios::binary);
    //ifstream file_input("/Users/zyzdiana/GitHub/AC297r-Volume-Registration/Thesis/testVol1.dat",std::ios::binary);
    ifstream file_input(filename,std::ios::binary);
    // Put some values in
    vector<double> buffer;
    int length = 32*32*32;
    buffer.reserve(length);
    copy(istream_iterator<double>(file_input),
        istream_iterator<double>(),
        back_inserter(buffer));
    file_input.close();

    vector3D Vol;
    
    Vol.resize(ZZ);
    for (int i = 0; i < ZZ; ++i) {
        Vol[i].resize(YY);
        
        for (int j = 0; j < YY; ++j)
            Vol[i][j].resize(XX);
    }
    for (int i = 0; i < ZZ; ++i) {
        for (int j = 0; j < YY; ++j){
            for (int k = 0; k < XX; ++k)
                Vol[i][j][k] = buffer[32*32*i+32*j+k];
        }
    }
    return Vol;
}

int clip(int input, int min, int max){
    if(input > max)
        return max;
    else if (input < min)
        return min;
    else
        return input;
}

double trilinear_interp(vector3D volume, double x, double y, double z){
    int shape = volume.size();

    int x0 = int(floor(x));
    int x1 = x0 + 1;
    int y0 = int(floor(y));
    int y1 = y0 + 1;
    int z0 = int(floor(z));
    int z1 = z0 + 1;

    x0 = clip(x0, 0, shape);
    x1 = clip(x1, 0, shape);
    y0 = clip(y0, 0, shape);
    y1 = clip(y1, 0, shape);
    z0 = clip(z0, 0, shape);
    z1 = clip(z1, 0, shape);

    //define some coefficients
    double xd = x - x0;
    double yd = y - y0;
    double zd = z - z0;

    //set up for the bilinear interpolation
    double C00 = volume[y0][x0][z0]*(1-xd) + volume[y0][x1][z0]*xd;
    double C10 = volume[y1][x0][z0]*(1-xd) + volume[y1][x1][z0]*xd;
    
    double C01 = volume[y0][x0][z1]*(1-xd) + volume[y0][x1][z1]*xd;
    double C11 = volume[y1][x0][z1]*(1-xd) + volume[y1][x1][z1]*xd;
    
    double C0 = C00*(1-yd) + C10*yd;
    double C1 = C01*(1-yd) + C11*yd;
    
    double C = C0*(1-zd) + C1*zd;
    return C; 
}

const int X_inv[64][64] ={
{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 9,-9,-9, 9, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-6, 6, 6,-6, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-6, 6, 6,-6, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 4,-4,-4, 4, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
{-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 9,-9, 0, 0,-9, 9, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-6, 6, 0, 0, 6,-6, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9, 0, 0,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0},
{ 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
{-27,27,27,-27,27,-27,-27,27,-18,-9,18, 9,18, 9,-18,-9,-18,18,-9, 9,18,-18, 9,-9,-18,18,18,-18,-9, 9, 9,-9,-12,-6,-6,-3,12, 6, 6, 3,-12,-6,12, 6,-6,-3, 6, 3,-12,12,-6, 6,-6, 6,-3, 3,-8,-4,-4,-2,-4,-2,-2,-1},
{18,-18,-18,18,-18,18,18,-18, 9, 9,-9,-9,-9,-9, 9, 9,12,-12, 6,-6,-12,12,-6, 6,12,-12,-12,12, 6,-6,-6, 6, 6, 6, 3, 3,-6,-6,-3,-3, 6, 6,-6,-6, 3, 3,-3,-3, 8,-8, 4,-4, 4,-4, 2,-2, 4, 4, 2, 2, 2, 2, 1, 1},
{-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0},
{18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1},
{-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1},
{ 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},
{-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0},
{18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1},
{-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1},
{ 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
{-12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1},
{ 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1}
};

vector1D X_inv_dot(vector1D Y){
    vector1D result;
    int shape = Y.size();
    result.resize(shape);
    for (int i=0; i<shape;++i){
        double tmp = 0;
        for (int j=0; j<shape; ++j){
           tmp += X_inv[i][j]*Y[j];
        }
        result[i] = tmp;
    }
    return result;
}

double dot(vector1D Y, vector1D A){
    int shape = Y.size();
    double result = 0;
    for (int i=0; i<shape;++i){
        result += Y[i]*A[i];
    }
    return result;
}

vector1D get_target_Y(double x, double y, double z){
    vector1D Y;
    Y.resize(64);
    Y[0] = 1.;
    Y[1] = x;
    Y[2] = x*x;
    Y[3] = x*x*x;
    Y[4] = y;
    Y[5] = x*y;
    Y[6] = x*x*y;
    Y[7] = x*x*x*y;
    Y[8] = y*y;
    Y[9] = x*y*y;
    Y[10] = x*x*y*y;
    Y[11] = x*x*x*y*y;
    Y[12] = y*y*y;
    Y[13] = x*y*y*y;
    Y[14] = x*x*y*y*y;
    Y[15] = x*x*x*y*y*y;
    

    Y[16] = z;
    Y[17] = x*z;
    Y[18] = x*x*z;
    Y[19] = x*x*x*z;
    Y[20] = y*z;
    Y[21] = x*y*z;
    Y[22] = x*x*y*z;
    Y[23] = x*x*x*y*z;
    Y[24] = y*y*z;
    Y[25] = x*y*y*z;
    Y[26] = x*x*y*y*z;
    Y[27] = x*x*x*y*y*z;
    Y[28] = y*y*y*z;
    Y[29] = x*y*y*y*z;
    Y[30] = x*x*y*y*y*z;
    Y[31] = x*x*x*y*y*y*z;

    Y[32] = z*z;
    Y[33] = x*z*z;
    Y[34] = x*x*z*z;
    Y[35] = x*x*x*z*z;
    Y[36] = y*z*z;
    Y[37] = x*y*z*z;
    Y[38] = x*x*y*z*z;
    Y[39] = x*x*x*y*z*z;
    Y[40] = y*y*z*z;
    Y[41] = x*y*y*z*z;
    Y[42] = x*x*y*y*z*z;
    Y[43] = x*x*x*y*y*z*z;
    Y[44] = y*y*y*z*z;
    Y[45] = x*y*y*y*z*z;
    Y[46] = x*x*y*y*y*z*z;
    Y[47] = x*x*x*y*y*y*z*z;

    Y[48] = z*z*z;
    Y[49] = x*z*z*z;
    Y[50] = x*x*z*z*z;
    Y[51] = x*x*x*z*z*z;
    Y[52] = y*z*z*z;
    Y[53] = x*y*z*z*z;
    Y[54] = x*x*y*z*z*z;
    Y[55] = x*x*x*y*z*z*z;
    Y[56] = y*y*z*z*z;
    Y[57] = x*y*y*z*z*z;
    Y[58] = x*x*y*y*z*z*z;
    Y[59] = x*x*x*y*y*z*z*z;
    Y[60] = y*y*y*z*z*z;
    Y[61] = x*y*y*y*z*z*z;
    Y[62] = x*x*y*y*y*z*z*z;
    Y[63] = x*x*x*y*y*y*z*z*z;
    return Y;
}

vector4D tricubic_derivatives(vector3D volume){
    int shape = volume.size();
    vector4D derivatives;
    derivatives.resize(shape+30);
    for (int i = 0; i < shape+30; ++i) {
        derivatives[i].resize(shape+30); 
        for (int j = 0; j < shape+30; ++j){
            derivatives[i][j].resize(shape+30);
            for (int k = 0; k < shape+30; ++k)
                derivatives[i][j][k].resize(64);
        }
    }

    for(int i=-15; i<shape+15; ++i){
        for(int j=-15; j<shape+15; ++j){
            for(int k=-15; k<shape+15; ++k){
                int x1 = i;
                int y1 = j;
                int z1 = k;
                
                int x0 = x1 - 1;
                int x2 = x1 + 1;
                int x3 = x2 + 1;
                int y0 = y1 - 1;
                int y2 = y1 + 1;
                int y3 = y2 + 1;
                int z0 = z1 - 1;
                int z2 = z1 + 1;
                int z3 = z2 + 1;

                //Wrap Around
                x0 = (x0 + shape) % shape;
                x1 = (x1 + shape) % shape;
                x2 = (x2 + shape) % shape;
                x3 = (x3 + shape) % shape;

                y0 = (y0 + shape) % shape;
                y1 = (y1 + shape) % shape;
                y2 = (y2 + shape) % shape;
                y3 = (y3 + shape) % shape;

                z0 = (z0 + shape) % shape;
                z1 = (z1 + shape) % shape;
                z2 = (z2 + shape) % shape;
                z3 = (z3 + shape) % shape;

                //Compute vector Y from known points
                vector1D Y;
                Y.resize(64);

                //values of f(x,y,z) at each corner.
                Y[0] = volume[y1][x1][z1];
                Y[1] = volume[y1][x1][z2];
                Y[2] = volume[y1][x2][z1];
                Y[3] = volume[y1][x2][z2];
                Y[4] = volume[y2][x1][z1];
                Y[5] = volume[y2][x1][z2];
                Y[6] = volume[y2][x2][z1];
                Y[7] = volume[y2][x2][z2];

                //values of df/dx at each corner.
                Y[8] = ((volume[y1][x1][z2]-volume[y1][x1][z0])/2.);
                Y[9] = ((volume[y1][x1][z3]-volume[y1][x1][z1])/2.);
                Y[10] = ((volume[y1][x2][z2]-volume[y1][x2][z0])/2.);
                Y[11] = ((volume[y1][x2][z3]-volume[y1][x2][z1])/2.);
                Y[12] = ((volume[y2][x1][z2]-volume[y2][x1][z0])/2.);
                Y[13] = ((volume[y2][x1][z3]-volume[y2][x1][z1])/2.);
                Y[14] = ((volume[y2][x2][z2]-volume[y2][x2][z0])/2.);
                Y[15] = ((volume[y2][x2][z3]-volume[y2][x2][z1])/2.);

                //values of df/dy at each corner.
                Y[16] = ((volume[y1][x2][z1]-volume[y1][x0][z1])/2.);
                Y[17] = ((volume[y1][x2][z2]-volume[y1][x0][z2])/2.);
                Y[18] = ((volume[y1][x3][z1]-volume[y1][x1][z1])/2.);
                Y[19] = ((volume[y1][x3][z2]-volume[y1][x1][z2])/2.);
                Y[20] = ((volume[y2][x2][z1]-volume[y2][x0][z1])/2.);
                Y[21] = ((volume[y2][x2][z2]-volume[y2][x0][z2])/2.);
                Y[22] = ((volume[y2][x3][z1]-volume[y2][x1][z1])/2.);
                Y[23] = ((volume[y2][x3][z2]-volume[y2][x1][z2])/2.);

                //values of df/dz at each corner.
                Y[24] = ((volume[y2][x1][z1]-volume[y0][x1][z1])/2.);
                Y[25] = ((volume[y2][x1][z2]-volume[y0][x1][z2])/2.);
                Y[26] = ((volume[y2][x2][z1]-volume[y0][x2][z1])/2.);
                Y[27] = ((volume[y2][x2][z2]-volume[y0][x2][z2])/2.);
                Y[28] = ((volume[y3][x1][z1]-volume[y1][x1][z1])/2.);
                Y[29] = ((volume[y3][x1][z2]-volume[y1][x1][z2])/2.);
                Y[30] = ((volume[y3][x2][z1]-volume[y1][x2][z1])/2.);
                Y[31] = ((volume[y3][x2][z2]-volume[y1][x2][z2])/2.);

                //values of d2f/dxdy at each corner.
                Y[32] = ((volume[y1][x2][z2]-volume[y1][x0][z2]-volume[y1][x2][z0]+volume[y1][x0][z0])/4.);
                Y[33] = ((volume[y1][x2][z3]-volume[y1][x0][z3]-volume[y1][x2][z1]+volume[y1][x0][z1])/4.);
                Y[34] = ((volume[y1][x3][z2]-volume[y1][x1][z2]-volume[y1][x3][z0]+volume[y1][x1][z0])/4.);
                Y[35] = ((volume[y1][x3][z3]-volume[y1][x1][z3]-volume[y1][x3][z1]+volume[y1][x1][z1])/4.);
                Y[36] = ((volume[y2][x2][z2]-volume[y2][x0][z2]-volume[y2][x2][z0]+volume[y2][x0][z0])/4.);
                Y[37] = ((volume[y2][x2][z3]-volume[y2][x0][z3]-volume[y2][x2][z1]+volume[y2][x0][z1])/4.);
                Y[38] = ((volume[y2][x3][z2]-volume[y2][x1][z2]-volume[y2][x3][z0]+volume[y2][x1][z0])/4.);
                Y[39] = ((volume[y2][x3][z3]-volume[y2][x1][z3]-volume[y2][x3][z1]+volume[y2][x1][z1])/4.);

                //values of d2f/dxdz at each corner.
                Y[40] = ((volume[y2][x1][z2]-volume[y0][x1][z2]-volume[y2][x1][z0]+volume[y0][x1][z0])/4.);
                Y[41] = ((volume[y2][x1][z3]-volume[y0][x1][z3]-volume[y2][x1][z1]+volume[y0][x1][z1])/4.);
                Y[42] = ((volume[y2][x2][z2]-volume[y0][x2][z2]-volume[y2][x2][z0]+volume[y0][x2][z0])/4.);
                Y[43] = ((volume[y2][x2][z3]-volume[y0][x2][z3]-volume[y2][x2][z1]+volume[y0][x2][z1])/4.);
                Y[44] = ((volume[y3][x1][z2]-volume[y1][x1][z2]-volume[y3][x1][z0]+volume[y1][x1][z0])/4.);
                Y[45] = ((volume[y3][x1][z3]-volume[y1][x1][z3]-volume[y3][x1][z1]+volume[y1][x1][z1])/4.);
                Y[46] = ((volume[y3][x2][z2]-volume[y1][x2][z2]-volume[y3][x2][z0]+volume[y1][x2][z0])/4.);
                Y[47] = ((volume[y3][x2][z3]-volume[y1][x2][z3]-volume[y3][x2][z1]+volume[y1][x2][z1])/4.);

                //values of d2f/dydz at each corner.
                Y[48] = ((volume[y2][x2][z1]-volume[y2][x0][z1]-volume[y0][x2][z1]+volume[y0][x0][z1])/4.);
                Y[49] = ((volume[y2][x2][z2]-volume[y2][x0][z2]-volume[y0][x2][z2]+volume[y0][x0][z2])/4.);
                Y[50] = ((volume[y2][x3][z1]-volume[y2][x1][z1]-volume[y0][x3][z1]+volume[y0][x1][z1])/4.);
                Y[51] = ((volume[y2][x3][z2]-volume[y2][x1][z2]-volume[y0][x3][z2]+volume[y0][x1][z2])/4.);
                Y[52] = ((volume[y3][x2][z1]-volume[y3][x0][z1]-volume[y1][x2][z1]+volume[y1][x0][z1])/4.);
                Y[53] = ((volume[y3][x2][z2]-volume[y3][x0][z2]-volume[y1][x2][z2]+volume[y1][x0][z2])/4.);
                Y[54] = ((volume[y3][x3][z1]-volume[y3][x1][z1]-volume[y1][x3][z1]+volume[y1][x1][z1])/4.);
                Y[55] = ((volume[y3][x3][z2]-volume[y3][x1][z2]-volume[y1][x3][z2]+volume[y1][x1][z2])/4.);

                //values of d3f/dxdydz at each corner.
                Y[56] = ((volume[y2][x2][z2]-volume[y2][x0][z2]-volume[y2][x2][z0]+volume[y2][x0][z0])
                          -(volume[y0][x2][z2]-volume[y0][x0][z2]-volume[y0][x2][z0]+volume[y0][x0][z0]))/8.;
                Y[57] = ((volume[y2][x2][z3]-volume[y2][x0][z3]-volume[y2][x2][z1]+volume[y2][x0][z1])
                          -(volume[y0][x2][z3]-volume[y0][x0][z3]-volume[y0][x2][z1]+volume[y0][x0][z1]))/8.;
                Y[58] = ((volume[y2][x3][z2]-volume[y2][x1][z2]-volume[y2][x3][z0]+volume[y2][x1][z0])
                          -(volume[y0][x3][z2]-volume[y0][x1][z2]-volume[y0][x3][z0]+volume[y0][x1][z0]))/8.;
                Y[59] = ((volume[y2][x3][z3]-volume[y2][x1][z3]-volume[y2][x3][z1]+volume[y2][x1][z1])
                          -(volume[y0][x3][z3]-volume[y0][x1][z3]-volume[y0][x3][z1]+volume[y0][x1][z1]))/8.;

                Y[60] = ((volume[y3][x2][z2]-volume[y3][x0][z2]-volume[y3][x2][z0]+volume[y3][x0][z0])
                          -(volume[y1][x2][z2]-volume[y1][x0][z2]-volume[y1][x2][z0]+volume[y1][x0][z0]))/8.;
                Y[61] = ((volume[y3][x2][z3]-volume[y3][x0][z3]-volume[y3][x2][z1]+volume[y3][x0][z1])
                          -(volume[y1][x2][z3]-volume[y1][x0][z3]-volume[y1][x2][z1]+volume[y1][x0][z1]))/8.;
                Y[62] = ((volume[y3][x3][z2]-volume[y3][x1][z2]-volume[y3][x3][z0]+volume[y3][x1][z0])
                          -(volume[y1][x3][z2]-volume[y1][x1][z2]-volume[y1][x3][z0]+volume[y1][x1][z0]))/8.;
                Y[63] = ((volume[y3][x3][z3]-volume[y3][x1][z3]-volume[y3][x3][z1]+volume[y3][x1][z1])
                          -(volume[y1][x3][z3]-volume[y1][x1][z3]-volume[y1][x3][z1]+volume[y1][x1][z1]))/8.;

                derivatives[j+15][i+15][k+15] = X_inv_dot(Y);
            }
        }
    }
    return derivatives;
}

double tricubic_interp(vector4D derivatives, double x, double y, double z){

    int x1 = int(floor(x));
    int y1 = int(floor(y));
    int z1 = int(floor(z));

    vector1D A = derivatives[y1+15][x1+15][z1+15];
    vector1D target_Y = get_target_Y(y-y1, x-x1, z-z1);
    return dot(target_Y,A);
}

double to_radian(double angle){
    return angle*PI/180.0;
}

vector1D rotate_coords_3d(double x, double y, double z, double theta, double wx, double wy, double wz, double ox, double oy, double oz){
    theta = to_radian(theta);
    // make sure w is a unit vetor:
    double tmp = wx*wx + wy*wy + wz*wz;
    if ((tmp != 1) && (tmp != 0)){
        double norm = sqrt(tmp);
        wx = wx/norm;
        wy = wy/norm;
        wz = wz/norm;
    }
    double s = sin(theta);
    double c = cos(theta);
    x = x - ox;
    y = y - oy;
    z = z - oz;

    vector1D result;
    result.resize(3);
    result[0] = c*x+s*(wy*z-wz*y)+(1-c)*(wx*x+wy*y+wz*z)*wx + ox;
    result[1] = c*y+s*(wz*x-wx*z)+(1-c)*(wx*x+wy*y+wz*z)*wy + oy;
    result[2] = c*z+s*(wx*y-wy*x)+(1-c)*(wx*x+wy*y+wz*z)*wz + oz;
    return result;
}
int main(){
    vector3D vol1 = read_file("/Users/zyzdiana/GitHub/AC297r-Volume-Registration/Thesis/testVol1.dat");
    // double res = trilinear_interp(vol1,-0.5,-1.5,2.5);
    //cout << res << endl;

    // vector4D derivs = tricubic_derivatives(vol1);

    // double interp = tricubic_interp(derivs,0.5,1.5,2.5);
    // cout << interp <<endl;


    // vector1D rot = rotate_coords_3d(10,14,18,15,0,0,1,15.5,15.5,15.5);
    // cout<< rot[0]<<endl;
    // cout<< rot[1]<<endl;
    // cout<< rot[2]<<endl;
    // ofstream output_file("vol1_derivs.dat",ios::binary);
    // ostream_iterator<double> output_iterator(output_file, "\n");
    // int shape = 32;

    // for (vector4D::iterator it_x = derivs.begin(); it_x != derivs.end(); ++it_x){
    //     for (vector3D::iterator it_y = it_x->begin(); it_y != it_x->end(); ++it_y){
    //         for (vector2D::iterator it_z = it_y->begin(); it_z != it_y->end(); ++it_z){
    //             copy(it_z->begin(), it_z->end(), output_iterator);
    //         }
    //     }
    // }

    // output_file.close();
    // for(int i=0; i<64;++i){
    //     cout << ' ' << Y[i];
    // }
    // cout << Y[0] << endl;
    // cout << Y[1] << endl;
    // cout << Y[2] << endl;
    // cout << Y[3] << endl;
    return 0;
}