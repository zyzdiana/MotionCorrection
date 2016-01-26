#ifndef Utils_h
#define Utils_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <complex>
#include <vector>
#include <math.h>
#include <string>


using namespace std;
template <typename T>
class Utils{
  public:
    const double PI = 3.141592653589793238463;
    typedef vector<vector<vector<vector<T> > > > vector4D;
    typedef vector<vector<vector<T> > > vector3D;
    typedef vector<vector<T> > vector2D;
    typedef vector<T> vector1D;

    int clip(int input, int min, int max){
        if(input > max)
            return max;
        else if (input < min)
            return min;
        else
            return input;
    }
};