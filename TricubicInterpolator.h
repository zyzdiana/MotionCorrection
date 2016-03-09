#ifndef TricubicInterpolator_h
#define TricubicInterpolator_h

#include "Interpolator3D.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
using namespace Eigen;
using namespace std;

template <
    typename VolumeT,
    typename coordT>
class TricubicInterpolator :
    public Interpolator3D<VolumeT, coordT> {

  public:
    typedef typename VolumeT::T T;
    typedef Matrix< T, Dynamic, Dynamic >  MatrixXT;
    typedef Matrix< T, 64, 1 > Vector64T;
    typedef Matrix< T, 1, 64 > RowVector64T;
    typedef Matrix< T, 3, 1 > CoordT;
    typedef Matrix< T, 3, Dynamic >  Matrix3X;

    const int derives_shape;
    const int cubeSize;
    const CoordT cubeCenter;
    MatrixXT generate_X_inv(){
        MatrixXT X_inv(64,64);
        X_inv << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                -3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 9,-9,-9, 9, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                -6, 6, 6,-6, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                -6, 6, 6,-6, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 4,-4,-4, 4, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
                -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 9,-9, 0, 0,-9, 9, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                -6, 6, 0, 0, 6,-6, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9, 0, 0,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0,
                 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0,
                -27,27,27,-27,27,-27,-27,27,-18,-9,18, 9,18, 9,-18,-9,-18,18,-9, 9,18,-18, 9,-9,-18,18,18,-18,-9, 9, 9,-9,-12,-6,-6,-3,12, 6, 6, 3,-12,-6,12, 6,-6,-3, 6, 3,-12,12,-6, 6,-6, 6,-3, 3,-8,-4,-4,-2,-4,-2,-2,-1,
                18,-18,-18,18,-18,18,18,-18, 9, 9,-9,-9,-9,-9, 9, 9,12,-12, 6,-6,-12,12,-6, 6,12,-12,-12,12, 6,-6,-6, 6, 6, 6, 3, 3,-6,-6,-3,-3, 6, 6,-6,-6, 3, 3,-3,-3, 8,-8, 4,-4, 4,-4, 2,-2, 4, 4, 2, 2, 2, 2, 1, 1,
                -6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0,
                18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1,
                -12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1,
                 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                -6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
                -6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0,
                18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1,
                -12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1,
                 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0,
                -12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1,
                 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1;
        return X_inv;
    }
    TricubicInterpolator(const VolumeT *volume) :
        Interpolator3D<VolumeT, coordT>(volume),
        X_inv(generate_X_inv()),
        derives_shape(volume->cubeSize+30),
        cubeSize(volume->cubeSize),
        cubeCenter(volume->cubeCenter),
        derivatives(derives_shape * derives_shape * derives_shape),
        target_YVec(target_YArr)
        {
            // a vector used to temporary store the polynomials of x, y, and z 
            Vector64T Y(64);

            for(int z = 0; z < derives_shape; ++z){
                for(int y = 0; y < derives_shape; ++y){
                    for(int x = 0; x < derives_shape; ++x){
                        int z1 = x - 15;
                        int y1 = y - 15;
                        int x1 = z - 15;
                        
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
                        x0 = (x0 + cubeSize) % cubeSize;
                        x1 = (x1 + cubeSize) % cubeSize;
                        x2 = (x2 + cubeSize) % cubeSize;
                        x3 = (x3 + cubeSize) % cubeSize;

                        y0 = (y0 + cubeSize) % cubeSize;
                        y1 = (y1 + cubeSize) % cubeSize;
                        y2 = (y2 + cubeSize) % cubeSize;
                        y3 = (y3 + cubeSize) % cubeSize;

                        z0 = (z0 + cubeSize) % cubeSize;
                        z1 = (z1 + cubeSize) % cubeSize;
                        z2 = (z2 + cubeSize) % cubeSize;
                        z3 = (z3 + cubeSize) % cubeSize;                    

                        //values of f(x,y,z) at each corner.
                        Y(0)= this->volume->at(z1, y1, x1);
                        Y(1)= this->volume->at(z2, y1, x1);
                        Y(2)= this->volume->at(z1, y1, x2);
                        Y(3)= this->volume->at(z2, y1, x2);
                        Y(4)= this->volume->at(z1, y2, x1);
                        Y(5)= this->volume->at(z2, y2, x1);
                        Y(6)= this->volume->at(z1, y2, x2);
                        Y(7)= this->volume->at(z2, y2, x2);

                        //values of df/dx
                        Y(8)= ((this->volume->at(z2, y1, x1)-this->volume->at(z0, y1, x1))/2.);
                        Y(9)= ((this->volume->at(z3, y1, x1)-this->volume->at(z1, y1, x1))/2.);
                        Y(10) = ((this->volume->at(z2, y1, x2)-this->volume->at(z0, y1, x2))/2.);
                        Y(11) = ((this->volume->at(z3, y1, x2)-this->volume->at(z1, y1, x2))/2.);
                        Y(12) = ((this->volume->at(z2, y2, x1)-this->volume->at(z0, y2, x1))/2.);
                        Y(13) = ((this->volume->at(z3, y2, x1)-this->volume->at(z1, y2, x1))/2.);
                        Y(14) = ((this->volume->at(z2, y2, x2)-this->volume->at(z0, y2, x2))/2.);
                        Y(15) = ((this->volume->at(z3, y2, x2)-this->volume->at(z1, y2, x2))/2.);

                        //values of df/dy
                        Y(16) = ((this->volume->at(z1, y1, x2)-this->volume->at(z1, y1, x0))/2.);
                        Y(17) = ((this->volume->at(z2, y1, x2)-this->volume->at(z2, y1, x0))/2.);
                        Y(18) = ((this->volume->at(z1, y1, x3)-this->volume->at(z1, y1, x1))/2.);
                        Y(19) = ((this->volume->at(z2, y1, x3)-this->volume->at(z2, y1, x1))/2.);
                        Y(20) = ((this->volume->at(z1, y2, x2)-this->volume->at(z1, y2, x0))/2.);
                        Y(21) = ((this->volume->at(z2, y2, x2)-this->volume->at(z2, y2, x0))/2.);
                        Y(22) = ((this->volume->at(z1, y2, x3)-this->volume->at(z1, y2, x1))/2.);
                        Y(23) = ((this->volume->at(z2, y2, x3)-this->volume->at(z2, y2, x1))/2.);

                        //values of df/df
                        Y(24) = ((this->volume->at(z1, y2, x1)-this->volume->at(z1, y0, x1))/2.);
                        Y(25) = ((this->volume->at(z2, y2, x1)-this->volume->at(z2, y0, x1))/2.);
                        Y(26) = ((this->volume->at(z1, y2, x2)-this->volume->at(z1, y0, x2))/2.);
                        Y(27) = ((this->volume->at(z2, y2, x2)-this->volume->at(z2, y0, x2))/2.);
                        Y(28) = ((this->volume->at(z1, y3, x1)-this->volume->at(z1, y1, x1))/2.);
                        Y(29) = ((this->volume->at(z2, y3, x1)-this->volume->at(z2, y1, x1))/2.);
                        Y(30) = ((this->volume->at(z1, y3, x2)-this->volume->at(z1, y1, x2))/2.);
                        Y(31) = ((this->volume->at(z2, y3, x2)-this->volume->at(z2, y1, x2))/2.);

                        //values of d2f/dxdy
                        Y(32) = ((this->volume->at(z2, y1, x2)-this->volume->at(z2, y1, x0)-this->volume->at(z0, y1, x2)+this->volume->at(z0, y1, x0))/4.);
                        Y(33) = ((this->volume->at(z3, y1, x2)-this->volume->at(z3, y1, x0)-this->volume->at(z1, y1, x2)+this->volume->at(z1, y1, x0))/4.);
                        Y(34) = ((this->volume->at(z2, y1, x3)-this->volume->at(z2, y1, x1)-this->volume->at(z0, y1, x3)+this->volume->at(z0, y1, x1))/4.);
                        Y(35) = ((this->volume->at(z3, y1, x3)-this->volume->at(z3, y1, x1)-this->volume->at(z1, y1, x3)+this->volume->at(z1, y1, x1))/4.);
                        Y(36) = ((this->volume->at(z2, y2, x2)-this->volume->at(z2, y2, x0)-this->volume->at(z0, y2, x2)+this->volume->at(z0, y2, x0))/4.);
                        Y(37) = ((this->volume->at(z3, y2, x2)-this->volume->at(z3, y2, x0)-this->volume->at(z1, y2, x2)+this->volume->at(z1, y2, x0))/4.);
                        Y(38) = ((this->volume->at(z2, y2, x3)-this->volume->at(z2, y2, x1)-this->volume->at(z0, y2, x3)+this->volume->at(z0, y2, x1))/4.);
                        Y(39) = ((this->volume->at(z3, y2, x3)-this->volume->at(z3, y2, x1)-this->volume->at(z1, y2, x3)+this->volume->at(z1, y2, x1))/4.);

                        //values of d2f/dxdf
                        Y(40) = ((this->volume->at(z2, y2, x1)-this->volume->at(z2, y0, x1)-this->volume->at(z0, y2, x1)+this->volume->at(z0, y0, x1))/4.);
                        Y(41) = ((this->volume->at(z3, y2, x1)-this->volume->at(z3, y0, x1)-this->volume->at(z1, y2, x1)+this->volume->at(z1, y0, x1))/4.);
                        Y(42) = ((this->volume->at(z2, y2, x2)-this->volume->at(z2, y0, x2)-this->volume->at(z0, y2, x2)+this->volume->at(z0, y0, x2))/4.);
                        Y(43) = ((this->volume->at(z3, y2, x2)-this->volume->at(z3, y0, x2)-this->volume->at(z1, y2, x2)+this->volume->at(z1, y0, x2))/4.);
                        Y(44) = ((this->volume->at(z2, y3, x1)-this->volume->at(z2, y1, x1)-this->volume->at(z0, y3, x1)+this->volume->at(z0, y1, x1))/4.);
                        Y(45) = ((this->volume->at(z3, y3, x1)-this->volume->at(z3, y1, x1)-this->volume->at(z1, y3, x1)+this->volume->at(z1, y1, x1))/4.);
                        Y(46) = ((this->volume->at(z2, y3, x2)-this->volume->at(z2, y1, x2)-this->volume->at(z0, y3, x2)+this->volume->at(z0, y1, x2))/4.);
                        Y(47) = ((this->volume->at(z3, y3, x2)-this->volume->at(z3, y1, x2)-this->volume->at(z1, y3, x2)+this->volume->at(z1, y1, x2))/4.);

                        //values of d2f/dydf
                        Y(48) = ((this->volume->at(z1, y2, x2)-this->volume->at(z1, y2, x0)-this->volume->at(z1, y0, x2)+this->volume->at(z1, y0, x0))/4.);
                        Y(49) = ((this->volume->at(z2, y2, x2)-this->volume->at(z2, y2, x0)-this->volume->at(z2, y0, x2)+this->volume->at(z2, y0, x0))/4.);
                        Y(50) = ((this->volume->at(z1, y2, x3)-this->volume->at(z1, y2, x1)-this->volume->at(z1, y0, x3)+this->volume->at(z1, y0, x1))/4.);
                        Y(51) = ((this->volume->at(z2, y2, x3)-this->volume->at(z2, y2, x1)-this->volume->at(z2, y0, x3)+this->volume->at(z2, y0, x1))/4.);
                        Y(52) = ((this->volume->at(z1, y3, x2)-this->volume->at(z1, y3, x0)-this->volume->at(z1, y1, x2)+this->volume->at(z1, y1, x0))/4.);
                        Y(53) = ((this->volume->at(z2, y3, x2)-this->volume->at(z2, y3, x0)-this->volume->at(z2, y1, x2)+this->volume->at(z2, y1, x0))/4.);
                        Y(54) = ((this->volume->at(z1, y3, x3)-this->volume->at(z1, y3, x1)-this->volume->at(z1, y1, x3)+this->volume->at(z1, y1, x1))/4.);
                        Y(55) = ((this->volume->at(z2, y3, x3)-this->volume->at(z2, y3, x1)-this->volume->at(z2, y1, x3)+this->volume->at(z2, y1, x1))/4.);

                        //values of d3f/dxdydf
                        Y(56) = ((this->volume->at(z2, y2, x2)-this->volume->at(z2, y2, x0)-this->volume->at(z0, y2, x2)+this->volume->at(z0, y2, x0))
                                  -(this->volume->at(z2, y0, x2)-this->volume->at(z2, y0, x0)-this->volume->at(z0, y0, x2)+this->volume->at(z0, y0, x0)))/8.;
                        Y(57) = ((this->volume->at(z3, y2, x2)-this->volume->at(z3, y2, x0)-this->volume->at(z1, y2, x2)+this->volume->at(z1, y2, x0))
                                  -(this->volume->at(z3, y0, x2)-this->volume->at(z3, y0, x0)-this->volume->at(z1, y0, x2)+this->volume->at(z1, y0, x0)))/8.;
                        Y(58) = ((this->volume->at(z2, y2, x3)-this->volume->at(z2, y2, x1)-this->volume->at(z0, y2, x3)+this->volume->at(z0, y2, x1))
                                  -(this->volume->at(z2, y0, x3)-this->volume->at(z2, y0, x1)-this->volume->at(z0, y0, x3)+this->volume->at(z0, y0, x1)))/8.;
                        Y(59) = ((this->volume->at(z3, y2, x3)-this->volume->at(z3, y2, x1)-this->volume->at(z1, y2, x3)+this->volume->at(z1, y2, x1))
                                  -(this->volume->at(z3, y0, x3)-this->volume->at(z3, y0, x1)-this->volume->at(z1, y0, x3)+this->volume->at(z1, y0, x1)))/8.;

                        Y(60) = ((this->volume->at(z2, y3, x2)-this->volume->at(z2, y3, x0)-this->volume->at(z0, y3, x2)+this->volume->at(z0, y3, x0))
                                  -(this->volume->at(z2, y1, x2)-this->volume->at(z2, y1, x0)-this->volume->at(z0, y1, x2)+this->volume->at(z0, y1, x0)))/8.;
                        Y(61) = ((this->volume->at(z3, y3, x2)-this->volume->at(z3, y3, x0)-this->volume->at(z1, y3, x2)+this->volume->at(z1, y3, x0))
                                  -(this->volume->at(z3, y1, x2)-this->volume->at(z3, y1, x0)-this->volume->at(z1, y1, x2)+this->volume->at(z1, y1, x0)))/8.;
                        Y(62) = ((this->volume->at(z2, y3, x3)-this->volume->at(z2, y3, x1)-this->volume->at(z0, y3, x3)+this->volume->at(z0, y3, x1))
                                  -(this->volume->at(z2, y1, x3)-this->volume->at(z2, y1, x1)-this->volume->at(z0, y1, x3)+this->volume->at(z0, y1, x1)))/8.;
                        Y(63) = ((this->volume->at(z3, y3, x3)-this->volume->at(z3, y3, x1)-this->volume->at(z1, y3, x3)+this->volume->at(z1, y3, x1))
                                  -(this->volume->at(z3, y1, x3)-this->volume->at(z3, y1, x1)-this->volume->at(z1, y1, x3)+this->volume->at(z1, y1, x1)))/8.;

                        // compute the index to the derivatives vector from loop indicies i, j and k
                        int idx = derives_shape * (derives_shape * y + x) + z;
                        // cout << "idx: " << idx << endl;
                        // store the derivatives
                        derivatives[idx] = X_inv*Y;
                    }
                }
            }
        }

    void fill_target_Y(const T z, const T y, const T x) const {
        T xSq = x*x;

        target_YArr[0] = (T) 1.0;
        target_YArr[1] = x;
        target_YArr[2] = xSq;
        target_YArr[3] = xSq * x; // x*x*x
        target_YArr[4] = target_YArr[0] * y; // y
        target_YArr[5] = target_YArr[1] * y; // x*y
        target_YArr[6] = target_YArr[2] * y; // x*x*y
        target_YArr[7] = target_YArr[3] * y; // x*x*x*y
        target_YArr[8] = target_YArr[4] * y; // y*y
        target_YArr[9] = target_YArr[5] * y; // x*y*y
        target_YArr[10] = target_YArr[6] * y; // x*x*y*y
        target_YArr[11] = target_YArr[7] * y; // x*x*x*y*y
        target_YArr[12] = target_YArr[8] * y; // y*y*y
        target_YArr[13] = target_YArr[9] * y; // x*y*y*y
        target_YArr[14] = target_YArr[10] * y; // x*x*y*y*y
        target_YArr[15] = target_YArr[11] * y; // x*x*x*y*y*y
        

        target_YArr[16] = target_YArr[0] * z; // z
        target_YArr[17] = target_YArr[1] * z; // x*z
        target_YArr[18] = target_YArr[2] * z; // x*x*z
        target_YArr[19] = target_YArr[3] * z; // x*x*x*z
        target_YArr[20] = target_YArr[4] * z; // y*z
        target_YArr[21] = target_YArr[5] * z; // x*y*z
        target_YArr[22] = target_YArr[6] * z; // x*x*y*z
        target_YArr[23] = target_YArr[7] * z; // x*x*x*y*z
        target_YArr[24] = target_YArr[8] * z; // y*y*z
        target_YArr[25] = target_YArr[9] * z; // x*y*y*z
        target_YArr[26] = target_YArr[10] * z; // x*x*y*y*z
        target_YArr[27] = target_YArr[11] * z; // x*x*x*y*y*z
        target_YArr[28] = target_YArr[12] * z; // y*y*y*z
        target_YArr[29] = target_YArr[13] * z; // x*y*y*y*z
        target_YArr[30] = target_YArr[14] * z; // x*x*y*y*y*z
        target_YArr[31] = target_YArr[15] * z; // x*x*x*y*y*y*z

        target_YArr[32] = target_YArr[16] * z; // z*z
        target_YArr[33] = target_YArr[17] * z; // x*z*z
        target_YArr[34] = target_YArr[18] * z; // x*x*z*z
        target_YArr[35] = target_YArr[19] * z; // x*x*x*z*z
        target_YArr[36] = target_YArr[20] * z; // y*z*z
        target_YArr[37] = target_YArr[21] * z; // x*y*z*z
        target_YArr[38] = target_YArr[22] * z; // x*x*y*z*z
        target_YArr[39] = target_YArr[23] * z; // x*x*x*y*z*z
        target_YArr[40] = target_YArr[24] * z; // y*y*z*z
        target_YArr[41] = target_YArr[25] * z; // x*y*y*z*z
        target_YArr[42] = target_YArr[26] * z; // x*x*y*y*z*z
        target_YArr[43] = target_YArr[27] * z; // x*x*x*y*y*z*z
        target_YArr[44] = target_YArr[28] * z; // y*y*y*z*z
        target_YArr[45] = target_YArr[29] * z; // x*y*y*y*z*z
        target_YArr[46] = target_YArr[30] * z; // x*x*y*y*y*z*z
        target_YArr[47] = target_YArr[31] * z; // x*x*x*y*y*y*z*z

        target_YArr[48] = target_YArr[32] * z; // z*z*z;
        target_YArr[49] = target_YArr[33] * z; // x*z*z*z;
        target_YArr[50] = target_YArr[34] * z; // x*x*z*z*z;
        target_YArr[51] = target_YArr[35] * z; // x*x*x*z*z*z;
        target_YArr[52] = target_YArr[36] * z; // y*z*z*z;
        target_YArr[53] = target_YArr[37] * z; // x*y*z*z*z;
        target_YArr[54] = target_YArr[38] * z; // x*x*y*z*z*z;
        target_YArr[55] = target_YArr[39] * z; // x*x*x*y*z*z*z;
        target_YArr[56] = target_YArr[40] * z; // y*y*z*z*z;
        target_YArr[57] = target_YArr[41] * z; // x*y*y*z*z*z;
        target_YArr[58] = target_YArr[42] * z; // x*x*y*y*z*z*z;
        target_YArr[59] = target_YArr[43] * z; // x*x*x*y*y*z*z*z;
        target_YArr[60] = target_YArr[44] * z; // y*y*y*z*z*z;
        target_YArr[61] = target_YArr[45] * z; // x*y*y*y*z*z*z;
        target_YArr[62] = target_YArr[46] * z; // x*x*y*y*y*z*z*z;
        target_YArr[63] = target_YArr[47] * z; // x*x*x*y*y*y*z*z*z;
    }

    std::vector<Vector64T>* get_derivatives(){
        return &derivatives;
    }

    virtual T interp(
        const coordT z,
        const coordT y,
        const coordT x) const {

        int x1 = z;
        int y1 = y;
        int z1 = x;

        fill_target_Y(y-y1, x-z1, z-x1);
        return target_YVec.dot(derivatives[derives_shape * (derives_shape*(y1+15) + (x1+15)) + (z1+15)]);
    }
    virtual Matrix3X compute_axis_derivatives(){
            axis_derivatives.resize(3, cubeSize*cubeSize*cubeSize);
            int idx = 0;

            for(int z = 0; z < cubeSize; ++z){
                for(int y = 0; y < cubeSize; ++y){
                    for(int x = 0; x < cubeSize; ++x){
                        axis_derivatives.col(idx) << (derivatives[derives_shape * (derives_shape*(y+15) + (z+15)) + (x+15)])(4), 
                                            (derivatives[derives_shape * (derives_shape*(y+15) + (z+15)) + (x+15)])(16),
                                            (derivatives[derives_shape * (derives_shape*(y+15) + (z+15)) + (x+15)])(1);

                        idx++;
                    }
                }
            }
        return axis_derivatives;
    }

  protected:
    const MatrixXT X_inv;
    std::vector<Vector64T> derivatives;
    mutable T target_YArr[64];
    const Map<Vector64T> target_YVec;
    Matrix3X axis_derivatives;
};

#endif
