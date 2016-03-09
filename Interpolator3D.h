#ifndef Interpolator3D_h
#define Interpolator3D_h

#include "Volume.h"
using namespace Eigen;
template < 
    typename VolumeT,
    typename coordT>
class Interpolator3D {
  public:
    typedef typename VolumeT::T T;
    typedef Matrix< T, 3, 1> CoordT;
    typedef Matrix< T, 3, Dynamic >  Matrix3X;

    const int cubeSize;
    Interpolator3D(const VolumeT *volume) :
        volume(volume), 
        cubeSize(volume->cubeSize){}

    virtual T interp(
        const coordT z,
        const coordT y,
        const coordT x) const = 0;

    virtual Matrix3X compute_axis_derivatives(){
            axis_derivatives.resize(3, cubeSize*cubeSize*cubeSize);
            int idx = 0;

            for(int z = 0; z < cubeSize; ++z){
                for(int y = 0; y < cubeSize; ++y){
                    for(int x = 0; x < cubeSize; ++x){
                        int z1 = x;
                        int y1 = y;
                        int x1 = z;
                        
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

                        axis_derivatives.col(idx) << this->volume->at(z1, y2, x1),
                                                     (this->volume->at(z1, y1, x2)-this->volume->at(z1, y1, x0))/T(2),
                                                     this->volume->at(z2, y1, x1);

                        idx++;
                    }
                }
            }
        return axis_derivatives;
    }
  protected:
    const VolumeT *volume;
    Matrix3X axis_derivatives;

};

#endif
