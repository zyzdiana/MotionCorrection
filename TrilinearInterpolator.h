#ifndef TrilinearInterpolator_h
#define TrilinearInterpolator_h

#include "Interpolator3D.h"

template <
    typename VolumeT,
    typename coordT>
class TrilinearInterpolator :
    public Interpolator3D<VolumeT, coordT> {
  public:
    typedef typename VolumeT::T T;

    TrilinearInterpolator(const VolumeT *volume) :
        Interpolator3D<VolumeT, coordT>(volume) {}

    virtual T interp(
        const coordT z,
        const coordT y,
        const coordT x) const {
        
        // mdt Feb 17/16
        // casting a float to int always takes the floor of the float
        int x0 = x;
        int x1 = x0 + 1;
        int y0 = y;
        int y1 = y0 + 1;
        int z0 = z;
        int z1 = z0 + 1;

        // clip range
        x0 = this->volume->clipIndex(x0);
        x1 = this->volume->clipIndex(x1);
        y0 = this->volume->clipIndex(y0);
        y1 = this->volume->clipIndex(y1);
        z0 = this->volume->clipIndex(z0);
        z1 = this->volume->clipIndex(z1);
        
        //define some coefficients
        const coordT xd = x - (coordT) x0;
        const coordT yd = y - (coordT) y0;
        const coordT zd = z - (coordT) z0;

        const coordT oneMinusXd = ((coordT) 1.0) - xd;
        const coordT oneMinusYd = ((coordT) 1.0) - yd;
        const coordT oneMinusZd = ((coordT) 1.0) - zd;

        //set up for the bilinear interpolation
        const T C00 = this->volume->at(z0, y0, x0)*oneMinusXd + this->volume->at(z0, y0, x1)*xd;
        const T C10 = this->volume->at(z0, y1, x0)*oneMinusXd + this->volume->at(z0, y1, x1)*xd;

        const T C01 = this->volume->at(z1, y0, x0)*oneMinusXd + this->volume->at(z1, y0, x1)*xd;
        const T C11 = this->volume->at(z1, y1, x0)*oneMinusXd + this->volume->at(z1, y1, x1)*xd;
        
        const T C0 = C00*oneMinusYd + C10*yd;
        const T C1 = C01*oneMinusYd + C11*yd;
        
        const T C = C0*oneMinusZd + C1*zd;

        // return result
        return C; 
    }
};

#endif
