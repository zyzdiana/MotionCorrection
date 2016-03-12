#ifndef TrilinearInterpolator_h
#define TrilinearInterpolator_h

#include "Interpolator3D.h"
template <
    typename VolumeT,
    typename CoordT>
class TrilinearInterpolator :
    public Interpolator3D<VolumeT, CoordT> {
  public:
    typedef typename VolumeT::value_type T;

    TrilinearInterpolator(const VolumeT *volume) :
        Interpolator3D<VolumeT, CoordT>(volume) {}

    virtual T interp(
        const CoordT z,
        const CoordT y,
        const CoordT x) const {
        
        // mdt Feb 17/16
        // casting a float to int always takes the floor of the float
        int x0 = x;
        int x1 = x0 + 1;
        int y0 = y;
        int y1 = y0 + 1;
        int z0 = z;
        int z1 = z0 + 1;

        // mdt Mar 10/16
        // wrap range
        x0 = this->volume->wrapIndex(x0);
        x1 = this->volume->wrapIndex(x1);
        y0 = this->volume->wrapIndex(y0);
        y1 = this->volume->wrapIndex(y1);
        z0 = this->volume->wrapIndex(z0);
        z1 = this->volume->wrapIndex(z1);
        
        //define some coefficients
        const CoordT xd = x - (CoordT) x0;
        const CoordT yd = y - (CoordT) y0;
        const CoordT zd = z - (CoordT) z0;

        const CoordT oneMinusXd = ((CoordT) 1.0) - xd;
        const CoordT oneMinusYd = ((CoordT) 1.0) - yd;
        const CoordT oneMinusZd = ((CoordT) 1.0) - zd;

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
