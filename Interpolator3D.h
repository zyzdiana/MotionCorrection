#ifndef Interpolator3D_h
#define Interpolator3D_h

#include "Volume.h"

template < 
    typename VolumeT,
    typename coordT>
class Interpolator3D {
  public:
    typedef typename VolumeT::T T;

    Interpolator3D(const VolumeT *volume) :
        volume(volume) {}

    virtual T interp(
        const coordT z,
        const coordT y,
        const coordT x) const = 0;

  protected:
    const VolumeT *volume;

};

#endif
