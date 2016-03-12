#ifndef Interpolator3D_h
#define Interpolator3D_h

template < 
    typename _VolumeT,
    typename coordT>
class Interpolator3D {
  public:
    typedef _VolumeT VolumeT;
    typedef typename _VolumeT::value_type T;

    const int cubeSize;
    Interpolator3D(const VolumeT *volume) :
        volume(volume), 
        cubeSize(volume->cubeSize){}

    virtual T interp(
        const coordT z,
        const coordT y,
        const coordT x) const = 0;

  protected:
    const VolumeT *volume;

};

#endif
