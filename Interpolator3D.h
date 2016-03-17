#ifndef Interpolator3D_h
#define Interpolator3D_h

template < 
    typename _VolumeT,
    typename _CoordT>
class Interpolator3D {
  public:
    typedef _VolumeT VolumeT;
    typedef typename _VolumeT::value_type T;
    typedef _CoordT CoordT;

    const int cubeSize;
    Interpolator3D(const VolumeT *volume) :
        volume(volume), 
        cubeSize(volume->cubeSize){}

    virtual T interp(
        const CoordT z,
        const CoordT y,
        const CoordT x) const = 0;

  protected:
    const VolumeT *volume;

};

#endif
