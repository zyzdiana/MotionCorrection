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

  protected:
    const VolumeT *volume;

};

#endif
