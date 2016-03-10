#ifndef Differentiator_h
#define Differentiator_h


template <typename VolumeT>
class Differentiator {
  public:
    virtual void xDerivative(VolumeT *dx) const = 0;
    virtual void yDerivative(VolumeT *dy) const = 0;
    virtual void zDerivative(VolumeT *dz) const = 0;
};

#endif
