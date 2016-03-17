#ifndef CentralDifferencesDifferentiator_h
#define CentralDifferencesDifferentiator_h

template <typename VolumeT>
class CentralDifferencesDifferentiator {
  public:
    typedef typename VolumeT::value_type T;

    CentralDifferencesDifferentiator(const VolumeT *volume) :
      volume(volume),
      cubeSize(volume->cubeSize) {}

    virtual void xDerivative(VolumeT *dx) const {
      for(size_t z = 0; z < cubeSize; z++) {
        for(size_t y = 0; y < cubeSize; y++) {
          for(size_t x = 0; x < cubeSize; x++) {
            size_t xMinus1 = volume->wrapIndex(x - 1);
            size_t xPlus1 = volume->wrapIndex(x + 1);

            dx->at(z, y, x) = (
                volume->at(z, y, xPlus1) -
                volume->at(z, y, xMinus1)
              ) * (T) 0.5;
          }
        }
      }
    }

    virtual void yDerivative(VolumeT *dy) const {
      for(size_t z = 0; z < cubeSize; z++) {
        for(size_t y = 0; y < cubeSize; y++) {
          size_t yMinus1 = volume->wrapIndex(y - 1);
          size_t yPlus1 = volume->wrapIndex(y + 1);

          for(size_t x = 0; x < cubeSize; x++) {
            dy->at(z, y, x) = (
                volume->at(z, yPlus1, x) -
                volume->at(z, yMinus1, x)
              ) * (T) 0.5;
          }
        }
      }
    }

    virtual void zDerivative(VolumeT *dz) const {
      for(size_t z = 0; z < cubeSize; z++) {
        size_t zMinus1 = volume->wrapIndex(z - 1);
        size_t zPlus1 = volume->wrapIndex(z + 1);

        for(size_t y = 0; y < cubeSize; y++) {
          for(size_t x = 0; x < cubeSize; x++) {
            dz->at(z, y, x) = (
                volume->at(zPlus1, y, x) -
                volume->at(zMinus1, y, x)
              ) * (T) 0.5;
          }
        }
      }
    }

  protected:
    const VolumeT *volume;
    const size_t cubeSize; 
};

#endif
