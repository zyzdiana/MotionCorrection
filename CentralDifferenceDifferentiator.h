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
      for(int z = 0; z < cubeSize; z++) {
        for(int y = 0; y < cubeSize; y++) {
          for(int x = 0; x < cubeSize; x++) {
            int xMinus1 = volume->wrapIndex(x - 1);
            int xPlus1 = volume->wrapIndex(x + 1);

            dx->at(z, y, x) = (
                volume->at(z, y, xPlus1) -
                volume->at(z, y, xMinus1)
              ) * (T) 0.5;
          }
        }
      }
    }

    virtual void yDerivative(VolumeT *dy) const {
      for(int z = 0; z < cubeSize; z++) {
        for(int y = 0; y < cubeSize; y++) {
          int yMinus1 = volume->wrapIndex(y - 1);
          int yPlus1 = volume->wrapIndex(y + 1);

          for(int x = 0; x < cubeSize; x++) {
            dy->at(z, y, x) = (
                volume->at(z, yPlus1, x) -
                volume->at(z, yMinus1, x)
              ) * (T) 0.5;
          }
        }
      }
    }

    virtual void zDerivative(VolumeT *dz) const {
      for(int z = 0; z < cubeSize; z++) {
        int zMinus1 = volume->wrapIndex(z - 1);
        int zPlus1 = volume->wrapIndex(z + 1);

        for(int y = 0; y < cubeSize; y++) {
          for(int x = 0; x < cubeSize; x++) {
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
    const int cubeSize; 
};

#endif
