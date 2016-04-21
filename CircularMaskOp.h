#ifndef CircularMaskOp_h
#define CircularMaskOp_h

#include <Eigen/Dense>

#include <cmath>

template <typename T>
class CircularMaskOpFunction {
  public:

    CircularMaskOpFunction(T cubeSize) :
      invRadius(((T) 2.0) / ((T) cubeSize)),
      eightOverThree(((T) 8.0) / ((T) 3.0))
      {}

    T maskValue(T z, T y, T x) {
      return maskValue( std::sqrt(z*z + y*y + x*x) );   
    }
    
    T maskValue(T radius) {
        T r = radius * invRadius;

        if(r < 0.75) {
          return ((T) 1.0); 
        }
        else if(r > 1) {
          return ((T) 0.0); 
        }
        
        return wcos(r * eightOverThree - ((T) 2.0));
    }

  protected:
    T wcos(T t) {
      if(t < -((T) 0.5) || t > ((T) 0.5)) {
        return 0;
      }
      else {
        return cos(M_PI * t); 
      }
    }

  protected:
    T invRadius;    
    const T eightOverThree;    
};

template < typename AtAddressableT, typename T >
class CircularMaskOp {
};


template < typename AtAddressableT, typename T>
class CircularMaskOp< SymmetricHalfVolumeAtAddressable<AtAddressableT>, T > {
  public:
    typedef typename AtAddressableT::value_type value_type;
    typedef SymmetricHalfVolumeAtAddressable<AtAddressableT>
      SymmetricHalfVolT;

    CircularMaskOp(const size_t cubeSize) :
      cubeSize(cubeSize),
      lastDim(SymmetricHalfVolT::getLastFourierDimension(cubeSize)),
      halfVolMask(cubeSize),
      halfVolMaskMap(halfVolMask.buffer, halfVolMask.totalPoints),
      maskFunction(cubeSize)
    {
      size_t offset = 0;
      for(int z = 0; z < cubeSize; z++) {
        int zIndex = z;
        if(zIndex >= lastDim) {
          zIndex -= cubeSize; 
        }

        for(int y = 0; y < cubeSize; y++) {
          int yIndex = y;
          if(yIndex >= lastDim) {
            yIndex -= cubeSize; 
          }
          
          for(int x = 0; x < lastDim; x++, offset++) {
            int xIndex = x;
            if(xIndex >= lastDim) {
              xIndex -= cubeSize; 
            }

            halfVolMask.at(offset) =
              maskFunction.maskValue(zIndex, yIndex, xIndex);
          }
        }
      }
    }

    void applyMask(SymmetricHalfVolT *halfVol) const {
      HalfVolMapT halfVolMap(halfVol->buffer, halfVol->totalPoints);

      halfVolMap = halfVolMap.array() * halfVolMaskMap.array();
    }


  protected:
    typedef Eigen::Map< Eigen::Matrix<value_type, Eigen::Dynamic, 1> >
      HalfVolMapT;

    size_t cubeSize;
    size_t lastDim;
    SymmetricHalfVolT halfVolMask;
    HalfVolMapT halfVolMaskMap;
    CircularMaskOpFunction<T> maskFunction;
};

template < typename AtAddressableT, typename T>
class CircularMaskOp< VolumeAtAddressable<AtAddressableT>, T > {
  public:
    typedef typename AtAddressableT::value_type value_type;
    typedef VolumeAtAddressable<AtAddressableT> VolumeT;

    CircularMaskOp(const size_t cubeSize) :
      cubeSize(cubeSize),
      volMask(cubeSize),
      volMaskMap(volMask.buffer, volMask.totalPoints),
      maskFunction(cubeSize)
    {
      size_t offset = 0;

      T startIndex = - ((T) (cubeSize - 1)) / ((T) 2.0);
      T endIndex = startIndex + cubeSize; 

      for(T z = startIndex; z < endIndex; z += ((T) 1.0)) {

        for(T y = startIndex; y < endIndex; y += ((T) 1.0)) {
          
          for(T x = startIndex; x < endIndex; x += ((T) 1.0), offset++) {
            volMask.at(offset) =
              maskFunction.maskValue(z, y, x);
          }
        }
      }
    }

    void applyMask(VolumeT *vol) const {
      VolumeMapT volMap(vol->buffer, vol->totalPoints);

      volMap = volMap.array() * volMaskMap.array();
    }

  protected:
    typedef Eigen::Map< Eigen::Matrix<value_type, Eigen::Dynamic, 1> >
      VolumeMapT;

    size_t cubeSize;
    size_t lastDim;
    VolumeT volMask;
    VolumeMapT volMaskMap;
    CircularMaskOpFunction<T> maskFunction;
};

#endif
