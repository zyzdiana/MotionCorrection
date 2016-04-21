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

template < typename VolumeT, typename value_type >
class CircularMaskOp_Base {
  protected:
    typedef Eigen::Map< Eigen::Matrix<value_type, Eigen::Dynamic, 1> >
      VolumeMapT;

  public:
    CircularMaskOp_Base(const size_t cubeSize) :
      cubeSize(cubeSize),
      volMask(cubeSize),
      volMaskMap(volMask.buffer, volMask.totalPoints)
      {}

    void applyMask(VolumeT *vol) const {
      VolumeMapT volMap(vol->buffer, vol->totalPoints);

      volMap = volMap.array() * volMaskMap.array();
    }
    
    void applyMask(const VolumeT *inVol, const VolumeT *outVol) const {
      VolumeMapT inVolMap(inVol->buffer, inVol->totalPoints);
      VolumeMapT outVolMap(outVol->buffer, outVol->totalPoints);

      outVolMap.noalias() = inVolMap.array() * volMaskMap.array();
    }

  protected:
    size_t cubeSize;
    size_t lastDim;
    VolumeT volMask;
    VolumeMapT volMaskMap;
};

template < typename AtAddressableT, typename T >
class CircularMaskOp {
};


template < typename AtAddressableT, typename T>
class CircularMaskOp< SymmetricHalfVolumeAtAddressable<AtAddressableT>, T > :
  public CircularMaskOp_Base<
    SymmetricHalfVolumeAtAddressable<AtAddressableT>,
    typename AtAddressableT::value_type
    >
    {
  public:
    typedef typename AtAddressableT::value_type value_type;
    typedef SymmetricHalfVolumeAtAddressable<AtAddressableT>
      SymmetricHalfVolT;

    typedef CircularMaskOp_Base<
      SymmetricHalfVolumeAtAddressable<AtAddressableT>,
      typename AtAddressableT::value_type
      > Parent;

    CircularMaskOp(const size_t cubeSize) :
      Parent(cubeSize),
      maskFunction(cubeSize)
    {
      const size_t lastDim = SymmetricHalfVolT::getLastFourierDimension(cubeSize);
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

            this->volMask.at(offset) =
              maskFunction.maskValue(zIndex, yIndex, xIndex);
          }
        }
      }
    }

  protected:
    CircularMaskOpFunction<T> maskFunction;
};

template < typename AtAddressableT, typename T>
class CircularMaskOp< VolumeAtAddressable<AtAddressableT>, T > :
  public CircularMaskOp_Base<
    VolumeAtAddressable<AtAddressableT>,
    typename AtAddressableT::value_type
    > {
  public:
    typedef typename AtAddressableT::value_type value_type;
    typedef VolumeAtAddressable<AtAddressableT> VolumeT;
    typedef CircularMaskOp_Base<
      VolumeAtAddressable<AtAddressableT>,
      typename AtAddressableT::value_type
      > Parent;

    CircularMaskOp(const size_t cubeSize) :
      Parent(cubeSize),
      maskFunction(cubeSize)
    {
      size_t offset = 0;

      T startIndex = - ((T) (cubeSize - 1)) / ((T) 2.0);
      T endIndex = startIndex + cubeSize; 

      for(T z = startIndex; z < endIndex; z += ((T) 1.0)) {

        for(T y = startIndex; y < endIndex; y += ((T) 1.0)) {
          
          for(T x = startIndex; x < endIndex; x += ((T) 1.0), offset++) {
            this->volMask.at(offset) =
              maskFunction.maskValue(z, y, x);
          }
        }
      }
    }

  protected:
    CircularMaskOpFunction<T> maskFunction;
};

#endif
