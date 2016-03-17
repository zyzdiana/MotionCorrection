#ifndef VolumeAtAddressable_h
#define VolumeAtAddressable_h

#include <ctgmath>

template <
  typename AtAddressableT,
  typename CoordT
  >
class VolumeAtAddressable : public AtAddressableT {
  public:
    typedef typename AtAddressableT::value_type value_type;

    VolumeAtAddressable(
      const size_t cubeSize
    ) :
    cubeSize(cubeSize),
    totalPoints(cubeSize * cubeSize * cubeSize),
    cubeCenter(cubeCenterFromCubeSize(cubeSize)),
    maxCubeIndex(cubeSize - 1)
    {}
    
    VolumeAtAddressable(
      const size_t cubeSize,
      const size_t atAddressableSize
    ) :
    AtAddressableT(atAddressableSize),
    cubeSize(cubeSize),
    totalPoints(cubeSize * cubeSize * cubeSize),
    cubeCenter(cubeCenterFromCubeSize(cubeSize)),
    maxCubeIndex(cubeSize - 1)
    {}
    
    VolumeAtAddressable(
      const size_t cubeSize,
      const AtAddressableT atAddressable
    ) :
    AtAddressableT(atAddressable),
    cubeSize(cubeSize),
    totalPoints(cubeSize * cubeSize * cubeSize),
    cubeCenter(cubeCenterFromCubeSize(cubeSize)),
    maxCubeIndex(cubeSize - 1)
    {}

    const value_type& at(const size_t z, const size_t y, const size_t x) const {
        return asAtAddressable().at((z * cubeSize + y) * cubeSize + x);
    }
    
    value_type& at(const size_t z, const size_t y, const size_t x) {
        return asAtAddressable().at((z * cubeSize + y) * cubeSize + x);
    } 
    
    const value_type& at(const size_t offset) const {
        return asAtAddressable().at(offset);
    }
    
    value_type& at(const size_t offset) {
        return asAtAddressable().at(offset);
    } 
    
    int clipIndex(const int index) const {
        if(index > maxCubeIndex) {
            return maxCubeIndex;
        }
        
        if (index < 0) {
            return 0;
        }
        
        return index;
    }
   
    int wrapIndex(const int index) const {
        return (index + cubeSize) % cubeSize;
    }
    
    size_t wrapIndex(const size_t index) const {
        return (index + cubeSize) % cubeSize;
    }
    
    float wrapIndex(const float index) const {
        return fmod(index + (float) cubeSize, (float) cubeSize);
    }
    
    double wrapIndex(const double index) const {
        return fmod(index + (double) cubeSize, (double) cubeSize);
    }
    
    template <typename TripleT>
    TripleT cubeCenterAsTriple() const {
      return TripleT(cubeCenter, cubeCenter, cubeCenter); 
    }

  protected:
    AtAddressableT& asAtAddressable() {  
      return static_cast<AtAddressableT&>(*this);
    }
    
    const AtAddressableT& asAtAddressable() const {  
      return static_cast<const AtAddressableT&>(*this);
    }

    static CoordT cubeCenterFromCubeSize(const float cubeSize) {
        return ((CoordT) cubeSize)/(CoordT)2.0 - (CoordT)0.5;
    }

  public:
    const size_t cubeSize; 
    const size_t totalPoints;
    const size_t maxCubeIndex; 
    const CoordT cubeCenter;
};

#endif
