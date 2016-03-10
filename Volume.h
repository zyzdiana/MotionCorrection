#ifndef Volume_h
#define Volume_h

template <
    typename _T,
    typename _StorageT,
    typename _CoordT>
class Volume {
  public:
    typedef _T T;
    typedef _StorageT StorageT;
    typedef _CoordT CoordT;

    Volume(
        const StorageT data,
        const size_t cubeSize) :
        cubeSize(cubeSize),
        cubeCenter(cubeCenterFromCubeSize(cubeSize)),
        data(data),
        maxCubeIndex(cubeSize - 1) {}

    const T& at(const size_t z, const size_t y, const size_t x) const {
        return at((z * cubeSize + y) * cubeSize + x);
    }
    
    T& at(const size_t z, const size_t y, const size_t x) {
        return at((z * cubeSize + y) * cubeSize + x);
    }

    const T& at(const size_t index) const {
        return data[index]; 
    }
    
    T& at(const size_t index) {
        return data[index]; 
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

  protected:
    template <typename TripleT>
    TripleT cubeCenterAsTriple() const {
      return TripleT(cubeCenter, cubeCenter, cubeCenter); 
    }

    static CoordT cubeCenterFromCubeSize(const float cubeSize) {
        return ((CoordT) cubeSize)/(CoordT)2.0 - (CoordT)0.5;
    }


  public:
    const size_t cubeSize;
    const CoordT cubeCenter;

  protected:
    const int maxCubeIndex;
    StorageT data;
};

#endif
