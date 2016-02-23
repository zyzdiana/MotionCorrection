#ifndef Volume_h
#define Volume_h

#include <Eigen/Dense>

template <
    typename _T,
    typename _StorageT >
class Volume {
  public:
    typedef _T T;
    typedef _StorageT StorageT;

    Volume(
        const StorageT data,
        const size_t cubeSize) :
        cubeSize(cubeSize),
        cubeCenter((cubeSize/2 - 0.5),(cubeSize/2 - 0.5),(cubeSize/2 - 0.5)),
        data(data),
        maxCubeIndex(cubeSize - 1) {}

    T at(const size_t z, const size_t y, const size_t x) const {
        return at((z * cubeSize + y) * cubeSize + x);
    }

    T at(const size_t index) const {
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



  public:
    const size_t cubeSize;
    const Eigen::Vector3f cubeCenter;

  protected:
    const int maxCubeIndex;
    const StorageT data;
};

#endif
