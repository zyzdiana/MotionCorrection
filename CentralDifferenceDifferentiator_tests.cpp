#include "catch.hpp"

#include "CentralDifferenceDifferentiator.h"

#include "Volume.h"

#include <vector>
#include <complex>

TEST_CASE("the derivative of a constant volume is zero everywhere") {
    typedef std::complex<float> dataT;
    typedef Volume<dataT, std::vector<dataT>, float > VolumeT; 

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = (dataT) 1.0; 
    }

    VolumeT volume(initialData, cubeSize);

    const dataT dataZero = 0;

    CentralDifferencesDifferentiator<VolumeT> differ(&volume);


    SECTION("in the z direction") {
      VolumeT dz(initialData, cubeSize);
      
      differ.zDerivative(&dz);

      for(size_t i = 0; i < cubeVectorLength; i++) {
          REQUIRE(dataZero == dz.at(i)); 
      }
    }
    
    SECTION("in the y direction") {
      VolumeT dy(initialData, cubeSize);
      
      differ.yDerivative(&dy);

      for(size_t i = 0; i < cubeVectorLength; i++) {
          REQUIRE(dataZero == dy.at(i)); 
      }
    }
    
    SECTION("in the x direction") {
      VolumeT dx(initialData, cubeSize);
      
      differ.xDerivative(&dx);

      for(size_t i = 0; i < cubeVectorLength; i++) {
          REQUIRE(dataZero == dx.at(i)); 
      }
    }
}
