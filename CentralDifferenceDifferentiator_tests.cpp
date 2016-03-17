#include "catch.hpp"

#include "CentralDifferenceDifferentiator.h"

#include "VolumeAtAddressable.h"

#include <vector>
#include <complex>

TEST_CASE("the derivative of a constant volume is zero everywhere") {
    typedef std::complex<float> dataT;
    typedef VolumeAtAddressable< std::vector<dataT>, float> VolumeT; 

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = (dataT) 1.0; 
    }

    VolumeT volume(cubeSize, initialData);

    const dataT dataZero = 0;

    CentralDifferencesDifferentiator<VolumeT> differ(&volume);


    SECTION("in the z direction") {
      VolumeT dz(cubeSize, cubeVectorLength);
      
      differ.zDerivative(&dz);

      for(size_t i = 0; i < cubeVectorLength; i++) {
          REQUIRE(dataZero == dz.at(i)); 
      }
    }
    
    SECTION("in the y direction") {
      VolumeT dy(cubeSize, cubeVectorLength);
      
      differ.yDerivative(&dy);

      for(size_t i = 0; i < cubeVectorLength; i++) {
          REQUIRE(dataZero == dy.at(i)); 
      }
    }
    
    SECTION("in the x direction") {
      VolumeT dx(cubeSize, cubeVectorLength);
      
      differ.xDerivative(&dx);

      for(size_t i = 0; i < cubeVectorLength; i++) {
          REQUIRE(dataZero == dx.at(i)); 
      }
    }
}
