#include "catch.hpp"

#include <complex>

#include "FFTOp.h"
#include "FFTWBuffer.h"
#include "VolumeAtAddressable.h"
#include "SymmetricHalfVolumeAtAddressable.h"

#include "CircularMaskOp.h"

#include "BinaryFile.h"

TEST_CASE("A circular mask can be instantiated") {
    typedef float dataT;
    typedef std::complex<float> complexT;
    typedef FFTOp<dataT> DataFFTOpT;
    typedef DataFFTOpT::spatialVolumeT DataVolumeT; 
    typedef DataFFTOpT::fourierVolumeT ComplexVolumeT; 
    typedef CircularMaskOp< DataVolumeT, dataT> DataCircularMaskOpT;
    typedef CircularMaskOp< ComplexVolumeT, dataT> ComplexDataCircularMaskOpT;

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    DataVolumeT initialData(cubeSize);

    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<DataVolumeT>::read(&initialData,
                "CircularMaskOp_tests/volInput.dat"));


    SECTION("and masking in the Fourier domain gives the precomputed answers") {
      ComplexDataCircularMaskOpT fourierMaskOp(cubeSize);
      ComplexVolumeT fourierData(cubeSize);
      DataFFTOpT fftOp(cubeSize);

      fftOp.forward(&initialData, &fourierData);
      
      fourierMaskOp.applyMask(&fourierData); 
       
      DataVolumeT resultData(cubeSize);
      
      fftOp.backward(&fourierData, &resultData);
    
      size_t pointsLength = resultData.totalPoints;
      std::vector<dataT> solutionData(pointsLength);
    
      REQUIRE(pointsLength * sizeof(dataT)
            == BinaryFile< std::vector<dataT> >::read(&solutionData,
                "CircularMaskOp_tests/CircularMaskOpFourierOutput.dat"));

      for(size_t i = 0; i < pointsLength; i++) {
        INFO("i: " << i);
        REQUIRE(resultData.at(i) == Approx(solutionData[i]).epsilon(0.005));
      }
    }
    
    SECTION("and masking in the image domain gives the precomputed answers") {
      DataCircularMaskOpT imageMaskOp(cubeSize);
      
      imageMaskOp.applyMask(&initialData); 
    
      size_t pointsLength = initialData.totalPoints;
      std::vector<dataT> solutionData(pointsLength);
    
      REQUIRE(pointsLength * sizeof(dataT)
            == BinaryFile< std::vector<dataT> >::read(&solutionData,
                "CircularMaskOp_tests/CircularMaskOpImageOutput.dat"));

      for(size_t i = 0; i < pointsLength; i++) {
        INFO("i: " << i);
        REQUIRE(initialData.at(i) == Approx(solutionData[i]));
      }
    }
}
