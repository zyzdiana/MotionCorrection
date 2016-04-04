#include "catch.hpp"

#include <complex>

#include "FFTOp.h"
#include "FFTWBuffer.h"
#include "VolumeAtAddressable.h"
#include "SymmetricHalfVolumeAtAddressable.h"

#include "BinaryFile.h"

TEST_CASE("Volumes can be allocated for the real-to-complex transform") {
    typedef float dataT;
    typedef std::complex<float> complexT;
    typedef FFTOp<dataT> DataFFTOpT;
    typedef DataFFTOpT::spatialVolumeT DataVolumeT; 
    typedef DataFFTOpT::fourierVolumeT ComplexVolumeT; 

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    DataVolumeT initialData(cubeSize);
    ComplexVolumeT fourierData(cubeSize);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData.at(i) = i; 
    }

    DataFFTOpT fftOp(cubeSize);

    fftOp.forward(&initialData, &fourierData);

    SECTION("and the FFT output gives the precomputed answers") {
      size_t pointsLength = fourierData.totalPoints;
      std::vector<complexT> solutionData(pointsLength);
    
      REQUIRE(pointsLength * sizeof(complexT)
            == BinaryFile< std::vector<complexT> >::read(&solutionData,
                "FFTOp_tests/FFTOpTestIndexArrayOutput.dat"));

      for(size_t i = 0; i < pointsLength; i++) {
        REQUIRE(fourierData.at(i).real() == Approx(solutionData[i].real()));
        REQUIRE(fourierData.at(i).imag() == Approx(solutionData[i].imag()));
      }
    }
    
    DataVolumeT backwardData(cubeSize);
    
    fftOp.backward(&fourierData, &backwardData);

    dataT invTotalRealPoints = 1.0/backwardData.totalPoints;

    SECTION("and the iFFT returns the original data scaled by ") {
      for(size_t i = 0; i < cubeVectorLength; i++) {
        REQUIRE(initialData.at(i) ==
          Approx(backwardData.at(i) * invTotalRealPoints));
      }
    }
}
