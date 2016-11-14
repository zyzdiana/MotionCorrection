#include "catch.hpp"

#include "CubicBSplineInterpolator.h"

#include "VolumeAtAddressable.h"

#include "Interpolator3D_tests.h"

#include "BinaryFile.h"

#include <vector>
#include <complex>

TEST_CASE("a cubic B-spline interpolator can be created from a volume") {
    typedef float dataT;
    typedef VolumeAtAddressable< std::vector<dataT> > VolumeT; 
    typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT;

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = i; 
    }

    VolumeT volume(cubeSize, initialData);

    InterpolatorT interpolator(&volume);

    InterpolatorTests<InterpolatorT>::approx_identity_tests(
      &interpolator, &volume, 0.0001);

}


TEST_CASE(
  "a cubic B-spline interpolator can be created from a constant-valued volume") {
    typedef float dataT;
    typedef VolumeAtAddressable< std::vector<dataT> > VolumeT; 
    typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT;

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    VolumeT volume(cubeSize);
    
    const dataT constValue = 1.0;

    for(size_t i = 0; i < cubeVectorLength; i++) {
        volume.at(i) = constValue; 
    }

    InterpolatorT interpolator(&volume);

    InterpolatorTests<InterpolatorT>::approx_identity_tests(
      &interpolator, &volume, 0.0001); 
    InterpolatorTests<InterpolatorT>::approx_constant_tests(
      &interpolator, &volume, constValue, 0.0001); 
}


TEST_CASE(
  "a cubic B-spline interpolator can be created from an image volume") {
    typedef float dataT;
    typedef VolumeAtAddressable< std::vector<dataT> > VolumeT; 
    typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT;

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    VolumeT volume(cubeSize);
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&volume,
                "CubicBSplineInterpolator_tests/testVolInput.dat"));

    InterpolatorT interpolator(&volume);

    InterpolatorTests<InterpolatorT>::approx_identity_tests(
      &interpolator, &volume, 0.0001); 

    SECTION(
      "and interpolating at 0.25-voxel intervals gives the precomputed answers") {

      size_t pointsLength = 2146689;
      std::vector<dataT> solutionData(pointsLength);
    
      REQUIRE(pointsLength * sizeof(dataT)
            == BinaryFile< std::vector<dataT> >::read(&solutionData,
                "CubicBSplineInterpolator_tests/testPointsOutput.dat"));
   
      size_t offset = 0;
      for(dataT z = 0; z <= cubeSize; z += 0.25) {
        dataT zWrapped = volume.wrapIndex(z);
        INFO("z: " << z);
        INFO("zWrapped: " << zWrapped);

        for(dataT y = 0; y <= cubeSize; y += 0.25) {
          dataT yWrapped = volume.wrapIndex(y);
          INFO("y: " << y);
          INFO("yWrapped: " << yWrapped);

          for(dataT x = 0; x <= cubeSize; x+= 0.25, offset += 1) {
            dataT xWrapped = volume.wrapIndex(x);
            INFO("x: " << x);
            INFO("xWrapped: " << xWrapped);
            INFO("offset: " << offset);
            // note that interp assumes points are already wrapped back
            // into the volume, so we need to fix this here
            REQUIRE(
              interpolator.interp(zWrapped, yWrapped, xWrapped) ==
              Approx(solutionData[offset]));
          }
        }
      }
    }
}

