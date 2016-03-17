#include "catch.hpp"

#include "TricubicInterpolator.h"
#include "CentralDifferenceDifferentiator.h"

#include "VolumeAtAddressable.h"

#include "Interpolator3D_tests.h"

#include "BinaryFile.h"

#include <vector>
#include <complex>

TEST_CASE("a tricubic interpolator can be created from a volume") {
    typedef float dataT;
    typedef VolumeAtAddressable< std::vector<dataT>, float> VolumeT; 
    typedef TricubicInterpolator<VolumeT, float> InterpolatorT;

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = i; 
    }

    VolumeT volume(cubeSize, initialData);

    CentralDifferencesDifferentiator<VolumeT> volDiffer(&volume);
    VolumeT dx(cubeSize, cubeVectorLength);
    volDiffer.xDerivative(&dx);

    VolumeT dy(cubeSize, cubeVectorLength);
    volDiffer.yDerivative(&dy);

    VolumeT dz(cubeSize, cubeVectorLength);
    volDiffer.zDerivative(&dz);

    CentralDifferencesDifferentiator<VolumeT> dxDiffer(&dx);
   
    VolumeT dxy(cubeSize, cubeVectorLength);
    dxDiffer.yDerivative(&dxy);
    
    VolumeT dxz(cubeSize, cubeVectorLength);
    dxDiffer.zDerivative(&dxz);

    CentralDifferencesDifferentiator<VolumeT> dyDiffer(&dy);

    VolumeT dyz(cubeSize, cubeVectorLength);
    dyDiffer.zDerivative(&dyz);

    CentralDifferencesDifferentiator<VolumeT> dxyDiffer(&dxy);
    
    VolumeT dxyz(cubeSize, cubeVectorLength);
    dxyDiffer.zDerivative(&dxyz);

    InterpolatorT interpolator(&volume, &dx, &dy, &dz, &dxy, &dxz, &dyz, &dxyz);

    InterpolatorTests<InterpolatorT>::identity_tests(&interpolator, &volume);

}


TEST_CASE(
  "a tricubic interpolator can be created from a constant-valued volume") {
    typedef float dataT;
    typedef VolumeAtAddressable< std::vector<dataT>, float> VolumeT; 
    typedef TricubicInterpolator<VolumeT, float> InterpolatorT;

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);
    
    const dataT constValue = 1.0;

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = constValue; 
    }

    VolumeT volume(cubeSize, initialData);

    CentralDifferencesDifferentiator<VolumeT> volDiffer(&volume);
    VolumeT dx(cubeSize, cubeVectorLength);
    volDiffer.xDerivative(&dx);

    VolumeT dy(cubeSize, cubeVectorLength);
    volDiffer.yDerivative(&dy);

    VolumeT dz(cubeSize, cubeVectorLength);
    volDiffer.zDerivative(&dz);

    CentralDifferencesDifferentiator<VolumeT> dxDiffer(&dx);
   
    VolumeT dxy(cubeSize, cubeVectorLength);
    dxDiffer.yDerivative(&dxy);
    
    VolumeT dxz(cubeSize, cubeVectorLength);
    dxDiffer.zDerivative(&dxz);

    CentralDifferencesDifferentiator<VolumeT> dyDiffer(&dy);

    VolumeT dyz(cubeSize, cubeVectorLength);
    dyDiffer.zDerivative(&dyz);

    CentralDifferencesDifferentiator<VolumeT> dxyDiffer(&dxy);
    
    VolumeT dxyz(cubeSize, cubeVectorLength);
    dxyDiffer.zDerivative(&dxyz);

    InterpolatorT interpolator(&volume, &dx, &dy, &dz, &dxy, &dxz, &dyz, &dxyz);

    InterpolatorTests<InterpolatorT>::identity_tests(&interpolator, &volume); 
    InterpolatorTests<InterpolatorT>::constant_tests(
      &interpolator, &volume, constValue); 
}

TEST_CASE(
  "a tricubic interpolator can be created from an image volume") {
    typedef float dataT;
    typedef VolumeAtAddressable< std::vector<dataT>, float> VolumeT; 
    typedef TricubicInterpolator<VolumeT, float> InterpolatorT;

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<dataT>::read(&initialData,
                "TricubicInterpolator_tests/testVolInput.dat"));

    VolumeT volume(cubeSize, initialData);

    CentralDifferencesDifferentiator<VolumeT> volDiffer(&volume);
    VolumeT dx(cubeSize, cubeVectorLength);
    volDiffer.xDerivative(&dx);

    VolumeT dy(cubeSize, cubeVectorLength);
    volDiffer.yDerivative(&dy);

    VolumeT dz(cubeSize, cubeVectorLength);
    volDiffer.zDerivative(&dz);

    CentralDifferencesDifferentiator<VolumeT> dxDiffer(&dx);
   
    VolumeT dxy(cubeSize, cubeVectorLength);
    dxDiffer.yDerivative(&dxy);
    
    VolumeT dxz(cubeSize, cubeVectorLength);
    dxDiffer.zDerivative(&dxz);

    CentralDifferencesDifferentiator<VolumeT> dyDiffer(&dy);

    VolumeT dyz(cubeSize, cubeVectorLength);
    dyDiffer.zDerivative(&dyz);

    CentralDifferencesDifferentiator<VolumeT> dxyDiffer(&dxy);
    
    VolumeT dxyz(cubeSize, cubeVectorLength);
    dxyDiffer.zDerivative(&dxyz);

    InterpolatorT interpolator(&volume, &dx, &dy, &dz, &dxy, &dxz, &dyz, &dxyz);

    InterpolatorTests<InterpolatorT>::identity_tests(&interpolator, &volume); 

    SECTION(
      "and interpolating at 0.25-voxel intervals gives the precomputed answers") {

      size_t pointsLength = 2146689;
      std::vector<dataT> solutionData(pointsLength);
    
      REQUIRE(pointsLength * sizeof(dataT)
            == BinaryFile<dataT>::read(&solutionData,
                "TricubicInterpolator_tests/testPointsOutput.dat"));
   
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
