#include "catch.hpp"

#include "TricubicInterpolator.h"
#include "CentralDifferenceDifferentiator.h"

#include "VolumeAtAddressable.h"

#include "Interpolator3D_tests.h"

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
    dyDiffer.zDerivative(&dxz);

    CentralDifferencesDifferentiator<VolumeT> dxyDiffer(&dxy);
    
    VolumeT dxyz(cubeSize, cubeVectorLength);
    dxyDiffer.zDerivative(&dxz);

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
    dyDiffer.zDerivative(&dxz);

    CentralDifferencesDifferentiator<VolumeT> dxyDiffer(&dxy);
    
    VolumeT dxyz(cubeSize, cubeVectorLength);
    dxyDiffer.zDerivative(&dxz);

    InterpolatorT interpolator(&volume, &dx, &dy, &dz, &dxy, &dxz, &dyz, &dxyz);

    InterpolatorTests<InterpolatorT>::identity_tests(&interpolator, &volume); 
    InterpolatorTests<InterpolatorT>::constant_tests(
      &interpolator, &volume, constValue); 
}
