#include "catch.hpp"

#include "Gauss_Newton.h"

#include "CubicBSplineInterpolator.h"
#include "VolumeAtAddressable.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTWBuffer.h"

TEST_CASE("a Gauss-Newton minimizer can be instantiated") {
    typedef float dataT;
    typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
    typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT; 
    typedef Gauss_Newton<InterpolatorT> MinimizerT; 
    typedef MinimizerT::ParamT ParamT;

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    VolumeT volume(cubeSize);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        volume.at(i) = ((dataT) i) / ((dataT) cubeVectorLength); 
    }

    InterpolatorT interpolator(&volume);

    CentralDifferencesDifferentiator<VolumeT> volDiffer(&volume);
    VolumeT dx(cubeSize, cubeVectorLength);
    volDiffer.xDerivative(&dx);

    VolumeT dy(cubeSize, cubeVectorLength);
    volDiffer.yDerivative(&dy);

    VolumeT dz(cubeSize, cubeVectorLength);
    volDiffer.zDerivative(&dz);
   
    double gradientAndHessianComputeTime;

    MinimizerT minimizer(&interpolator, &dz, &dy, &dx,
      &gradientAndHessianComputeTime);
      
    WARN("elapsed time computing gradient and Hessian: "
      << gradientAndHessianComputeTime << " ms");

    SECTION("and registering an image with itself produces 0 transformation") {
      ParamT initialParam;
      initialParam << 0, 0, 0, 0, 0, 0;
  
      ParamT finalParam;
 
      size_t maxSteps = 20;

      dataT paramUpdateNormLimit = 1e-10;
  
      double elapsedTime;
      size_t elapsedSteps;
  
      minimizer.minimize(&volume, &initialParam, &finalParam,
        maxSteps, paramUpdateNormLimit, &elapsedSteps, &elapsedTime);
 
        for(int i = 0; i < 6; i++) {
          REQUIRE(0 == Approx(finalParam(i)));
        }

      WARN("elapsed time: " << elapsedTime << " ms");
      WARN("elapsed steps: " << elapsedSteps);
      WARN("finalParam: " << finalParam.transpose());
    }
}

