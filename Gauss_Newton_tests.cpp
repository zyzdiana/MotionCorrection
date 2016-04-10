#include "catch.hpp"

#include "Gauss_Newton.h"

#include "CubicBSplineInterpolator.h"
#include "VolumeAtAddressable.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTWBuffer.h"

#include "BinaryFile.h"

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

      const dataT paramUpdate2NormLimit = 1e-6;
      const dataT paramUpdateInfinityNormLimit = 0;
  
      double elapsedTime;
      size_t elapsedSteps;
  
      minimizer.minimize(&volume, &initialParam, &finalParam,
        maxSteps, paramUpdate2NormLimit, paramUpdateInfinityNormLimit,
        &elapsedSteps, &elapsedTime);
 
        for(int i = 0; i < 6; i++) {
          REQUIRE(0 == Approx(finalParam(i)));
        }

      WARN("elapsed time: " << elapsedTime << " ms");
      WARN("elapsed steps: " << elapsedSteps);
      WARN("finalParam: " << finalParam.transpose());
    }
}


TEST_CASE("a Gauss-Newton minimizer can be instantiated from image data") {
    typedef float dataT;
    typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
    typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT; 
    typedef Gauss_Newton<InterpolatorT> MinimizerT; 
    typedef MinimizerT::ParamT ParamT;

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    VolumeT refVolume(cubeSize);
    VolumeT newVolume(cubeSize);

    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&refVolume,
                "Gauss_Newton_tests/refVolInput.dat"));
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&newVolume,
                "Gauss_Newton_tests/newVolInput.dat"));

    InterpolatorT interpolator(&refVolume);

    CentralDifferencesDifferentiator<VolumeT> volDiffer(&refVolume);
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

    SECTION("and registering two images returns identical result to Mathematica") {
      ParamT initialParam;
      initialParam << 0, 0, 0, 0, 0, 0;
  
      ParamT finalParam;

      // We force this to go exactly 20 steps, so that we get to the same
      // point as the Mathematica code
      size_t maxSteps = 20;

      const dataT paramUpdate2NormLimit = 0;
      const dataT paramUpdateInfinityNormLimit = 0;
  
      double elapsedTime;
      size_t elapsedSteps;
  
      minimizer.minimize(&newVolume, &initialParam, &finalParam,
        maxSteps, paramUpdate2NormLimit, paramUpdateInfinityNormLimit,
        &elapsedSteps, &elapsedTime);

      std::vector<dataT> paramSolution(6);

      REQUIRE(6 * sizeof(dataT)
            == BinaryFile< std::vector<dataT> >::read(&paramSolution,
                "Gauss_Newton_tests/parameterOutput.dat"));

      for(int i = 0; i < 6; i++) {
        REQUIRE(paramSolution[i] == Approx(finalParam(i)));
      }

      WARN("elapsed time: " << elapsedTime << " ms");
      WARN("elapsed steps: " << elapsedSteps);
    }
}

