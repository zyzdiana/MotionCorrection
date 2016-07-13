#include "catch.hpp"

#include "Weighted_Gauss_Newton_New_Grad.h"

#include "CubicBSplineInterpolator.h"
#include "VolumeAtAddressable.h"

#include "TwoNormConvergenceTest.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTOp.h"
#include "FFTWBuffer.h"

#include "CircularMaskOp.h"

#include "BinaryFile.h"

TEST_CASE("a weighted Gauss-Newton minimizer using new-image gradients can be instantiated") {
  typedef float dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT; 

  const size_t cubeSize = 32;
  const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;
 
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > MaskVolumeT; 
  typedef CircularMaskOp< VolumeT, MaskVolumeT> CircularMaskOpT;
  typedef TwoNormConvergenceTest<dataT> ConvergenceTestT;
  typedef Weighted_Gauss_Newton_New_Grad<
    InterpolatorT, ConvergenceTestT, CircularMaskOpT > MinimizerT; 
  typedef MinimizerT::ParamT ParamT;
 
  VolumeT volume(cubeSize);
  
  for(size_t i = 0; i < cubeVectorLength; i++) {
      volume.at(i) = ((dataT) i) / ((dataT) cubeVectorLength); 
  }

  InterpolatorT interpolator(&volume);
  
  double gradientAndHessianComputeTime;

  CentralDifferencesDifferentiator<VolumeT> volDiffer(&volume);
  VolumeT dx(cubeSize, cubeVectorLength);
  volDiffer.xDerivative(&dx);

  VolumeT dy(cubeSize, cubeVectorLength);
  volDiffer.yDerivative(&dy);

  VolumeT dz(cubeSize, cubeVectorLength);
  volDiffer.zDerivative(&dz);

  CircularMaskOpT imageMaskOp(cubeSize);
  
  MinimizerT minimizer(&interpolator, cubeSize, &imageMaskOp);

  SECTION("and registering an image with itself produces 0 transformation") {
    ParamT initialParam;
    initialParam << 0, 0, 0, 0, 0, 0;

    ParamT finalParam;

    size_t maxSteps = 20;
    const dataT stepSizeScale = 0.25;
    const dataT stepSizeLimit = 1e-5;

    const dataT paramUpdate2NormLimit = 1e-6;
    ConvergenceTestT convergenceTest(paramUpdate2NormLimit);

    double elapsedTime;
    size_t elapsedSteps;

    minimizer.minimize(&volume, &dz, &dy, &dx,
      &initialParam, &finalParam,
      maxSteps, stepSizeScale, stepSizeLimit,
      &convergenceTest, 
      &elapsedSteps, &elapsedTime, &gradientAndHessianComputeTime);
    
    WARN("elapsed time computing gradient and Hessian: "
      << gradientAndHessianComputeTime << " ms");

    WARN("elapsed time: " << elapsedTime << " ms");
    WARN("elapsed steps: " << elapsedSteps);
    WARN("finalParam: " << finalParam.transpose());

    for(int i = 0; i < 6; i++) {
      REQUIRE(0 == Approx(finalParam(i)));
    }
  } 
}


TEST_CASE("a weighted Gauss-Newton minimizer using new-image gradients can be instantiated from image data") {
  typedef double dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef CubicBSplineInterpolator<VolumeT, dataT> InterpolatorT; 

  const size_t cubeSize = 32;
  const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;


  typedef std::complex<dataT> complexT;
  typedef FFTOp<dataT> DataFFTOpT;
  typedef DataFFTOpT::spatialVolumeT DataVolumeT; 
  typedef DataFFTOpT::fourierVolumeT ComplexVolumeT; 
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > ImageMaskVolumeT; 
  typedef CircularMaskOp< DataVolumeT, ImageMaskVolumeT> DataCircularMaskOpT;
  typedef TwoNormConvergenceTest<dataT> ConvergenceTestT;
  typedef SymmetricHalfVolumeAtAddressable< FFTWBuffer<dataT> >
    FourierMaskVolumeT; 
  typedef CircularMaskOp< ComplexVolumeT, FourierMaskVolumeT>
    ComplexDataCircularMaskOpT;
  typedef Weighted_Gauss_Newton_New_Grad<
    InterpolatorT, ConvergenceTestT, DataCircularMaskOpT > MinimizerT; 
  typedef MinimizerT::ParamT ParamT;
  
  ComplexDataCircularMaskOpT fourierMaskOp(cubeSize);
  DataVolumeT maskedRefVolume(cubeSize);
  DataVolumeT maskedNewVolume(cubeSize);

  {
    VolumeT refVolume(cubeSize);
    VolumeT newVolume(cubeSize);
  
    ComplexVolumeT fourierData(cubeSize);
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&refVolume,
                "Weighted_Gauss_Newton_New_Grad_tests/refVolInput.dat"));
    
    DataFFTOpT fftOp(cubeSize);
  
    fftOp.forward(&refVolume, &fourierData);
    
    fourierMaskOp.applyMask(&fourierData); 
        
    fftOp.backward(&fourierData, &maskedRefVolume);
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&newVolume,
                "Weighted_Gauss_Newton_New_Grad_tests/newVolInput.dat"));
    
    fftOp.forward(&newVolume, &fourierData);
    
    fourierMaskOp.applyMask(&fourierData); 
        
    fftOp.backward(&fourierData, &maskedNewVolume);
  }

  InterpolatorT interpolator(&maskedRefVolume);

  CentralDifferencesDifferentiator<VolumeT> volDiffer(&maskedNewVolume);
  VolumeT dx(cubeSize, cubeVectorLength);
  volDiffer.xDerivative(&dx);

  VolumeT dy(cubeSize, cubeVectorLength);
  volDiffer.yDerivative(&dy);

  VolumeT dz(cubeSize, cubeVectorLength);
  volDiffer.zDerivative(&dz);
  
  double gradientAndHessianComputeTime;
  
  DataCircularMaskOpT imageMaskOp(cubeSize);

  MinimizerT minimizer(&interpolator, cubeSize, &imageMaskOp);
        
  WARN("elapsed time computing gradient and Hessian: "
    << gradientAndHessianComputeTime << " ms");

  SECTION("and registering two images returns identical result to Mathematica") {
    ParamT initialParam;
    initialParam << 0, 0, 0, 0, 0, 0;
  
    ParamT finalParam;

    // We force this to stop as soon as we take a bad step, so that we get to
    // the same point as the Mathematica code
    size_t maxSteps = 20;
    const dataT stepSizeScale = 0.25;
    const dataT stepSizeLimit = 1.0;
 
    double elapsedTime;
    size_t elapsedSteps;
  
    minimizer.minimize(&maskedNewVolume, &dz, &dy, &dx,
      &initialParam, &finalParam,
      maxSteps, stepSizeScale, stepSizeLimit,
      NULL, 
      &elapsedSteps, &elapsedTime);

    std::vector<dataT> paramSolution(6);

    WARN("elapsed time: " << elapsedTime << " ms");
    WARN("elapsed steps: " << elapsedSteps);
    WARN("finalParam: " << finalParam.transpose());

    REQUIRE(6 * sizeof(dataT)
          == BinaryFile< std::vector<dataT> >::read(&paramSolution,
              "Weighted_Gauss_Newton_New_Grad_tests/parameterOutput.dat"));

    for(int i = 0; i < 6; i++) {
      REQUIRE(paramSolution[i] == Approx(finalParam(i)).epsilon(0.001));
    }
  }
}

