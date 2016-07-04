#include "Weighted_Gauss_Newton_Ref_Grad.h"
#include "Weighted_Gauss_Newton_New_Grad.h"

#include "TricubicInterpolator.h"
#include "CubicBSplineInterpolator.h"
#include "VolumeAtAddressable.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTOp.h"
#include "FFTWBuffer.h"

#include "CircularMaskOp.h"

#include "BinaryFile.h"

#include <iostream>
#include <fstream>

#include <string>
#include <sstream>

#include <tclap/CmdLine.h>

int main(int argc, char* argv[]) {
  typedef float dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef TricubicInterpolator<VolumeT, float> TricubicInterpolatorT; 
  typedef CubicBSplineInterpolator<VolumeT, float> CubicBSplineInterpolatorT; 
  
  size_t cubeSize;
  std::string basePath;
  std::string outputPath;
  dataT translationScaleMM;
  dataT rotationScaleMM;

  try {
    TCLAP::CmdLine cmd("Registering vNav volumes", ' ', "dev");

    TCLAP::ValueArg<std::string> inputFileArg("i", "input", "Path to a volume data, including file name, with the final \"xx.dat\" removed.", true, "", "path", cmd);
    TCLAP::ValueArg<std::string> outputFileArg("o", "output", "Path to output file that will be written.", true, "", "path", cmd);
    TCLAP::ValueArg<unsigned int> widthArg("", "width", "Number of voxels along the side of the vNav", true, 32, "integer", cmd);
    TCLAP::ValueArg<dataT> resolutionArg("", "res", "Size of a voxel's side in mm", true, 8, "mm", cmd);

    cmd.parse(argc, argv);

    basePath = inputFileArg.getValue();
    outputPath = outputFileArg.getValue();

    cubeSize = widthArg.getValue();

    translationScaleMM = resolutionArg.getValue();
    rotationScaleMM = 100.0;

  }   catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(1);
  } 

  const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

  typedef std::complex<float> complexT;
  typedef FFTOp<dataT> DataFFTOpT;
  typedef DataFFTOpT::spatialVolumeT DataVolumeT; 
  typedef DataFFTOpT::fourierVolumeT ComplexVolumeT; 
  typedef CircularMaskOp< DataVolumeT, DataVolumeT> DataCircularMaskOpT;
  typedef CircularMaskOp< ComplexVolumeT, 
    SymmetricHalfVolumeAtAddressable< FFTWBuffer<dataT> > >
    ComplexDataCircularMaskOpT;
  typedef Weighted_Gauss_Newton_Ref_Grad<
    TricubicInterpolatorT, DataCircularMaskOpT > TricubicRefGradMinimizerT; 
  typedef Weighted_Gauss_Newton_Ref_Grad<
    CubicBSplineInterpolatorT, DataCircularMaskOpT > CubicBSplineRefGradMinimizerT; 
  typedef Weighted_Gauss_Newton_New_Grad<
    TricubicInterpolatorT, DataCircularMaskOpT > TricubicNewGradMinimizerT; 
  typedef Weighted_Gauss_Newton_New_Grad<
    CubicBSplineInterpolatorT, DataCircularMaskOpT > CubicBSplineNewGradMinimizerT; 
  typedef CubicBSplineRefGradMinimizerT::ParamT ParamT;
  
  ComplexDataCircularMaskOpT fourierMaskOp(cubeSize);
  DataVolumeT maskedRefVolume(cubeSize);
  ComplexVolumeT fourierData(cubeSize);

  
  DataFFTOpT fftOp(cubeSize);
  std::ofstream outputFile(outputPath);

  size_t step = 0;

  for(int baseIndex = 0; baseIndex <= 36; baseIndex += 36) {

    {
      VolumeT refVolume(cubeSize);
  
      std::stringstream refPath;
      refPath << basePath << baseIndex << ".dat";
      
      if(cubeVectorLength * sizeof(dataT)
              != BinaryFile<VolumeT>::read(&refVolume, refPath.str())) {
        std::cerr << "read incorrect length from file: " << refPath.str()
          << std::endl;
        exit(1);
      }
       
      fftOp.forward(&refVolume, &fourierData);

      fourierMaskOp.applyMask(&fourierData); 
          
      fftOp.backward(&fourierData, &maskedRefVolume);
    }

    BinaryFile<VolumeT>::write(&maskedRefVolume,
      "debug_filtered_ref_volume.dat");

    CubicBSplineInterpolatorT cubicBSplineInterpolator(&maskedRefVolume);
  
    CentralDifferencesDifferentiator<VolumeT> volDiffer(&maskedRefVolume);
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
     
    TricubicInterpolatorT tricubicInterpolator(&maskedRefVolume,
      &dx, &dy, &dz, &dxy, &dxz, &dyz, &dxyz);

    DataCircularMaskOpT imageMaskOp(cubeSize);
  
    CubicBSplineRefGradMinimizerT cubicBSplineRefGradMinimizer(
      &cubicBSplineInterpolator, &imageMaskOp,
      &dz, &dy, &dx, NULL);
    
    TricubicRefGradMinimizerT tricubicRefGradMinimizer(
      &tricubicInterpolator, &imageMaskOp,
      &dz, &dy, &dx, NULL);
    
    CubicBSplineNewGradMinimizerT cubicBSplineNewGradMinimizer(
      &cubicBSplineInterpolator, cubeSize, &imageMaskOp);
    
    TricubicNewGradMinimizerT tricubicNewGradMinimizer(
      &tricubicInterpolator, cubeSize, &imageMaskOp);
          
    DataVolumeT maskedNewVolume(cubeSize);
    VolumeT newVolume(cubeSize);
    ParamT initialParam;
    ParamT finalParam;
 
    for(unsigned int i = baseIndex + 1; i < baseIndex + 36; i++, step++) { 
      std::stringstream newPath;
      newPath << basePath << i << ".dat";
      
      if(cubeVectorLength * sizeof(dataT)
              != BinaryFile<VolumeT>::read(&newVolume, newPath.str())) {
        std::cerr << "read incorrect length from file: " << newPath.str()
          << std::endl;
        exit(1);
      }
      
      fftOp.forward(&newVolume, &fourierData);
      
      fourierMaskOp.applyMask(&fourierData); 
          
      fftOp.backward(&fourierData, &maskedNewVolume);
      
      initialParam << 0, 0, 0, 0, 0, 0;
      
      const size_t maxSteps = 50;
      const dataT stepSizeScale = 0.25;
      const dataT stepSizeLimit = 1e-5;
  
      const dataT paramUpdate2NormLimit = 0;
      //const dataT paramUpdateInfinityNormLimit = 1e-6;
      const dataT paramUpdateInfinityNormLimit = 0; 
      const dataT paramUpdateMMLimit = 0.01;
      
      std::cout << "----------" << std::endl;
      std::cout << "step " << step << std::endl;
    
      { 
        double tricubicRefGradElapsedTime;
        size_t tricubicRefGradElapsedSteps;
        
        tricubicRefGradMinimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          paramUpdate2NormLimit, paramUpdateInfinityNormLimit,
          paramUpdateMMLimit, translationScaleMM, rotationScaleMM,
          &tricubicRefGradElapsedSteps, &tricubicRefGradElapsedTime);
        
        outputFile << tricubicRefGradElapsedTime << " " << tricubicRefGradElapsedSteps
          << " " << finalParam.transpose() << std::endl;
      
        std::cout << "tricubic ref grad elapsed time: " << tricubicRefGradElapsedTime << " ms" << std::endl;
        std::cout << "tricubic ref grad elapsed steps: " << tricubicRefGradElapsedSteps << std::endl;
      }
      
      {
        double cubicBSplineRefGradElapsedTime;
        size_t cubicBSplineRefGradElapsedSteps;
        
        cubicBSplineRefGradMinimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          paramUpdate2NormLimit, paramUpdateInfinityNormLimit,
          paramUpdateMMLimit, translationScaleMM, rotationScaleMM,
          &cubicBSplineRefGradElapsedSteps, &cubicBSplineRefGradElapsedTime); 
        
        outputFile << cubicBSplineRefGradElapsedTime << " " << cubicBSplineRefGradElapsedSteps
          << " " << finalParam.transpose() << std::endl;

        std::cout << "cubic B-spline ref grad elapsed time: " << cubicBSplineRefGradElapsedTime << " ms" << std::endl;
        std::cout << "cubic B-spline ref grad elapsed steps: " << cubicBSplineRefGradElapsedSteps << std::endl;
      }
      
      CentralDifferencesDifferentiator<VolumeT> newVolDiffer(&maskedNewVolume);
      VolumeT newVolDx(cubeSize, cubeVectorLength);
      newVolDiffer.xDerivative(&newVolDx);
  
      VolumeT newVolDy(cubeSize, cubeVectorLength);
      newVolDiffer.yDerivative(&newVolDy);
  
      VolumeT newVolDz(cubeSize, cubeVectorLength);
      newVolDiffer.zDerivative(&newVolDz);
      
      { 
        double tricubicNewGradElapsedTime;
        size_t tricubicNewGradElapsedSteps;
        
        tricubicNewGradMinimizer.minimize(&maskedNewVolume,
          &newVolDz, &newVolDy, &newVolDx, 
          &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          paramUpdate2NormLimit, paramUpdateInfinityNormLimit,
          paramUpdateMMLimit, translationScaleMM, rotationScaleMM,
          &tricubicNewGradElapsedSteps, &tricubicNewGradElapsedTime);
        
        outputFile << tricubicNewGradElapsedTime << " " << tricubicNewGradElapsedSteps
          << " " << finalParam.transpose() << std::endl;

        std::cout << "tricubic new grad elapsed time: " << tricubicNewGradElapsedTime << " ms" << std::endl;
        std::cout << "tricubic new grad elapsed steps: " << tricubicNewGradElapsedSteps << std::endl;
      }
      
      {
        double cubicBSplineNewGradElapsedTime;
        size_t cubicBSplineNewGradElapsedSteps;
        
        cubicBSplineNewGradMinimizer.minimize(&maskedNewVolume,
          &newVolDz, &newVolDy, &newVolDx, 
          &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          paramUpdate2NormLimit, paramUpdateInfinityNormLimit,
          paramUpdateMMLimit, translationScaleMM, rotationScaleMM,
          &cubicBSplineNewGradElapsedSteps, &cubicBSplineNewGradElapsedTime); 
        
        outputFile << cubicBSplineNewGradElapsedTime << " " << cubicBSplineNewGradElapsedSteps
          << " " << finalParam.transpose() << std::endl;

        std::cout << "cubic B-spline new grad elapsed time: " << cubicBSplineNewGradElapsedTime << " ms" << std::endl;
        std::cout << "cubic B-spline new grad elapsed steps: " << cubicBSplineNewGradElapsedSteps << std::endl;
      }

    }
  }

  return 0;
}
