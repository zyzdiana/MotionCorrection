#include "Weighted_Gauss_Newton_Ref_Grad.h"

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
  typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT; 
  
  size_t cubeSize;
  std::string basePath;
  std::string outputPath;
  
  try {
    TCLAP::CmdLine cmd("Registering vNav volumes", ' ', "dev");

    TCLAP::ValueArg<std::string> inputFileArg("i", "input", "Path to a volume data, including file name, with the final \"xx.dat\" removed.", true, "", "path", cmd);
    TCLAP::ValueArg<std::string> outputFileArg("o", "output", "Path to output file that will be written.", true, "", "path", cmd);
    TCLAP::ValueArg<unsigned int> widthArg("", "width", "Number of voxels along the side of the vNav", true, 32, "integer", cmd);

    cmd.parse(argc, argv);

    basePath = inputFileArg.getValue();
    outputPath = outputFileArg.getValue();

    cubeSize = widthArg.getValue();

  }   catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(1);
  } 

  const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

  typedef std::complex<float> complexT;
  typedef FFTOp<dataT> DataFFTOpT;
  typedef DataFFTOpT::spatialVolumeT DataVolumeT; 
  typedef DataFFTOpT::fourierVolumeT ComplexVolumeT; 
  typedef CircularMaskOp< DataVolumeT, dataT> DataCircularMaskOpT;
  typedef CircularMaskOp< ComplexVolumeT, dataT> ComplexDataCircularMaskOpT;
  typedef Weighted_Gauss_Newton_Ref_Grad<
    InterpolatorT, DataCircularMaskOpT > MinimizerT; 
  typedef Weighted_Gauss_Newton_Ref_Grad<
    InterpolatorT, DataCircularMaskOpT> MinimizerT; 
  typedef MinimizerT::ParamT ParamT;
  
  ComplexDataCircularMaskOpT fourierMaskOp(cubeSize);
  DataVolumeT maskedRefVolume(cubeSize);
  ComplexVolumeT fourierData(cubeSize);

  
  DataFFTOpT fftOp(cubeSize);
  std::ofstream outputFile(outputPath);

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


    InterpolatorT interpolator(&maskedRefVolume);
  
    CentralDifferencesDifferentiator<VolumeT> volDiffer(&maskedRefVolume);
    VolumeT dx(cubeSize, cubeVectorLength);
    volDiffer.xDerivative(&dx);
  
    VolumeT dy(cubeSize, cubeVectorLength);
    volDiffer.yDerivative(&dy);
  
    VolumeT dz(cubeSize, cubeVectorLength);
    volDiffer.zDerivative(&dz);
    
    double gradientAndHessianComputeTime;
    
    DataCircularMaskOpT imageMaskOp(cubeSize);
  
    MinimizerT minimizer(&interpolator, &imageMaskOp,
      &dz, &dy, &dx,
      &gradientAndHessianComputeTime);
          
    DataVolumeT maskedNewVolume(cubeSize);
    VolumeT newVolume(cubeSize);
    ParamT initialParam;
    ParamT finalParam;
    
    for(unsigned int i = baseIndex + 1; i < baseIndex + 36; i++) { 
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
      
      size_t maxSteps = 20;
  
      const dataT paramUpdate2NormLimit = 0;
      const dataT paramUpdateInfinityNormLimit = 1e-6;
     
      double elapsedTime;
      size_t elapsedSteps;
      
      minimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
        maxSteps, paramUpdate2NormLimit, paramUpdateInfinityNormLimit,
        &elapsedSteps, &elapsedTime);
     
      std::cout << "----------" << std::endl;
      std::cout << "step " << i << std::endl;
      std::cout << "elapsed time: " << elapsedTime << " ms" << std::endl;
      std::cout << "elapsed steps: " << elapsedSteps << std::endl;
      std::cout << "finalParam: " << finalParam.transpose() << std::endl;
      outputFile << elapsedTime << " " << elapsedSteps << " " << finalParam.transpose() << std::endl;
    }
  }

  return 0;
}
