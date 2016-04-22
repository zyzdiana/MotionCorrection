#include "Weighted_Gauss_Newton_Ref_Grad.h"

#include "CubicBSplineInterpolator.h"
#include "VolumeAtAddressable.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTOp.h"
#include "FFTWBuffer.h"

#include "CircularMaskOp.h"

#include "BinaryFile.h"

#include <iostream>

#include <string>
#include <sstream>

int main(int argc, char* argv[]) {
  typedef float dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT; 

  if(argc != 2) {
    std::cerr << "wrong number of arguments." << std::endl; 
    exit(1);
  }

  const size_t cubeSize = 32;
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

  
  std::string basePath(argv[1]);
  DataFFTOpT fftOp(cubeSize);

  {
    VolumeT refVolume(cubeSize);
 
    std::stringstream refPath;
    refPath << basePath << "0.dat";
    
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
  
  for(unsigned int i = 1; i < 72; i++) {
 
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
  }

  return 0;
}
