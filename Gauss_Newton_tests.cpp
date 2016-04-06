#include "catch.hpp"

#include "Gauss_Newton.h"

#include "CubicBSplineInterpolator.h"
#include "VolumeAtAddressable.h"

#include "CentralDifferenceDifferentiator.h"


TEST_CASE("a Gauss-Newton minimizer can be instantiated") {
    typedef float dataT;
    typedef VolumeAtAddressable< std::vector<dataT> > VolumeT; 
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
    
    MinimizerT minimizer(&interpolator, &dz, &dy, &dx);

    SECTION("and registering an image with itself produces 0 transformation") {
      ParamT initialParam;
      initialParam << 0, 0, 0, 0, 0, 0;
  
      ParamT finalParam;
  
      dataT paramUpdateNormLimit = 1e-6;
  
      double elapsedTime;
      size_t elapsedSteps;
  
      minimizer.minimize(&volume, &initialParam, &finalParam,
        paramUpdateNormLimit, &elapsedSteps, &elapsedTime);
 
        for(int i = 0; i < 6; i++) {
          REQUIRE(0 == Approx(finalParam(i)));
        }

      WARN("elapsed time: " << elapsedTime << " ms");
      WARN("elapsed steps: " << elapsedSteps);
      WARN("finalParam: " << finalParam.transpose());
    }
}

