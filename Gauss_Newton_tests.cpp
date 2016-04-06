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

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    VolumeT volume(cubeSize);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        volume.at(i) = i; 
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

    ParamT initialParam;
    initialParam << 0, 0, 0, 0, 0, 0;

    ParamT finalParam;

    minimizer.minimize(&volume, &initialParam, &finalParam);

/*
    SECTION("and all the values can be read back linearly") {
        for(size_t i = 0; i < cubeVectorLength; i++) {
            REQUIRE(initialData[i] == volume.at(i)); 
        }
    }
    
    SECTION("and all the values can be read back via coordinates") {
        size_t i = 0;

        for(size_t z = 0; z < cubeSize; z++) {
            for(size_t y = 0; y < cubeSize; y++) {
                for(size_t x = 0; x < cubeSize; x++) {
                    REQUIRE(initialData[i] == volume.at(z, y, x));
                    i++;
                }
            }
        }
    }
*/
}

