#include "catch.hpp"

#include "TricubicInterpolator.h"

#include "Interpolator3D_tests.h"

#include <vector>
#include <complex>

TEST_CASE("a tricubic interpolator can be created from a volume") {
    //typedef std::complex<float> make dataT;
    typedef float dataT;
    typedef Volume<dataT, std::vector<dataT>, float > VolumeT;
    typedef TricubicInterpolator<VolumeT, float> InterpolatorT;

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = i; 
    }

    VolumeT volume(initialData, cubeSize);

    VolumeT dx(cubeSize);
    VolumeT dy(cubeSize);
    VolumeT dz(cubeSize);
    VolumeT dxy(cubeSize);
    VolumeT dxz(cubeSize);
    VolumeT dyz(cubeSize);
    VolumeT dxyz(cubeSize);

    InterpolatorT interpolator(&volume, &dx, &dy, &dz, &dxy, &dxz, &dyz, &dxyz);

    InterpolatorTests<InterpolatorT>::tests(&interpolator, &volume); 
}
