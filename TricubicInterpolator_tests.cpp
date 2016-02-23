#include "catch.hpp"

#include "TricubicInterpolator.h"

#include <vector>
#include <complex>

TEST_CASE("a tricubic interpolator can be created from a volume") {
    //typedef std::complex<float> make dataT;
    typedef float dataT;

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = i; 
    }

    typedef Volume<dataT, std::vector<dataT> > VolumeT;
    VolumeT volume(initialData, cubeSize);

    TricubicInterpolator<VolumeT, float> interpolator(&volume);

    SECTION("and all the points are returned when interpolated exactly") { 
        for(size_t z = 0; z < cubeSize; z++) {
            for(size_t y = 0; y < cubeSize; y++) {
                for(size_t x = 0; x < cubeSize; x++) {
                    REQUIRE(interpolator.interp(z, y, x) == volume.at(z, y, x));
                }
            }
        }
    }
}
