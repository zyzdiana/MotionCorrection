#include "catch.hpp"

#include "TrilinearInterpolator.h"

#include <vector>
#include <complex>

TEST_CASE("a trilinear interpolator can be created from a volume") {
    typedef std::complex<float> dataT;

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = i; 
    }

    typedef Volume<dataT, std::vector<dataT> > VolumeT;
    VolumeT volume(initialData, cubeSize);

    TrilinearInterpolator<VolumeT, float> interpolator(&volume);

    SECTION("and all the points are returned when interpolated exactly") { 
        for(size_t z = 0; z < cubeSize; z++) {
            for(size_t y = 0; y < cubeSize; y++) {
                for(size_t x = 0; x < cubeSize; x++) {
                    REQUIRE(interpolator.interp(z, y, x) == volume.at(z, y, x));
                }
            }
        }
    }
    
    SECTION("and all the points are averaged when interpolated in x") { 
        for(size_t z = 0; z < cubeSize; z++) {
            for(size_t y = 0; y < cubeSize; y++) {
                for(size_t x = 0; x < cubeSize - 1; x++) {
                    dataT avg =
                        (volume.at(z, y, x) + volume.at(z, y, x + 1))
                        * (float) 0.5;
                    REQUIRE(interpolator.interp(z, y, x + 0.5) == avg);
                }
            }
        }
    }
    
    SECTION("and all the points are averaged when interpolated in y") { 
        for(size_t z = 0; z < cubeSize; z++) {
            for(size_t y = 0; y < cubeSize - 1; y++) {
                for(size_t x = 0; x < cubeSize; x++) {
                    dataT avg =
                        (volume.at(z, y, x) + volume.at(z, y + 1, x))
                        * (float) 0.5;
                    REQUIRE(interpolator.interp(z, y + 0.5, x) == avg);
                }
            }
        }
    }
    
    SECTION("and all the points are averaged when interpolated in z") { 
        for(size_t z = 0; z < cubeSize - 1; z++) {
            for(size_t y = 0; y < cubeSize; y++) {
                for(size_t x = 0; x < cubeSize; x++) {
                    dataT avg =
                        (volume.at(z, y, x) + volume.at(z + 1, y, x))
                        * (float) 0.5;
                    REQUIRE(interpolator.interp(z + 0.5, y, x) == avg);
                }
            }
        }
    }
}
