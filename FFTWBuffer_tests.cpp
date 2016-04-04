#include "catch.hpp"

#include <complex>

#include "FFTWBuffer.h"

TEST_CASE("a volume of floats can be created from a vector") {
    typedef float dataT;
    typedef std::complex<float> complexT;
    typedef FFTWBuffer< dataT> FFTWBufferT; 

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    FFTWBufferT initialData(cubeVectorLength);
    FFTWBufferT fourierData(cubeSize * cubeSize * (cubeSize / 2 + 1));

/*
    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = i; 
    }

    VolumeT volume(cubeSize, initialData);

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
