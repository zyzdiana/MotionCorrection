#include "catch.hpp"

#include "Volume.h"

#include <vector>
#include <complex>

TEST_CASE("a volume of complex floats can be created from a vector") {
    typedef std::complex<float> dataT;
   
    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = i; 
    }

    Volume<dataT, std::vector<dataT> > volume(initialData, cubeSize);

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
}
