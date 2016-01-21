#include "catch.hpp"

#include "BinaryFile.h"

#include <vector>
#include <complex>

TEST_CASE("a vector of complex floats can be written") {
    typedef std::complex<float> dataT;

    std::vector< dataT > initialData(100);
    
    for(size_t i = 0; i < 100; i++) {
        initialData[i] = dataT(i, -i); 
    }

    REQUIRE(100 * sizeof(dataT)
        == BinaryFile<dataT>::write(&initialData, "temp.dat"));

    SECTION("and read back into a new vector") {
        std::vector< dataT > readData(100);

        REQUIRE(100 * sizeof(dataT)
            == BinaryFile<dataT>::read(&readData, "temp.dat"));

        SECTION("and the data is identical") {
            for(size_t i = 0; i < 100; i++) {
                REQUIRE(initialData[i] == readData[i]); 
            }
        }
    }
}
