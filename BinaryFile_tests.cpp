#include "catch.hpp"

#include "ReadFile.h"
#include "BinaryFile.h"
#include "Interpolation_3D.h"

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

TEST_CASE("rotation of zero degrees returns original volume for trilinear interpolation") {

    typedef Interpolation_3D<float> interp;
    typedef interp::vector4D vector4D;
    typedef interp::vector3D vector3D;
    typedef interp::vector2D vector2D;
    typedef interp::vector1D vector1D;
    typedef std::complex<float> dataT;

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;
    float ox,oy,oz = float(num_slice-1)/2;
    vector3D vol1 = ReadFile<float>::read_volume(path, num_slice);

    float theta = 0.0;
    for(int wx = 0; wx < 2; ++wx){
        for(int wy = 0; wy < 2; ++wy){
            for(int wz = 0; wz < 2; ++wz){
                vector3D dest = interp::rotate_volume_trilinear(vol1,theta,wx,wy,wz);
                SECTION("the rotation with trilinear is identical") {
                    for(int i = 0; i < num_slice; ++i){
                        for(int j = 0; j < num_slice; ++j){
                            for(int k = 0; k < num_slice; ++k){
                                REQUIRE(dest[i][j][k] == vol1[i][j][k]);
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("rotation of 90 degrees returns original volume for trilinear interpolation") {

    typedef Interpolation_3D<float> interp;
    typedef interp::vector4D vector4D;
    typedef interp::vector3D vector3D;
    typedef interp::vector2D vector2D;
    typedef interp::vector1D vector1D;
    typedef std::complex<float> dataT;

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;
    float ox,oy,oz = float(num_slice-1)/2;
    vector3D vol1 = ReadFile<float>::read_volume(path, num_slice);

    float theta = 90.0;
    SECTION("rotation around x axis") {
        float wx = 1;
        float wy = 0;
        float wz = 0;
        vector3D dest = interp::rotate_volume_trilinear(vol1,theta,wx,wy,wz);
        string ref_path = "/Users/zyzdiana/GitHub/AC297r-Volume-Registration/Thesis/rot90_x.dat";
        vector3D ref = ReadFile<float>::read_reference_volume(ref_path, num_slice);
        for(int i = 0; i < num_slice; ++i){
            for(int j = 0; j < num_slice; ++j){
                for(int k = 0; k < num_slice; ++k){
                    REQUIRE(dest[i][j][k] == vol1[i][j][k]);
                }
            }
        }
    }
    SECTION("rotation around y axis") {
        float wx = 0;
        float wy = 1;
        float wz = 0;
        vector3D dest = interp::rotate_volume_trilinear(vol1,theta,wx,wy,wz);
        string ref_path = "/Users/zyzdiana/GitHub/AC297r-Volume-Registration/Thesis/rot90_y.dat";
        vector3D ref = ReadFile<float>::read_reference_volume(ref_path, num_slice);
        for(int i = 0; i < num_slice; ++i){
            for(int j = 0; j < num_slice; ++j){
                for(int k = 0; k < num_slice; ++k){
                    REQUIRE(dest[i][j][k] == vol1[i][j][k]);
                }
            }
        }
    }
    SECTION("rotation around z axis") {
        float wx = 0;
        float wy = 0;
        float wz = 1;
        vector3D dest = interp::rotate_volume_trilinear(vol1,theta,wx,wy,wz);
        string ref_path = "/Users/zyzdiana/GitHub/AC297r-Volume-Registration/Thesis/rot90_z.dat";
        vector3D ref = ReadFile<float>::read_reference_volume(ref_path, num_slice);
        for(int i = 0; i < num_slice; ++i){
            for(int j = 0; j < num_slice; ++j){
                for(int k = 0; k < num_slice; ++k){
                    REQUIRE(dest[i][j][k] == vol1[i][j][k]);
                }
            }
        }
    }
}

// TEST_CASE("rotation of zero degrees returns original volume for tricubic interpolation") {

//     typedef Interpolation_3D<float> interp;
//     typedef interp::vector4D vector4D;
//     typedef interp::vector3D vector3D;
//     typedef interp::vector2D vector2D;
//     typedef interp::vector1D vector1D;
//     typedef std::complex<float> dataT;

//     string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
//     int num_slice = 32;
//     float ox,oy,oz = float(num_slice-1)/2;
//     vector3D vol1 = ReadFile<float>::read_volume(path, num_slice);

//     float theta = 0.0;
//     for(int wx = 0; wx < 2; ++wx){
//         for(int wy = 0; wy < 2; ++wy){
//             for(int wz = 0; wz < 2; ++wz){
//                 vector3D dest = interp::rotate_volume_tricubic(vol1,theta,wx,wy,wz);
//                 SECTION("the rotation with trilinear is identical") {
//                     for(int i = 0; i < num_slice; ++i){
//                         for(int j = 0; j < num_slice; ++j){
//                             for(int k = 0; k < num_slice; ++k){
//                                 REQUIRE(dest[i][j][k] == vol1[i][j][k]);
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }


