#include "catch.hpp"
#include "ReadFile.h"
#include "BinaryFile.h"
#include "TrilinearInterpolator.h"
#include "TricubicInterpolator.h"
#include "Utils.h"

#include <math.h>
#include <vector>
#include <complex>

/* Test cases for Trilinear Interpolation */
TEST_CASE("rotation of zero degrees using trilinear interpolation returns itself") {

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<float> dataVector(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector, path, cubeSize);

    typedef Volume<float, std::vector<float> > VolumeT;
    VolumeT volume(dataVector, cubeSize);

    TrilinearInterpolator<VolumeT, float> interpolator(&volume);


    float theta = 0.0;
    for(int wx = 0; wx < 2; ++wx){
        for(int wy = 0; wy < 2; ++wy){
            for(int wz = 0; wz < 2; ++wz){

                for(size_t z = 0; z < cubeSize; z++) {
                    for(size_t y = 0; y < cubeSize; y++) {
                        for(size_t x = 0; x < cubeSize; x++) {
                            REQUIRE(interpolator.interp(z, y, x) == volume.at(z, y, x));
                        }
                    }
                }
                
            }
        }
    }
}

TEST_CASE("rotation of zero degrees returns original volume for tricubic interpolation") {

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<float> dataVector(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector, path, cubeSize);

    typedef Volume<float, std::vector<float> > VolumeT;
    VolumeT volume(dataVector, cubeSize);

    TricubicInterpolator<VolumeT, float> interpolator(&volume);


    float theta = 0.0;
    for(int wx = 0; wx < 2; ++wx){
        for(int wy = 0; wy < 2; ++wy){
            for(int wz = 0; wz < 2; ++wz){

                for(size_t z = 0; z < cubeSize; z++) {
                    for(size_t y = 0; y < cubeSize; y++) {
                        for(size_t x = 0; x < cubeSize; x++) {
                            REQUIRE(interpolator.interp(z, y, x) == volume.at(z, y, x));
                        }
                    }
                }
                
            }
        }
    }
}

TEST_CASE("Testing Rotation of 90 and 180 degrees"){

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<float> dataVector(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector, path, cubeSize);

    typedef Volume<float, std::vector<float> > VolumeT;
    VolumeT volume(dataVector, cubeSize);
    Vector3f cubeCenter = volume.cubeCenter;

    TrilinearInterpolator<VolumeT, float> linear_interpolator(&volume);
    TricubicInterpolator<VolumeT, float> cubic_interpolator(&volume);


    for(int wx = 0; wx < 2; ++wx){
        std::stringstream sswx;
        sswx << wx;

        for(int wy = 0; wy < 2; ++wy){
            std::stringstream sswy;
            sswy << wy;

            for(int wz = 0; wz < 2; ++wz){
                std::stringstream sswz;
                sswz << wz;

                SECTION("Rotation of 90 degrees using trilinear interpolation returns same volume as python"){
                    float theta = 90.0;
                    Matrix3f R;
                    R = Utils::get_rotation_matrix(theta, wx, wy, wz);
                    Vector3f dest_coords;
                
                    std::stringstream sstheta;
                    sstheta << (int) theta;

                    string path_ref =
                      "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/linear_Vol1_" +
                      sstheta.str() +
                      sswx.str() +
                      sswy.str() + 
                      sswz.str() + ".dat";
                    std::vector<float> ref_vector = ReadFile<float>::read_reference_volume(path_ref, cubeVectorLength);
                    VolumeT vol_ref(ref_vector, cubeSize);

                    for(size_t z = 0; z < cubeSize; z++) {
                        for(size_t y = 0; y < cubeSize; y++) {
                            for(size_t x = 0; x < cubeSize; x++) {

                                dest_coords = R*(Vector3f(y,x,z)-cubeCenter)+cubeCenter;
                                //float tmp = linear_interpolator.interp(dest_coords(1), dest_coords(0), dest_coords(2));
                                REQUIRE(abs(linear_interpolator.interp(dest_coords(1), dest_coords(0), dest_coords(2)) - vol_ref.at(x, y, z)) <= 1e-5);
                            }
                        }
                    }
                }

                SECTION("Rotation of 180 degrees using trilinear interpolation returns same volume as python"){
                    float theta = 180.0;
                    Matrix3f R;
                    R = Utils::get_rotation_matrix(theta, wx, wy, wz);
                    Vector3f dest_coords;

                    std::stringstream sstheta;
                    sstheta << (int) theta;

                    string path_ref =
                      "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/linear_Vol1_" +
                      sstheta.str() +
                      sswx.str() +
                      sswy.str() + 
                      sswz.str() + ".dat";
                    std::vector<float> ref_vector = ReadFile<float>::read_reference_volume(path_ref, cubeVectorLength);
                    VolumeT vol_ref(ref_vector, cubeSize);

                    for(size_t z = 0; z < cubeSize; z++) {
                        for(size_t y = 0; y < cubeSize; y++) {
                            for(size_t x = 0; x < cubeSize; x++) {

                                dest_coords = R*(Vector3f(y,x,z)-cubeCenter)+cubeCenter;
                                //float tmp = linear_interpolator.interp(dest_coords(1), dest_coords(0), dest_coords(2));
                                REQUIRE(abs(linear_interpolator.interp(dest_coords(1), dest_coords(0), dest_coords(2)) - vol_ref.at(x, y, z)) <= 1e-5);
                            }
                        }
                    } 
                }

                SECTION("Rotation of 90 degrees using tricubic interpolation returns same volume as python"){
                    float theta = 90.0;
                    Matrix3f R;
                    R = Utils::get_rotation_matrix(theta, wx, wy, wz);
                    Vector3f dest_coords;

                    std::stringstream sstheta;
                    sstheta << (int) theta;

                    string path_ref =
                      "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/linear_Vol1_" +
                      sstheta.str() +
                      sswx.str() +
                      sswy.str() + 
                      sswz.str() + ".dat";
                    std::vector<float> ref_vector = ReadFile<float>::read_reference_volume(path_ref, cubeVectorLength);
                    VolumeT vol_ref(ref_vector, cubeSize);

                    for(size_t z = 0; z < cubeSize; z++) {
                        for(size_t y = 0; y < cubeSize; y++) {
                            for(size_t x = 0; x < cubeSize; x++) {
                                dest_coords = R*(Vector3f(y,x,z)-cubeCenter)+cubeCenter;
                                //float tmp = linear_interpolator.interp(dest_coords(1), dest_coords(0), dest_coords(2));
                                REQUIRE(abs(cubic_interpolator.interp(dest_coords(1), dest_coords(0), dest_coords(2)) - vol_ref.at(x, y, z)) <= 1e-5);
                            }
                        }
                    } 

                }
                SECTION("Rotation of 180 degrees using tricubic interpolation returns same volume as python"){
                    float theta = 180.0;
                    Matrix3f R;
                    R = Utils::get_rotation_matrix(theta, wx, wy, wz);
                    Vector3f dest_coords;

                    std::stringstream sstheta;
                    sstheta << (int) theta;

                    string path_ref =
                      "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/linear_Vol1_" +
                      sstheta.str() +
                      sswx.str() +
                      sswy.str() + 
                      sswz.str() + ".dat";
                    std::vector<float> ref_vector = ReadFile<float>::read_reference_volume(path_ref, cubeVectorLength);
                    VolumeT vol_ref(ref_vector, cubeSize);

                    for(size_t z = 0; z < cubeSize; z++) {
                        for(size_t y = 0; y < cubeSize; y++) {
                            for(size_t x = 0; x < cubeSize; x++) {
                                dest_coords = R*(Vector3f(y,x,z)-cubeCenter)+cubeCenter;
                                REQUIRE(abs(cubic_interpolator.interp(dest_coords(1), dest_coords(0), dest_coords(2)) - vol_ref.at(x, y, z)) <= 1e-5);
                            }
                        }
                    } 

                }
            }
        }
    }
}

// TEST_CASE("rotation of 90 degrees using trilinear interpolation returns the same volume as python") {

//     typedef Interpolation_3D<float> interp;
//     typedef interp::vector1D vector1D;
//     typedef interp::pointer1D pointer1D;
//     typedef std::complex<float> dataT;

//     string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
//     int num_slice = 32;
//     int cubeVectorLength = num_slice*num_slice*num_slice;
//     vector1D vol1(cubeVectorLength); 
//     ReadFile<float>::read_volume(&vol1, path, num_slice);

//     float theta = 90.0;
//     for(int wx = 0; wx < 2; ++wx){
//         for(int wy = 0; wy < 2; ++wy){
//             for(int wz = 0; wz < 2; ++wz){
//                 vector1D dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);
//                 string path_ref = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/Vol1_" + std::to_string(int(theta)) + std::to_string(int(wx)) + std::to_string(int(wy)) + std::to_string(int(wz)) + ".dat";
//                 //printf(path_ref, int(theta),wx, wy, wz);
//                 vector1D vol_ref = ReadFile<float>::read_reference_volume(path_ref, cubeVectorLength);
//                 SECTION("the interpolated volume is the same as python output") {
//                     for(int i = 0; i < vol1.size(); ++i){
//                         REQUIRE(dest[i] == vol_ref[i]);
//                     }
//                 }
//             }
//         }
//     }
// }

