#include "catch.hpp"
#include "ReadFile.h"
#include "BinaryFile.h"
#include "Interpolation_3D.h"

#include <math.h>
#include <vector>
#include <complex>

/* Test cases for Trilinear Interpolation */
TEST_CASE("rotation of zero degrees using trilinear interpolation returns itself") {

    typedef Interpolation_3D<float> interp;
    typedef interp::vector1D vector1D;
    typedef interp::pointer1D pointer1D;
    typedef std::complex<float> dataT;

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;
    vector1D vol1(num_slice*num_slice*num_slice); 
    ReadFile<float>::read_volume(&vol1, path, num_slice);

    float theta = 0.0;
    for(int wx = 0; wx < 2; ++wx){
        for(int wy = 0; wy < 2; ++wy){
            for(int wz = 0; wz < 2; ++wz){
                vector1D dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);
                SECTION("the rotation with trilinear is identical with itself") {
                    for(int i = 0; i < vol1.size(); ++i){
                        REQUIRE(dest[i] == vol1[i]);
                    }
                }
            }
        }
    }
}

TEST_CASE("rotation of zero degrees returns original volume for tricubic interpolation") {

    typedef Interpolation_3D<float> interp;
    typedef interp::vector1D vector1D;
    typedef interp::pointer1D pointer1D;
    typedef std::complex<float> dataT;

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;
    vector1D vol1(num_slice*num_slice*num_slice); 
    ReadFile<float>::read_volume(&vol1, path, num_slice);

    // pre compute derivatives
    int derivs_shape = num_slice + 30;
    pointer1D derivs(derivs_shape*derivs_shape*derivs_shape);
    interp::tricubic_derivatives(&vol1, &derivs);

    float theta = 0.0;
    for(int wx = 0; wx < 2; ++wx){
        for(int wy = 0; wy < 2; ++wy){
            for(int wz = 0; wz < 2; ++wz){
                vector1D dest_cub = interp::rotate_volume_tricubic(&derivs,theta,wx,wy,wz);
                SECTION("the rotation using tricubic is identical") {
                    for(int i = 0; i < vol1.size(); ++i){
                        REQUIRE(dest_cub[i] == vol1[i]);
                    }
                }
            }
        }
    }
}

TEST_CASE("Testing Rotation of 90 and 180 degrees"){
    typedef Interpolation_3D<float> interp;
    typedef interp::vector1D vector1D;
    typedef interp::pointer1D pointer1D;
    typedef std::complex<float> dataT;

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;
    int vector_length = num_slice*num_slice*num_slice;
    vector1D vol1(num_slice*num_slice*num_slice); 
    ReadFile<float>::read_volume(&vol1, path, num_slice); 

    // pre compute derivatives
    int derivs_shape = num_slice + 30;
    pointer1D derivs(derivs_shape*derivs_shape*derivs_shape);
    interp::tricubic_derivatives(&vol1, &derivs);
    
    // int wx = 1;
    // int wy = 0;
    // int wz = 0;
    for(int wx = 0; wx < 2; ++wx){
        for(int wy = 0; wy < 2; ++wy){
            for(int wz = 0; wz < 2; ++wz){
                SECTION("Rotation of 90 degrees using trilinear interpolation returns same volume as python"){
                    float theta = 90.0;
                    vector1D dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);  
                    string path_ref = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/linear_Vol1_" + std::to_string(int(theta)) + std::to_string(int(wx)) + std::to_string(int(wy)) + std::to_string(int(wz)) + ".dat";
                    vector1D vol1_ref = ReadFile<float>::read_reference_volume(path_ref, vector_length);
                    for(int i = 0; i < vol1.size(); ++i){
                        REQUIRE(abs(dest[i] - vol1_ref[i]) <= 1e-5);
                    }     
                }

                SECTION("Rotation of 180 degrees using trilinear interpolation returns same volume as python"){
                    float theta = 180.0;
                    vector1D dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);  
                    string path_ref = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/linear_Vol1_" + std::to_string(int(theta)) + std::to_string(int(wx)) + std::to_string(int(wy)) + std::to_string(int(wz)) + ".dat";
                    vector1D vol1_ref = ReadFile<float>::read_reference_volume(path_ref, vector_length);
                    for(int i = 0; i < vol1.size(); ++i){
                        REQUIRE(abs(dest[i] - vol1_ref[i]) <= 1e-5);
                    }  
                }

                SECTION("Rotation of 90 degrees using tricubic interpolation returns same volume as python"){
                    float theta = 90.0;
                    vector1D dest_cub = interp::rotate_volume_tricubic(&derivs,theta,wx,wy,wz);   
                    string path_ref = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/cubic_Vol1_" + std::to_string(int(theta)) + std::to_string(int(wx)) + std::to_string(int(wy)) + std::to_string(int(wz)) + ".dat";
                    vector1D vol1_ref = ReadFile<float>::read_reference_volume(path_ref, vector_length);
                    for(int i = 0; i < vol1.size(); ++i){
                        REQUIRE(abs(dest_cub[i] - vol1_ref[i]) <= 1e-5);
                    }     
                }
                SECTION("Rotation of 180 degrees using tricubic interpolation returns same volume as python"){
                    float theta = 180.0;
                    vector1D dest_cub = interp::rotate_volume_tricubic(&derivs,theta,wx,wy,wz);   
                    string path_ref = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/cubic_Vol1_" + std::to_string(int(theta)) + std::to_string(int(wx)) + std::to_string(int(wy)) + std::to_string(int(wz)) + ".dat";
                    vector1D vol1_ref = ReadFile<float>::read_reference_volume(path_ref, vector_length);
                    for(int i = 0; i < vol1.size(); ++i){
                        REQUIRE(abs(dest_cub[i] - vol1_ref[i]) <= 1e-5);
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
//     int vector_length = num_slice*num_slice*num_slice;
//     vector1D vol1(vector_length); 
//     ReadFile<float>::read_volume(&vol1, path, num_slice);

//     float theta = 90.0;
//     for(int wx = 0; wx < 2; ++wx){
//         for(int wy = 0; wy < 2; ++wy){
//             for(int wz = 0; wz < 2; ++wz){
//                 vector1D dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);
//                 string path_ref = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/Vol1_" + std::to_string(int(theta)) + std::to_string(int(wx)) + std::to_string(int(wy)) + std::to_string(int(wz)) + ".dat";
//                 //printf(path_ref, int(theta),wx, wy, wz);
//                 vector1D vol1_ref = ReadFile<float>::read_reference_volume(path_ref, vector_length);
//                 SECTION("the interpolated volume is the same as python output") {
//                     for(int i = 0; i < vol1.size(); ++i){
//                         REQUIRE(dest[i] == vol1_ref[i]);
//                     }
//                 }
//             }
//         }
//     }
// }

