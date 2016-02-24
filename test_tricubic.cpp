#include "ReadFile.h"
#include "BinaryFile.h"
#include "TrilinearInterpolator.h"
#include "TricubicInterpolator.h"
#include "Utils.h"

#include <vector>
#include <complex>
#include <string>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <iostream>


using namespace std;

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

int main(){
    typedef float dataT;
    typedef Volume<dataT, std::vector<dataT> > VolumeT;

    // const size_t cubeSize = 10;
    // const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    // std::vector<dataT> initialData(cubeVectorLength);

    // for(size_t i = 0; i < cubeVectorLength; i++) {
    //     initialData[i] = i; 
    // }

    // typedef Volume<dataT, std::vector<dataT> > VolumeT;
    // VolumeT volume(initialData, cubeSize);


    // TricubicInterpolator<VolumeT, float> interpolator(&volume);

    // for(size_t z = 0; z < cubeSize; z++) {
    //     for(size_t y = 0; y < cubeSize; y++) {
    //         for(size_t x = 0; x < cubeSize; x++) {
    //             // cout << "index z, y, x: " << z << " " <<  y << " " << x<< endl;
    //             // float tmp = interpolator.interp(z, y, x);
    //             // cout << "interp result " << tmp << " volume(zyx): " <<  volume.at(z, y, x) << endl;

    //             if(interpolator.interp(z, y, x) != volume.at(z, y, x)){
    //                 cout << "index z, y, x: " << z << " " <<  y << " " << x<< endl;
    //             }
    //         }
    //     }
    // }

    // double wall0 = get_wall_time();
    // double cpu0  = get_cpu_time();

    // string path = "test_data/8mm_iso_x_rot_0_5_to_2_5_deg_z_trans_rep_0_slice_0.dat";
    // const size_t cubeSize = 32;
    // const size_t cubeVectorLength = cubeSize*cubeSize*cubeSize;

    // std::vector<dataT> initialData(cubeVectorLength);
    // ReadFile<float>::read_volume(&initialData, path, cubeSize);
    // VolumeT volume(initialData, cubeSize);

    // double wall1 = get_wall_time();
    // double cpu1  = get_cpu_time();
    // cout << "Wall Time ReadVolume = " << wall1 - wall0 << endl;
    // cout << "CPU Time ReadVolume  = " << cpu1  - cpu0  << endl;

    // wall0 = get_wall_time();
    // cpu0  = get_cpu_time();   

    // TricubicInterpolator<VolumeT, float> interpolator(&volume);

    // for(size_t z = 0; z < cubeSize; z++) {
    //     for(size_t y = 0; y < cubeSize; y++) {
    //         for(size_t x = 0; x < cubeSize; x++) {
    //             // cout << "index z, y, x: " << z << " " <<  y << " " << x<< endl;
    //             // float tmp = interpolator.interp(z, y, x);
    //             // cout << "interp result " << tmp << " volume(zyx): " <<  volume.at(z, y, x) << endl;

    //             if(interpolator.interp(z, y, x) != volume.at(z, y, x)){
    //                 cout << "index z, y, x: " << z << " " <<  y << " " << x<< endl;
    //             }
    //         }
    //     }
    // }

    // wall1 = get_wall_time();
    // cpu1  = get_cpu_time();
    // cout << "Wall Time Cubic Interplation = " << wall1 - wall0 << endl;
    // cout << "CPU Time Cubic Interplation = " << cpu1  - cpu0  << endl;


    string path = "test_data/8mm_iso_x_rot_0_5_to_2_5_deg_z_trans_rep_0";

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<float> dataVector(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector, path, cubeSize);

    typedef Volume<float, std::vector<float> > VolumeT;
    VolumeT volume(dataVector, cubeSize);

    TrilinearInterpolator<VolumeT, float> linear_interpolator(&volume);
    TricubicInterpolator<VolumeT, float> cubic_interpolator(&volume);

    int wx = 1;
    int wy = 0;
    int wz = 0;

    float theta = 90.0;
//
//    string path_ref = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/linear_Vol1_" + std::to_string(int(theta)) + std::to_string(int(wx)) + std::to_string(int(wy)) + std::to_string(int(wz)) + ".dat";
//    std::vector<float> ref_vector = ReadFile<float>::read_reference_volume(path_ref, cubeVectorLength);
//    VolumeT vol_ref(ref_vector, cubeSize);
    Matrix3f R;
    R = Utils::get_rotation_matrix(theta, wx, wy, wz);
    Vector3f dest_coords;

    while(true) {
    for(size_t z = 0; z < cubeSize; z++) {
        for(size_t y = 0; y < cubeSize; y++) {
            for(size_t x = 0; x < cubeSize; x++) {

                dest_coords.noalias() = R*(Vector3f(y,x,z)-volume.cubeCenter)+volume.cubeCenter;

                float tmp = cubic_interpolator.interp(dest_coords(1), dest_coords(0), dest_coords(2));
                //float tmp = linear_interpolator.interp(z,y,x);
                //if(abs(tmp - vol_ref.at(x, y, z)) > 1e-5){
                //    cout << "idx(z,y,x)" << z << "," << y << "," << x << " values: " << tmp << ", " << vol_ref.at(x, y, z) << endl;
                //}
            }
        }
    }
    }

}
