#include "ReadFile.h"
#include "BinaryFile.h"
#include "Interpolation_3D.h"

#include <vector>
#include <complex>
#include <string>
#include <time.h>
#include <sys/time.h>

using namespace std;

typedef Interpolation_3D<float> interp;
typedef interp::vector1D vector1D;
typedef interp::pointer1D pointer1D;
typedef complex<float> dataT;

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
    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;
    int vector_length = num_slice*num_slice*num_slice;

    vector1D vol1(32*32*32); 

    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

    ReadFile<float>::read_volume(&vol1, path, num_slice);

    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    cout << "Wall Time ReadFile = " << wall1 - wall0 << endl;
    cout << "CPU Time ReadFile  = " << cpu1  - cpu0  << endl;

    float theta = 90.0;
    float wx = 1;
    float wy = 0;
    float wz = 0;

    // /* Test Trilinear Interpolation   */
    // wall0 = get_wall_time();
    // cpu0  = get_cpu_time();

    // vector1D dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);

    // wall1 = get_wall_time();
    // cpu1  = get_cpu_time();
    // cout << "Wall Time Linear = " << wall1 - wall0 << endl;
    // cout << "CPU Time Linear  = " << cpu1  - cpu0  << endl;

    // string path_ref_linear = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/linear_Vol1_" + std::to_string(int(theta)) + std::to_string(int(wx)) + std::to_string(int(wy)) + std::to_string(int(wz)) + ".dat";
    // vector1D vol1_ref_linear = ReadFile<float>::read_reference_volume(path_ref_linear, vector_length);
    // for (int i = 0; i < 32*32*32; ++i){
    //     if(abs(dest[i] - vol1_ref_linear[i]) >= 1e-5)
    //         cout << i << " " << dest[i] << " " << vol1_ref_linear[i] << endl;
    // }

    /*  Test Tricubic Interpolation   */
    wall0 = get_wall_time();
    cpu0  = get_cpu_time();   

    int derivs_shape = num_slice + 30;
    pointer1D derivs(derivs_shape*derivs_shape*derivs_shape);
    interp::tricubic_derivatives(&vol1, &derivs);

    wall1 = get_wall_time();
    cpu1  = get_cpu_time();
    cout << "Wall Time Cubic Derivatives = " << wall1 - wall0 << endl;
    cout << "CPU Time Cubic Derivatives = " << cpu1  - cpu0  << endl;


    wall0 = get_wall_time();
    cpu0  = get_cpu_time();   

    vector1D dest_cub = interp::rotate_volume_tricubic(&derivs,theta,wx,wy,wz);

    wall1 = get_wall_time();
    cpu1  = get_cpu_time();
    cout << "Wall Time Cubic Interpolation = " << wall1 - wall0 << endl;
    cout << "CPU Time Cubic Interpolation = " << cpu1  - cpu0  << endl;

    string path_ref = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/cubic_Vol1_" + std::to_string(int(theta)) + std::to_string(int(wx)) + std::to_string(int(wy)) + std::to_string(int(wz)) + ".dat";
    cout << path_ref << endl;
    vector1D vol1_ref = ReadFile<float>::read_reference_volume(path_ref, vector_length);
    for (int i = 0; i < 32*32*32; ++i){
        if(abs(dest_cub[i] - vol1_ref[i]) >= 1e-7)
            cout << i << " " << dest_cub[i]- vol1_ref[i] << endl;
            //cout << i << " " << abs(dest_cub[i] - vol1_ref[i]) <<endl;
    }

    // for (int i = 0; i < 32; ++i){
    //         cout << i << " " << vol1_ref[i] << endl;
    //         //cout << i << " " << abs(dest_cub[i] - vol1_ref[i]) <<endl;
    // }

    // cout << "size of trilinear interp result: " << sizeof(dest_cub[0]) << endl;
    // cout << "size of loaded interp result: " << sizeof(vol1_ref[0]) << endl;
    cout << endl;

    // cout << endl;

    // vector1D Y(64);
    // interp::get_target_Y(&Y, 0, 0, 27);

    // for (int i = 0; i < 64; ++i){
    //     cout << (*(derivs[27+15]+i)) << " " << Y[i] << " ";
    // }
    //cout << derivs[27+15] << " " << derivs[28+15] << " " << derivs[29+15] << " " << derivs[30+15] << " " << derivs[31+15] << endl;
    return 0;

}
