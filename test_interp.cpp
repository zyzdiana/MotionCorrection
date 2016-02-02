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

    vector1D vol1(32*32*32); 

    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

    ReadFile<float>::read_volume(&vol1, path, num_slice);

    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    cout << "Wall Time ReadFile = " << wall1 - wall0 << endl;
    cout << "CPU Time ReadFile  = " << cpu1  - cpu0  << endl;

    // for (vector1D::iterator it_x = vol1.begin(); it_x != vol1.end(); ++it_x){
    //     cout << " " << *it_x; 
    // }


    float theta = 0.0;
    float wx = 0;
    float wy = 1;
    float wz = 0;

    /* Test Trilinear Interpolation   */
    wall0 = get_wall_time();
    cpu0  = get_cpu_time();

    vector1D dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);

    wall1 = get_wall_time();
    cpu1  = get_cpu_time();
    cout << "Wall Time Linear = " << wall1 - wall0 << endl;
    cout << "CPU Time Linear  = " << cpu1  - cpu0  << endl;

    for (int i = 0; i < 32*32*32; ++i){
        if(dest[i] != vol1[i])
            cout << dest[i] << " ";
    }

    // while(wx!=0){
    //     dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);
    // } 
    // int derivs_shape = num_slice+30;
    // pointer1D derivs(derivs_shape*derivs_shape*derivs_shape);
    // for (auto it_x = derivs.begin(); it_x != derivs.end(); ++it_x){
    //     cout << " " << *(*it_x); 
    // }

    /*  Test Tricubic Derivatives   */
    // interp::tricubic_derivatives(&vol1, &derivs);
    // int i=0;
    // int j=0;
    // int k=0;
    // float* p = derivs[derivs_shape*derivs_shape*(j+15) + derivs_shape*(i+15) + (k+15)];
    // cout << *(p) << endl;
    // cout << *(p+1) << endl;
    // cout << *(p+2) << endl;


    /*  Test Tricubic Interpolation   */
    // wall0 = get_wall_time();
    // cpu0  = get_cpu_time();   

    // vector1D dest_cub = interp::rotate_volume_tricubic(&vol1,theta,wx,wy,wz);

    // wall1 = get_wall_time();
    // cpu1  = get_cpu_time();
    // cout << "Wall Time Cubic = " << wall1 - wall0 << endl;
    // cout << "CPU Time Cubic  = " << cpu1  - cpu0  << endl;


    // for (int i = 0; i < 32*32*32; ++i){
    //     if(dest_cub[i] != vol1[i])
    //         cout << i << " ";
    // }
    return 0;
}
