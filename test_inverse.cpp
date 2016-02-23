#include "ReadFile.h"
#include "BinaryFile.h"
#include "Interpolation_3D.h"
#include "Gauss_Newton.h"

#include <vector>
#include <complex>
#include <string>
#include <time.h>
#include <math.h>
#include <sys/time.h>

using namespace std;

typedef Interpolation_3D<float> interp;
typedef interp::vector1D vector1D;
typedef interp::pointer1D pointer1D;
typedef complex<float> dataT;
typedef Gauss_Newton<float> GN;

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
    cout << endl;
    /************** Test matrix inverse  ****************/

    // float rotMatrix[3][3] = {{1,2,1},{3,2,1},{1,2,3}};
    // float inverse[3][3];
    // GN::inverse3(rotMatrix,inverse);
    // for (int i = 0; i < 3; ++i){
    //     for (int j = 0; j < 3; ++j){
    //         cout << inverse[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // float rotMatrix1[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    // GN::inverse3(rotMatrix1,inverse);
    // for (int i = 0; i < 3; ++i){
    //     for (int j = 0; j < 3; ++j){
    //         cout << inverse[i][j] << " ";
    //     }
    //     cout << endl;
    // }



    /************** Test rotate coords  ****************/
    // float result[3];
    // float params[6] = {2,-1,1,0,0,0};
    // int shape = 32;
    // GN::rotate_coords_transformation(result, 1, 1, 1, params,15.5,15.5,15.5, 16, shape);
    // for (int i = 0; i < 3; ++i){
    //     cout << result[i] << " ";
    // }
    // cout << endl;

    // GN::rotate_coords_transformation(result, 1.5, 1.5, 1.5, params,15.5,15.5,15.5, 16, shape);
    // for (int i = 0; i < 3; ++i){
    //     cout << result[i] << " ";
    // }
    // cout << endl;
    // cout << endl;


    /************** Read in a volume and compute derivatives  ****************/
    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    string path0 = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_1";
    int num_slice = 32;


    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

    vector1D vol1(32*32*32); 
    vector1D vol0(32*32*32);

    ReadFile<float>::read_volume(&vol1, path, num_slice);
    ReadFile<float>::read_volume(&vol0, path0, num_slice);

    int derivs_shape = num_slice + 30;
    pointer1D derivs(derivs_shape*derivs_shape*derivs_shape);
    interp::tricubic_derivatives(&vol1, &derivs);

    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    cout << "Wall Time ReadFile and compute derivs = " << wall1 - wall0 << endl;
    cout << "CPU Time ReadFile and compute derivs = " << cpu1  - cpu0  << endl;


    wall0 = get_wall_time();
    cpu0  = get_cpu_time();   
    pointer1D gradP(num_slice*num_slice*num_slice);
    pointer1D gradP_T(6);
    // for (int i = 0; i < 32*32; ++i){
    //     for (int j = 0; j < 3; ++j){
    //         cout << gradP[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    //cout << gradP[0] << " " << endl;

    GN::get_gradient_P(&gradP, &derivs, 1.0 , false);
    GN::transpose(&gradP_T, &gradP);
    wall1 = get_wall_time();
    cpu1  = get_cpu_time();
    cout << "Wall Time compute gradP = " << wall1 - wall0 << endl;
    cout << "CPU Time compute gradP = " << cpu1  - cpu0  << endl;

    // for (int j = 0; j < 32*32; ++j){
    //     for (int i = 0; i < 6; ++i){
    //         cout << j << " " << gradP[j][i] << " ";
    //     }
    //     cout << endl;
    // }

    // pointer1D gradP(num_slice*num_slice*num_slice);
    // int counter = 1;
    // while(counter != 0){
    //     GN::get_gradient_P(&gradP, &derivs, 1.0 , false);
    // }

    // pointer1D gradP_T(6);
    // GN::transpose(&gradP_T,&gradP);
    // for (int j =0; j < gradP_T.size(); ++j){
    //     for (int i = 0; i < gradP.size(); ++i){
    //         if(gradP_T[j][i] != gradP[i][j]){
    //             cout << i << " " << j << endl;
    //         }
    //     }
    // }
    /************** test mask  ****************/
    // vector1D mask_weights(32*32*32); 
    // //int num_slice = 32;
    // int vector_length = num_slice*num_slice*num_slice;

    // GN::get_mast_weights(&mask_weights, 16, 0.4);

    // string path_ref = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/mask_weight.dat";
    // vector1D mask_ref = ReadFile<float>::read_reference_volume(path_ref, vector_length);
    // for (int i = 0; i < vector_length; ++i){
    //     if(abs(mask_ref[i] - mask_weights[i]) > 10e-6)
    //         cout << i << " " << mask_ref[i] - mask_weights[i]<< endl;
    // }
    // cout << mask_ref[1391] << " " << mask_weights[1391] << endl;

    // for (vector1D::iterator it = mask_weights.begin(); it != mask_weights.end(); ++it){
    //     cout  << *it << " ";
    // }

    /************** test mask  ****************/
    float P_initial[6] = {0,0,0,0,0,0};
    float *P;
    float divide_factor = 1.0;
    float alpha = 1.0;n
    float decrease_factor = 0.25;
    int max_iter = 10;

    wall0 = get_wall_time();
    cpu0  = get_cpu_time();  

    P = GN::gauss_newton(&vol0, &vol1, &derivs, P_initial, divide_factor, alpha, decrease_factor, max_iter, false);

    wall1 = get_wall_time();
    cpu1  = get_cpu_time();
    cout << "Wall Time Gauss Newton = " << wall1 - wall0 << endl;
    cout << "CPU Time Gauss Newton  = " << cpu1  - cpu0  << endl;

    cout << endl;
    cout << "Parameters at min error:  ";
    for (int i = 0; i < 6; ++i) cout << P[i] << " ";
    cout << endl;
    cout << "Translations (mm): " ; 
    for (int i = 0; i < 3; ++i) cout << P[i]*8.0 << " ";
    cout << endl;
    cout << "Rotations (degree): " ; 
    for (int i = 3; i < 6; ++i) cout << P[i]*180.0/PI << " ";
    cout << endl;

    cout << endl;
    return 0;
}
