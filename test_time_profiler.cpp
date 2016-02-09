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

int main(){
    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;

    vector1D vol1(32*32*32); 

    ReadFile<float>::read_volume(&vol1, path, num_slice);

    float theta = 0.0;
    float wx = 1;
    float wy = 0;
    float wz = 0;

    int derivs_shape = num_slice + 30;
    pointer1D derivs(derivs_shape*derivs_shape*derivs_shape);
    interp::tricubic_derivatives(&vol1, &derivs);

    vector1D dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);

    // while(wx!=0){
    //     dest = interp::rotate_volume_trilinear(&vol1,theta,wx,wy,wz);
    // }
    
    vector1D dest_cub;
    for (int i=0; i < 10; ++i){
        vector1D dest_cub = interp::rotate_volume_tricubic(&derivs,theta,wx,wy,wz);
    }

    return 0;
}
