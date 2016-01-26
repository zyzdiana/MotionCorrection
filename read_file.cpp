#include "ReadFile.h"
#include "BinaryFile.h"
#include "Interpolation_3D.h"

#include <vector>
#include <complex>
#include <string>

using namespace std;

typedef Interpolation_3D<float> interp;
typedef interp::vector4D vector4D;
typedef interp::vector3D vector3D;
typedef interp::vector2D vector2D;
typedef interp::vector1D vector1D;
typedef complex<float> dataT;

int main(){
    //string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    //vector3D vol1 = ReadFile<float>::read_volume(path, 32);

    // for (vector3D::iterator it_x = vol1.begin(); it_x != vol1.end(); ++it_x){
    //     for (vector2D::iterator it_y = it_x->begin(); it_y != it_x->end(); ++it_y){
    //         for (vector1D::iterator it_z = it_y->begin(); it_z != it_y->end(); ++it_z){
    //             cout << " " << *it_z; 
    //         }
    //     }
    // }

    //float x = , float y, float z, float theta, float wx, float wy, float wz, float ox, float oy, float oz
    // vector1D rotated_coords = interp::rotate_coords_3d(0.0,0.0,0.0,0.0,0.0,0.0,1.0,15.5,15.5,15.5);
    // for (vector1D::iterator it = rotated_coords.begin(); it != rotated_coords.end(); ++it){
    //     cout << " " << *it; 
    // }
    // cout << " " <<  endl;

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;
    vector3D vol1 = ReadFile<float>::read_volume(path, num_slice);

    string ref_path = "/Users/zyzdiana/GitHub/AC297r-Volume-Registration/Thesis/rot90_y.dat";
    vector3D ref = ReadFile<float>::read_reference_volume(ref_path, num_slice);
    
    float theta = 90.0;
    float wx = 1;
    float wy = 0;
    float wz = 0;

    vector3D dest = interp::rotate_volume_trilinear(vol1,theta,wx,wy,wz);
    //vector3D dest = interp::rotate_volume_tricubic(vol1,theta,wx,wy,wz);
    
    for(int i = 0; i < num_slice; ++i){
        for(int j = 0; j < num_slice; ++j){
            for(int k = 0; k < num_slice; ++k){
                cout << dest[i][j][k] << " ";
                // if(dest[i][j][k] != ref[j][i][k]){
                //     cout <<i<<","<<j<<","<<k<<endl;
                // }
            }
        }
    }
    //double interp = interp::trilinear_interp(vol1,0,31,0);
    //cout << interp << endl;
    // float ox,oy,oz = float(num_slice-1)/2;
    // vector1D point(3,15);
    // vector1D rotated_point = interp::rotate_coords_3d(point[0],point[1],point[2], theta, wx, wy, wz, ox, oy, oz);
    // for (vector1D::iterator it = rotated_point.begin(); it != rotated_point.end(); ++it){
    //     cout << " " << *it; 
    // }
    // cout << " " <<  endl;
    return 0;
}
