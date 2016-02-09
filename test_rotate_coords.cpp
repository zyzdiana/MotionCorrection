#include "ReadFile.h"
#include "BinaryFile.h"
#include "Interpolation_3D.h"

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

int main(){
    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;

    vector1D vol1(32*32*32); 

    ReadFile<float>::read_volume(&vol1, path, num_slice);

    float theta = 0.0;
    float wx = 0;
    float wy = 1;
    float wz = 0;

    vector1D dest_coords(3);
    float i = 0;
    float j = 0;
    float k = 0;

    float ox = 15.5;
    float oy = 15.5;
    float oz = 15.5;
    // for (i = 0; i < 32; ++ i){
    //     interp::rotate_coords_3d(&dest_coords, i, j, k, theta, wx, wy, wz, ox, oy, oz); 
    //     cout << dest_coords[0] << " " << dest_coords[1] << " " << dest_coords[2] << " " << endl;     
    // }

    for (i = 0; i < 32; ++ i){
        for (j = 0; j < 32; ++ j){
            for (k = 0; k < 32; ++ k){
                interp::rotate_coords_3d(&dest_coords, i, j, k, theta, wx, wy, wz, ox, oy, oz); 
                if (dest_coords[0] != i){
                    cout << "i " << i << endl; 
                }
                if (dest_coords[1] != j){
                    cout << "j " << j << endl; 
                }
                if (dest_coords[2] != k){
                    cout << "k " << k << endl; 
                }
                //cout << dest_coords[0] << " " << dest_coords[1] << " " << dest_coords[2] << " " << endl;  
            }
        }   
    }

    //interp::rotate_coords_3d(&dest_coords, i, j, k, theta, wx, wy, wz, ox, oy, oz);
    
    //cout << dest_coords[0] << " " << dest_coords[1] << " " << dest_coords[2] << " " << endl;
    //interp::rotate_coords_3d(vector1D *result, float x, float y, float z, float theta, float wx, float wy, float wz, float ox, float oy, float oz)
    return 0;
}
