#include "ReadFile.h"
#include "BinaryFile.h"
#include "Interpolation_3D.h"

#include <vector>
#include <complex>
#include <string>

using namespace std;

typedef Interpolation_3D<float> interp;
typedef interp::vector1D vector1D;
typedef complex<float> dataT;

int main(){
    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;
    vector1D vol1(num_slice*num_slice*num_slice); 
    ReadFile<float>::read_volume(&vol1, path, num_slice);

    /* print volume */
    // for (vector1D::iterator it_x = vol1.begin(); it_x != vol1.end(); ++it_x){
    //     cout << " " << *it_x; 
    // }

    string path1 = "/Users/zyzdiana/Dropbox/THESIS/C++_Test_Code/testVol11111.dat";
    vector1D vol1_ref = ReadFile<float>::read_reference_volume(path1, 32*32*32);  
    for (int i = 0; i < vol1_ref.size(); ++i){
        // if (vol1[i] != vol1_ref[i]){
        //     cout << i << endl;
        // }
        cout << vol1_ref[i]<< " ";
    }
    cout << endl;
    return 0;
}
