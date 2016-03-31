#include "ReadFile.h"
#include "BinaryFile.h"
#include "TrilinearInterpolator.h"
#include "TricubicInterpolator.h"
#include "Utils.h"
#include "Gauss_Newton.h"

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
    typedef VolumeT::T T;
    typedef Matrix< T, 64, 1 >     Vector64T;

    typedef Matrix< float, 6, 1 >     Vector6f;
    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    string path1 = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_1";

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;
    // Moving volume
    std::vector<float> dataVector(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector, path, cubeSize);

    typedef Volume<float, std::vector<float> > VolumeT;
    VolumeT volume(dataVector, cubeSize);

    //Target volume
    std::vector<float> dataVector1(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector1, path1, cubeSize);

    VolumeT volume1(dataVector1, cubeSize);

    TrilinearInterpolator<VolumeT, float> linear_interpolator(&volume);
    TricubicInterpolator<VolumeT, float> cubic_interpolator(&volume);


    Gauss_Newton<VolumeT, float> gn(&cubic_interpolator);
 
    while(true){
        gn.gauss_newton(&volume1, 10);
    }


}
