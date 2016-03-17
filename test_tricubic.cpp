#include "ReadFile.h"
#include "BinaryFile.h"
#include "TricubicInterpolator.h"
#include "Utils.h"

#include "CentralDifferenceDifferentiator.h"
#include "VolumeAtAddressable.h"

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
    string path = "test_data/8mm_iso_x_rot_0_5_to_2_5_deg_z_trans_rep_0";

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<float> dataVector(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector, path, cubeSize);

    typedef VolumeAtAddressable< std::vector<float>, float> VolumeT; 

    VolumeT volume(cubeSize, dataVector);

    CentralDifferencesDifferentiator<VolumeT> volDiffer(&volume);
    VolumeT dx(cubeSize, cubeVectorLength);
    volDiffer.xDerivative(&dx);

    VolumeT dy(cubeSize, cubeVectorLength);
    volDiffer.yDerivative(&dy);

    VolumeT dz(cubeSize, cubeVectorLength);
    volDiffer.zDerivative(&dz);

    CentralDifferencesDifferentiator<VolumeT> dxDiffer(&dx);
   
    VolumeT dxy(cubeSize, cubeVectorLength);
    dxDiffer.yDerivative(&dxy);
    
    VolumeT dxz(cubeSize, cubeVectorLength);
    dxDiffer.zDerivative(&dxz);

    CentralDifferencesDifferentiator<VolumeT> dyDiffer(&dy);

    VolumeT dyz(cubeSize, cubeVectorLength);
    dyDiffer.zDerivative(&dxz);

    CentralDifferencesDifferentiator<VolumeT> dxyDiffer(&dxy);
    
    VolumeT dxyz(cubeSize, cubeVectorLength);
    dxyDiffer.zDerivative(&dxz);
    
    TricubicInterpolator<VolumeT, float>
      cubic_interpolator(&volume, &dx, &dy, &dz, &dxy, &dxz, &dyz, &dxyz);

    int wx = 1;
    int wy = 0;
    int wz = 0;

    float theta = 90.0;
    Matrix3f R;
    R = Utils::get_rotation_matrix(theta, wx, wy, wz);
    Vector3f dest_coords;
    Vector3f cubeCenter;
    cubeCenter << volume.cubeCenter, volume.cubeCenter, volume.cubeCenter;

    while(true) {
    for(size_t z = 0; z < cubeSize; z++) {
        for(size_t y = 0; y < cubeSize; y++) {
            for(size_t x = 0; x < cubeSize; x++) {

                dest_coords.noalias() =
                  R * (Vector3f(y,x,z) - cubeCenter) + cubeCenter;

                float tmp = cubic_interpolator.interp(dest_coords(1), dest_coords(0), dest_coords(2));
            }
        }
    }
    }

}
