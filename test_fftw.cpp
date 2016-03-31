#include <stdlib.h>
#include <fftw3.h>

#include "ReadFile.h"
#include "BinaryFile.h"
#include "TrilinearInterpolator.h"
#include "TricubicInterpolator.h"
#include "Utils.h"
#include "Gauss_Newton.h"
#include "FFTOp.h"

#include <vector>
#include <complex>
#include <string>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <iostream>

using namespace std;
using namespace Eigen;

int main(void) {

    typedef float dataT;
    typedef Volume<dataT, std::vector<dataT>, float > VolumeT; 
    typedef VolumeT::T T;
    typedef Matrix< T, 64, 1 >    Vector64T;
    typedef Matrix< float, 6, 1 > Vector6f;
    typedef Matrix< float, 3, 1 > coordT;

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";

    const size_t cubeSize = 32;
    const float radius = cubeSize/2;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    // Moving volume
    std::vector<float> dataVector(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector, path, cubeSize);

    VolumeT volume(dataVector, cubeSize);
    coordT cubeCenter = volume.cubeCenterAsTriple<coordT>();

    FFTOp fft = FFTOp(cubeSize);
    std::vector<float> transformed(cubeVectorLength);
    fft.transform(&dataVector, &transformed);


  return 0;
}
