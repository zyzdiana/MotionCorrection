#include <stdlib.h>
#include <fftw3.h>

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

    static float window(float n, float radius, float d = 0.4){
        float tmp = n/radius - 0.75;
        if (tmp < 0){
            return 1.0;
        }
        else{
            if((tmp/d > -0.5) and (tmp/d < 0.5)){
                return cos(M_PI*(tmp/d));
            }
            else{
                return 0.0;
            }
        }
    }
int main(void) {

    typedef float dataT;
    typedef Volume<dataT, std::vector<dataT>, float > VolumeT; 
    typedef VolumeT::T T;
    typedef Matrix< T, 64, 1 >    Vector64T;
    typedef Matrix< float, 6, 1 > Vector6f;
    typedef Matrix< float, 3, 1 > coordT;

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    string path1 = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_1";

    const size_t cubeSize = 32;
    const float radius = cubeSize/2;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;
    // Moving volume
    std::vector<float> dataVector(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector, path, cubeSize);

    VolumeT volume(dataVector, cubeSize);
    coordT cubeCenter = volume.cubeCenterAsTriple<coordT>();

    //Target volume
    std::vector<float> dataVector1(cubeVectorLength);
    ReadFile<float>::read_volume(&dataVector1, path1, cubeSize);

    VolumeT volume1(dataVector1, cubeSize);



  float input[cubeSize][cubeSize][cubeSize], input_back[cubeSize][cubeSize][cubeSize];
  for(int z = 0; z < cubeSize; ++z){
      for(int y = 0; y < cubeSize; ++y){
          for(int x = 0; x < cubeSize; ++x){
            input[z][y][x] = volume.at(z,y,z);
          }
      }
  }
  fftwf_complex out[cubeSize][cubeSize][cubeSize];
  fftwf_plan p, q;

  /* forward Fourier transform, save the result in 'out' */
  p = fftwf_plan_dft_r2c_3d(cubeSize,cubeSize,cubeSize, input, out, FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);

  // int idx = 0;
  // for(int z = 0; z < cubeSize; ++z){
  //     for(int y = 0; y < cubeSize; ++y){
  //         for(int x = 0; x < cubeSize; ++x){
  //           T n = (coordT(z,y,x)-cubeCenter).norm();
  //           out[idx][0] *= window(n, radius);
  //           out[idx][1] *= window(n, radius);
  //           idx += 1;
  //         }
  //     }
  // }
  //cout << "idx: " << idx << " length: " << cubeVectorLength/2+1 << endl;
  // for (int i = 0; i < cubeVectorLength/2+1; i++)
  //   cout <<"freq: " << i << " " << out[i][0] << " " << out[i][1]<<endl;

  // /* backward Fourier transform, save the result in 'in2' */
  // q = fftwf_plan_dft_c2r_1d(cubeVectorLength, out, &transformed[0], FFTW_ESTIMATE);
  // fftwf_execute(q);
  // /* normalize */
  // for (int i = 0; i < cubeVectorLength; i++) {
  //   transformed[i] *= 1./cubeVectorLength;
  // }

  // // for (int i = 0; i < cubeVectorLength; i++)
  // //   printf("recover: %3d %+9.5f vs. %+9.5f \n",
  // //       i, dataVector[i], transformed[i]);

  // VolumeT transformed_vol(cubeSize);
  // fftwf_destroy_plan(q);

  // fftwf_cleanup();


  return 0;
}
