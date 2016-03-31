#ifndef FFTOp_h
#define FFTOp_h

#include "Interpolation_3D.h"
#include "TricubicInterpolator.h"
#include "CentralDifferenceDifferentiator.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>

using namespace std;
using namespace Eigen;

class FFTOp{
  public:
    //const size_t cubeSize;
    size_t cubeVectorLength;

    FFTOp(const size_t cubeSize){
      cubeVectorLength = cubeSize*cubeSize*cubeSize;
	  size_t N_out =  cubeVectorLength/2 + 1;
	  dataVector.resize(cubeVectorLength);

	  p_forward = fftwf_plan_dft_r2c_1d(cubeVectorLength, &dataVector[0], output, FFTW_ESTIMATE);
	  p_backward = fftwf_plan_dft_c2r_1d(cubeVectorLength, output, &dataVector[0], FFTW_ESTIMATE);
    }

 	~FFTOp(){
 		fftwf_destroy_plan(p_forward);
 		fftwf_destroy_plan(p_backward);
 		fftwf_cleanup();
 	}

    void transform(std::vector<float>* volume, std::vector<float>* transformed){
		fftwf_execute_dft_r2c(p_forward, &(volume->at(0)), output);
		// int idx = 0;
		// for(int z = 0; z < cubeSize; ++z){
		//   for(int y = 0; y < cubeSize; ++y){
		//       for(int x = 0; x < cubeSize; ++x){
		//         T n = (coordT(z,y,x)-cubeCenter).norm();
		//         out[idx][0] *= window(n, radius);
		//         out[idx][1] *= window(n, radius);
		//         idx += 1;
		//       }
		//   }
		// }
		fftwf_execute_dft_c2r(p_backward, output, &(transformed->at(0)));
    }
  protected:
  	fftwf_plan p_forward;
  	fftwf_plan p_backward;
  	std::vector<float> dataVector;
  	fftwf_complex *output;
};
#endif