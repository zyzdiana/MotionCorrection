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


int main(void) {
  int N = 32;
  int N_out =  N/2 + 1;
  float in[N];
  fftwf_complex out[N_out];
  float in_back[N];
  fftwf_plan p, q;
  int i;

  /* prepare a cosine wave */
  for (i = 0; i < N; i++) {
    in[i] = cos(M_PI*i/N);
  }
  for (i = 0; i < N; i++)
    printf("input: %3d %+9.5f \n", i, in[i]);

  /* forward Fourier transform, save the result in 'out' */
  p = fftwf_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

  fftwf_execute(p);
  for (i = 0; i < N_out; i++)
    printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][0], out[i][1]);
  fftwf_destroy_plan(p);

  /* backward Fourier transform, save the result in 'in2' */
  printf("\nInverse transform:\n");
  q = fftwf_plan_dft_c2r_1d(N, out, in_back,FFTW_ESTIMATE);
  fftwf_execute(q);
  /* normalize */
  for (i = 0; i < N; i++) {
    in_back[i] *= 1./N;
  }
  for (i = 0; i < N; i++)
    printf("recover: %3d %+9.5f vs. %+9.5f \n",
        i, in[i], in_back[i]);
  fftwf_destroy_plan(q);

  fftwf_cleanup();
  return 0;
}
