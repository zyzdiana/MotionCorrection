#include "FFTOp.h"

#include <complex>

template <>
typename FFTWPlanType<float>::PlanType
FFTOp<float>::createForwardPlan(size_t cubeSize) {
    spatialVolumeT spatial(cubeSize);
    fourierVolumeT fourier(cubeSize);

	  return fftwf_plan_dft_r2c_3d(
      cubeSize, cubeSize, cubeSize,
      spatial.buffer, (fftwf_complex*) fourier.buffer, FFTW_ESTIMATE);
}

template <>
typename FFTWPlanType<float>::PlanType
FFTOp<float>::createBackwardPlan(size_t cubeSize) {
    fourierVolumeT fourier(cubeSize);
    spatialVolumeT spatial(cubeSize);

	  return fftwf_plan_dft_c2r_3d(
      cubeSize, cubeSize, cubeSize,
      (fftwf_complex*) fourier.buffer, spatial.buffer, FFTW_ESTIMATE);
}

template <>
FFTOp<float>::FFTOp(size_t cubeSize):
  cubeSize(cubeSize),
  forwardPlan(createForwardPlan(cubeSize)),
  backwardPlan(createBackwardPlan(cubeSize)) {
}

template <>
FFTOp<float>::~FFTOp() {
 		  fftwf_destroy_plan(forwardPlan);
 		  fftwf_destroy_plan(backwardPlan);
}


template <>
void FFTOp<float>::forward (
      spatialVolumeT *spatial,
      fourierVolumeT *fourier) const {
 		  fftwf_execute_dft_r2c(
        forwardPlan, spatial->buffer, (fftwf_complex*) fourier->buffer);
}

template <>
void FFTOp<float>::backward (
      fourierVolumeT *fourier,
      spatialVolumeT *spatial) const {
 		  fftwf_execute_dft_c2r(
        backwardPlan, (fftwf_complex*) fourier->buffer, spatial->buffer);
}
