#include "FFTWBuffer.h"
#include <complex>

#include "fftw3.h"

template <>
FFTWBuffer<float>::FFTWBuffer(size_t numElements) :
  numElements(numElements),
  buffer((float*) fftwf_malloc(numElements * sizeof(float))) {
}

template <>
FFTWBuffer<float>::~FFTWBuffer() {
  fftwf_free(buffer);
}

template <>
FFTWBuffer< std::complex<float> >::FFTWBuffer(size_t numElements) :
  numElements(numElements),
  buffer((std::complex<float>*) fftwf_malloc(numElements * sizeof(std::complex<float>))) {
}

template <>
FFTWBuffer< std::complex<float> >::~FFTWBuffer() {
  fftwf_free(buffer);
}


template <>
FFTWBuffer<double>::FFTWBuffer(size_t numElements) :
  numElements(numElements),
  buffer((double*) fftw_malloc(numElements * sizeof(double))) {
}

template <>
FFTWBuffer<double>::~FFTWBuffer() {
  fftw_free(buffer);
}

template <>
FFTWBuffer< std::complex<double> >::FFTWBuffer(size_t numElements) :
  numElements(numElements),
  buffer((std::complex<double>*) fftw_malloc(numElements * sizeof(std::complex<double>))) {
}

template <>
FFTWBuffer< std::complex<double> >::~FFTWBuffer() {
  fftw_free(buffer);
}
