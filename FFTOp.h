#ifndef FFTOp_h
#define FFTOp_h

#include "VolumeAtAddressable.h"
#include "SymmetricHalfVolumeAtAddressable.h"
#include "FFTWBuffer.h"

#include "fftw3.h"

#include <complex>

template <typename T>
class FFTWPlanType {
};

template <>
class FFTWPlanType<float> {
  public:
    typedef fftwf_plan PlanType;
};

template <>
class FFTWPlanType< std::complex<float> > {
  public:
    typedef fftwf_plan PlanType;
};

template <>
class FFTWPlanType<double> {
  public:
    typedef fftw_plan PlanType;
};

template <>
class FFTWPlanType< std::complex<double> > {
  public:
    typedef fftw_plan PlanType;
};

template <typename T>
class FFTFourierType {
};

template <>
class FFTFourierType<float> {
  public:
    typedef std::complex<float> FourierType;
    typedef
      SymmetricHalfVolumeAtAddressable< FFTWBuffer< FourierType > >
      FourierVolumeType;
};


template <typename _spatialT>
class FFTOp {
  public:
    typedef _spatialT spatialT;
    typedef typename FFTFourierType<spatialT>::FourierType fourierT;
    typedef VolumeAtAddressable< FFTWBuffer<spatialT> > spatialVolumeT;
    typedef
      typename FFTFourierType<spatialT>::FourierVolumeType
      fourierVolumeT;

  public:
    FFTOp(size_t cubeSize);

 	  ~FFTOp();


    void forward( 
      spatialVolumeT *spatial,
      fourierVolumeT *fourier) const;
    
    void backward( 
      fourierVolumeT *fourier,
      spatialVolumeT *spatial) const;

  protected:
    static typename FFTWPlanType<spatialT>::PlanType
      createForwardPlan(size_t cubeSize);
    
    static typename FFTWPlanType<spatialT>::PlanType
      createBackwardPlan(size_t cubeSize);

  protected:
    const typename FFTWPlanType<spatialT>::PlanType forwardPlan;
    const typename FFTWPlanType<spatialT>::PlanType backwardPlan;
    const size_t cubeSize;
};


#endif
