#ifndef Gauss_Newton_Ref_Grad_h
#define Gauss_Newton_Ref_Grad_h

#include "Gauss_Newton_Base.h"

template <
  typename _InterpolatorT,
  typename _ConvergenceTestT 
  >
class Gauss_Newton_Ref_Grad :
  Gauss_Newton_Base<_InterpolatorT, _ConvergenceTestT>{
  public:
    typedef Gauss_Newton_Base<_InterpolatorT, _ConvergenceTestT> Parent;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::ConvergenceTestT ConvergenceTestT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;

    Gauss_Newton_Ref_Grad(
      const InterpolatorT *interpRef, 
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      double *gradientAndHessianComputeTime
      ) :
      Parent(interpRef, refdz->cubeSize) {
      
      this->generateResidualGradientAndApproxHessian(
        &(this->residualGradient), &(this->approxResidualHessian),
        &(this->pointList),
        refdz, refdy, refdx, this->cubeSize, this->cubeCenter,
        gradientAndHessianComputeTime);
     
//      std::cout << "approxResidualHessian:" << std::endl <<
//        approxResidualHessian << std::endl;

      this->residualHessianLDL.compute(this->approxResidualHessian);
    }
    
    
    void minimize(
      const VolumeT *newVolume,
      const ParamT *initialParam,
      ParamT *finalParam,
      const size_t maxSteps = 20,
      const T stepSizeScale = 0.25,
      const T stepSizeLimit = 0,
      const ConvergenceTestT *convergenceTest = NULL, 
      size_t *elapsedSteps = NULL,
      double *elapsedTime = NULL 
      ) {

      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }

      Parent::minimize(newVolume, initialParam, finalParam,
        maxSteps, stepSizeScale, stepSizeLimit,
        convergenceTest, 
        elapsedSteps, NULL);
      
      if(NULL != elapsedTime) { 
        gettimeofday(&timeAfter, NULL);
  
        *elapsedTime =
          ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
          ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
      }
    }

  protected: 
    typedef typename Parent::NewVolVecT NewVolVecT;

};


#endif
