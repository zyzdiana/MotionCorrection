#ifndef Moving_Weighted_Gauss_Newton_New_Grad_h
#define Moving_Weighted_Gauss_Newton_New_Grad_h

#include "Gauss_Newton.h"

#include "StaticWeightedResidualOp.h"

#include "MovingWeightedResidualGradientAndHessian.h"

#include <fcntl.h>
#include <unistd.h>
#ifdef LINUX
#include <stdlib.h>
#endif 

#include <iostream>

#include <cfloat>

template <
  typename _InterpolatorT,
  typename _ParamAccumulatorT,
  typename _WeightFuncT,
  typename _WeightGradientFuncT,
  typename _ConvergenceTestT = void
  >
class Moving_Weighted_Gauss_Newton_New_Grad : 
  Gauss_Newton <
    StaticWeightedResidualOp<_WeightFuncT, _InterpolatorT>,
    MovingWeightedResidualGradientAndHessian<
      typename _InterpolatorT::VolumeT,
      typename _InterpolatorT::CoordT,
      _WeightFuncT,
      _WeightGradientFuncT,
      StaticWeightedResidualOp<_WeightFuncT, _InterpolatorT>
      >,
    _InterpolatorT, 
    _ParamAccumulatorT,
    _ConvergenceTestT >{
  public:
    typedef MovingWeightedResidualGradientAndHessian<
      typename _InterpolatorT::VolumeT,
      typename _InterpolatorT::CoordT,
      _WeightFuncT,
      _WeightGradientFuncT,
      StaticWeightedResidualOp<_WeightFuncT, _InterpolatorT>
      > ResidualGradientAndHessianT;
    typedef Gauss_Newton <
      StaticWeightedResidualOp<_WeightFuncT, _InterpolatorT>,
      ResidualGradientAndHessianT, 
      _InterpolatorT,
      _ParamAccumulatorT,
      _ConvergenceTestT > Parent;
    typedef typename Parent::ResidualOpT ResidualOpT;
    typedef _WeightFuncT WeightFuncT;
    typedef _WeightGradientFuncT WeightGradientFuncT;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::ConvergenceTestT ConvergenceTestT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;

    Moving_Weighted_Gauss_Newton_New_Grad(
      const InterpolatorT *interpRef,
      const size_t cubeSize,
      WeightFuncT *weightFunc,
      WeightGradientFuncT *weightGradientFunc
      ) :
      Parent(
        interpRef,
        new ResidualOpT(cubeSize, weightFunc, interpRef),
        new ResidualGradientAndHessianT(
          cubeSize, weightFunc, weightGradientFunc, NULL, NULL),
        cubeSize) {
        this->residualGradientAndHessian->setResidualOp(this->residualOp);
        this->residualGradientAndHessian->setInitialPointList(
          &(this->pointList));
        }

  protected:
    typedef typename Parent::NewVolVecT NewVolVecT;
    typedef typename Parent::PointListT PointListT;

  public:
    
    void minimize(
      const VolumeT *newVolume,
      const VolumeT *newdz,
      const VolumeT *newdy,
      const VolumeT *newdx,
      const ParamT *initialParam,
      ParamT *finalParam,
      const size_t maxSteps = 20,
      const T stepSizeScale = 0.25,
      ConvergenceTestT *convergenceTest = NULL, 
      size_t *elapsedSteps = NULL, 
      size_t *elapsedSearchSteps = NULL, 
      double *elapsedTime = NULL,
      double *gradientAndHessianComputeTime = NULL
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }
      
      this->computeResidual(newVolume, initialParam);

      this->residualGradientAndHessian->
        initializeResidualGradientAndApproxHessian(
          &(this->pointList),
          initialParam,
          newdz, newdy, newdx, 
          &(this->residualGradient), &(this->approxResidualHessian),
          &(this->residualHessianLDL),
          gradientAndHessianComputeTime);

      Parent::minimize(newVolume, newdz, newdy, newdx,
        initialParam, finalParam,
        maxSteps, stepSizeScale, 
        convergenceTest, 
        elapsedSteps, elapsedSearchSteps);

      if(NULL != elapsedTime) { 
        gettimeofday(&timeAfter, NULL);
  
        *elapsedTime =
          ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
          ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
      }
    }

};


#endif
