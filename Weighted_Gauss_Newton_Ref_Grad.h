#ifndef Weighted_Gauss_Newton_Ref_Grad_h
#define Weighted_Gauss_Newton_Ref_Grad_h

#include "Gauss_Newton_Base.h"

#include <fcntl.h>
#include <unistd.h>
#ifdef LINUX
#include <stdlib.h>
#endif 

#include <iostream>

#include <cfloat>

template <
  typename _InterpolatorT,
  typename _CircularMaskOpT
  >
class Weighted_Gauss_Newton_Ref_Grad : Gauss_Newton_Base<_InterpolatorT>{
  public:
    typedef Gauss_Newton_Base<_InterpolatorT> Parent;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;
    typedef _CircularMaskOpT CircularMaskOpT;

    Weighted_Gauss_Newton_Ref_Grad(
      const InterpolatorT *interpRef,
      const CircularMaskOpT *circularMaskOp,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      double *gradientAndHessianComputeTime
      ) :
      Parent(interpRef, refdz->cubeSize),
      circularMaskOp(circularMaskOp) {
     
      VolumeT weightedRefdz(this->cubeSize);
      VolumeT weightedRefdy(this->cubeSize);
      VolumeT weightedRefdx(this->cubeSize);

      circularMaskOp->applyMask(refdz, &weightedRefdz); 
      circularMaskOp->applyMask(refdy, &weightedRefdy); 
      circularMaskOp->applyMask(refdx, &weightedRefdx);

      this->generateResidualGradientAndApproxHessian(
        &(this->residualGradient), &(this->approxResidualHessian),
        &weightedRefdz, &weightedRefdy, &weightedRefdx,
        this->cubeSize, this->cubeCenter,
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
      const T paramUpdate2NormLimit = 0,
      const T paramUpdateInfinityNormLimit = 0,
      size_t *elapsedSteps = NULL, 
      double *elapsedTime = NULL 
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }

      ParamT curParam = *initialParam;
      ParamT prevParam = curParam;

      T prevResidualNorm = FLT_MAX;

      T stepSize = 1.0;

      NewVolVecT newVolVec(
        newVolume->buffer,
        this->cubeSize * this->cubeSize * this->cubeSize, 1);

      ParamT reducedResidual;
      size_t step = 0;

      for(; step < maxSteps; step++) {
//        std::cout << "-------" << std::endl; 
//        std::cout << "step " << step << std::endl; 
//        std::cout << "stepSize " << stepSize << std::endl; 
//        std::cout << "curParam: " << std::endl <<
//          curParam.transpose() << std::endl; 

        this->computeResidual(newVolume, &newVolVec, &curParam);

        circularMaskOp->applyMask(&this->residual);

//        std::cout << "residual[0]: " << residual(0) << std::endl;
//        std::cout << "residualGradient.col(0): " << std::endl <<
//          residualGradient.col(0) << std::endl;

        T newResidualNorm = this->residual.norm();
        
//        std::cout << "newResidualNorm " << newResidualNorm << std::endl; 

        // If the residual has become smaller since the last step, then proceed
        if(newResidualNorm <= prevResidualNorm) {
          prevResidualNorm = newResidualNorm;
          prevParam = curParam;

          reducedResidual.noalias() = this->residualGradient * this->residual;
      
  //        std::cout << "reducedResidual: " << std::endl <<
  //          reducedResidual << std::endl;
  
          // This equation solves the parameter update, but with signs negated
          ParamT negParamUpdate;
          negParamUpdate.noalias() = this->residualHessianLDL.solve(reducedResidual);
          
  //        std::cout << "negParamUpdate: " << std::endl <<
  //          negParamUpdate << std::endl;
  
          // Subtract the negated update (i.e., add the correct-sign update!)
          curParam -= negParamUpdate * stepSize;
  
  
          //checks for convergence
          if(paramUpdate2NormLimit > 0) {
            if(negParamUpdate.norm() < paramUpdate2NormLimit) {
              step++; 
              break; 
            }
          }
  
          if(paramUpdateInfinityNormLimit > 0) {
            bool largerThanLimit = false;
  
            for (int i = 0; (! largerThanLimit) && (i < 6); ++i) {
              largerThanLimit =  
                (std::abs(negParamUpdate(i)) >  paramUpdateInfinityNormLimit);
            }
            if (! largerThanLimit) {
              step++; 
              break;
            }
          }
        }
        // If the residual has become larger since the last update, step back
        // and reduce the step size
        else {
          curParam = prevParam;
          stepSize *= stepSizeScale;

          if(stepSizeLimit > 0) {
            if(stepSizeLimit > stepSize) {
              step++;
              break;
            }
          }
        }
      }

      *finalParam = curParam;

      if(NULL != elapsedSteps) {
        *elapsedSteps = step; 
      }

      if(NULL != elapsedTime) { 
        gettimeofday(&timeAfter, NULL);
  
        *elapsedTime =
          ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
          ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
      }
    }
    
  protected: 
    typedef typename Parent::NewVolVecT NewVolVecT;
    const CircularMaskOpT *circularMaskOp;

};


#endif
