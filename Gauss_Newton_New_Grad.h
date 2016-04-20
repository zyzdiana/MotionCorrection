#ifndef Gauss_Newton_New_Grad_h
#define Gauss_Newton_New_Grad_h

#include "Gauss_Newton_Base.h"

template <
  typename _InterpolatorT 
  >
class Gauss_Newton_New_Grad : Gauss_Newton_Base<_InterpolatorT>{
  public:
    typedef Gauss_Newton_Base<_InterpolatorT> Parent;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;

    Gauss_Newton_New_Grad(
      const InterpolatorT *interpRef, 
      const size_t cubeSize
      ) :
      Parent(interpRef, cubeSize) {}
    
    void minimize(
      const VolumeT *newVolume,
      const VolumeT *newdz,
      const VolumeT *newdy,
      const VolumeT *newdx,
      const ParamT *initialParam,
      ParamT *finalParam,
      const size_t maxSteps = 20,
      const T paramUpdate2NormLimit = 0,
      const T paramUpdateInfinityNormLimit = 0,
      size_t *elapsedSteps = NULL, 
      double *elapsedTime = NULL,
      double *gradientAndHessianComputeTime = NULL
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }
      
      this->generateResidualGradientAndApproxHessian(
        &(this->residualGradient), &(this->approxResidualHessian),
        newdz, newdy, newdx, this->cubeSize, this->cubeCenter,
        gradientAndHessianComputeTime);
     
      this->residualHessianLDL.compute(this->approxResidualHessian);

      ParamT curParam = *initialParam;

      NewVolVecT newVolVec(
        newVolume->buffer,
        this->cubeSize * this->cubeSize * this->cubeSize, 1);

      ParamT reducedResidual;
      size_t step = 0;
      
      for(; step < maxSteps; step++) {
//        std::cout << "-------" << std::endl; 
//        std::cout << "step " << step << std::endl; 
//        std::cout << "curParam: " << std::endl <<
//          curParam.transpose() << std::endl; 

        this->computeResidual(newVolume, &newVolVec, &curParam);

//        std::cout << "residual[0]: " << residual(0) << std::endl;
//        std::cout << "residualGradient.col(0): " << std::endl <<
//          residualGradient.col(0) << std::endl;

        // at this point we could check if this new residual is better
        // than the previous residual, and if not, we could respond
        // by taking some other search step. This test can be expensive,
        // though, so for right now we skip it and hope things are improving

        reducedResidual.noalias() = this->residualGradient * this->residual;
     
//        std::cout << "reducedResidual: " << std::endl <<
//          reducedResidual << std::endl;

        // This equation solves the parameter update, but with signs negated
        ParamT negParamUpdate;
        negParamUpdate.noalias() = this->residualHessianLDL.solve(reducedResidual);

        // Subtract the negated update (i.e., add the correct-sign update!)
        curParam -= negParamUpdate;


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

};


#endif
