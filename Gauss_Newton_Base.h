#ifndef Gauss_Newton_Base_h
#define Gauss_Newton_Base_h

//#include <iostream>

#include <sys/time.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

template <
  typename _InterpolatorT 
  >
class Gauss_Newton_Base{
  public:
    typedef _InterpolatorT InterpolatorT;
    typedef typename InterpolatorT::VolumeT VolumeT;
    typedef typename InterpolatorT::CoordT CoordT;
    typedef typename VolumeT::value_type T;
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    Gauss_Newton_Base(
      const InterpolatorT *interpRef, 
      const size_t cubeSize 
      ) :
      interpRef(interpRef),
      cubeSize(cubeSize),
      cubeCenter(cubeCenterFromCubeSize(cubeSize)),
      residualGradient(6, cubeSize * cubeSize * cubeSize),
      pointList(3, cubeSize * cubeSize * cubeSize),
      transformedPointList(3, cubeSize * cubeSize * cubeSize),
      interpPoints(cubeSize * cubeSize * cubeSize, 1),
      residual(cubeSize * cubeSize * cubeSize, 1) {
      
      generatePointList(&pointList, cubeSize, cubeCenter);
    }
    

  protected:
    typedef Eigen::Matrix< T, 6, Eigen::Dynamic > ResidualGradientT;
    typedef Eigen::Matrix< T, 6, 6 > ResidualHessianT;
    typedef Eigen::LDLT< ResidualHessianT, Eigen::Upper > ResidualHessianLDLT;
    typedef Eigen::Matrix< T, 3, Eigen::Dynamic > PointListT;
    typedef Eigen::Matrix< T, 3, 1 > PointT;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> ResidualT;
    typedef Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, 1 > > NewVolVecT;
    
    static CoordT cubeCenterFromCubeSize(const size_t cubeSize) {
        return ((CoordT) cubeSize)/(CoordT)2.0 - (CoordT)0.5;
    }

    static void generateResidualGradientAndApproxHessian(
      ResidualGradientT *residualGradient,
      ResidualHessianT *approxResidualHessian,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      const size_t cubeSize,
      const CoordT cubeCenter,
      double *elapsedTime = NULL 
      ) {

     
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }

      size_t offset = 0;

      Eigen::Matrix< T, 3, 6 > MMatrix;
      Eigen::Matrix< T, 1, 3> tempDerivVector;
      Eigen::Matrix< T, 1, 6> tempGradVector;

      approxResidualHessian->setZero(6, 6);
//      std::cout << "initial approxResidualHessian: " << *approxResidualHessian << std::endl;

      for(size_t z = 0; z < cubeSize; z++) {
        CoordT zIndex = ((CoordT) z) - cubeCenter;

        for(size_t y = 0; y < cubeSize; y++) {
          CoordT yIndex = ((CoordT) y) - cubeCenter;

          for(size_t x = 0; x < cubeSize; x++, offset++) {
            CoordT xIndex = ((CoordT) x) - cubeCenter;

            // Note we define the MMatrix differently than in the
            // Mathematica code, because we're going to use the ordering
            // dz, dy, dx instead of dx, dy, dz
            MMatrix <<
              -1,  0,  0,       0, -xIndex,  yIndex,
               0, -1,  0,  xIndex,       0, -zIndex,
               0,  0, -1, -yIndex,  zIndex,       0;

            tempDerivVector <<
              refdz->at(offset), 
              refdy->at(offset), 
              refdx->at(offset);

            
            tempGradVector.noalias() = tempDerivVector * MMatrix; 

            approxResidualHessian->template selfadjointView<Eigen::Upper>().rankUpdate(
              tempGradVector.transpose());
            
            residualGradient->col(offset) = tempGradVector;
          }
        }
      }
      
      if(NULL != elapsedTime) { 
        gettimeofday(&timeAfter, NULL);
  
        *elapsedTime =
          ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
          ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
      }
    }

    static void generatePointList(
      PointListT *pointList, const size_t cubeSize, const CoordT cubeCenter) {

      size_t offset = 0;

      for(size_t z = 0; z < cubeSize; z++) {
        CoordT zCoord = ((CoordT) z) - cubeCenter; 
        
        for(size_t y = 0; y < cubeSize; y++) {
          CoordT yCoord = ((CoordT) y) - cubeCenter; 
          
          for(size_t x = 0; x < cubeSize; x++, offset++) {
            CoordT xCoord = ((CoordT) x) - cubeCenter;

            pointList->col(offset) = PointT(zCoord, yCoord, xCoord);
          }
        }
      }
    }

    static void transformPointListWithParam(
      const ParamT *param,
      const PointListT *pointList,
      PointListT *transformedPointList) {

      typedef Eigen::AngleAxis<T> RotationT;
      typedef Eigen::Translation<T, 3> TranslationT;

      const Eigen::Matrix<T, 3, 1> rotVec = param->tail(3);

      const T angle = rotVec.norm();

      Eigen::Matrix<T, 3, 1> rotAxis;
      if(0 == angle) {
        rotAxis << 1, 0, 0;
      }
      else {
        rotAxis = rotVec.normalized();
      }
      
      transformedPointList->noalias() =
        ( RotationT(-angle, rotAxis) * TranslationT(-param->head(3)) ) *
        (*pointList);
    }

    void computeResidual(
      const VolumeT *newVol,
      const NewVolVecT *newVolVec,
      const ParamT *param) {

      transformPointListWithParam(
        param, &pointList, &transformedPointList);

      size_t pointListLength = cubeSize * cubeSize * cubeSize;

      PointT cubeCenterPoint(cubeCenter, cubeCenter, cubeCenter);

      for(size_t offset = 0; offset < pointListLength; offset++) {
        PointT curPoint =
          transformedPointList.col(offset) + cubeCenterPoint;

        interpPoints(offset, 0) =
          interpRef->interp(
            newVol->wrapIndex(curPoint(0)),
            newVol->wrapIndex(curPoint(1)),
            newVol->wrapIndex(curPoint(2))
          ); 
      }
 
//      std::cout << "interpPoints[0]: " << interpPoints(0) << std::endl;
//      std::cout << "newVolVec[0]: " << (*newVolVec)(0) << std::endl;

      residual.noalias() = interpPoints - (*newVolVec);
    }

  protected:
    const InterpolatorT *interpRef;
    const size_t cubeSize;
    const CoordT cubeCenter;
    ResidualGradientT residualGradient;
    ResidualHessianT approxResidualHessian;
    ResidualHessianLDLT residualHessianLDL;
    PointListT pointList;
    PointListT transformedPointList;
    ResidualT interpPoints;
    ResidualT residual;
};


#endif
