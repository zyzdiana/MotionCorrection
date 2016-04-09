#ifndef Gauss_Newton_h
#define Gauss_Newton_h

#include <iostream>

#include <sys/time.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

template <
  typename _InterpolatorT 
  >
class Gauss_Newton{
  public:
    typedef _InterpolatorT InterpolatorT;
    typedef typename InterpolatorT::VolumeT VolumeT;
    typedef typename InterpolatorT::CoordT CoordT;
    typedef typename VolumeT::value_type T;
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    Gauss_Newton(
      const InterpolatorT *interpRef, 
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      double *gradientAndHessianComputeTime
      ) :
      interpRef(interpRef),
      cubeSize(refdz->cubeSize),
      cubeCenter(cubeCenterFromCubeSize(cubeSize)),
      residualGradient(6, cubeSize * cubeSize * cubeSize),
      pointList(3, cubeSize * cubeSize * cubeSize),
      transformedPointList(3, cubeSize * cubeSize * cubeSize),
      interpPoints(cubeSize * cubeSize * cubeSize, 1),
      residual(cubeSize * cubeSize * cubeSize, 1) {
      
      generateResidualGradientAndApproxHessian(
        &residualGradient, &approxResidualHessian,
        refdz, refdy, refdx, cubeSize, cubeCenter,
        gradientAndHessianComputeTime);
      
      residualHessianLDL.compute(approxResidualHessian);

      generatePointList(&pointList, cubeSize, cubeCenter);
    }
    
    void minimize(
      const VolumeT *newVolume,
      const ParamT *initialParam,
      ParamT *finalParam,
      const size_t maxSteps = 20,
      const T paramUpdateNormLimit = 1e-6,
      size_t *elapsedSteps = NULL, 
      double *elapsedTime = NULL 
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }

      ParamT curParam = *initialParam;

      NewVolVecT newVolVec(
        newVolume->buffer, cubeSize * cubeSize * cubeSize, 1);

      ParamT reducedResidual;
      size_t step = 0;
      
      for(; step < maxSteps; step++) {
//        std::cout << "-------" << std::endl; 
//        std::cout << "step " << step << std::endl; 
//        std::cout << "curParam: " << std::endl <<
//          curParam.transpose() << std::endl; 

        computeResidual(newVolume, &newVolVec, &curParam);

        // at this point we could check if this new residual is better
        // than the previous residual, and if not, we could respond
        // by taking some other search step. This test can be expensive,
        // though, so for right now we skip it and hope things are improving

        reducedResidual.noalias() = residualGradient * residual;
      
//        std::cout << "reducedResidual: " << std::endl <<
//          reducedResidual.transpose() << std::endl; 

        ParamT paramUpdate;
        paramUpdate.noalias() = residualHessianLDL.solve(reducedResidual);

        curParam -= paramUpdate;


        //check for convergence
        if(paramUpdate.norm() < paramUpdateNormLimit) {
          break; 
        }
      }

      *finalParam = curParam;

      if(NULL != elapsedSteps) {
        *elapsedSteps = step + 1; 
      }

      if(NULL != elapsedTime) { 
        gettimeofday(&timeAfter, NULL);
  
        *elapsedTime =
          ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
          ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
      }
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

      
//      std::cout << "param:" << std::endl <<
//        *param << std::endl;
      
//      std::cout << "pointList(0):" << std::endl <<
//        pointList.col(0) << std::endl;
      
//      std::cout << "transformedPointList(0):" << std::endl <<
//        transformedPointList.col(0) << std::endl;
      
//      std::cout << "interp(0,0,0):" << std::endl <<
//        interpRef->interp(0,0,0) << std::endl;

//      std::cout << "interpPoints:" << std::endl <<
//        interpPoints.head(10) << std::endl;
      
//      std::cout << "newVolVec:" << std::endl <<
//        newVolVec->head(10) << std::endl;

      residual.noalias() = (*newVolVec) - interpPoints;
      
//      std::cout << "residual:" << std::endl <<
//        residual.head(10) << std::endl;
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

/*
    typedef typename InterpolatorT::T T;
    typedef typename InterpolatorT::VolumeT VolumeT;
    typedef typename InterpolatorT::CoordT CoordT;

    // define some alias
    typedef Eigen::Matrix< T, 6, Dynamic >  Matrix_gradP;
    typedef Eigen::Matrix< T, Dynamic, 1 >  Vector_flatR;
    typedef Eigen::Matrix< T, 3, Dynamic >  Matrix3X;
    typedef Eigen::Matrix< T, 4, Dynamic >  Matrix4X;
    typedef Eigen::Matrix< T, 6, 6 >  Matrix66T;
    typedef Eigen::Matrix< T, 6, 3 >  Matrix_M;
    typedef Eigen::Matrix< T, 3, 3 > Matrix3T;
    typedef Eigen::Matrix< T, 3, 1 > coordT;
    typedef Eigen::Matrix< T, 4, 1 > PointT;
    typedef Eigen::Matrix< T, 64, 1 > Vector64T;
    typedef Eigen::Matrix< T, 6, 1 > Vector6T;

    typedef Eigen::Translation< T, 3> Translation3T;
    typedef Eigen::AngleAxis< T > AngleAxisT;


  public:
    Gauss_Newton(
        const InterpolatorT  *interpRef, 
        const VolumeT *refdx,
        const VolumeT *refdy,
        const VolumeT *refdz) :
        interpRef(interpRef),
        refdx(refdx),
        refdy(refdy),
        refdz(refdz),
        cubeSize(refdx->cubeSize),
        cubeCenter(cubeCenterFromCubeSize(cubeSize)){
            // initialization
            gradP.resize(6, cubeSize*cubeSize*cubeSize);
            grid_points.resize(4, cubeSize*cubeSize*cubeSize);
            weights.resize(cubeSize*cubeSize*cubeSize);

            // get axis derivatives from the interpolator
            Matrix_M M;
            coordT axis_derivatives;

            VolumeT dz(cubeSize);
            VolumeT dy(cubeSize);
            VolumeT dx(cubeSize);
            differentiator->zDerivative(&dz);
            differentiator->yDerivative(&dy);
            differentiator->xDerivative(&dx);
            // compute radius for masking
            T radius = cubeSize/2;

            int idx = 0;
            for(int z = 0; z < cubeSize; ++z){
                for(int y = 0; y < cubeSize; ++y){
                    for(int x = 0; x < cubeSize; ++x){
                        grid_points.col(idx) = PointT(z-cubeCenter(2),y-cubeCenter(1),x-cubeCenter(0),1);
                        // compute mask
                        T n = (coordT(z,y,x)-cubeCenter).norm();
                        weights[idx] = window(n, radius);
                        //weights[idx] = 1;
                        // compute the gradient using axis derivatives from the interpolator
                        M = get_M(coordT(z,y,x)-cubeCenter);

                        axis_derivatives << dx.at(idx) , dy.at(idx) , dz.at(idx);
                        gradP.col(idx) = M * axis_derivatives * weights[idx];
                        idx++;
                    }
                }
            }

        }

  protected:
    const size_t cubeSize;
    const CoordT cubeCenter;
    const std::vector<T> weights;
    Matrix4X grid_points;



    static T window(T n, T radius, T d = 0.4){
        T tmp = n/radius - 0.75;
        if (tmp < 0){
            return 1.0;
        }
        else{
            if((tmp/d > -0.5) and (tmp/d < 0.5)){
                return cos(M_PI*(tmp/d));
            }
            else{
                return 0.0;
            }
        }
    }



    Matrix4X rotate_coords_transformation(const Vector6T params){
        //compute vector norm
        T theta = params.tail(3).norm();

        // apply translation
        //coordT trans_vec(-params(1),-params(0),-params(2));
        Transform<T,3,Affine> translation(Translation3T(-params.head(3)));
        //Transform<T,3,Affine> translation( ( Translation3T(trans_vec) ) );
        Matrix4X transformed_coords = translation.matrix() * grid_points;

        //if rotation is zero
        if (theta == 0){
            Transform<T,3,Affine> t_back( ( Translation3T(cubeCenter) ) );
            transformed_coords = t_back.matrix() * transformed_coords;
            return transformed_coords;
        }
        Transform<T,3,Affine> rotation = Translation3T(cubeCenter)  * AngleAxisT(theta, params.tail(3)/theta).inverse();
        transformed_coords =  rotation * transformed_coords;

        for (int idx = 0; idx < cubeSize*cubeSize*cubeSize; ++idx){
            for (int i = 0; i < 3; ++i){
                transformed_coords(i,idx) = fmod(transformed_coords(i,idx) + cubeSize, cubeSize);
            }
        }
        return transformed_coords;
    }


    Vector6T gauss_newton(VolumeT *vol_target, int max_iter){

        Vector6T P_old;
        P_old.fill(0);
        Vector6T P_cur;
        P_cur.fill(0);
        Vector6T P_new;
        P_new.fill(0);

        std::vector<T> errors;
        errors.push_back(10);

        // local variables for the loop
        Vector_flatR flatR(cubeSize*cubeSize*cubeSize);
        Vector6T deltaP;
        Matrix66T iJrTJr;


        int counter = 0;
        Matrix4X dest_coords;
        T error = 0;
        while(counter < max_iter){
            //cout << "Counter: " << counter << endl;

            P_old = P_cur;
            P_cur = P_new;

            int idx = 0;
            error = 0;
            dest_coords = rotate_coords_transformation(P_cur);
            for(int i = 0; i < cubeSize; ++i){
                for(int j = 0; j < cubeSize; ++j){
                    for(int k = 0; k < cubeSize; ++k){
                        flatR(idx) = vol_target->at(i,j,k) - interpolator->interp(dest_coords(0,idx),dest_coords(1,idx),dest_coords(2,idx)) * weights[idx];
                        error += flatR(idx)*flatR(idx);
                        ++idx;
                    }
                }
            }
            //if error increases, step back
            if(error > errors.back()){
                P_old = P_cur;
            }
            else{
                errors.push_back(error);
                deltaP = (gradP * gradP.transpose()).inverse() * (-gradP*flatR);
                P_new = P_cur - deltaP;

                //check for convergence
                int count = 0;
                for (int i = 0; i < 6; ++i){
                    if (abs(deltaP(i)) <  1e-5){
                        count += 1;
                    }
                }
                if (count == 6){//converged
                    cout << "Converged in " << counter+1 << " iterations!" << endl;
                    break;
                }
            }
            ++counter;
        }
        return P_new;
    }
  protected:
    const Interpolator3D<VolumeT, CoordT> *interpolator;
    const VolumeT *volume;
    const Matrix_gradP gradP;
    Matrix3X axis_derivatives;
};
*/

#endif
