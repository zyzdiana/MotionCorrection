#ifndef Gauss_Newton_h
#define Gauss_Newton_h

#include "Interpolation_3D.h"
#include "TricubicInterpolator.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <complex>
#include <vector>
#include <string>

using namespace std;
using namespace Eigen;

template <
    typename VolumeT,
    typename coordT>
class Gauss_Newton{
  public:
    typedef typename VolumeT::T T;

    // define some alias
    typedef Matrix< T, 6, Dynamic >  Matrix_gradP;
    typedef Matrix< T, Dynamic, 1 >  Vector_flatR;
    typedef Matrix< T, 3, Dynamic >  Matrix3X;
    typedef Matrix< T, 4, Dynamic >  Matrix4X;
    typedef Matrix< T, 6, 6 >  Matrix66T;
    typedef Matrix< T, 6, 3 >  Matrix_M;
    typedef Matrix< T, 3, 3 > Matrix3T;
    typedef Matrix< T, 3, 1 > CoordT;
    typedef Matrix< T, 4, 1 > PointT;
    typedef Matrix< T, 64, 1 > Vector64T;
    typedef Matrix< T, 6, 1 > Vector6T;

    typedef Translation< T, 3> Translation3T;
    typedef AngleAxis< T > AngleAxisT;

    const int derives_shape; 
    const int cubeSize;
    const CoordT cubeCenter;
    std::vector<T> weights;
    Matrix4X grid_points;


    Gauss_Newton(TricubicInterpolator<VolumeT, coordT>  *interpolator) :
        interpolator(interpolator),
        derives_shape(interpolator->derives_shape),
        cubeSize(interpolator->cubeSize),
        cubeCenter(interpolator->cubeCenter) {
            // initialization
            gradP.resize(6, cubeSize*cubeSize*cubeSize);
            grid_points.resize(4, cubeSize*cubeSize*cubeSize);
            weights.resize(cubeSize*cubeSize*cubeSize);

            // get axis derivatives from the interpolator
            Matrix_M M;
            axis_derivatives = interpolator->compute_axis_derivatives();
            // compute radius for masking
            T radius = cubeSize/2;

            int idx = 0;
            for(int z = 0; z < cubeSize; ++z){
                for(int y = 0; y < cubeSize; ++y){
                    for(int x = 0; x < cubeSize; ++x){
                        grid_points.col(idx) = PointT(z-cubeCenter(2),y-cubeCenter(1),x-cubeCenter(0),1);
                        // compute mask
                        T n = (CoordT(z,y,x)-cubeCenter).norm();
                        weights[idx] = window(n, radius);
                        // compute the gradient using axis derivatives from the interpolator
                        M = get_M(CoordT(z,y,x)-cubeCenter);
                        gradP.col(idx) = M * axis_derivatives.col(idx) * weights[idx];
                        idx++;
                    }
                }
            }
        }

    static T window(T n, T radius, T d = 0.4){
        T tmp = n/radius - 0.75;
        if (tmp < 0){
            return 1.0;
        }
        else{
            if((tmp/d > -0.5) and (tmp/d < 0.5)){
                return cos(PI*(tmp/d));
            }
            else{
                return 0.0;
            }
        }
    }


    // CoordT rotate_coords_transformation(T x, T y, T z, const Vector6T params){
    //     //compute vector norm
    //     T theta = params.tail(3).norm();
    //     //if rotation is zero
    //     if (theta == 0){
    //         return CoordT(y-params(0),x-params(1),z-params(2));
    //     }

    //     Matrix3T R;

    //     R = AngleAxisf(theta, params.tail(3)/theta).inverse();
    //     CoordT dest_coord =  R*(CoordT(y-params(0),x-params(1),z-params(2)) - cubeCenter) +cubeCenter;
    //     for (int i = 0; i < 3; ++i){
    //         dest_coord(i) = fmod(dest_coord(i) + cubeSize, cubeSize);
    //     }
    //     return dest_coord;
    // }

    // CoordT rotate_coords_transformation(T x, T y, T z, const Vector6T params){
    //     //compute vector norm
    //     T theta = params.tail(3).norm();
    //     //if rotation is zero
    //     if (theta == 0){
    //         return CoordT(y,x,z) - params.head(3);
    //     }
    //     Transform<T,3,Affine> transform = Translation3T(cubeCenter)  * AngleAxisT(theta, params.tail(3)/theta).inverse() *  Translation3T(-cubeCenter) * Translation3T(-params.head(3));

    //     CoordT dest_coord =  transform * CoordT(y,x,z);

    //     for (int i = 0; i < 3; ++i){
    //         dest_coord(i) = fmod(dest_coord(i) + cubeSize, cubeSize);
    //     }
    //     return dest_coord;
    // }

    Matrix4X rotate_coords_transformation(const Vector6T params){
        //compute vector norm
        T theta = params.tail(3).norm();

        // apply translation
        Transform<T,3,Affine> translation(Translation3T(-params.head(3)));
        Matrix4X transformed_coords = translation.matrix() * grid_points;
        //if rotation is zero
        if (theta == 0){
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

    Matrix_M get_M(CoordT pt){
            Matrix_M M;
            M(0, 0) = 0;
            M(1, 0) = 0;
            M(2, 0) = -1;
            M(3, 0) = -pt(1);
            M(4, 0) = pt(0);
            M(5, 0) = 0;

            M(0, 1) = 0;
            M(1, 1) = -1;
            M(2, 1) = 0;
            M(3, 1) = pt(2);
            M(4, 1) = 0;
            M(5, 1) = -pt(0);

            M(0, 2) = -1;
            M(1, 2) = 0;
            M(2, 2) = 0;
            M(3, 2) = 0;
            M(4, 2) = -pt(2);
            M(5, 2) = pt(1);
            return M;
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

            P_old = P_cur;
            P_cur = P_new;

            int idx = 0;
            error = 0;
            dest_coords = rotate_coords_transformation(P_cur);
            for(int i = 0; i < cubeSize; ++i){
                for(int j = 0; j < cubeSize; ++j){
                    for(int k = 0; k < cubeSize; ++k){
                        //dest_coords = rotate_coords_transformation(j,i,k, P_cur);
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
                    if (abs(P_new(i) - P_cur(i)) <  1e-5){
                        count += 1;
                    }
                }
                if (count == 6){//converged
                    // cout << "Converged in " << counter+1 << " iterations!" << endl;
                    break;
                }
            }
            ++counter;
        }
        return P_new;
    }
  protected:
    const Interpolator3D<VolumeT, coordT> *interpolator;
    Matrix_gradP gradP;
    Matrix3X axis_derivatives;
};


#endif