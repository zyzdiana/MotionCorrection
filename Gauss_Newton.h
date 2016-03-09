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
    typedef Matrix< T, Dynamic, Dynamic >  Matrix_gradP;
    typedef Matrix< T, Dynamic, 1 >  Vector_flatR;
    typedef Matrix< T, 3, Dynamic >  Matrix3X;
    typedef Matrix< T, 6, 6 >  Matrix66T;
    typedef Matrix< T, 6, 3 >  Matrix_M;
    typedef Matrix< T, 3, 3 > Matrix3T;
    typedef Matrix< T, 3, 1 > CoordT;
    typedef Matrix< T, 64, 1 > Vector64T;
    typedef Matrix< T, 6, 1 > Vector6T;

    const int derives_shape; 
    const int cubeSize;
    const CoordT cubeCenter;

    Gauss_Newton(TricubicInterpolator<VolumeT, coordT>  *interpolator) :
        interpolator(interpolator),
        derives_shape(interpolator->derives_shape),
        cubeSize(interpolator->cubeSize),
        cubeCenter(interpolator->cubeCenter) {
            // compute the gradient using axis derivatives from the interpolator
            gradP.resize(6, cubeSize*cubeSize*cubeSize);
            Matrix_M M;
            axis_derivatives = interpolator->compute_axis_derivatives();
            int idx = 0;

            for(int z = 0; z < cubeSize; ++z){
                for(int y = 0; y < cubeSize; ++y){
                    for(int x = 0; x < cubeSize; ++x){
                        M = get_M(z-cubeCenter(2), y-cubeCenter(1), x-cubeCenter(0));
                        gradP.col(idx) = M * axis_derivatives.col(idx);
                        idx++;
                    }
                }
            }
        }


    CoordT rotate_coords_transformation(T x, T y, T z, const Vector6T params){
        //compute vector norm
        T theta = sqrt(pow(params(3), 2) + pow(params(4), 2) + pow(params(5), 2) );
        //if rotation is zero
        if (theta == 0){
            return CoordT(y-params(0),x-params(1),z-params(2));
        }

        Matrix3T R;

        R = AngleAxisf(theta, params.tail(3)/params.tail(3).norm()).inverse();
        CoordT dest_coord =  R*(CoordT(y-params(0),x-params(1),z-params(2)) - cubeCenter) +cubeCenter;
        for (int i = 0; i < 3; ++i){
            dest_coord(i) = fmod(dest_coord(i) + cubeSize, cubeSize);
        }
        return dest_coord;
    }

    Matrix_M get_M(T x1,T x2,T x3){
            Matrix_M M;
            M(0, 0) = 0;
            M(1, 0) = 0;
            M(2, 0) = -1;
            M(3, 0) = -x2;
            M(4, 0) = x1;
            M(5, 0) = 0;

            M(0, 1) = 0;
            M(1, 1) = -1;
            M(2, 1) = 0;
            M(3, 1) = x3;
            M(4, 1) = 0;
            M(5, 1) = -x1;

            M(0, 2) = -1;
            M(1, 2) = 0;
            M(2, 2) = 0;
            M(3, 2) = 0;
            M(4, 2) = -x3;
            M(5, 2) = x2;
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
        CoordT dest_coords;
        T error = 0;
        while(counter < max_iter){
            cout << counter << " ";

            P_old = P_cur;
            P_cur = P_new;

            int idx = 0;
            error = 0;
            for(int i = 0; i < cubeSize; ++i){
                for(int j = 0; j < cubeSize; ++j){
                    for(int k = 0; k < cubeSize; ++k){
                        dest_coords = rotate_coords_transformation(j,i,k, P_cur);
                        flatR(idx) = vol_target->at(i,j,k) - interpolator->interp(dest_coords(0),dest_coords(1),dest_coords(2));
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
                    cout << "Converged in " << counter+1 << " iterations!" << endl;
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