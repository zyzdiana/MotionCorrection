#ifndef Gauss_Newton_h
#define Gauss_Newton_h

#include "Interpolation_3D.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <complex>
#include <vector>
#include <string>

using namespace std;
static const int MAX_GRADP_SIZE = 40*40*40;

template <typename T>
class Gauss_Newton{
  public:
	typedef Interpolation_3D<T> interp;
	typedef vector<T> vector1D;
	typedef vector<T*> pointer1D; //vector of pointers for storing the coefficients.


    static T window(T n, T radius, T d){
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


    static void get_mast_weights(vector1D *mask_weights, T radius, T d){
        int shape = cbrt(mask_weights->size());

        T ox = (T(shape)-1)/2;
        T oy = (T(shape)-1)/2;
        T oz = (T(shape)-1)/2;

        for(int i = 0; i < shape; ++i){
            for(int j = 0; j < shape; ++j){
                for(int k = 0; k < shape; ++k){
                    T n = sqrt((i-ox)*(i-ox) + (j-oy)*(j-oy) + (k-oz)*(k-oz));
                    mask_weights->at(j+shape*i+shape*shape*k) = window(n, radius, d);
                }
            }
        }
    }

    //function to compute the inverse of a 3 by 3 matrix
    static void inverse3(const T rot_Matrix[3][3], T inverse[3][3]){
        T det = 0;
        // compute determinant
        for (int i = 0; i < 3; ++i){
            det += (rot_Matrix[0][i]*(rot_Matrix[1][(i+1)%3] * rot_Matrix[2][(i+2)%3] - rot_Matrix[1][(i+2)%3]*rot_Matrix[2][(i+1)%3]));
        }
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                inverse[j][i] = ((rot_Matrix[(i+1)%3][(j+1)%3] * rot_Matrix[(i+2)%3][(j+2)%3]) - (rot_Matrix[(i+1)%3][(j+2)%3] * rot_Matrix[(i+2)%3][(j+1)%3]) )/det;
            }
        }       
    }


    static void rotate_coords_transformation(T result[3], T x, T y, T z, const T params[6], T ox, T oy, T oz, T k, int shape){
        //compute vector norm
        T l = sqrt(pow(params[3], 2) + pow(params[4], 2) + pow(params[5], 2) )/k;
        
        //if rotation is zero
        if (l == 0){
            result[0] = fmod((x-params[1]+shape), shape);
            result[1] = fmod((y-params[0]+shape), shape);
            result[2] = fmod((z-params[2]+shape), shape); 
            return;
        }


        T s = sin(l/2.0);
        T alpha = cos(l/2.0);
        T beta = s*params[3]/l;
        T gamma = s*params[4]/l;
        T delta = s*params[5]/l;

        T rotMatrix[3][3];
        T R[3][3];

        rotMatrix[0][0] = alpha*alpha+beta*beta-gamma*gamma-delta*delta;
        rotMatrix[0][1] = 2*(beta*gamma-alpha*delta);
        rotMatrix[0][2] = 2*(beta*delta+alpha*gamma);

        rotMatrix[1][0] = 2*(beta*gamma+alpha*delta);
        rotMatrix[1][1] = alpha*alpha-beta*beta+gamma*gamma-delta*delta;
        rotMatrix[1][2] = 2*(gamma*delta-alpha*beta);

        rotMatrix[2][0] = 2*(beta*delta-alpha*gamma);
        rotMatrix[2][1] = 2*(gamma*delta+alpha*beta);
        rotMatrix[2][2] = alpha*alpha-beta*beta-gamma*gamma+delta*delta;

        inverse3(rotMatrix, R);

        T yy = x - ox - params[1];
        T xx = y - oy - params[0];
        T zz = z - oz - params[2];

        result[1] = fmod((R[0][0]*xx + R[0][1]*yy + R[0][2]*zz + oy + shape), shape);
        result[0] = fmod((R[1][0]*xx + R[1][1]*yy + R[1][2]*zz + ox + shape), shape);
        result[2] = fmod((R[2][0]*xx + R[2][1]*yy + R[2][2]*zz + oz + shape), shape);
    }


    static void get_M(T M[6][3], T x1_org,T x2_org,T x3_org, T k){
            T x1 = x1_org/k;
            T x2 = x2_org/k;
            T x3 = x3_org/k;
            M[0][0] = 0;
            M[1][0] = 0;
            M[2][0] = -1;
            M[3][0] = -x2;
            M[4][0] = x1;
            M[5][0] = 0;

            M[0][1] = 0;
            M[1][1] = -1;
            M[2][1] = 0;
            M[3][1] = x3;
            M[4][1] = 0;
            M[5][1] = -x1;

            M[0][2] = -1;
            M[1][2] = 0;
            M[2][2] = 0;
            M[3][2] = 0;
            M[4][2] = -x3;
            M[5][2] = x2;
    }


    // Function for computing the dot product of Mtensor and derivatives
    static T* dot3(T grad_P[6] ,const T M[6][3], const T* derivs){ 
        int indicies[3] = {1,4,16};

        for (int i = 0; i < 6; ++i){
            T tmp = 0;
            for (int j = 0; j < 3; ++ j){
                tmp += M[i][j] * derivs[indicies[j]];
            }
            grad_P[i] = tmp;
        }
        return grad_P;
    }

    static T* dot3_mask(T grad_P[6] ,const T M[6][3], const T* derivs, T mask_weight){ 
        int indicies[3] = {1,4,16};

        for (int i = 0; i < 6; ++i){
            T tmp = 0;
            for (int j = 0; j < 3; ++ j){
                tmp += M[i][j] * derivs[indicies[j]];
            }
            grad_P[i] = tmp * mask_weight;
        }
        return grad_P;
    }

    static void get_gradient_P(pointer1D *gradientP, const pointer1D *derivatives, T divide_factor, bool mask){
        int derives_shape = cbrt(derivatives->size());
        int shape = derives_shape - 30;

        static T gradP[MAX_GRADP_SIZE][6];

        //find center of the volume
        T ox = (float(shape)-1)/2;
        T oy = (float(shape)-1)/2;
        T oz = (float(shape)-1)/2;

        T M[6][3];
        if(mask){
            vector1D mask_weights(shape*shape*shape);
            T rad = float(shape)/2.0;
            get_mast_weights(&mask_weights, rad, 0.4);
            int idx = 0;
            for(int i = 0; i < shape; ++i){
                for(int j = 0; j < shape; ++j){
                    for(int k = 0; k < shape; ++k){
                        get_M(M, i-ox, j-oy, k-oz, divide_factor);
                        gradientP->at(idx) = dot3_mask(gradP[idx], M, 
                                                      derivatives->at(derives_shape*derives_shape*(i+15) + derives_shape*(j+15) + (k+15)), 
                                                      mask_weights[idx]);
                        ++idx;
                    }
                }
            }           
        }
        else{
            int idx = 0;
            for(int i = 0; i < shape; ++i){
                for(int j = 0; j < shape; ++j){
                    for(int k = 0; k < shape; ++k){
                        get_M(M, i-ox, j-oy, k-oz, divide_factor);
                        // gradientP->at(idx) = dot3(gradP[idx], M, derivatives->at(derives_shape*derives_shape*(i+15) + derives_shape*(j+15) + (k+15)));
                        gradientP->at(idx) = dot3(gradP[idx], M, derivatives->at(derives_shape*derives_shape*(i+15) + derives_shape*(j+15) + (k+15)));
                        idx++;
                    }
                }
            }
        }
    }

    static void transpose(pointer1D* matrix_T, const pointer1D* matrix){
        static T gradP_T[6][MAX_GRADP_SIZE];
        for (int j =0; j < matrix_T->size(); ++j){
            for (int i = 0; i < matrix->size(); ++i){
                gradP_T[j][i] = (matrix->at(i))[j];
            }
            matrix_T->at(j) = gradP_T[j];
        }
    }

    static void dotn(vector1D *result, pointer1D *gradP_T, vector1D *residual){
        for (int i = 0; i < gradP_T->size(); ++i){
            T tmp = 0;
            for (int j = 0; j < residual->size(); ++j){
                tmp +=  (gradP_T->at(i))[j] * residual->at(j);
            }
            result->at(i) = -tmp;
        }
    }

    static void dot(T result[6], T matrix[36], vector1D * vector){
        for (int i = 0; i < 6; ++i){
            T tmp = 0;
            for (int j = 0; j < 6; ++j){
                tmp += matrix[i*6+j] * vector->at(j);
            }
            result[i] = tmp;
        }
    }

    static void dot_matrix(T res[6][6], pointer1D *gradP_T, pointer1D *gradP){
        for (int i = 0; i < 6; ++i){
            for (int j = 0; j < 6; ++j){
                T tmp = 0;
                for (int idx = 0; idx< gradP->size(); ++idx) 
                    tmp +=  (gradP_T->at(i))[idx] * (gradP->at(idx))[j];
                res[i][j] = tmp;
            }
        }
    }

    static void inverse6(T input[6][6], T matrix[36]){
        int size = 6;
        for (int i = 0; i < size*size; ++i){
            matrix[i] = input[i/6][i%6];
        }
        for (int i = 1; i < size; ++i) 
            matrix[i] /= matrix[0]; // normalize row 0
        for (int i = 1; i < size; ++i) { 
            for (int j = i; j < size; ++j){ // do a column of L
                T sum = 0.0;
                for (int k = 0; k < i; ++k) 
                    sum += matrix[j*size+k] * matrix[k*size+i];
                matrix[j*size+i] -= sum;
            }

            if (i == size-1) continue;
            for (int j=i+1; j < size; ++j){  // do a row of U
                T sum = 0.0;
                for (int k = 0; k < i; ++k)
                    sum += matrix[i*size+k]*matrix[k*size+j];
                matrix[i*size+j] = (matrix[i*size+j]-sum) / matrix[i*size+i];
            }
        }
        for ( int i = 0; i < size; ++i)  // invert L
            for ( int j = i; j < size; ++j ){
                T x = 1.0;
                if ( i != j ) {
                x = 0.0;
                for ( int k = i; k < j; ++k ) 
                    x -= matrix[j*size+k]*matrix[k*size+i];
                }
                matrix[j*size+i] = x / matrix[j*size+j];
            }
        for ( int i = 0; i < size; ++i )   // invert U
            for ( int j = i; j < size; ++j ){
                if ( i == j ) continue;
                T sum = 0.0;
                for ( int k = i; k < j; ++k )
                    sum += matrix[k*size+j]*( (i==k) ? 1.0 : matrix[i*size+k] );
                matrix[i*size+j] = -sum;
            }

        for ( int i = 0; i < size; ++i )   // final inversion
            for ( int j = 0; j < size; ++j ){
                T sum = 0.0;
                for ( int k = ((i>j)?i:j); k < size; ++k )  
                    sum += ((j==k)?1.0:matrix[j*size+k])*matrix[k*size+i];
                matrix[j*size+i] = sum;
            }
    }

    static T* gauss_newton(vector1D *Vol1, vector1D *Vol2, pointer1D *Vol2_derivatives, T P_initial[6],
                                T divide_factor, T alpha, T decrease_factor, int max_iter, bool mask){
        static T P[6];
        int shape = cbrt(Vol1->size());
        T rad = float(shape)/2.0;

        //find center of the volume
        T ox = (float(shape)-1)/2;
        T oy = (float(shape)-1)/2;
        T oz = (float(shape)-1)/2;  

        T P_old[6];
        T P_cur[6];
        T P_new[6];

        copy(P_initial, P_initial+6, P_old);
        copy(P_initial, P_initial+6, P_cur);
        copy(P_initial, P_initial+6, P_new);

        // cout << "P old start " ; for (int i = 0; i < 6; ++i) cout << P_old[i] << " "; cout << endl;
        // cout << "P cur start " ; for (int i = 0; i < 6; ++i) cout << P_cur[i] << " "; cout << endl;
        // cout << "P new start " ; for (int i = 0; i < 6; ++i) cout << P_new[i] << " "; cout << endl;

        vector1D errors;
        errors.push_back(10);

        pointer1D Vol2_grad_P(shape*shape*shape);
        pointer1D Vol2_grad_P_T(6);
        get_gradient_P(&Vol2_grad_P, Vol2_derivatives, divide_factor, mask);
        transpose(&Vol2_grad_P_T, &Vol2_grad_P);

        // local variables for the loop
        vector1D dest(shape*shape*shape);
        vector1D flatR(shape*shape*shape);
        vector1D Jr_rp(6);
        T JrTJr[6][6];
        T iJrTJr[36];
        T tmp_res[6];

        int counter = 0;
        T dest_coords[3];
        T error = 0;
        while(counter < max_iter){
            cout << counter << " ";

            copy(P_cur, P_cur+6, P_old);
            copy(P_new, P_new+6, P_cur);

            if(mask){
                //compute mask weights
                vector1D mask_weights(shape*shape*shape);
                T rad = T(shape)/2.0;
                get_mast_weights(&mask_weights, rad, 0.4); 

                for(int i = 0; i < shape; ++i){
                    for(int j = 0; j < shape; ++j){
                        for(int k = 0; k < shape; ++k){
                            rotate_coords_transformation(dest_coords, i, j, k, P_cur, ox, oy, oz, divide_factor, shape);
                            int idx = shape*shape*j+shape*i+k;
                            dest[idx] = interp::tricubic_interp(Vol2_derivatives,dest_coords[0],dest_coords[1],dest_coords[2])*mask_weights[idx];
                        }
                    }
                }
            }
            else{
                for(int i = 0; i < shape; ++i){
                    for(int j = 0; j < shape; ++j){
                        for(int k = 0; k < shape; ++k){
                            rotate_coords_transformation(dest_coords, i, j, k, P_cur, ox, oy, oz, divide_factor, shape);
                            dest[shape*shape*j+shape*i+k] = interp::tricubic_interp(Vol2_derivatives,dest_coords[0],dest_coords[1],dest_coords[2]);
                        }
                    }
                }
            }

            // compute the residuals and ssd
            error = 0;
            for (int i = 0; i < flatR.size(); ++i){
                flatR[i] = Vol1->at(i) - dest[i];
                error += flatR[i] * flatR[i];
            }

            //if error increases, step back and decrease alpha
            if(error > errors.back()){
                alpha *= decrease_factor;
                copy(P_cur, P_cur+6, P_old);
            }
            else{
                errors.push_back(error);
                dotn(&Jr_rp, &Vol2_grad_P_T, &flatR);
                dot_matrix(JrTJr, &Vol2_grad_P_T, &Vol2_grad_P);
                inverse6(JrTJr,iJrTJr);
                dot(tmp_res, iJrTJr, &Jr_rp);

                for (int i = 0; i < 6; ++i){
                    P_new[i] = P_cur[i] - alpha * tmp_res[i];
                }
                //cout << "P new after computation " ; for (int i = 0; i < 6; ++i) cout << P_new[i] << " "; cout << endl;
                int count = 0;
                for (int i = 0; i < 6; ++i){
                    if (abs(P_new[i] - P_cur[i]) <  1e-5){
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

        copy(P_old, P_old+6, P);
        return P;
    }
};


#endif