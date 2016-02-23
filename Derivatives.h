#ifndef Derivatives_h
#define Derivatives_h

#include "Volume.h"

template <typename VolumeT>

class Derivatives {
  public:
    typedef typename VolumeT::T T;
    int derives_shape;
    std::vector<T*> derivatives;

    Derivatives(
        const VolumeT data,
        const size_t derivSize) :
        derivsSize(derivsSize),
        data(data) {}

    derives_shape(const VolumeT *volume){
        derives_shape = volume->cubeSize + 30;
        derivatives.resize(derivSize * derivSize * derivSize);

        // allocate this big static array to store all the coefficients
        static T arrayCoeffs[derivSize][64];

        // a vector used to temporary store the polynomials of x, y, and z 
        std::vector<T> Y(64);

        for(int z = 0; z < derivSize; ++z){
            for(int y = 0; y < derivSize; ++y){
                for(int x = 0; x < derivSize; ++x){
                    int z1 = z - 15;
                    int y1 = y - 15;
                    int x1 = x - 15;
                    
                    int x0 = x1 - 1;
                    int x2 = x1 + 1;
                    int x3 = x2 + 1;
                    int y0 = y1 - 1;
                    int y2 = y1 + 1;
                    int y3 = y2 + 1;
                    int z0 = z1 - 1;
                    int z2 = z1 + 1;
                    int z3 = z2 + 1;

                    //Wrap Around
                    x0 = (x0 + shape) % shape;
                    x1 = (x1 + shape) % shape;
                    x2 = (x2 + shape) % shape;
                    x3 = (x3 + shape) % shape;

                    y0 = (y0 + shape) % shape;
                    y1 = (y1 + shape) % shape;
                    y2 = (y2 + shape) % shape;
                    y3 = (y3 + shape) % shape;

                    z0 = (z0 + shape) % shape;
                    z1 = (z1 + shape) % shape;
                    z2 = (z2 + shape) % shape;
                    z3 = (z3 + shape) % shape;                    

                    //values of f(x,y,z) at each corner.
                    Y[0] = this->volume->at(y1, x1, z1);
                    Y[1] = this->volume->at(y1, x1, z2);
                    Y[2] = this->volume->at(y1, x2, z1);
                    Y[3] = this->volume->at(y1, x2, z2);
                    Y[4] = this->volume->at(y2, x1, z1);
                    Y[5] = this->volume->at(y2, x1, z2);
                    Y[6] = this->volume->at(y2, x2, z1);
                    Y[7] = this->volume->at(y2, x2, z2);

                    //values of df/dx
                    Y[8] = ((this->volume->at(y1, x1, z2)-this->volume->at(y1, x1, z0))/2.);
                    Y[9] = ((this->volume->at(y1, x1, z3)-this->volume->at(y1, x1, z1))/2.);
                    Y[10] = ((this->volume->at(y1, x2, z2)-this->volume->at(y1, x2, z0))/2.);
                    Y[11] = ((this->volume->at(y1, x2, z3)-this->volume->at(y1, x2, z1))/2.);
                    Y[12] = ((this->volume->at(y2, x1, z2)-this->volume->at(y2, x1, z0))/2.);
                    Y[13] = ((this->volume->at(y2, x1, z3)-this->volume->at(y2, x1, z1))/2.);
                    Y[14] = ((this->volume->at(y2, x2, z2)-this->volume->at(y2, x2, z0))/2.);
                    Y[15] = ((this->volume->at(y2, x2, z3)-this->volume->at(y2, x2, z1))/2.);

                    //values of df/dy
                    Y[16] = ((this->volume->at(y1, x2, z1)-this->volume->at(y1, x0, z1))/2.);
                    Y[17] = ((this->volume->at(y1, x2, z2)-this->volume->at(y1, x0, z2))/2.);
                    Y[18] = ((this->volume->at(y1, x3, z1)-this->volume->at(y1, x1, z1))/2.);
                    Y[19] = ((this->volume->at(y1, x3, z2)-this->volume->at(y1, x1, z2))/2.);
                    Y[20] = ((this->volume->at(y2, x2, z1)-this->volume->at(y2, x0, z1))/2.);
                    Y[21] = ((this->volume->at(y2, x2, z2)-this->volume->at(y2, x0, z2))/2.);
                    Y[22] = ((this->volume->at(y2, x3, z1)-this->volume->at(y2, x1, z1))/2.);
                    Y[23] = ((this->volume->at(y2, x3, z2)-this->volume->at(y2, x1, z2))/2.);

                    //values of df/dz
                    Y[24] = ((this->volume->at(y2, x1, z1)-this->volume->at(y0, x1, z1))/2.);
                    Y[25] = ((this->volume->at(y2, x1, z2)-this->volume->at(y0, x1, z2))/2.);
                    Y[26] = ((this->volume->at(y2, x2, z1)-this->volume->at(y0, x2, z1))/2.);
                    Y[27] = ((this->volume->at(y2, x2, z2)-this->volume->at(y0, x2, z2))/2.);
                    Y[28] = ((this->volume->at(y3, x1, z1)-this->volume->at(y1, x1, z1))/2.);
                    Y[29] = ((this->volume->at(y3, x1, z2)-this->volume->at(y1, x1, z2))/2.);
                    Y[30] = ((this->volume->at(y3, x2, z1)-this->volume->at(y1, x2, z1))/2.);
                    Y[31] = ((this->volume->at(y3, x2, z2)-this->volume->at(y1, x2, z2))/2.);

                    //values of d2f/dxdy
                    Y[32] = ((this->volume->at(y1, x2, z2)-this->volume->at(y1, x0, z2)-this->volume->at(y1, x2, z0)+this->volume->at(y1, x0, z0))/4.);
                    Y[33] = ((this->volume->at(y1, x2, z3)-this->volume->at(y1, x0, z3)-this->volume->at(y1, x2, z1)+this->volume->at(y1, x0, z1))/4.);
                    Y[34] = ((this->volume->at(y1, x3, z2)-this->volume->at(y1, x1, z2)-this->volume->at(y1, x3, z0)+this->volume->at(y1, x1, z0))/4.);
                    Y[35] = ((this->volume->at(y1, x3, z3)-this->volume->at(y1, x1, z3)-this->volume->at(y1, x3, z1)+this->volume->at(y1, x1, z1))/4.);
                    Y[36] = ((this->volume->at(y2, x2, z2)-this->volume->at(y2, x0, z2)-this->volume->at(y2, x2, z0)+this->volume->at(y2, x0, z0))/4.);
                    Y[37] = ((this->volume->at(y2, x2, z3)-this->volume->at(y2, x0, z3)-this->volume->at(y2, x2, z1)+this->volume->at(y2, x0, z1))/4.);
                    Y[38] = ((this->volume->at(y2, x3, z2)-this->volume->at(y2, x1, z2)-this->volume->at(y2, x3, z0)+this->volume->at(y2, x1, z0))/4.);
                    Y[39] = ((this->volume->at(y2, x3, z3)-this->volume->at(y2, x1, z3)-this->volume->at(y2, x3, z1)+this->volume->at(y2, x1, z1))/4.);

                    //values of d2f/dxdz
                    Y[40] = ((this->volume->at(y2, x1, z2)-this->volume->at(y0, x1, z2)-this->volume->at(y2, x1, z0)+this->volume->at(y0, x1, z0))/4.);
                    Y[41] = ((this->volume->at(y2, x1, z3)-this->volume->at(y0, x1, z3)-this->volume->at(y2, x1, z1)+this->volume->at(y0, x1, z1))/4.);
                    Y[42] = ((this->volume->at(y2, x2, z2)-this->volume->at(y0, x2, z2)-this->volume->at(y2, x2, z0)+this->volume->at(y0, x2, z0))/4.);
                    Y[43] = ((this->volume->at(y2, x2, z3)-this->volume->at(y0, x2, z3)-this->volume->at(y2, x2, z1)+this->volume->at(y0, x2, z1))/4.);
                    Y[44] = ((this->volume->at(y3, x1, z2)-this->volume->at(y1, x1, z2)-this->volume->at(y3, x1, z0)+this->volume->at(y1, x1, z0))/4.);
                    Y[45] = ((this->volume->at(y3, x1, z3)-this->volume->at(y1, x1, z3)-this->volume->at(y3, x1, z1)+this->volume->at(y1, x1, z1))/4.);
                    Y[46] = ((this->volume->at(y3, x2, z2)-this->volume->at(y1, x2, z2)-this->volume->at(y3, x2, z0)+this->volume->at(y1, x2, z0))/4.);
                    Y[47] = ((this->volume->at(y3, x2, z3)-this->volume->at(y1, x2, z3)-this->volume->at(y3, x2, z1)+this->volume->at(y1, x2, z1))/4.);

                    //values of d2f/dydz
                    Y[48] = ((this->volume->at(y2, x2, z1)-this->volume->at(y2, x0, z1)-this->volume->at(y0, x2, z1)+this->volume->at(y0, x0, z1))/4.);
                    Y[49] = ((this->volume->at(y2, x2, z2)-this->volume->at(y2, x0, z2)-this->volume->at(y0, x2, z2)+this->volume->at(y0, x0, z2))/4.);
                    Y[50] = ((this->volume->at(y2, x3, z1)-this->volume->at(y2, x1, z1)-this->volume->at(y0, x3, z1)+this->volume->at(y0, x1, z1))/4.);
                    Y[51] = ((this->volume->at(y2, x3, z2)-this->volume->at(y2, x1, z2)-this->volume->at(y0, x3, z2)+this->volume->at(y0, x1, z2))/4.);
                    Y[52] = ((this->volume->at(y3, x2, z1)-this->volume->at(y3, x0, z1)-this->volume->at(y1, x2, z1)+this->volume->at(y1, x0, z1))/4.);
                    Y[53] = ((this->volume->at(y3, x2, z2)-this->volume->at(y3, x0, z2)-this->volume->at(y1, x2, z2)+this->volume->at(y1, x0, z2))/4.);
                    Y[54] = ((this->volume->at(y3, x3, z1)-this->volume->at(y3, x1, z1)-this->volume->at(y1, x3, z1)+this->volume->at(y1, x1, z1))/4.);
                    Y[55] = ((this->volume->at(y3, x3, z2)-this->volume->at(y3, x1, z2)-this->volume->at(y1, x3, z2)+this->volume->at(y1, x1, z2))/4.);

                    //values of d3f/dxdydz
                    Y[56] = ((this->volume->at(y2, x2, z2)-this->volume->at(y2, x0, z2)-this->volume->at(y2, x2, z0)+this->volume->at(y2, x0, z0))
                              -(this->volume->at(y0, x2, z2)-this->volume->at(y0, x0, z2)-this->volume->at(y0, x2, z0)+this->volume->at(y0, x0, z0)))/8.;
                    Y[57] = ((this->volume->at(y2, x2, z3)-this->volume->at(y2, x0, z3)-this->volume->at(y2, x2, z1)+this->volume->at(y2, x0, z1))
                              -(this->volume->at(y0, x2, z3)-this->volume->at(y0, x0, z3)-this->volume->at(y0, x2, z1)+this->volume->at(y0, x0, z1)))/8.;
                    Y[58] = ((this->volume->at(y2, x3, z2)-this->volume->at(y2, x1, z2)-this->volume->at(y2, x3, z0)+this->volume->at(y2, x1, z0))
                              -(this->volume->at(y0, x3, z2)-this->volume->at(y0, x1, z2)-this->volume->at(y0, x3, z0)+this->volume->at(y0, x1, z0)))/8.;
                    Y[59] = ((this->volume->at(y2, x3, z3)-this->volume->at(y2, x1, z3)-this->volume->at(y2, x3, z1)+this->volume->at(y2, x1, z1))
                              -(this->volume->at(y0, x3, z3)-this->volume->at(y0, x1, z3)-this->volume->at(y0, x3, z1)+this->volume->at(y0, x1, z1)))/8.;

                    Y[60] = ((this->volume->at(y3, x2, z2)-this->volume->at(y3, x0, z2)-this->volume->at(y3, x2, z0)+this->volume->at(y3, x0, z0))
                              -(this->volume->at(y1, x2, z2)-this->volume->at(y1, x0, z2)-this->volume->at(y1, x2, z0)+this->volume->at(y1, x0, z0)))/8.;
                    Y[61] = ((this->volume->at(y3, x2, z3)-this->volume->at(y3, x0, z3)-this->volume->at(y3, x2, z1)+this->volume->at(y3, x0, z1))
                              -(this->volume->at(y1, x2, z3)-this->volume->at(y1, x0, z3)-this->volume->at(y1, x2, z1)+this->volume->at(y1, x0, z1)))/8.;
                    Y[62] = ((this->volume->at(y3, x3, z2)-this->volume->at(y3, x1, z2)-this->volume->at(y3, x3, z0)+this->volume->at(y3, x1, z0))
                              -(this->volume->at(y1, x3, z2)-this->volume->at(y1, x1, z2)-this->volume->at(y1, x3, z0)+this->volume->at(y1, x1, z0)))/8.;
                    Y[63] = ((this->volume->at(y3, x3, z3)-this->volume->at(y3, x1, z3)-this->volume->at(y3, x3, z1)+this->volume->at(y3, x1, z1))
                              -(this->volume->at(y1, x3, z3)-this->volume->at(y1, x1, z3)-this->volume->at(y1, x3, z1)+this->volume->at(y1, x1, z1)))/8.;

                    // compute the index to the derivatives vector from loop indicies i, j and k
                    int idx = derives_shape * (derives_shape*(y+15) + (x+15)) + (z+15);

                    // store the derivatives
                    derivatives->at(idx) = X_inv_dot(&Y, arrayCoeffs[idx]);
                }
            }
        }
    }

    T at(const size_t z, const size_t y, const size_t x) const {
        return at((z * derivsSize + y) * derivsSize + x);
    }

    T at(const size_t index) const {
        return data[index]; 
    }

  public:
    const size_t derivsSize;

  protected:
    const StorageT data;
};

#endif
