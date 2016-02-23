#ifndef Utils_h
#define Utils_h

#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace Eigen;
using namespace std;

static const float PI = 3.14159265358979323846264;

class Utils{
  public:
	static float to_radian(float angle){
        return angle*PI/180.0;
    }


    static Matrix3f get_rotation_matrix(float theta, float wx, float wy, float wz){
        theta = to_radian(theta);

        Matrix3f m;

        float tmp = wx*wx + wy*wy + wz*wz;
        if ((tmp != 1) && (tmp != 0)){
            float norm = sqrt(tmp);
            wx = wx/norm;
            wy = wy/norm;
            wz = wz/norm;
        }

        m = AngleAxisf(theta, Vector3f(wx,wy,wz));
        return m;
    }
    
};



#endif