#include "catch.hpp"

#include "ReadFile.h"
#include "BinaryFile.h"
#include "Interpolation_3D.h"

#include <vector>
#include <complex>

TEST_CASE("rotation of zero degrees") {

    typedef Interpolation_3D<float> interp;
    typedef interp::vector4D vector4D;
    typedef interp::vector3D vector3D;
    typedef interp::vector2D vector2D;
    typedef interp::vector1D vector1D;
    typedef std::complex<float> dataT;

    string path = "/Users/zyzdiana/Dropbox/THESIS/Oct_13_navs/8mm_iso_x_rot_3_0_to_5_0_deg_z_trans_rep_0";
    int num_slice = 32;
    vector3D vol1 = ReadFile<float>::read_volume(path, num_slice);

    float theta = 0.0;
    float wx = 0;
    float wy = 0;
    float wz = 0;
    float ox,oy,oz = float(num_slice-1)/2;
    vector1D point(3,15);
    vector1D rotated_point = interp::rotate_coords_3d(point[0],point[1],point[2], theta, wx, wy, wz, ox, oy, oz);
    REQUIRE(point) ==  rotated_point;
}
