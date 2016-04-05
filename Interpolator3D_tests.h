#ifndef Interpolator3D_tests_h
#define Interpolator3D_tests_h

template <typename InterpolatorT>
class InterpolatorTests {
  public:
    typedef typename InterpolatorT::VolumeT VolumeT;
    typedef typename InterpolatorT::CoordT CoordT;
    typedef typename InterpolatorT::T T;
  
    static void identity_tests(
      const InterpolatorT *interpolator,
      const VolumeT *volume) {

      const size_t cubeSize = volume->cubeSize;

      SECTION("and all the points are returned when interpolated exactly") { 
        for(size_t z = 0; z < cubeSize; z++) {
          for(size_t y = 0; y < cubeSize; y++) {
            for(size_t x = 0; x < cubeSize; x++) {
              INFO("(z,y,x): (" << z << ", " << y << ", " << x << ")"); 
              REQUIRE(interpolator->interp(z, y, x) ==
                volume->at(z, y, x));
            }
          }
        }
      }
    }
    
    static void approx_identity_tests(
      const InterpolatorT *interpolator,
      const VolumeT *volume,
      const T tolerance) {

      const size_t cubeSize = volume->cubeSize;

      SECTION("and all the points are returned when interpolated exactly") { 
        for(size_t z = 0; z < cubeSize; z++) {
          for(size_t y = 0; y < cubeSize; y++) {
            for(size_t x = 0; x < cubeSize; x++) {
              INFO("(z,y,x): (" << z << ", " << y << ", " << x << ")"); 
              REQUIRE(interpolator->interp(z, y, x) ==
                Approx(volume->at(z, y, x)).epsilon(tolerance));
            }
          }
        }
      }
    }
   
    static void constant_tests(
      const InterpolatorT *interpolator,
      const VolumeT *volume,
      const T constValue) {

      const CoordT cubeSize = volume->cubeSize;

      const CoordT stepSize = 0.25;

      SECTION("and interpolating any point returns the constant value") { 
          for(CoordT z = 0; z < cubeSize; z += stepSize) {
              for(CoordT y = 0; y < cubeSize; y += stepSize) {
                  for(CoordT x = 0; x < cubeSize; x += stepSize) {
                      INFO("(z,y,x): (" << z << ", " << y << ", " << x << ")"); 
                      REQUIRE(constValue == interpolator->interp(z, y, x));
                  }
              }
          }
      }
    }
    
    static void approx_constant_tests(
      const InterpolatorT *interpolator,
      const VolumeT *volume,
      const T constValue,
      const T tolerance) {

      const CoordT cubeSize = volume->cubeSize;

      const CoordT stepSize = 0.25;

      SECTION("and interpolating any point returns the constant value") { 
        for(CoordT z = 0; z < cubeSize; z += stepSize) {
          for(CoordT y = 0; y < cubeSize; y += stepSize) {
            for(CoordT x = 0; x < cubeSize; x += stepSize) {
              INFO("(z,y,x): (" << z << ", " << y << ", " << x << ")"); 
              INFO("error: " << constValue - interpolator->interp(z, y, x)); 
              REQUIRE(constValue ==
                Approx(interpolator->interp(z, y, x)).epsilon(tolerance));
            }
          }
        }
      }
    }
};

#endif
