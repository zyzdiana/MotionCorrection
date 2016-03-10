#ifndef Interpolator3D_tests_h
#define Interpolator3D_tests_h

template <typename InterpolatorT>
class InterpolatorTests {
  public:
    typedef typename InterpolatorT::VolumeT VolumeT;
    typedef typename InterpolatorT::T T;
  
    static void tests(
      const InterpolatorT *interpolator,
      const VolumeT *volume) {

      const size_t cubeSize = volume->cubeSize;

      SECTION("and all the points are returned when interpolated exactly") { 
          for(size_t z = 0; z < cubeSize; z++) {
              for(size_t y = 0; y < cubeSize; y++) {
                  for(size_t x = 0; x < cubeSize; x++) {
                      REQUIRE(interpolator->interp(z, y, x) ==
                        volume->at(z, y, x));
                  }
              }
          }
      }
    
    }
};

#endif
