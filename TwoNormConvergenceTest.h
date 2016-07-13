#ifndef TwoNormConvergenceTest_h
#define TwoNormConvergenceTest_h


template< typename T >
class TwoNormConvergenceTest {
  public:
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    TwoNormConvergenceTest(const T paramUpdate2NormLimit) :
      paramUpdate2NormLimit(paramUpdate2NormLimit) {}

    bool operator()(ParamT *paramUpdate) const {
      return (paramUpdate->norm() < paramUpdate2NormLimit);
    }

  protected:
    T paramUpdate2NormLimit;
};

#endif
