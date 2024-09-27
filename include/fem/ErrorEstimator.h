#ifndef _ERROR_ESTIMATOR_H_
#define _ERROR ESTIMATOR_H_

#include "fem/FEFunction.h"
#include "fem/FEExpression.h"

namespace chemfem{
  namespace fem{
    class ErrorEstimator
    {
    public:

      ErrorEstimator(const FEFunction&);

      void AddVolumeResidualTerm(ScalarFunction);
      void AddEdgeJumpTerm();

      Vector Assemble();
    private:
      
      std::vector<FEExpression> Terms;
      FEFunction u;
    };
  };
};

#endif
