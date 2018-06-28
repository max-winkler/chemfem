#ifndef _FE_FUNCTION_H_
#define _FE_FUNCTION_H_

#include <iostream>

#include "linalg/Vector.h"
#include "fem/FESpace.h"

namespace chemfem{
  namespace fem{

    /**
     * Each instance of this class represents a function from some finite element space.
     */
    class FEFunction{
    public:
      /**
       * Copy constructor. Usually called when an FE function is returned 
       * from another function.
       */
      FEFunction(const FEFunction&);

      /**
       * Creates the zero-function in the given finite element space.
       */
      FEFunction(const FESpace&);

      /**
       * Returns a reference to the finite element space.
       */
      const FESpace& GetFESpace() const;

      /**
       * Returns a constant reference to the value of the i-th degree of freedom.
       */
      const double& operator[](size_t) const;
      
    private:
      const FESpace* Space;
      Vector Data;
      
    };
  };
};

#endif
