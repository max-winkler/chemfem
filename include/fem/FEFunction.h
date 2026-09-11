#ifndef _FE_FUNCTION_H_
#define _FE_FUNCTION_H_

#include <iostream>

#include "linalg/Vector.h"
#include "fem/FESpace.h"
#include "fem/PointValues.h"

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

      FEFunction& operator=(const FEFunction&) = default;

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

      /**
       * Build an FE function from the vector of the values of the free degrees of freedom.
       * The degrees of freedom at Dirichlet nodes are taken from the FE space.
       */
      void CreateFunction(const Vector&);

      /**
       * Initializes the FE function by the coefficients of the degrees of freedom.
       * Essential boundary conditions are ignored.
       */
      void SetCoefficients(const Vector&);

      /**
       * Value, gradient and Hessian in a quadrature point, e.g. of the solution of the
       * previous time step inside the integrand of a form
       */
      PointValues Evaluate(const QuadPoint&) const;

      void WriteVtk(const std::string&) const;
      
    private:
      const FESpace* Space;
      Vector Data;

      /// The last point Evaluate() was called for and its result. An integrand is called
      /// once per basis function in the same point, so this saves the repeated work.
      mutable QuadPoint CachedPoint;
      mutable PointValues CachedValues;
      mutable bool CacheValid = false;
      
    };
  };
};

#endif
