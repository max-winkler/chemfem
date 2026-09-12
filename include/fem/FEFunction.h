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

      /**
       * Function value in a quadrature point. Considerably cheaper than Evaluate, which also
       * computes the gradient and the Hessian.
       */
      double Value(const QuadPoint&) const;

      /// Integral of the function over the domain, divided by the area of the domain
      double Mean() const;

      /**
       * Subtracts the mean value, so that the integral of the function vanishes. Used for a
       * pressure that is only determined up to a constant. Requires an element that
       * represents a constant by equal coefficients, as Lagrange and Crouzeix-Raviart do.
       */
      void SubtractMean();

      void WriteVtk(const std::string&) const;
      
    private:
      const FESpace* Space;
      Vector Data;

      /// What the cache below holds: nothing, the function value only, or all of PointValues
      enum CacheContent {NOTHING, VALUE_ONLY, EVERYTHING};

      /// The last point that was evaluated and its result. An integrand is called once per
      /// basis function in the same point, so this saves the repeated work.
      mutable QuadPoint CachedPoint;
      mutable PointValues CachedValues;
      mutable CacheContent CacheValid = NOTHING;
      
    };
  };
};

#endif
