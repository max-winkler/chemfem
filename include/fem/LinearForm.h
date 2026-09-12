#ifndef _LINEAR_FORM_
#define _LINEAR_FORM_

#include <functional>
#include <vector>

#include "linalg/Vector.h"
#include "fem/FESpace.h"
#include "fem/FEExpression.h"
#include "fem/PointValues.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace fem{

    /**
     * This class represents a linear form which stores and assembles the vector for the right-hand
     * side or Neumann boundary conditions.
     *
     * Terms beyond the predefined ones are given as an integrand, which gets the values of
     * the test function v in a quadrature point, see also BilinearForm.
     */
    class LinearForm
    {
    public:

      /// Integrand for the test function v
      typedef std::function<double(const PointValues& v)> Integrand;

      /// Integrand of a volume term that depends on the quadrature point as well
      typedef std::function<double(const QuadPoint& p, const PointValues& v)> PointIntegrand;

      /// Integrand of a boundary term that depends on the quadrature point or the edge
      typedef std::function<double(const QuadPoint& p, const EdgeGeometry& edge,
                                   const PointValues& v)> BoundaryIntegrand;

      /**
       * Constructor which associates the linear form with a function space.
       */
      LinearForm(const FESpace&);

      /**
       * Adds a volume force. This is a function handle to the function defining the right-hand
       * side of the partial differential equation.
       */
      void AddVolumeForce(ScalarFunction);

      /**
       * Adds a Neumann boundary condition. Requires a function handle to the function definiting the
       * boundary condition.
       */
      void AddNeumannBC(ScalarFunction);

      /// Adds the integral of the integrand over all cells
      void AddVolumeTerm(Integrand);

      /// Adds the integral of the integrand over all cells
      void AddVolumeTerm(PointIntegrand);

      /**
       * Adds the integral of the integrand over the part of the boundary where the
       * indicator is true, over the whole boundary if it is omitted
       */
      void AddBoundaryTerm(Integrand, BoundaryIndicator = nullptr);

      /**
       * Adds the integral of the integrand over the part of the boundary where the
       * indicator is true, over the whole boundary if it is omitted
       */
      void AddBoundaryTerm(BoundaryIntegrand, BoundaryIndicator = nullptr);

      /**
       * Assembles the load vector. Before calling this routine all terms that are required should be
       * added to the linear form.
       */
      void Assemble();

      /**
       * Returns a reference to the load vector. Before calling this routine the assemble function
       * has to be invoked.
       */
      Vector& LoadVector();

      const FESpace& GetTestSpace() const;

    private:
      /// A boundary term, given by one of the two kinds of integrands
      struct BoundaryTerm
      {
        Integrand integrand;
        BoundaryIntegrand point_integrand;
        BoundaryIndicator part;
      };

      /**
       * Stores the terms added to the linear form in a vector.
       */
      std::vector<FEExpression> Terms;

      std::vector<Integrand> VolumeTerms;
      std::vector<PointIntegrand> PointVolumeTerms;
      std::vector<BoundaryTerm> BoundaryTerms;

      /**
       * Stores a reference to the test space.
       */
      const FESpace& TestSpace;

      /**
       * The load vector created in the assemble routine.
       */
      Vector Vec;
    };

  };
};

#endif
