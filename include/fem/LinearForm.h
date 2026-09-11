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
     * Terms beyond the predefined ones are given as functors, which get the values of
     * the test function v in a quadrature point, see also BilinearForm.
     */
    class LinearForm
    {
    public:

      /// Integrand of a volume term for the test function v
      typedef std::function<double(const QuadPoint& p, const CellGeometry& cell,
                                   const PointValues& v)> VolumeIntegrand;

      /// Integrand of a boundary term, evaluated on a boundary edge of the cell
      typedef std::function<double(const QuadPoint& p, const CellGeometry& cell,
                                   const EdgeGeometry& edge,
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

      /// Adds the integral of the functor over all cells
      void AddVolumeTerm(VolumeIntegrand);

      /**
       * Adds the integral of the functor over the part of the boundary where the
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
      struct BoundaryTerm
      {
        BoundaryIntegrand integrand;
        BoundaryIndicator part;
      };

      /**
       * Stores the terms added to the linear form in a vector.
       */
      std::vector<FEExpression> Terms;

      std::vector<VolumeIntegrand> VolumeTerms;
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
