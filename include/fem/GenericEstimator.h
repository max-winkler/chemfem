#ifndef _GENERIC_ESTIMATOR_H_
#define _GENERIC_ESTIMATOR_H_

#include <functional>
#include <vector>

#include "fem/FEExpression.h"
#include "fem/FEFunction.h"

#include "linalg/Coordinate.h"
#include "linalg/Vector2D.h"

namespace chemfem{
  namespace fem{

    /// State of the discrete solution at a point
    struct SolutionState
    {
      /// Value of the discrete solution
      double value;
      /// Gradient of the discrete solution, in physical coordinates
      chemfem::linalg::Vector2D gradient;
      /**
       * Laplacian of the discrete solution. Exact (and zero) for P1. For higher
       * order elements this is reported as zero and a warning is issued, because
       * Element only exposes Value() and Gradient().
       */
      double laplacian;
    };

    /// Geometry of the cell an indicator is accumulated on
    struct CellGeometry
    {
      /// Diameter of the cell, h_T
      double h;
      /// Area of the cell
      double area;
      /// Index of the cell in the mesh
      size_t index;
    };

    /// Geometry of one edge of that cell
    struct EdgeGeometry
    {
      /// Length of the edge, h_E
      double h;
      /// Unit normal, pointing out of the cell the indicator is accumulated on
      chemfem::linalg::Vector2D normal;
      /// Local index of the edge within that cell, 0..2
      int local_index;
      /// True if the edge has no neighbor
      bool boundary;
    };

    /// Jump of the function value across an edge
    double Jump(const SolutionState& u, const SolutionState& u_out);

    /// Jump of the normal derivative across an edge, [du/dn]
    double NormalJump(const SolutionState& u, const SolutionState& u_out,
                      const EdgeGeometry& edge);

    /**
     * Integrand of a volume term, evaluated at a quadrature point inside a cell.
     */
    typedef std::function<double(const chemfem::linalg::Coordinate& pos,
                                 const CellGeometry& cell,
                                 const SolutionState& u)> VolumeIntegrand;

    /**
     * Integrand of an edge term, evaluated at a quadrature point on one edge of a
     * cell. Everything is seen from the cell the indicator is accumulated on, so an
     * interior edge is visited once from either side. On a boundary edge u_out
     * repeats u.
     */
    typedef std::function<double(const chemfem::linalg::Coordinate& pos,
                                 const CellGeometry& cell,
                                 const EdgeGeometry& edge,
                                 const SolutionState& u,
                                 const SolutionState& u_out)> EdgeIntegrand;

    /// Which edges an edge term is evaluated on
    enum EdgeSelection {INTERIOR_EDGES, BOUNDARY_EDGES, ALL_EDGES};

    /**
     * Cell-wise error estimator assembled from user-defined integrands.
     *
     * The estimator stores only the terms, not the solution, so it is defined once
     * and evaluated for every solution of an adaptive loop. For the Poisson problem
     * the terms are predefined:
     *
     * \code
     *   GenericEstimator eta = GenericEstimator::Residual(f);
     *   Vector Indicators = eta.Assemble(Solution);
     * \endcode
     *
     * Any other estimator is composed from user-defined integrands:
     *
     * \code
     *   struct WeightedJump
     *   {
     *     double operator()(const Coordinate& pos, const CellGeometry&,
     *                       const EdgeGeometry& edge,
     *                       const SolutionState& u, const SolutionState& u_out) const
     *     {
     *       double j = a(pos) * NormalJump(u, u_out, edge);
     *       return 0.5 * edge.h / a(pos) * j*j;
     *     }
     *   };
     *
     *   GenericEstimator eta;
     *   eta.AddEdgeTerm(WeightedJump());
     * \endcode
     *
     * Every term contributes its integral over the cell resp. over one edge of the
     * cell, so the integrands return the squared residual where that is wanted. The
     * result is the sum of all terms per cell, i.e. eta_T^2 for the example.
     */
    class GenericEstimator
    {
    public:

      /// Volume residual of the Poisson problem, \f$h_T^2 |f + \Delta u_h|^2\f$
      struct VolumeResidual
      {
        ScalarFunction f;

        explicit VolumeResidual(ScalarFunction f) : f(f) {}

        double operator()(const chemfem::linalg::Coordinate& pos,
                          const CellGeometry& cell,
                          const SolutionState& u) const;
      };

      /**
       * Jump of the normal derivative across an edge,
       * \f$\tfrac12 h_E |[\partial u_h/\partial n]|^2\f$. The factor one half
       * distributes the edge over the two cells sharing it.
       */
      struct EdgeJump
      {
        double operator()(const chemfem::linalg::Coordinate& pos,
                          const CellGeometry& cell,
                          const EdgeGeometry& edge,
                          const SolutionState& u,
                          const SolutionState& u_out) const;
      };

      /**
       * The standard residual based a posteriori error estimator for
       * \f$-\Delta u = f\f$,
       *
       * \f[ \eta_T^2 = h_T^2 \|f + \Delta u_h\|_{L^2(T)}^2
       *              + \tfrac12 \sum_E h_E \|[\partial u_h/\partial n]\|_{L^2(E)}^2 \f]
       *
       * which bounds the error in the H1 seminorm. It is reliable and efficient, so
       * the ratio of estimator to true error settles at a constant.
       */
      static GenericEstimator Residual(ScalarFunction f);

      /// Adds \f$\int_T g\,dx\f$ to the indicator of every cell T
      void AddVolumeTerm(VolumeIntegrand);

      /// Adds \f$\int_E g\,ds\f$ to the indicator of the cell(s) adjacent to E
      void AddEdgeTerm(EdgeIntegrand, EdgeSelection = INTERIOR_EDGES);

      /**
       * Evaluates all terms for the given discrete solution. The returned vector
       * holds one value per cell, in the cell ordering of the mesh the solution
       * lives on.
       */
      chemfem::linalg::Vector Assemble(const FEFunction&) const;

    private:
      std::vector<VolumeIntegrand> VolumeTerms;
      std::vector<EdgeIntegrand> EdgeTerms;
      std::vector<EdgeSelection> EdgeTermSelection;
    };

  };
};

#endif
