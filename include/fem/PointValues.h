#ifndef _POINT_VALUES_H_
#define _POINT_VALUES_H_

#include <cstddef>
#include <vector>

#include "fem/Element.h"

#include "linalg/Coordinate.h"
#include "linalg/Matrix2D.h"
#include "linalg/Vector.h"
#include "linalg/Vector2D.h"

namespace chemfem{
  namespace fem{

    /**
     * Values of a function in a point, in physical coordinates. The error estimator gets
     * them for the discrete solution, the integrands of the forms for the trial and the
     * test function.
     */
    struct PointValues
    {
      double value;
      chemfem::linalg::Vector2D gradient;
      chemfem::linalg::Matrix2D hessian;
      /// The trace of the Hessian
      double laplacian;
    };

    /**
     * Values of a vector valued function in a point, in physical coordinates. The integrands
     * of a form on a vector valued space get them for the trial and the test function.
     */
    struct VectorValues
    {
      chemfem::linalg::Vector2D value;
      /// Row i holds the gradient of the component i
      chemfem::linalg::Matrix2D gradient;
      /// Trace of the gradient
      double divergence;
      /// The scalar rotation, du_1/dx - du_0/dy
      double curl;
    };

    /// A quadrature point, in physical coordinates and on the reference element of its cell
    struct QuadPoint
    {
      chemfem::linalg::Coordinate x;
      size_t cell;
      double xi, eta;
    };

    /// Geometry of the cell a term is integrated over
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
      /// Unit normal, pointing out of the cell
      chemfem::linalg::Vector2D normal;
      /// Local index of the edge within the cell, 0..2
      int local_index;
      /// True if the edge has no neighbor
      bool boundary;
    };

    /// Values of the basis function k of the element in (xi,eta) on the reference element
    PointValues ReferenceValues(const Element&, int k, double xi, double eta);

    /// Reference values of all basis functions in all points, stored as [q*NrDof + k]
    std::vector<PointValues> TabulateReference(const Element&, const chemfem::linalg::Vector& Xi,
                                               const chemfem::linalg::Vector& Eta);

    /**
     * Maps values on the reference element to the cell with the inverse transposed
     * Jacobian B^-T of the affine reference map: the gradient becomes B^-T g, the Hessian
     * B^-T H B^-1.
     */
    PointValues MapFromReference(const PointValues&, const chemfem::linalg::Matrix2D& InvJacT);

    /**
     * The same for a basis function of a product element: it is the given scalar function in
     * the given component and zero in the others.
     */
    VectorValues MapFromReference(const PointValues&, const chemfem::linalg::Matrix2D& InvJacT,
                                  int component);

    /**
     * Maps a vector valued reference basis function to a cell with the contravariant Piola
     * transform, s = (1/det J) J s_ref. The divergence transforms as div s = (1/det J) div
     * s_ref, which is why the flux of s through an edge equals the one of s_ref through the
     * corresponding reference edge. The gradient is not provided.
     */
    VectorValues MapFromReference(const VectorRefValues&, const chemfem::linalg::Matrix2D& Jac,
                                  double det);

  };
};

#endif
