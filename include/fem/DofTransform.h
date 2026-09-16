#ifndef _DOF_TRANSFORM_H_
#define _DOF_TRANSFORM_H_

#include "fem/FESpace.h"

#include "linalg/DenseMatrix.h"
#include "linalg/Matrix2D.h"
#include "linalg/Vector.h"

namespace chemfem{
  namespace fem{

    /**
     * The local DOFs of a cell and the global ones are related by sigma_local = C sigma_global,
     * where the row i of C holds the coefficients of the local DOF i in terms of the global
     * ones. C follows from what each DOF measures, which Element::Dof reports, together with
     * the convention for the global DOFs, so no element implements it itself.
     *
     * A DOF that is a point value contributes the unit vector, and C is the identity for every
     * element built from those. A derivative along an edge, as a Hermite element has, becomes a
     * combination of the two global derivatives at its vertex, with the edge vector as the
     * coefficients:
     *
     *   d v / d (v_m - v_n) = (v_m - v_n)_1 dv/dx + (v_m - v_n)_2 dv/dy
     *
     * The local matrix and load vector of a cell then have to be transformed before they are
     * inserted, A <- C_test^T A C_trial and b <- C_test^T b, and the coefficients when the
     * function is evaluated, c_local <- C c_global.
     */
    /**
     * Whether the space needs the transformation below, which the caller has to ask once per
     * assembly. It also reports the DOF types that are not implemented yet and returns false
     * for them, so an unusable element is named on the console instead of quietly producing a
     * wrong system.
     */
    bool NeedsDofTransform(const FESpace&);

    /// The test side of a local matrix, A <- C^T A. Does nothing without a transformation.
    void TransformLocalRows(chemfem::linalg::DenseMatrix&, int NrColumns, const FESpace&,
                            size_t cell, const chemfem::linalg::Matrix2D& Jac);

    /// The trial side of a local matrix, A <- A C
    void TransformLocalColumns(chemfem::linalg::DenseMatrix&, int NrRows, const FESpace&,
                               size_t cell, const chemfem::linalg::Matrix2D& Jac);

    /// A local load vector, b <- C^T b
    void TransformLocalVector(chemfem::linalg::Vector&, const FESpace&, size_t cell,
                              const chemfem::linalg::Matrix2D& Jac);

  };
};

#endif
