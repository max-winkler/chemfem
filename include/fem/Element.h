#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "linalg/Coordinate.h"
#include "linalg/Matrix2D.h"
#include "linalg/Vector2D.h"

using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector2D;

namespace chemfem{
  namespace fem{

    /// Finite element type.
    enum FEType {Lagrange, CrouzeixRaviart};

    /**
     * This class represents a single finite element. This is a virtual class and
     * one should use some child class.
     */
    class Element{

    public:

      /**
       * Initialize by type and degree.
       */
      Element(FEType, int);

      /**
       * Returns the number of local degrees of freedom.
       */
      int NrDof() const;

      /**
       * Number of components of the basis functions, 1 for a scalar element. A basis
       * function of an element with several components has exactly one component that does
       * not vanish, and Value, Gradient and Hessian describe that one.
       */
      int NrComponents() const;

      /// The component a local basis function belongs to
      int Component(int k) const { return k % nr_components; }

      /// Index of a local basis function within the scalar element behind it
      int ScalarIndex(int k) const { return k / nr_components; }

      /// Number of DOFs on each vertex of the cell
      int DofsPerVertex() const;

      /// Number of DOFs on each edge of the cell, without the vertices
      int DofsPerEdge() const;

      /// Number of DOFs in the interior of the cell
      int DofsInterior() const;

      /**
       * Returns the type of the finite element.
       */
      FEType Type() const;

      /**
       * Returns the degree of the finite element.
       */
      int Degree() const;

      /**
       * Return the function value of the trial functions
       */
      virtual double Value(int, double, double) const = 0;

      /**
       * Returns the gradient of the trial function.
       */
      virtual Vector2D Gradient(int, double, double) const = 0;

      /**
       * Returns the Hessian of the trial function on the reference element. It is
       * symmetric, so only one of the two off-diagonal entries carries information.
       */
      virtual Matrix2D Hessian(int, double, double) const = 0;

      /**
       * Reference coordinates of the point whose function value the local DOF is. Used
       * for the interpolation.
       */
      virtual chemfem::linalg::Coordinate NodalPoint(int) const = 0;

    protected:
      FEType type;
      int nr_dof;
      int degree;
      int dofs_per_vertex, dofs_per_edge, dofs_interior;

      /**
       * The components are interleaved: the local basis function k belongs to the component
       * k % nr_components of the scalar basis function k / nr_components. The blocks of the
       * vertices, edges and the interior are thereby preserved, each of them nr_components
       * times as long.
       */
      int nr_components;
    };

    /**
     * Reference coordinates of the point at parameter s in [0,1] on the local edge k,
     * which runs from vertex k to vertex k+1 of the reference triangle
     */
    inline void EdgeToRefCoords(int edge, double s, double& xi, double& eta)
    {
      switch(edge)
        {
        case 0: xi = s;      eta = 0.;      break;
        case 1: xi = 1.-s;   eta = s;       break;
        default: xi = 0.;    eta = 1.-s;    break;
        }
    }

  };
};

#endif
