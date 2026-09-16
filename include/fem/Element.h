#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "linalg/Coordinate.h"
#include "linalg/DenseMatrix.h"
#include "linalg/Matrix2D.h"
#include "linalg/Vector.h"
#include "linalg/Vector2D.h"

using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector2D;

namespace chemfem{
  namespace fem{

    /// Finite element type.
    enum FEType {Lagrange, CrouzeixRaviart, DG, RaviartThomas};

    /**
     * What the coefficient of a degree of freedom measures. Together with the convention for
     * the global DOFs this determines how a local DOF is expressed in terms of the global
     * ones, so an element describes its DOFs and does not implement that transformation
     * itself. A scoped enum because there is already a struct PointValues next door.
     */
    enum class DofType
    {
      /// The function value in NodalPoint(k)
      PointValue,
      /// The derivative along one edge, taken at one of its endpoints, as Hermite has
      EdgeDirectionalDerivative,
      /// The derivative normal to one edge, as Argyris and Morley have
      EdgeNormalDerivative,
      /// An integral over one edge, such as the flux of a Raviart-Thomas element
      EdgeMoment
    };

    /**
     * What a local DOF measures and which entities of the cell it belongs to. A directional
     * derivative needs both: the vertex fixes which global derivative pair it contributes to,
     * the edge fixes the direction. Entities that do not apply are -1, and index counts the
     * DOFs of the same kind on the same entity.
     */
    struct DofDescriptor
    {
      DofType type;
      int vertex;
      int edge;
      int index;
    };


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
       * Reference coordinates of the point the local DOF belongs to: a vertex, an edge
       * midpoint or the barycenter. For a DOF that is not a point value, such as the flux
       * of a Raviart-Thomas element, it is the point its edge is represented by.
       */
      virtual chemfem::linalg::Coordinate NodalPoint(int) const = 0;

      /**
       * What the local DOF k measures and where it sits. The default reports a point value on
       * the entity that the numbering of DofManager puts it on, which is correct for every
       * element whose coefficients are function values. An element with other functionals
       * overrides it.
       */
      virtual DofDescriptor Dof(int k) const
      {
        const int nv = 3*dofs_per_vertex;
        const int ne = 3*dofs_per_edge;

        if(k < nv)
          return DofDescriptor{DofType::PointValue, k/dofs_per_vertex, -1,
                               k%dofs_per_vertex};

        if(k < nv + ne)
          return DofDescriptor{DofType::PointValue, -1, (k - nv)/dofs_per_edge,
                               (k - nv)%dofs_per_edge};

        return DofDescriptor{DofType::PointValue, -1, -1, k - nv - ne};
      }

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
     * An element whose shape functions are scalar. A basis function of a vector valued
     * space built from such an element, as ProductElement does, is that scalar function
     * placed in the component Component(k), so the accessors below describe it completely.
     */
    class ScalarElement : public Element
    {
    public:
      ScalarElement(FEType type, int degree) : Element(type, degree) {}

      /// Value of the shape function k in the reference point
      virtual double Value(int, double, double) const = 0;

      /// Its gradient on the reference element
      virtual Vector2D Gradient(int, double, double) const = 0;

      /**
       * Its Hessian on the reference element. It is symmetric, so only one of the two
       * off-diagonal entries carries information.
       */
      virtual Matrix2D Hessian(int, double, double) const = 0;
    };

    /**
     * An element whose shape functions are vector fields already on the reference element,
     * such as Raviart-Thomas. They are not built from a scalar function, and a cell maps
     * them with a Piola transform rather than by leaving the values unchanged.
     */
    class VectorElement : public Element
    {
    public:
      VectorElement(FEType type, int degree) : Element(type, degree) {}

      /// Value of the shape function k in the reference point
      virtual Vector2D Value(int, double, double) const = 0;

      /// Its gradient on the reference element, row i holds the gradient of component i
      virtual Matrix2D Gradient(int, double, double) const = 0;

      /// Trace of the gradient
      double Divergence(int k, double xi, double eta) const
      {
        return Gradient(k, xi, eta).Trace();
      }
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
