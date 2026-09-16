#include "fem/ProductElement.h"

#include <iostream>

using chemfem::linalg::Coordinate;
using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector2D;

namespace chemfem{
  namespace fem{

    ProductElement::ProductElement(const ScalarElement& Scalar, int components)
      : ScalarElement(Scalar.Type(), Scalar.Degree()), Scalar(Scalar)
    {
      if(components < 1)
	std::cerr << "Error: A product element needs at least one component.\n";

      nr_components = components;

      nr_dof = Scalar.NrDof()*components;
      dofs_per_vertex = Scalar.DofsPerVertex()*components;
      dofs_per_edge = Scalar.DofsPerEdge()*components;
      dofs_interior = Scalar.DofsInterior()*components;
    }

    double ProductElement::Value(int k, double xi, double eta) const
    {
      return Scalar.Value(ScalarIndex(k), xi, eta);
    }

    Vector2D ProductElement::Gradient(int k, double xi, double eta) const
    {
      return Scalar.Gradient(ScalarIndex(k), xi, eta);
    }

    Matrix2D ProductElement::Hessian(int k, double xi, double eta) const
    {
      return Scalar.Hessian(ScalarIndex(k), xi, eta);
    }

    Coordinate ProductElement::NodalPoint(int k) const
    {
      return Scalar.NodalPoint(ScalarIndex(k));
    }

    DofDescriptor ProductElement::Dof(int k) const
    {
      // The kind and the entity are the ones of the scalar element, only the running index
      // within the entity is multiplied by the components
      DofDescriptor d = Scalar.Dof(ScalarIndex(k));

      d.index = d.index*nr_components + Component(k);

      return d;
    }

    const ScalarElement& ProductElement::ScalarPart() const
    {
      return Scalar;
    }

  }
}
