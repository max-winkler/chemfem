#include "fem/Element.h"

namespace chemfem{
  namespace fem{

    Element::Element(FEType type, int degree)
      : type(type), degree(degree), dofs_per_vertex(0), dofs_per_edge(0), dofs_interior(0),
	nr_components(1) {}

    int Element::NrDof() const
    {
      return nr_dof;
    }

    int Element::NrComponents() const
    {
      return nr_components;
    }

    int Element::DofsPerVertex() const
    {
      return dofs_per_vertex;
    }

    int Element::DofsPerEdge() const
    {
      return dofs_per_edge;
    }

    int Element::DofsInterior() const
    {
      return dofs_interior;
    }

    FEType Element::Type() const
    {
      return type;
    }

    int Element::Degree() const
    {
      return degree;
    }
  }
}
