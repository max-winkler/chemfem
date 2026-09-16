#include "fem/DofTransform.h"

#include <iostream>

using chemfem::linalg::DenseMatrix;
using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector;
using chemfem::linalg::Vector2D;

namespace chemfem{
  namespace fem{

    namespace {

      /**
       * The vector of the local edge e, pointing away from the vertex v it is taken at. The
       * local edge e runs from vertex e to vertex e+1, and the columns of the Jacobian are
       * the edges from vertex 0.
       */
      Vector2D EdgeVector(int e, int v, const Matrix2D& Jac)
      {
        Vector2D t;

        switch(e)
          {
          case 0:  t = Vector2D(Jac.a00, Jac.a10); break;                        // P1 - P0
          case 1:  t = Vector2D(Jac.a01 - Jac.a00, Jac.a11 - Jac.a10); break;    // P2 - P1
          default: t = Vector2D(-Jac.a01, -Jac.a11);                             // P0 - P2
          }

        // The edge runs from vertex e to vertex e+1, so it points away from the first and
        // towards the second
        if(v == (e+1)%3)
          t = Vector2D(-t.x, -t.y);

        return t;
      }

      /**
       * The two derivative DOFs of the vertex v and the rows of C that belong to them.
       * Returns false if the vertex has no such pair.
       *
       * The global derivative DOFs of a node measure the derivatives along its frame, so the
       * row of a local DOF holds the frame components of the edge vector it is taken along.
       * For the Cartesian frame of an interior node these are the components themselves.
       */
      bool VertexBlock(const FESpace& Space, size_t cell, int v, const Matrix2D& Jac,
                       int slot[2], Vector2D crow[2])
      {
        const Element& E = Space.RefElement();

        Vector2D dir[2];
        int found = 0;

        for(int k=0; k<E.NrDof() && found<2; ++k)
          {
            const DofDescriptor d = E.Dof(k);

            if(d.vertex != v)
              continue;

            if(d.type != DofType::EdgeDirectionalDerivative)
              continue;

            slot[found] = k;
            dir[found] = EdgeVector(d.edge, v, Jac);
            ++found;
          }

        if(found != 2)
          return false;

        const NodeFrame frame = Space.VertexFrame(cell, v);

        for(int i=0; i<2; ++i)
          crow[i] = Vector2D(dot(dir[i], frame.r1), dot(dir[i], frame.r2));

        return true;
      }

      /**
       * Whether the element needs the transformation. Reports the DOF types that are not
       * implemented yet only if asked to, so the places that run per cell stay quiet.
       */
      bool TransformNeeded(const Element& E, bool report)
      {
        bool directional = false, moment = false;

        for(int k=0; k<E.NrDof(); ++k)
          switch(E.Dof(k).type)
            {
            case DofType::PointValue:
              break;

            case DofType::EdgeDirectionalDerivative:
              directional = true;
              break;

            case DofType::EdgeNormalDerivative:
              if(report)
                std::cerr << "Error: DOFs that measure a normal derivative are not "
                          << "implemented yet, the assembled system would be wrong.\n";
              return false;

            case DofType::EdgeMoment:
              moment = true;
              break;
            }

        // An edge moment alone needs nothing here, it is oriented by the sign that
        // FESpace::LocalSign applies to the basis function itself. Next to a directional
        // derivative the two would have to be combined into one transformation.
        if(directional && moment)
          {
            if(report)
              std::cerr << "Error: Elements with both directional derivative and edge "
                        << "moment DOFs are not implemented yet, the assembled system "
                        << "would be wrong.\n";
            return false;
          }

        return directional;
      }
    }

    bool NeedsDofTransform(const FESpace& Space)
    {
      return TransformNeeded(Space.RefElement(), true);
    }

    void TransformLocalRows(DenseMatrix& LocalMatrix, int NrColumns, const FESpace& Space,
                            size_t cell, const Matrix2D& Jac)
    {
      for(int v=0; v<3; ++v)
        {
          int slot[2];
          Vector2D crow[2];

          if(!VertexBlock(Space, cell, v, Jac, slot, crow))
            continue;

          const int p = slot[0], q = slot[1];

          for(int l=0; l<NrColumns; ++l)
            {
              const double ap = LocalMatrix[p][l], aq = LocalMatrix[q][l];

              LocalMatrix[p][l] = crow[0].x*ap + crow[1].x*aq;
              LocalMatrix[q][l] = crow[0].y*ap + crow[1].y*aq;
            }
        }
    }

    void TransformLocalColumns(DenseMatrix& LocalMatrix, int NrRows, const FESpace& Space,
                               size_t cell, const Matrix2D& Jac)
    {
      for(int v=0; v<3; ++v)
        {
          int slot[2];
          Vector2D crow[2];

          if(!VertexBlock(Space, cell, v, Jac, slot, crow))
            continue;

          const int p = slot[0], q = slot[1];

          for(int k=0; k<NrRows; ++k)
            {
              const double ap = LocalMatrix[k][p], aq = LocalMatrix[k][q];

              LocalMatrix[k][p] = crow[0].x*ap + crow[1].x*aq;
              LocalMatrix[k][q] = crow[0].y*ap + crow[1].y*aq;
            }
        }
    }

    void TransformLocalVector(Vector& LocalVector, const FESpace& Space, size_t cell,
                              const Matrix2D& Jac)
    {
      for(int v=0; v<3; ++v)
        {
          int slot[2];
          Vector2D crow[2];

          if(!VertexBlock(Space, cell, v, Jac, slot, crow))
            continue;

          const int p = slot[0], q = slot[1];

          const double bp = LocalVector[p], bq = LocalVector[q];

          LocalVector[p] = crow[0].x*bp + crow[1].x*bq;
          LocalVector[q] = crow[0].y*bp + crow[1].y*bq;
        }
    }

    void GatherLocalCoefficients(const FESpace& Space, const Vector& Data, size_t cell,
                                 std::vector<double>& Local)
    {
      const size_t n = Space.NrLocalDof();

      if(Local.size() != n)
        Local.resize(n);

      for(size_t k=0; k<n; ++k)
        Local[k] = Data[Space.GetGlobalIndex(cell, k)];

      const Element& E = Space.RefElement();

      if(!TransformNeeded(E, false))
        return;

      const Matrix2D Jac = Space.GetMesh().Jacobian(cell);

      for(int v=0; v<3; ++v)
        {
          int slot[2];
          Vector2D crow[2];

          if(!VertexBlock(Space, cell, v, Jac, slot, crow))
            continue;

          const int p = slot[0], q = slot[1];
          const double cp = Local[p], cq = Local[q];

          // C, not its transpose: the local DOF p measures the derivative along dir[0],
          // which is that combination of the two global derivatives at the vertex
          Local[p] = crow[0].x*cp + crow[0].y*cq;
          Local[q] = crow[1].x*cp + crow[1].y*cq;
        }
    }

  }
}
