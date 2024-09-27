#include <iostream>
#include <cmath>

#include "mesh/CellInfo.h"
#include "fem/ErrorEstimator.h"
#include "quadrature/QuadFormula.h"

using chemfem::linalg::DenseMatrix;
using chemfem::mesh::Node;
using chemfem::mesh::CellInfo;
using chemfem::quadrature::QuadratureFormula;
using chemfem::quadrature::QUAD_FORMULA;

namespace chemfem{
  namespace fem{

    ErrorEstimator::ErrorEstimator(const FEFunction& u) : u(u) {}
    
    void ErrorEstimator::AddVolumeResidualTerm(ScalarFunction f)
    {
      Terms.push_back(FEExpression(ExpressionType::VOLUME_RESIDUAL, f));
    }

    void ErrorEstimator::AddEdgeJumpTerm()
    {
      Terms.push_back(FEExpression(ExpressionType::EDGE_JUMP));
    }

    Vector ErrorEstimator::Assemble() {

      std::cerr << "WARNING: Implementation of error estimation is not finished yet.\n";
      std::cerr << "         Only the volume residual for P1 elements works as far.\n";
      std::cerr << "         All other terms in the error estimator will be ignored.\n";
      
      const FESpace& Space = u.GetFESpace();
      const Mesh& mesh = Space.GetMesh();

      // Store local FE Function and their derivatives
      double *Value = new double[Space.NrLocalDof()];
      
      // Initialize quadrature formulas
      QuadratureFormula QuadFormula(QUAD_FORMULA::GAUSS_7);
      
      Vector Xi, Eta, Weights;
      QuadFormula.FormulaData(Weights, Xi, Eta);

      // Iterators for quadrature points
      Vector::const_iterator Wq, Xiq, Etaq;

      // Iterate over all cells
      auto Cells = mesh.GetCellList();
      size_t CellInd = 0;

      // Vector of cell estimators
      Vector Res(Cells.size());
      
      for(auto Cell = Cells.begin(); Cell != Cells.end(); ++Cell, ++CellInd)
        {
          // Used for transformation to reference element
          CellInfo Info = mesh.GetCellInfo(CellInd);
          double det = mesh.Determinant(CellInd);
                      
          const Node& x0 = mesh.Nodes[Cell->LocNode[0]];
          Vector b(2); b[0] = x0.getX(); b[1] = x0.getY();
                      
          DenseMatrix Jac = mesh.Jacobian(CellInd);

          double Val = 0.;
          
          // Iterate over all terms
          for(auto Term = Terms.begin(); Term != Terms.end(); ++Term)
            {              
              switch(Term->GetType())
                {
                case VOLUME_RESIDUAL:
                       
                  for(Wq = Weights.begin(), Xiq = Xi.begin(), Etaq = Eta.begin();
                      Wq != Weights.end(); ++Wq, ++Xiq, ++Etaq)
                    {
                      Vector XiEtaq(2);
                      XiEtaq[0] = *Xiq;
                      XiEtaq[1] = *Etaq;

                      Vector XYq = b + Jac*XiEtaq;

                      Val += (*Wq) * pow(Info.Diam() * Term->EvalCoeff(XYq[0], XYq[1]), 2.) * det;
                    }
                      
                  break;
                default:
                  if(CellInd == 0)
                    std::cerr << "WARNING: Expression type not valid or not implemented yet for error estimation.\n";
                  break;
                }
            }

          Res[CellInd] = Val;
          
        }

      delete[] Value;

      return Res;
    }
  }
}
