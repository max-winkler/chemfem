#include <iostream>
#include <cmath>
#include <cassert>

#include "fem/ErrorEstimator.h"
#include "quadrature/QuadFormula.h"
#include "mesh/Mesh.h"

using chemfem::linalg::DenseMatrix;
using chemfem::mesh::Mesh;
using chemfem::mesh::Node;
using chemfem::mesh::Edge;
using chemfem::mesh::Cell;
using chemfem::mesh::CellInfo;
using chemfem::mesh::EdgeType;
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
      
      const FESpace& Space = u.GetFESpace();
      const Mesh& mesh = Space.GetMesh();

      // Store local FE Function and their derivatives
      double *Value = new double[Space.NrLocalDof()];
      
      // Initialize quadrature formulas
      QuadratureFormula QuadFormula(QUAD_FORMULA::GAUSS_7);
      QuadratureFormula LineQuadFormula(QUAD_FORMULA::LINE_GAUSS_5);
      
      Vector Xi, Eta, Weights;
      QuadFormula.FormulaData(Weights, Xi, Eta);

      Vector XiLine, EtaLine, WeightsLine;
      QuadFormula.FormulaData(WeightsLine, XiLine, EtaLine);
      
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
          DenseMatrix InvJac(Jac.Transpose().Invert());
          
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

                      Val += (*Wq) * pow(Info.Diam() * Term->EvalCoeff(Coordinate{XYq[0], XYq[1]}), 2.) * det;
                    }
	        
	        if(CellInd == 0)
                    std::cerr << "WARNING: Volume residual estimator only implemented for P1 elements.\n";
                  break;
                case EDGE_JUMP:
                  for(int edge_idx=0; edge_idx<3; ++edge_idx)
                    {
                      // Normal vector
                      Vector n = Info.Normal(edge_idx);
		  double hE = Info.EdgeLength(edge_idx);

		  // Get neighbor cell
		  Edge edge = mesh.GetEdgeList()[Cell->LocEdge[edge_idx]];
		  
		  if(edge.Type() != EdgeType::INTERFACE_EDGE)
		    continue;
		  
		  int NeighCellInd = edge.GetNeighbor(CellInd);
		  
		  const mesh::Cell& NeighCell = mesh.GetCellList()[NeighCellInd];

		  DenseMatrix NeighJac = mesh.Jacobian(NeighCellInd);
		  DenseMatrix NeighInvJac(Jac.Transpose().Invert());		  

		  size_t neigh_edge_idx = NeighCell.EdgeIndex(Cell->LocEdge[edge_idx]);
		  
                      // Iterate over quadrature points
                      for(Wq = WeightsLine.begin(), Xiq = XiLine.begin(); Wq != WeightsLine.end(); ++Wq, ++Xiq)
                        {
                          // Compute gradient in quadrature point
                          Vector GradValue(2);
		      Vector NeighGradValue(2);
		      
                          for(size_t k=0; k<Space.NrLocalDof(); ++k)
                            {
			// \todo Add the following lines as a method for FEFunction
			
                              // Line coordinates to triangle coordinates
                              double xi, eta;
                              switch(edge_idx){
                              case 0: xi = *Xiq;    eta = 0.; break;
                              case 1: xi = 1.-*Xiq; eta = *Xiq; break;
                              case 2: xi = 0.;      eta = 1.-*Xiq; break;
                              }
			
                              GradValue += u[Space.GetGlobalIndex(CellInd,k)]
			  * Space.RefElement().Gradient(k, xi, eta);
			
			// Do the same for the neighbor (mind opposite edge orientation)
                              switch(neigh_edge_idx){
                              case 0: xi = 1.-*Xiq;    eta = 0.; break;
                              case 1: xi = *Xiq; eta = 1.-*Xiq; break;
                              case 2: xi = 0.;      eta = *Xiq; break;
                              }
			
                              NeighGradValue += u[Space.GetGlobalIndex(NeighCellInd,k)]
			  * Space.RefElement().Gradient(k, xi, eta);
			
                            }
		      
                          Val += 0.5*(*Wq) * pow(dot(InvJac*GradValue - NeighInvJac*NeighGradValue, n), 2.)
		        * hE * hE;
                        }
                    }	                          
                  break;
	      default:
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
