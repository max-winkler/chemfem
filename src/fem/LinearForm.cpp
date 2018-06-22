#include "fem/LinearForm.h"

#include "quadrature/QuadFormula.h"

using chemfem::linalg::Vector;
using chemfem::linalg::DenseMatrix;

using chemfem::mesh::Mesh;
using chemfem::mesh::Node;
using chemfem::mesh::Cell;

using chemfem::quadrature::QuadratureFormula;
using chemfem::quadrature::QUAD_FORMULA;

namespace chemfem{
  namespace fem{

    LinearForm::LinearForm(const FESpace& TestSpace) : TestSpace(TestSpace) {}

    void LinearForm::AddVolumeForce(double (*F)(double, double))
    {
      FEExpression Expression(VOLUME_FORCE, F);
      Terms.push_back(Expression);
    }

    void LinearForm::AddNeumannBC(double (*G)(double, double))
    {
      FEExpression Expression(VOLUME_FORCE, G);
      Terms.push_back(Expression);
    }

    Vector& LinearForm::LoadVector()
    {
      return Vec;
    }

    void LinearForm::Assemble()
    {
      Vec = Vector(TestSpace.NrFreeDof());

      // TODO: Select correct quadrature formula once it is implemented
      QuadratureFormula QuadFormula(QUAD_FORMULA::GAUSS_7);

      Vector Xi, Eta, Weights;
      QuadFormula.FormulaData(Weights, Xi, Eta);

      double *TestFuncValue = new double[TestSpace.DofPerCell];
      
      // Iterate over all cells
      int CellInd;
      std::vector<Cell>::const_iterator cell;
      for(cell = TestSpace.mesh.Cells.begin(), CellInd=0;
	  cell != TestSpace.mesh.Cells.end(); ++cell, ++CellInd)
	{
	  double det = cell->Determinant();
	  
	  Node* x0 = cell->LocNode[0];
	  Vector b(2); b[0] = x0->getX(); x0->getY();
	  
	  DenseMatrix Jac = cell->Jacobian();

	  Vector LocVec(TestSpace.DofPerCell);
	  
	  // Iterate over all quadrature points
	  Vector::const_iterator Wq, Xiq, Etaq;
	  
	  for(Wq = Weights.begin(), Xiq = Xi.begin(), Etaq = Eta.begin();
	      Wq != Weights.end(); ++Wq, ++Xiq, ++Etaq)
	    {
	      // Determine Quadrature points in world element
	      Vector XiEtaq(2);
	      XiEtaq[0] = *Xiq;
	      XiEtaq[1] = *Etaq;

	      Vector XYq = b + Jac*XiEtaq;

	      // Function value of test functions
	      for(int k=0; k<TestSpace.DofPerCell; ++k)
		TestFuncValue[k] = TestSpace.RefElement().Value(k, *Xiq, *Etaq);
	      
	      // Iterate over all terms
	      for(std::vector<FEExpression>::const_iterator Term = Terms.begin();
		  Term != Terms.end(); ++Term)
		{
		  double CoeffVal = Term->Coeff(XYq[0], XYq[1]);
			
		  switch(Term->Type)
		    {
		    case VOLUME_FORCE:
			  
		      for(int k=0; k<TestSpace.DofPerCell; ++k)
			LocVec[k] += 0.5 * (*Wq) * CoeffVal * TestFuncValue[k] * det;
		      break;
		      /*
			case NEUMANN_BC:
			    
			break;
		      */
		    default:
		      std::cerr << "Assembly of FE expressions of type " << Term->Type
				<< " not implemented yet.\n";
		      return;
		    }
		      
		} // loop over Terms
	    } // loop over quadrature points
	  for(int k=0; k<TestSpace.DofPerCell; ++k)
	    {
	      int GlobalIndex = TestSpace.GetGlobalIndex(CellInd, k);
	      if(TestSpace.DofType[GlobalIndex])
	
		Vec[TestSpace.DofIndex[GlobalIndex]] += LocVec[k];
	    }
		
	} // loop over cells

    }
  }
}
