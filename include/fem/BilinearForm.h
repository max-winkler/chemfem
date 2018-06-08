#ifndef _BILINEAR_FORM_H_
#define _BILINEAR_FORM_H_

#include "linalg/SparseMatrix.h"
#include "fem/FESpace.h"
#include "fem/FEExpression.h"

using chemfem::linalg::SparseMatrix;

namespace chemfem{
  namespace fem{

    /**
     * This class represents a finite element bilinear form. In  linear algebra context 
     * this corresponds to a matrix (e.g. stiffness or mass matrix).
     */
    class BilinearForm
    {
    public:

      /**
       * Constructor which initializes an empty bilinear form for a given 
       * trial and test space.
       */
      BilinearForm(const FESpace&, const FESpace&);

      /// Adds a Laplace term (\nabla u,\nabla v) to the bilinear form
      void AddLaplaceTerm();

      /**
       * Adds the diffusion term div(a\nabla u) to the bilinear form. 
       * The coefficient a:R->R passed to the function is the diffusion coefficient.
       * In case of a(x)=1 this corresponds to the Laplace operator. In this case
       * the function AddLaplaceTerm should be used.
       */
      void AddDiffusionTerm(double (*)(double, double));

      /**
       * Adds a reaction term (c u,v) to the bilinear form. The reaction parameter
       * is given as a function pointer.
       */      
      void AddReactionTerm(double (*)(double, double));

      /**
       * Assembles the finite element matrix. 
       */
      void Assemble();

      /**
       * Returns the matrix which corresponds to the bilinear form. Before calling this function
       * the assemble routine has to be invoked. Otherwise, an empty matrix is returned.
       */
      SparseMatrix& SystemMatrix();
      
    private:
      const FESpace& TrialSpace, TestSpace;
      SparseMatrix Matrix;
      std::vector<FEExpression> Terms;
    };
    
  };
};

#endif
