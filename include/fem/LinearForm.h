#ifndef _LINEAR_FORM_
#define _LINEAR_FORM_

#include <vector>

#include "linalg/Vector.h"
#include "fem/FESpace.h"
#include "fem/FEExpression.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace fem{

    /**
     * This class represents a linear form which stores and assembles the vector for the right-hand 
     * side or Neumann boundary conditions.
     */
    class LinearForm
    {
    public:
      /**
       * Constructor which associates the linear form with a function space.
       */
      LinearForm(const FESpace&);
      
      /**
       * Adds a volume force. This is a function handle to the function defining the right-hand 
       * side of the partial differential equation.
       */
      void AddVolumeForce(double (*)(double, double));

      /**
       * Adds a Neumann boundary condition. Requires a function handle to the function definiting the 
       * boundary condition.
       */
      void AddNeumannBC(double (*)(double, double));

      /**
       * Assembles the load vector. Before calling this routine all terms that are required should be 
       * added to the linear form.
       */
      void Assemble();

      /**
       * Returns a reference to the load vector. Before calling this routine the assemble function
       * has to be invoked.
       */
      Vector& LoadVector();
      
    private:
      /**
       * Stores the terms added to the linear form in a vector.
       */
      std::vector<FEExpression> Terms;

      /**
       * Stores a reference to the test space.
       */
      const FESpace& TestSpace;

      /**
       * The load vector created in the assemble routine.
       */
      Vector Vec;
    };
    
  };
};

#endif
