#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <iostream>

namespace chemfem{
  namespace linalg{

    /**
     * This class represents a vector and is used for the storage of finite element functions.
     */
    class Vector
    {
    public:

      /**
       * Create an empty vector with 0 components.
       */
      Vector();
      
      /**
       * Constructor initalizing a vector with n components. The initial value is zero.
       */
      Vector(const size_t);

      /**
       * Destructor, deletes used memory.
       */
      ~Vector();
      
      /**
       * Returns the length of the vector.
       */
      size_t size() const;

      /**
       * Access a single element of a vector.
       */
      double& operator[](const size_t);

      /** 
       * Copy an instance of the class Vector.
       */
      Vector& operator=(const Vector&);

      /**
       * Copy the values of an array to the vector
       */
      Vector& operator=(const double*);
    private:
      /// Dimension of the vector
      int n;
      // Data vector
      double* data;
    };    
    
  };
};
  

#endif
