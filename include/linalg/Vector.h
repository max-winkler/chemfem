#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <vector>

namespace chemfe{
  namespace linalg{

    /**
     * This class represents a vector and is used for the storage of finite element functions.
     */
    class Vector
    {
    public:
      /**
       * Constructor initalizing a vector with n components. The initial value is zero.
       */
      Vector(const size_t);

      /**
       * Returns the length of the vector.
       */
      size_t size();

      /**
       * Access a single element of a vector.
       */
      double& operator[](const size_t);
      
    private:
      /// Dimension of the vector
      const int n;
      std::vector<double> data;
    };    
    
  };
};
  

#endif
