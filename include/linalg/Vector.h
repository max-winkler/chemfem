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
       * Constant iterator class.
       */
      class const_iterator{
      public:
	const_iterator();
	const_iterator(const const_iterator&);
	const_iterator(const double*);
	const_iterator& operator=(const const_iterator& it);
	bool operator==(const const_iterator&);
	bool operator!=(const const_iterator&);
	const_iterator& operator++();
	const double& operator*() const;
      private:
	const double *cur;
      };
      
      /**
       * Create an empty vector with 0 components.
       */
      Vector();
      
      /**
       * Constructor initalizing a vector with n components. The initial value is zero.
       */
      Vector(const size_t);

      /**
       * Copy constructor.
       */
      Vector(const Vector&);
      
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
       * Access a single element of a constant vector.
       */
      const double& operator[](const size_t) const;

      /** 
       * Copy an instance of the class Vector.
       */
      Vector& operator=(const Vector&);

      /**
       * Copy the values of an array to the vector
       */
      Vector& operator=(const double*);

      /**
       * Computes the sum of 2 vectors.
       */
      Vector operator+(const Vector&) const;

      /**
       * Computes the difference of 2 vectors.
       */
      Vector operator-(const Vector&) const;

      /**
       * Returns the scalar vector product a*v.
       */
      friend Vector operator*(const double&, const& Vector);
      
      /**
       * Computes expressions like a*x+y for real values a and vectors x and y.
       */
      void axpy(double a, const Vector&, Vector&) const;
      
      /**
       * For console or file output.
       */
      friend std::ostream& operator<<(std::ostream&, const Vector&);
      
      /**
       * Computes the p-norm of the vector. If the argument if neglected, the two-norm is returned.
       */
      double Norm(double b = 2.) const;
      
      /**
       * Returns an iterator pointing to the beginning of the vector.
       */
      const_iterator begin() const;

      /**
       * Returns an iterator pointing to the end of the vector.
       */
      const_iterator end() const;
      
    private:
      /// Dimension of the vector
      int n;
      // Data vector
      double* data;

    };

    // Free methods
    double dot(const Vector&, const Vector&);
  };
};
  

#endif
