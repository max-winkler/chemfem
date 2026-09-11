#include "fem/PointValues.h"

using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector;

namespace chemfem{
  namespace fem{

    PointValues ReferenceValues(const Element& E, int k, double xi, double eta)
    {
      PointValues r;
      r.value = E.Value(k, xi, eta);
      r.gradient = E.Gradient(k, xi, eta);
      r.hessian = E.Hessian(k, xi, eta);
      r.laplacian = r.hessian.Trace();
      return r;
    }

    std::vector<PointValues> TabulateReference(const Element& E, const Vector& Xi,
                                               const Vector& Eta)
    {
      const size_t n = E.NrDof();
      std::vector<PointValues> table(Xi.size()*n);

      for(size_t q=0; q<Xi.size(); ++q)
        for(size_t k=0; k<n; ++k)
          table[q*n + k] = ReferenceValues(E, k, Xi[q], Eta[q]);

      return table;
    }

    PointValues MapFromReference(const PointValues& ref, const Matrix2D& InvJacT)
    {
      PointValues v;
      v.value = ref.value;
      v.gradient = InvJacT * ref.gradient;
      v.hessian = InvJacT * ref.hessian * InvJacT.Transpose();
      v.laplacian = v.hessian.Trace();
      return v;
    }

  }
}
