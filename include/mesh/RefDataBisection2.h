#ifndef _REF_TYPE_BISECTION_2_H_
#define _REF_TYPE_BISECTION_2_H_

#include "mesh/RefData.h"

namespace chemfem{
  namespace mesh{

    /**
     * This class is derived from RefData and implements a regular refinement of a single cell.
     */
    class RefDataBisection2 : public RefData
    {
    public:
      /// Constructor which sets the local member variables.
      RefDataBisection2();
    };
    
  };
};

#endif
