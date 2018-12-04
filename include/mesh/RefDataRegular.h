#ifndef _REF_TYPE_REGULAR_H_
#define _REF_TYPE_REGULAR_H_

#include "mesh/RefData.h"

namespace chemfem{
  namespace mesh{

    /**
     * This class is derived from RefData and implements a regular refinement of a single cell.
     */
    class RefDataRegular : public RefData
    {
    public:
      /// Constructor which sets the local member variables.
      RefDataRegular();
    };

  };
};

#endif
