#ifndef _COORDINATE_H_
#define _COORDINATE_H_

namespace chemfem{
  namespace fem{

    /**
     * A point in the plane. Nothing but the location, so that the same type serves
     * as the argument of a coefficient, a right hand side or an exact solution.
     */
    struct Coordinate
    {
      double x, y;
    };

  };
};

#endif
