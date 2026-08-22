#ifndef _COORDINATE_H_
#define _COORDINATE_H_

namespace chemfem{
  namespace linalg{

    /**
     * A point in the plane. A place, not a direction: see Vector2D for the latter.
     */
    struct Coordinate
    {
      double x, y;
    };

  };
};

#endif
