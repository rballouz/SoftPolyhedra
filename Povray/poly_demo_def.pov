#include "colors.inc"
#include "shapes.inc"
#include "polyhedron.inc"

camera {
    location <3, -3, 3>
    look_at <0, 0, 0>
}
  
light_source { <20, 20, 20> color White }

object {polyhedron1
        rotate z*360
        }
