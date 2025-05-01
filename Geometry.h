#pragma once
#include "Vector.h"
#include "Ray.h"

namespace rt {
struct Hit { double t; Vec3 P,N; };
class Geometry {
public:
    Vec3  albedo {1,1,1};
    bool  isMirror= false;
    bool  isTransparent = false;
    virtual bool intersect(const Ray&,Hit&) const = 0;
    virtual ~Geometry() =default;
};

} 
