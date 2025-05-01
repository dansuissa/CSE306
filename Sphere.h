#pragma once
#include "Geometry.h"    
namespace rt{
class Sphere : public Geometry{
public:
    Vec3  C;               
    double R;              
    bool   invertNormal;   
    Sphere(const Vec3& centre, double radius,
           const Vec3& alb,
           bool mirror= false,
           bool transparent=false,
           bool invNorm=false)
      : C(centre)
      , R(radius)
      , invertNormal(invNorm)
    {
        albedo= alb;
        isMirror= mirror;
        isTransparent=transparent;
    }
    bool intersect(const Ray& r, Hit& h) const override;
};

} 
