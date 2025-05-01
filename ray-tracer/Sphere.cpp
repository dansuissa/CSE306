#include "Sphere.h"
#include <cmath>

namespace rt{

bool Sphere::intersect(const Ray& r, Hit& h) const{
    Vec3 OC = r.O - C;
    double b = dot(r.u, OC);
    double c= OC.norm2() - R*R;
    double disc = b*b - c;
    if (disc < 0) return false;
    double sq = std::sqrt(disc);
    double t1 = -b-sq;
    double t2 = -b +sq;
    double t= (t1 > 1e-5 ? t1
                           : (t2 >1e-5 ? t2 : 1e30));
    if (t>= 1e30) return false;
    h.t=t;
    h.P= r.O+r.u*t;
    h.N = (h.P - C);
    h.N.normalize();
    if (invertNormal)
        h.N =-h.N;
    return true;
}

} 
