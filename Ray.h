#pragma once
#include "Vector.h"
namespace rt {
struct Ray {Vec3 O, u;Ray(const Vec3& o,const Vec3& d):O(o),u(d){}};
} 
