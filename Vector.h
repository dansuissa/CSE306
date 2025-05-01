#pragma once
#include <algorithm>
#include "MathUtils.h"
namespace rt{
struct Vec3{
    double x{}, y{}, z{};
    constexpr Vec3() = default;
    constexpr Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3  operator+(const Vec3& o) const {return {x+o.x, y + o.y, z+o.z};}
    Vec3  operator-(const Vec3& o) const {return {x- o.x, y-o.y, z-o.z};}
    Vec3  operator-() const {return {-x,-y,-z};}
    Vec3  operator*(double s) const {return {x*s, y* s, z*s};}
    Vec3  operator/(double s) const {return {x /s, y/s, z/s};}
    Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this;}

    /* geometry ----------------------------------------------------------- */
    double norm2() const {return x*x + y*y + z*z;}
    double norm() const { return std::sqrt(norm2());}
    void normalize(){double n = norm(); x/=n; y/=n; z/=n;}
};
inline Vec3 operator*(double s,const Vec3& v){return v*s;}
inline double dot(const Vec3& a,const Vec3& b){return a.x*b.x + a.y*b.y + a.z*b.z;}
inline Vec3  cross(const Vec3& a,const Vec3& b){
    return { a.y*b.z - a.z*b.y,
             a.z*b.x - a.x*b.z,
             a.x*b.y - a.y*b.x };
}
} 
