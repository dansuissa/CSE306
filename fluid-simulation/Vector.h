#pragma once
#include <cmath>
#include <ostream>

struct Vec2
{
    double x{}, y{};
    Vec2() = default;
    Vec2(double X, double Y) : x(X), y(Y) {}

    Vec2  operator+(const Vec2& b) const { return {x + b.x, y + b.y}; }
    Vec2  operator-(const Vec2& b) const { return {x - b.x, y - b.y}; }
    Vec2  operator*(double s) const { return {x * s, y * s};}
    Vec2  operator/(double s) const { return {x / s, y / s};}
    Vec2& operator+=(const Vec2& b) { x+=b.x; y+=b.y; return *this; }

    double norm2() const { return x*x + y*y; }
    double norm() const { return std::sqrt(norm2()); }
};

inline double dot(const Vec2& a, const Vec2& b) { return a.x*b.x + a.y*b.y; }
inline Vec2 operator*(double s, const Vec2& v)  { return {v.x*s, v.y*s};   }

inline std::ostream& operator<<(std::ostream& o,const Vec2& p)
{ return o << '(' << p.x << ',' << p.y << ')'; }
