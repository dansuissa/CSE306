#include "PowerDiagram.h"
#include <cmath>
#include <cstddef>
static Polygon clip_power(const Polygon& C,
                          const Vec2& pi, const Vec2& pj,
                          double wi, double wj)
{
    Polygon out;
    if (C.v.empty()) return out;
    Vec2   d{ pj.x - pi.x, pj.y - pi.y };
    double d2 = d.x*d.x + d.y*d.y;
    Vec2   m{ (pi.x+pj.x)*0.5 + d.x*(wi-wj)/(2*d2),
              (pi.y+pj.y)*0.5 + d.y*(wi-wj)/(2*d2) };
    auto side = [&](const Vec2& P)
    {
        return (P.x - m.x)*d.x + (P.y - m.y)*d.y;
    };
    const auto& v = C.v;
    for (std::size_t k = 0; k < v.size(); ++k)
    {
        std::size_t j = (k + 1) % v.size();
        const Vec2& A = v[k];
        const Vec2& B = v[j];
        double sa = side(A);
        double sb = side(B);
        if (sa * sb < 0.0)
        {
            double t = sa / (sa - sb);
            out.v.push_back({ A.x + t*(B.x - A.x),
                              A.y + t*(B.y - A.y) });
        }
        if (sb <= 0.0)
            out.v.push_back(B);
    }
    return out;
}
PowerDiagram::PowerDiagram(const std::vector<Vec2>& sites,
                           const std::vector<double>& weights)
: p(sites), w(weights) {}
void PowerDiagram::compute()
{
    Polygon box;
    box.v = { {0,0}, {1,0}, {1,1}, {0,1} };

    cells.assign(p.size(), {});

    for (std::size_t i = 0; i < p.size(); ++i)
    {
        Polygon cell = box;
        for (std::size_t j = 0; j < p.size(); ++j)
            if (i != j)
                cell = clip_power(cell, p[i], p[j], w[i], w[j]);
        cells[i] = std::move(cell);
    }
}
double PowerDiagram::area(std::size_t idx) const
{
    const auto& v = cells[idx].v;
    if (v.size() < 3) return 0.0;
    double a = 0.0;
    for (std::size_t k = 0; k < v.size(); ++k)
    {
        std::size_t j = (k + 1) % v.size();
        a += v[k].x * v[j].y - v[j].x * v[k].y;
    }
    return std::fabs(a) * 0.5;
}
