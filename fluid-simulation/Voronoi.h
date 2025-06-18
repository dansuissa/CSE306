#pragma once
#include <vector>
#include "Vector.h"
#include "Polygon.h"
#include "nanoflann.hpp"

struct PointCloudAdapter
{
    const std::vector<Vec2>& m_points;
    PointCloudAdapter(const std::vector<Vec2>& points) : m_points(points) {}
    inline size_t kdtree_get_point_count() const { return m_points.size(); }
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if (dim == 0) return m_points[idx].x;
        else return m_points[idx].y;
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};
class Voronoi
{
public:
    Voronoi() = default;
    explicit Voronoi(const std::vector<Vec2>& pts) : P(pts) {}
    void compute();
    std::vector<Vec2>    P;     
    std::vector<Polygon> cells; 
};
