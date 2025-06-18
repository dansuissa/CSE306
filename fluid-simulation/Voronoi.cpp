#include "Voronoi.h"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <nanoflann.hpp>
namespace
{
    inline bool inside(const Vec2& X, const Vec2& Pi, const Vec2& Pj)
    {
        Vec2 M = 0.5*(Pi + Pj);      
        return dot(X - M , Pj - Pi) < 0.0;
    }
    inline Vec2 intersect(const Vec2& A, const Vec2& B,
                          const Vec2& Pi, const Vec2& Pj)
    {
        Vec2 M = 0.5*(Pi + Pj);
        Vec2 N = Pj - Pi;            

        double denom = dot(B - A, N);
        if (std::abs(denom) < 1e-12) return A; 
        double t = dot(M - A, N) / denom;
        return A + t * (B - A);
    }

    Polygon clip_by_bisector(const Polygon& poly,
                             const Vec2& Pi, const Vec2& Pj)
    {
        Polygon out;
        const size_t n = poly.v.size();
        if (n == 0) return out;

        for (size_t cur = 0; cur < n; ++cur)
        {
            size_t prev = (cur + n - 1) % n;
            const Vec2& A = poly.v[prev];
            const Vec2& B = poly.v[cur];
            bool Ain = inside(A, Pi, Pj);
            bool Bin = inside(B, Pi, Pj);
            if (Ain && Bin)           
            {
                out.v.push_back(B);
            }
            else if (Ain && !Bin)     
            {
                out.v.push_back(intersect(A,B,Pi,Pj));
            }
            else if (!Ain && Bin)      
            {
                out.v.push_back(intersect(A,B,Pi,Pj));
                out.v.push_back(B);
            }
        }
        return out;
    }
}

void Voronoi::compute()
{
    cells.assign(P.size(), {});
    PointCloudAdapter pc_adapter(P);
    using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PointCloudAdapter>,
        PointCloudAdapter, 2 /* dim */>;
    my_kd_tree_t index(2, pc_adapter, {10 /* max leaf */});
    Polygon big_box;
    big_box.v = { {-1,-1}, {2,-1}, {2,2}, {-1,2} };
#pragma omp parallel for schedule(dynamic,1)
    for (ptrdiff_t i = 0; i < (ptrdiff_t)P.size(); ++i)
    {
        size_t k = 20; 
        while (true) {
            size_t num_neighbors_to_find = std::min(k, P.size() > 1 ? P.size() - 1 : 0);
            if (num_neighbors_to_find == 0) {
                cells[i] = big_box;
                break;
            }
            Polygon cell = big_box;
            std::vector<size_t> ret_index(num_neighbors_to_find);
            std::vector<double> out_dist_sqr(num_neighbors_to_find);
            nanoflann::KNNResultSet<double> resultSet(num_neighbors_to_find);
            resultSet.init(&ret_index[0], &out_dist_sqr[0]);
            index.findNeighbors(resultSet, &P[i].x, nanoflann::SearchParameters());

            for (size_t neighbor_idx = 0; neighbor_idx < num_neighbors_to_find; ++neighbor_idx) {
                size_t j = ret_index[neighbor_idx];
                if (i == (ptrdiff_t)j) continue;
                cell = clip_by_bisector(cell, P[i], P[j]);
                if (cell.v.empty()) break;
            }

            if (cell.v.empty()) {
                k = P.size(); 
                continue;
            }

            double max_dist_sq = 0.0;
            for (const auto& v : cell.v) {
                max_dist_sq = std::max(max_dist_sq, (v - P[i]).norm2());
            }

            double furthest_neighbor_dist_sq = out_dist_sqr[num_neighbors_to_find - 1];
            if (furthest_neighbor_dist_sq >= 4.0 * max_dist_sq || num_neighbors_to_find >= P.size() - 1) {
                cells[i] = std::move(cell);
                break;
            }

            k *= 2; 
        }
    }
}
