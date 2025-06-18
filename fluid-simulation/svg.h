#pragma once
#include <cstdio>
#include <string>
#include "Polygon.h"
#include "Vector.h"

inline void save_svg(const std::vector<Polygon>& polys,
                     const std::string&          filename,
                     const std::vector<Vec2>*    red_points = nullptr,
                     const std::vector<Vec2>*    black_points = nullptr)
{
    FILE* f = std::fopen(filename.c_str(), "w");
    if (!f) return;
    std::fprintf(f,
      "<svg xmlns='http://www.w3.org/2000/svg' width='1000' height='1000'>\n");
    for (const auto& poly : polys) {
        if (poly.v.empty()) continue;
        std::fprintf(f,"<polygon points='");
        for (const auto& p : poly.v)
            std::fprintf(f,"%f,%f ", p.x*1000.0, 1000.0 - p.y*1000.0);
        std::fprintf(f,"' fill='none' stroke='black' stroke-width='1'/>\n");
    }

    if (red_points) {
        for (const auto& p : *red_points)
            std::fprintf(f,
                "<circle cx='%f' cy='%f' r='2' fill='red'/>\n",
                p.x*1000.0, 1000.0 - p.y*1000.0);
    }
    if (black_points) {
        for (const auto& p : *black_points)
            std::fprintf(f,
                "<circle cx='%f' cy='%f' r='2' fill='black'/>\n",
                p.x*1000.0, 1000.0 - p.y*1000.0);
    }
    std::fprintf(f,"</svg>\n");
    std::fclose(f);
}
