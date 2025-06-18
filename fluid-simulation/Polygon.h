#pragma once
#include <vector>
#include "Vector.h"

struct Polygon
{
	std::vector<Vec2> v;     

	Vec2 centroid() const {
		double cx = 0.0, cy = 0.0, area = 0.0;
		const size_t n = v.size();
		if (n < 3) return {0, 0};
		for (size_t i = 0; i < n; ++i) {
			size_t j = (i + 1) % n;
			double cross = v[i].x * v[j].y - v[j].x * v[i].y;
			cx += (v[i].x + v[j].x) * cross;
			cy += (v[i].y + v[j].y) * cross;
			area += cross;
		}
		area *= 0.5;
		if (std::abs(area) < 1e-12) return {0, 0};
		cx /= (6.0 * area);
		cy /= (6.0 * area);
		return {cx, cy};
	}
};	