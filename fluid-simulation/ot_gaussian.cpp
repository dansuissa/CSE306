#include <random>
#include "OT.hpp"

int main()
{
    Polygon box; box.v = { {0,0},{1,0},{1,1},{0,1} };

    const int N = 3000;
    std::vector<Vec2> sites(N);
    std::default_random_engine eng(1);
    std::uniform_real_distribution<double> U(0.,1.);

    for (auto& s : sites) { s.x = U(eng); s.y = U(eng); }

    OT solver(box, sites, 400);
    solver.run();
    solver.saveSVG("power_uniform.svg");
    return 0;
}