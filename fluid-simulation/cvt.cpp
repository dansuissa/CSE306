#include <vector>
#include <random>
#include <cstdio>
#include <chrono>
#include "Voronoi.h"
#include "svg.h"

constexpr int N = 20000;        
constexpr int ITER = 1;          
constexpr int REDDOTS = 1000;    

int main() {
    std::vector<Vec2> sites(N);
    std::mt19937 gen(42);
    std::uniform_real_distribution<> U(0, 1);
    for (auto& s : sites) {
        s.x = U(gen);
        s.y = U(gen);
    }
    printf("Starting optimized Voronoi computation for %d sites...\n", N);
    auto start = std::chrono::high_resolution_clock::now();
    Voronoi vor(sites);
    vor.compute();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    printf("Computation finished in %lld ms.\n", (long long)duration.count());
    std::vector<Vec2> samples(REDDOTS);
    for (auto& s : samples) {
        s.x = U(gen);
        s.y = U(gen);
    }
    save_svg(vor.cells, "cvt_uniform_optimized.svg", &samples, &sites);
    std::puts("cvt_uniform_optimized.svg  – ready.");
    return 0;
} 