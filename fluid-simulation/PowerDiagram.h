#pragma once
#include <vector>
#include "Vector.h"          
#include "Polygon.h"         

class PowerDiagram
{
public:
    PowerDiagram(const std::vector<Vec2>& sites,
                 const std::vector<double>& weights);

    void   compute();                
    double area(std::size_t i) const; 

    std::vector<Polygon> cells;      

private:
    std::vector<Vec2>    p;          
    std::vector<double>  w;          
};