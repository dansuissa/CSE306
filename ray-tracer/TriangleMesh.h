#pragma once
#include <vector>
#include <stack>
#include <cstdio>
#include "Vector.h"
#include "Ray.h"
#include "AABB.h"
#include "Geometry.h"

namespace rt {
struct TriangleIndices{
    int v0,v1,v2, n0,n1,n2;
    TriangleIndices(int a,int b,int c,int na,int nb,int nc)
     : v0(a),v1(b),v2(c), n0(na),n1(nb),n2(nc){}
};

class TriangleMesh : public Geometry{
public:
    std::vector<Vec3> vertices, normals;
    std::vector<TriangleIndices> indices;
    void readOBJ (const char* file);
    void buildBVH();
    bool  intersect(const Ray& r, Hit& h) const override;

private:
    struct Node {
        AABB bbox;
        int  start,end; 
        int  left=-1,right=-1;
    };
    std::vector<Node> bvh;
    int root = -1;
    int buildNode(int start,int end);
    bool intersectTri(const Ray&,const TriangleIndices&,Hit&) const;
    bool bruteIntersect(const Ray&,Hit&) const;
};
} 
