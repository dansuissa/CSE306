#define _CRT_SECURE_NO_WARNINGS
#include "stb_image_write.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include "MathUtils.h"
#include "Vector.h"
#include "Ray.h"
#include "Sphere.h"
#include "TriangleMesh.h"
using namespace rt;

struct Scene {
    std::vector<Geometry*> objects;
    Vec3  lightPos;
    double lightI =1e5;
    void add(Geometry* g){ objects.push_back(g); }
    bool intersect(const Ray& r,Hit& h,int& id) const {
        bool hit=false; h.t=1e30;
        for(int i=0;i<(int)objects.size();i++){
            Hit tmp; if(objects[i]->intersect(r,tmp)&&tmp.t<h.t){
                h=tmp; id=i; hit=true;
            }
        }
        return hit;
    }
    Vec3 shade(const Ray& r,int depth=5) const {
        if(depth==0) return {0,0,0};
        Hit h; int id;
        if(!intersect(r,h,id)) return {0,0,0};
        const auto* obj = objects[id];
        if(obj->isMirror){
            Vec3 R = r.u-h.N*(2*dot(r.u,h.N));
            return shade(Ray(h.P+h.N*1e-4,R),depth-1);
        }
        Vec3 L = lightPos - h.P; double d2=L.norm2(); L.normalize();
        double lam = std::max(0.0,dot(L,h.N));
        Vec3 col = (lightI/(4*PI*d2))*(obj->albedo/PI)*lam;
        Ray shadow(h.P+h.N*1e-4,L);
        Hit hs; int sid;
        if(intersect(shadow,hs,sid) && hs.t*hs.t<d2) col={0,0,0};
        return col;
    }
};

int main(){
    const int W=512,H=512,spp=32;
    Scene scene; scene.lightPos={-10,20,40};
    scene.add(new Sphere({0,0,-1000},940,{0.4,0.8,0.7}));
    scene.add(new Sphere({0,0, 1000},940,{0.9,0.4,0.3}));
    scene.add(new Sphere({0, -1000,0},990,{0.3,0.4,0.7}));
    scene.add(new Sphere({-1000,0,0},940,{0.5,0.5,0.6}));
    scene.add(new Sphere({1000,0,0},940,{0.6,0.5,0.1}));
    scene.add(new Sphere({0,1000,0},940,{0.2,0.5,0.9}));
    auto* cat = new TriangleMesh();
    cat->readOBJ("cat.obj");
    double scale=0.6; Vec3 trans{0,-10,0};
    for(auto& v: cat->vertices) v = v*scale + trans;
    cat->albedo={0.95,0.95,0.98};
    cat->buildBVH();
    scene.add(cat);
    std::vector<unsigned char> img(W*H*3);
    #pragma omp parallel for schedule(dynamic,1)
    for(int y=0;y<H;y++){
        for(int x=0;x<W;x++){
            Vec3 acc{0,0,0};
            for(int s=0;s<spp;s++){
                double dx=0.5,dy=0.5;
                double d =-W/(2.0*std::tan(60.0*PI/180.0/2));
                Vec3 dir{ x-W/2.0+dx, H/2.0-y+dy, d};
                dir.normalize();
                acc += scene.shade(Ray({0,0,55},dir));
            }
            acc = acc / double(spp);
            acc.x = std::pow(acc.x,1/2.2);
            acc.y = std::pow(acc.y,1/2.2);
            acc.z = std::pow(acc.z,1/2.2);
            int idx = 3*(y*W+x);
            img[idx  ] = (unsigned char)(std::min(1.0,acc.x)*255);
            img[idx+1] = (unsigned char)(std::min(1.0,acc.y)*255);
            img[idx+2] = (unsigned char)(std::min(1.0,acc.z)*255);
        }
    }
    stbi_write_png("cat.png",W,H,3,img.data(),0);
    return 0;
}
