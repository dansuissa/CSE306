#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include "MathUtils.h"      
#include "Vector.h"          
#include "Ray.h"             
#include "Sphere.h"      
#include "stb_image_write.h"  
using namespace rt;
struct Scene{
    std::vector<Sphere*> objects;
    Vec3   lightPos = {-10,20, 40};
    double lightI= 1e5;
    void add(Sphere* s) { objects.push_back(s); }
    bool intersect(const Ray& r, Hit& rec, Sphere*& hitObj) const{
        bool hitAny =false;
        double bestT= 1e30;
        for (auto* obj : objects) {
            Hit h;
            if (obj->intersect(r, h) && h.t < bestT) {
                bestT = h.t;
                rec= h;
                hitObj= obj;
                hitAny= true;
            }
        }
        return hitAny;
    }

    Vec3 getColor(const Ray& r, int depth = 7) const {
        if (depth==0) return {0,0,0};
        Hit  rec; Sphere* obj = nullptr;
        if (!intersect(r, rec, obj)) return {0,0,0};
        Vec3 P=rec.P;
        Vec3 N =rec.N;
        Vec3 I =r.u;
        if (obj->isMirror) {
            Vec3 R = I-N *(2 *dot(I, N));
            return getColor(Ray(P +N *1e-4, R), depth - 1);
        }
        if (obj->isTransparent) {
            double n1 =1.0, n2 =1.8;
            Vec3 Ndir = N;
            bool inside=false;
            if (dot(I, Ndir) > 0) {
                std::swap(n1, n2);
                Ndir= -Ndir;
                inside = true;
            }
            double cosi = dot(I, Ndir);
            double eta= n1/n2;
            double k= 1 - eta*eta *(1-cosi*cosi);
            if (k < 0) {
                Vec3 R = I-Ndir * (2 * dot(I,Ndir));
                return getColor(Ray(P + Ndir*1e-4, R),depth -1);
            }
            Vec3 refr = I*eta - Ndir*(eta*cosi + std::sqrt(k));
            refr.normalize();
            double R0 = sqr((n1 - n2)/(n1 + n2));
            double Fr = R0 +(1-R0)*std::pow(1 - std::fabs(cosi),5);
            Fr = std::min(1.0, Fr *(inside ? 1.2 :1.0));
            Vec3 R = I - Ndir *(2 * dot(I, Ndir));
            Vec3 colR = getColor(Ray(P + Ndir*1e-4, R), depth - 1);
            Vec3 colT = getColor(Ray(P - Ndir*1e-4, refr), depth - 1);
            if (inside) {
                colT = Vec3(colT.x*0.95, colT.y*0.95, colT.z*0.98);
            }
            return colR*Fr + colT*(1-Fr);
        }
        Vec3 L = lightPos - P;
        double d2 = L.norm2();
        L.normalize();
        double lam = std::max(0.0, dot(L, N));
        Vec3 base = (lightI/(4*PI*d2)) * (obj->albedo/PI) * lam;
        Ray shadow(P + N*1e-4, L);
        Hit sh; Sphere* shObj = nullptr;
        if (intersect(shadow, sh, shObj) && sh.t*sh.t < d2) {
            if (shObj && shObj->isTransparent)
                base = base * 0.8;   
            else
                base = {0,0,0};
        }
        return base;
    }
};

int main(){
    const int W = 512, H = 512, spp = 10;
    Scene scene;

    scene.add(new Sphere({0,0,-1000}, 940, {0.4,0.8,0.7}));
    scene.add(new Sphere({0,0,1000}, 940, {0.9,0.4,0.3}));
    scene.add(new Sphere({0,-1000,0}, 990, {0.3,0.4,0.7}));
    scene.add(new Sphere({-1000,0,0}, 940, {0.5,0.5,0.6}));
    scene.add(new Sphere({ 1000,0,0}, 940, {0.6,0.5,0.1}));
    scene.add(new Sphere({ 0,1000,0}, 940, {0.2,0.5,0.9}));
    scene.add(new Sphere({-20,0, 0}, 10,   {0.5,0.5,0.5}, true,  false));      
    scene.add(new Sphere({0,0, 0},10,   {0.5,0.5,0.5}, false, true ));      
    scene.add(new Sphere({20, 0,0},10,   {0.5,0.5,0.5}, false, true ));      
    scene.add(new Sphere({ 20,0, 0},9.9,  {0.5,0.5,0.5}, false, true, true));  
    std::vector<unsigned char> img(W*H*3);
    #pragma omp parallel for schedule(dynamic,1)
    for(int y=0; y<H; y++){
      for(int x=0; x<W; x++){
        Vec3 col{0,0,0};
        for(int s=0; s<spp; s++){
          double dx = 0.5, dy = 0.5;
          double d = -W/(2.0*std::tan(60.0*PI/180.0/2));
          Vec3 dir{ x - W/2.0+dx, H/2.0 -y+dy,d};
          dir.normalize();
          col += scene.getColor(Ray({0,0,55},dir),7);
        }
        col = col / double(spp);
        col.x = std::pow(col.x,1/2.2);
        col.y = std::pow(col.y,1/2.2);
        col.z = std::pow(col.z,1/2.2);
        int idx = 3*(y*W + x);
        img[idx  ] = (unsigned char)(std::min(1.0, col.x)*255);
        img[idx+1] = (unsigned char)(std::min(1.0, col.y)*255);
        img[idx+2] = (unsigned char)(std::min(1.0, col.z)*255);
      }
    }
    stbi_write_png("spheres.png", W,H,3, img.data(),0);
    return 0;
}
