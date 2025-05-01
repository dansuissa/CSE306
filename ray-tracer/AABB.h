#pragma once
#include "Vector.h"
#include "Ray.h"
#include <limits>

namespace rt{

struct AABB{
    Vec3 bmin, bmax;
    void reset() {
        const double inf = 1e30;
        bmin ={ inf, inf,inf};
        bmax = {-inf, -inf, -inf};
    }
    void expand(const Vec3& v){
        bmin.x = std::min(bmin.x,v.x);
        bmin.y = std::min(bmin.y,v.y);
        bmin.z=std::min(bmin.z,v.z);
        bmax.x = std::max(bmax.x,v.x);
        bmax.y =std::max(bmax.y,v.y);
        bmax.z = std::max(bmax.z,v.z);
    }
    bool intersect(const Ray& r,double& tNear) const{
        double t0 = 0, t1 = std::numeric_limits<double>::max();
        for(int i=0;i<3;i++){
            double orig =i==0? r.O.x : i==1? r.O.y : r.O.z;
            double dir  = i==0? r.u.x : i==1? r.u.y : r.u.z;
            double minb = i==0? bmin.x: i==1? bmin.y: bmin.z;
            double maxb = i==0? bmax.x: i==1? bmax.y: bmax.z;
            double inv= 1.0/dir;
            double tA =(minb-orig)*inv;
            double tB = (maxb-orig)*inv;
            if(tA>tB) std::swap(tA,tB);
            t0 = std::max(t0,tA);
            t1 = std::min(t1,tB);
            if(t0>t1) return false;
        }
        tNear=t0;
        return true;
    }
};

} 
