#include "TriangleMesh.h"
#include <limits>
#include <cstring>  
namespace rt{
void TriangleMesh::readOBJ(const char* fname)
{
    FILE* f = std::fopen(fname,"r");
    if(!f){ std::fprintf(stderr,"Cannot open %s\n",fname); return; }
    char line[512];
    while(std::fgets(line,512,f)){
        if(std::strncmp(line,"v ",2)==0){
            Vec3 v; std::sscanf(line,"v %lf %lf %lf",&v.x,&v.y,&v.z);
            vertices.push_back(v);
        }
        else if(std::strncmp(line,"vn",2)==0){
            Vec3 n; std::sscanf(line,"vn %lf %lf %lf",&n.x,&n.y,&n.z);
            normals.push_back(n);
        }
        else if(std::strncmp(line,"f ",2)==0){
            int i0,i1,i2,j0,k0,j1,k1,j2,k2, n;
            n = std::sscanf(line+2,"%d/%d/%d %d/%d/%d %d/%d/%d",
                            &i0,&j0,&k0,&i1,&j1,&k1,&i2,&j2,&k2);
            if(n==9){ indices.emplace_back(i0-1,i1-1,i2-1,
                                           k0-1,k1-1,k2-1); continue;}
            n = std::sscanf(line+2,"%d//%d %d//%d %d//%d",
                            &i0,&k0,&i1,&k1,&i2,&k2);
            if(n==6){ indices.emplace_back(i0-1,i1-1,i2-1,
                                           k0-1,k1-1,k2-1); continue;}
            n = std::sscanf(line+2,"%d %d %d",&i0,&i1,&i2);
            if(n==3){ indices.emplace_back(i0-1,i1-1,i2-1,-1,-1,-1);}
        }
    }
    std::fclose(f);
    std::fprintf(stderr,"Loaded: %zu verts  %zu tris\n",
                 vertices.size(),indices.size());
}
int TriangleMesh::buildNode(int start,int end)
{
    Node node;
    node.start=start; node.end=end; node.bbox.reset();
    for(int i=start;i<end;i++){
        node.bbox.expand(vertices[indices[i].v0]);
        node.bbox.expand(vertices[indices[i].v1]);
        node.bbox.expand(vertices[indices[i].v2]);
    }
    int idx = (int)bvh.size(); bvh.push_back(node);

    const int leafMax=8;
    if(end-start <= leafMax) return idx;

    Vec3 diag= node.bbox.bmax-node.bbox.bmin;
    int axis  = diag.x>diag.y ? (diag.x>diag.z?0:2)
                              : (diag.y>diag.z?1:2);
    double mid=
        axis==0 ? (node.bbox.bmin.x + node.bbox.bmax.x)*0.5 :
        axis==1 ? (node.bbox.bmin.y + node.bbox.bmax.y)*0.5 :
                  (node.bbox.bmin.z + node.bbox.bmax.z)*0.5 ;
    int pivot=start;
    for(int i=start;i<end;i++){
        Vec3 c = (vertices[indices[i].v0] +
                  vertices[indices[i].v1] +
                  vertices[indices[i].v2]) / 3.0;
        double val = axis==0? c.x : axis==1? c.y : c.z;
        if(val<mid){std::swap(indices[i],indices[pivot]); ++pivot;}
    }
    if(pivot<=start+ leafMax || pivot>=end-leafMax) return idx;
    bvh[idx].left =buildNode(start,pivot);
    bvh[idx].right = buildNode(pivot,end);
    return idx;
}
void TriangleMesh::buildBVH(){
    bvh.clear();
    root = buildNode(0,(int)indices.size());
}
bool TriangleMesh::intersectTri(const Ray& r,
                                const TriangleIndices& t,
                                Hit& h) const
{
    const Vec3& A =vertices[t.v0];
    const Vec3& B= vertices[t.v1];
    const Vec3& C= vertices[t.v2];
    Vec3 e1= B -A, e2= C -A;
    Vec3 P = cross(r.u,e2);
    double det =dot(e1, P);
    if(std::fabs(det)<1e-8) return false;
    double inv = 1.0/ det;
    Vec3 S = r.O -A;
    double u = inv*dot(S,P); if(u<0||u>1) return false;
    Vec3 Q = cross(S,e1);
    double v = inv *dot(r.u,Q); if(v<0||u+v>1) return false;
    double tHit = inv*dot(e2,Q); if(tHit<1e-5) return false;
    h.t = tHit;
    h.P = r.O + r.u*tHit;
    if(!normals.empty() && t.n0>=0){
        Vec3 n0=normals[t.n0],
             n1=normals[t.n1],
             n2=normals[t.n2];
        h.N = n0*(1-u-v)+n1*u +n2*v;
    }else{
        h.N=cross(e1,e2);
    }
    h.N.normalize();
    return true;
}

bool TriangleMesh::bruteIntersect(const Ray& r,Hit& h) const
{
    bool hit=false;
    double bestT = std::numeric_limits<double>::max();
    for(const auto& tri: indices){
        Hit tmp;
        if(intersectTri(r,tri,tmp) && tmp.t<bestT){
            bestT=tmp.t; h=tmp; hit=true;
        }
    }
    return hit;
}

bool TriangleMesh::intersect(const Ray& r,Hit& h) const
{
    if(root<0) return bruteIntersect(r,h);

    bool hit=false;
    double bestT = std::numeric_limits<double>::max();
    std::stack<int> st; st.push(root);

    while(!st.empty()){
        int ni = st.top(); st.pop();
        const Node& n = bvh[ni];
        double tNear;
        if(!n.bbox.intersect(r,tNear) || tNear>bestT) continue;

        if(n.left<0 && n.right<0){
            for(int i=n.start;i<n.end;i++){
                Hit tmp;
                if(intersectTri(r,indices[i],tmp) && tmp.t<bestT){
                    bestT=tmp.t; h=tmp; hit=true;
                }
            }
        }else{
            if(n.left >=0) st.push(n.left);
            if(n.right>=0) st.push(n.right);
        }
    }
    return hit;
}
} 
