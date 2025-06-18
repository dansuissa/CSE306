#pragma once
#include <vector>
#include <cmath>
#include <cstdio>
#include "Vector.h"
#include "Polygon.h"
#include "PowerDiagram.h"
#include "lbfgs.h"
class OT
{
public:
    OT(const Polygon& clip,
       const std::vector<Vec2>& sites,
       int max_it = 400)
    : clipPoly(clip), p(sites), maxIter(max_it) {}

    int  run();                                  
    void saveSVG(const char* file);  

private:
    const Polygon             clipPoly;
    std::vector<Vec2>         p;             
    std::vector<double>       lambda;         
    int                       maxIter;
    std::vector<double>       w;              
    PowerDiagram              PD = PowerDiagram(p,{});
    static lbfgsfloatval_t _eval (void*,const lbfgsfloatval_t*,
                                  lbfgsfloatval_t*,int,lbfgsfloatval_t);
    static int            _prog (void*,const lbfgsfloatval_t*,
                                  const lbfgsfloatval_t*,
                                  lbfgsfloatval_t,lbfgsfloatval_t,
                                  lbfgsfloatval_t,lbfgsfloatval_t,
                                  int,int,int);

    lbfgsfloatval_t eval (const lbfgsfloatval_t* x, lbfgsfloatval_t* g,
                          int n, lbfgsfloatval_t step);
    int  prog (const lbfgsfloatval_t* x,const lbfgsfloatval_t* g,
               lbfgsfloatval_t fx,lbfgsfloatval_t xnorm,
               lbfgsfloatval_t gnorm,lbfgsfloatval_t step,
               int n,int k,int ls);
};
inline int OT::run()
{
    const int N = static_cast<int>(p.size());
    lambda.resize(N);
    Vec2 C{0.5,0.5};
    double tot = 0.0;
    for (int i=0;i<N;++i)
    {
        Vec2 d{C.x-p[i].x, C.y-p[i].y};
        double r2 = d.x*d.x + d.y*d.y;
        lambda[i] = std::exp(-r2/0.2);
        tot+= lambda[i];
    }
    for (double& v:lambda) v /= tot;
    w.assign(N, 0.0);
    lbfgs_parameter_t prm; lbfgs_parameter_init(&prm);
    prm.max_iterations = maxIter;

    lbfgsfloatval_t fx;
    int ret = lbfgs(N, w.data(), &fx, _eval, _prog, this, &prm);
    std::printf("LBFGS terminated (status %d)  f = %.8g\n",ret,fx);

    PD = PowerDiagram(p,w); PD.compute();
    return ret;
}
inline lbfgsfloatval_t OT::_eval(void*inst,const lbfgsfloatval_t*x,
                                 lbfgsfloatval_t*g,int n,lbfgsfloatval_t s)
{
    return static_cast<OT*>(inst)->eval(x,g,n,s);
}
inline int OT::_prog(void*inst,const lbfgsfloatval_t*x,const lbfgsfloatval_t*g,
                     lbfgsfloatval_t fx,lbfgsfloatval_t xn,lbfgsfloatval_t gn,
                     lbfgsfloatval_t st,int n,int k,int ls)
{
    return static_cast<OT*>(inst)->prog(x,g,fx,xn,gn,st,n,k,ls);
}
inline lbfgsfloatval_t OT::eval(const lbfgsfloatval_t* x,
                                lbfgsfloatval_t*       g,
                                int n, lbfgsfloatval_t)
{
    w.assign(x,x+n);
    double mean = 0.0;
    for (double wi:w) mean += wi;
    mean /= n;                           
    for (double& wi:w) wi -= mean;      
    PD = PowerDiagram(p,w);
    PD.compute();
    lbfgsfloatval_t f = 0.0;
    for (int i=0;i<n;++i)
    {
        double Ai =PD.area(i);
        g[i]= Ai - lambda[i];
        f+= -w[i]*Ai + lambda[i]*w[i];
        const auto& V= PD.cells[i].v;
        if (V.size()>2)
        {
            const Vec2& a0 = V[0];
            auto sq = [&](const Vec2& q){
                double dx=q.x-p[i].x, dy=q.y-p[i].y;
                return dx*dx + dy*dy;
            };
            for (std::size_t k=1;k+1<V.size();++k)
            {
                const Vec2& a1 = V[k];
                const Vec2& a2 = V[k+1];
                double tri = std::fabs((a1.x-a0.x)*(a2.y-a0.y) -
                                       (a2.x-a0.x)*(a1.y-a0.y))*0.5;
                f += tri/6.0 * ( sq(a0)+sq(a1)+sq(a2)
                               + sq(a0)+sq(a1)+sq(a2) );
            }
        }
    }
    double gmean = 0.0;
    for (int i=0;i<n;++i) gmean += g[i];
    gmean /= n;                         
    for (int i=0;i<n;++i) g[i] -= gmean; 
    return -f;                         
}
inline int OT::prog(const lbfgsfloatval_t*,const lbfgsfloatval_t*,
                    lbfgsfloatval_t fx,lbfgsfloatval_t,
                    lbfgsfloatval_t gn,lbfgsfloatval_t,
                    int n,int k,int)
{
    if (k%20==0)
        std::printf("iter %4d   f = %-12g   |∇f|∞ = %.3e\n",k,fx,gn);
    return 0;
}
inline void OT::saveSVG(const char* file)
{
    FILE* f=std::fopen(file,"w");
    std::fprintf(f,"<svg viewBox='0 0 1000 1000'"
                   " xmlns='http://www.w3.org/2000/svg'>\n");
    for (const auto& poly:PD.cells)
    {
        std::fprintf(f,"<polygon points='");
        for (auto&q:poly.v)
            std::fprintf(f,"%.3f,%.3f ",q.x*1000,1000-q.y*1000);
        std::fprintf(f,"' fill='none' stroke='black'/>\n");
    }
    for (auto&s:p)
        std::fprintf(f,"<circle cx='%.3f' cy='%.3f' r='2' fill='red'/>\n",
                     s.x*1000,1000-s.y*1000);
    std::fprintf(f,"</svg>\n");
    std::fclose(f);
}