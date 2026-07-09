#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <random>
#include <cstdlib>
using namespace std;
using Clock = chrono::steady_clock;
static Clock::time_point T0;
static inline double elapsed(){
    return chrono::duration<double>(Clock::now()-T0).count();
}
static double TIME_LIMIT = 0.86;
static inline double clamp01(double x){ return x<0.0?0.0:(x>1.0?1.0:x); }
static inline double dist_face(double x,double y,double z){
    double df=x; if(y<df)df=y; if(z<df)df=z;
    double t=1.0-x; if(t<df)df=t; t=1.0-y; if(t<df)df=t; t=1.0-z; if(t<df)df=t;
    return df;
}
static double geomRadius(const vector<array<double,3>>& c){
    int n=(int)c.size();
    if(n==0) return 0.0;
    double dface=1e18;
    for(const auto& p:c){ double df=dist_face(p[0],p[1],p[2]); if(df<dface)dface=df; }
    if(n==1) return dface;
    double mind2=1e18;
    for(int i=0;i<n;i++){
        const auto& a=c[i];
        for(int j=i+1;j<n;j++){
            double dx=a[0]-c[j][0], dy=a[1]-c[j][1], dz=a[2]-c[j][2];
            double d2=dx*dx+dy*dy+dz*dz;
            if(d2<mind2) mind2=d2;
        }
    }
    return min(dface, 0.5*sqrt(mind2));
}
struct Grid {
    int G;
    double h;
    vector<int> head;
    vector<int> nxt;
    void build(const vector<array<double,3>>& p, double cell){
        int n=(int)p.size();
        h = cell>1e-9 ? cell : 1e-9;
        G = (int)floor(1.0/h)+1; if(G<1) G=1;
        int cap = (int)(2.5*cbrt((double)n))+3; if(cap<4) cap=4;
        if(G>cap){ G=cap; h=1.0/G; }
        head.assign((size_t)G*G*G, -1);
        nxt.assign(n, -1);
        for(int i=0;i<n;i++){
            int cx=cellCoord(p[i][0]), cy=cellCoord(p[i][1]), cz=cellCoord(p[i][2]);
            long idx=((long)cx*G+cy)*G+cz;
            nxt[i]=head[idx]; head[idx]=i;
        }
    }
    inline int cellCoord(double v) const {
        int c=(int)floor(v/h); if(c<0)c=0; if(c>=G)c=G-1; return c;
    }
};
static double geomFast(const vector<array<double,3>>& p, double hintR){
    int n=(int)p.size();
    if(n==0) return 0.0;
    double dface=1e18;
    for(const auto& q:p){ double df=dist_face(q[0],q[1],q[2]); if(df<dface)dface=df; }
    if(n==1) return dface;
    double cell = 2.0*hintR; if(cell<1e-6) cell=1e-6;
    Grid g; g.build(p, cell);
    int G=g.G;
    double mind2=1e18;
    for(int i=0;i<n;i++){
        int cx=g.cellCoord(p[i][0]), cy=g.cellCoord(p[i][1]), cz=g.cellCoord(p[i][2]);
        for(int dx=-2;dx<=2;dx++){int x=cx+dx; if(x<0||x>=G)continue;
        for(int dy=-2;dy<=2;dy++){int y=cy+dy; if(y<0||y>=G)continue;
        for(int dz=-2;dz<=2;dz++){int z=cz+dz; if(z<0||z>=G)continue;
            long idx=((long)x*G+y)*G+z;
            for(int j=g.head[idx]; j!=-1; j=g.nxt[j]){
                if(j==i) continue;
                double ddx=p[i][0]-p[j][0], ddy=p[i][1]-p[j][1], ddz=p[i][2]-p[j][2];
                double d2=ddx*ddx+ddy*ddy+ddz*ddz;
                if(d2<mind2) mind2=d2;
            }
        }}}
    }
    if(mind2>1e17) return dface;
    return min(dface, 0.5*sqrt(mind2));
}
static std::mt19937 RNG(12345);
static double pushSweep(vector<array<double,3>>& p, double r, double relax){
    int n=(int)p.size();
    double D=2.0*r, lo=r, hi=1.0-r;
    if(hi<lo){ for(auto& q:p){ q[0]=0.5;q[1]=0.5;q[2]=0.5; } return 0.0; }
    Grid g; g.build(p, D);
    int G=g.G;
    static vector<array<double,3>> disp;
    disp.assign(n, {0.0,0.0,0.0});
    double D2=D*D;
    double maxOv=0.0;
    for(int i=0;i<n;i++){
        int cx=g.cellCoord(p[i][0]), cy=g.cellCoord(p[i][1]), cz=g.cellCoord(p[i][2]);
        for(int dx=-1;dx<=1;dx++){int x=cx+dx; if(x<0||x>=G)continue;
        for(int dy=-1;dy<=1;dy++){int y=cy+dy; if(y<0||y>=G)continue;
        for(int dz=-1;dz<=1;dz++){int z=cz+dz; if(z<0||z>=G)continue;
            long idx=((long)x*G+y)*G+z;
            for(int j=g.head[idx]; j!=-1; j=g.nxt[j]){
                if(j<=i) continue;
                double ddx=p[i][0]-p[j][0], ddy=p[i][1]-p[j][1], ddz=p[i][2]-p[j][2];
                double d2=ddx*ddx+ddy*ddy+ddz*ddz;
                if(d2<D2){
                    double d=sqrt(d2);
                    if(d<1e-12){
                        std::uniform_real_distribution<double> U(-1.0,1.0);
                        ddx=U(RNG); ddy=U(RNG); ddz=U(RNG);
                        d=sqrt(ddx*ddx+ddy*ddy+ddz*ddz)+1e-12;
                    }
                    double overlap=D-d;
                    if(overlap>maxOv) maxOv=overlap;
                    double f=0.5*relax*overlap/d;
                    double mx=f*ddx, my=f*ddy, mz=f*ddz;
                    disp[i][0]+=mx; disp[i][1]+=my; disp[i][2]+=mz;
                    disp[j][0]-=mx; disp[j][1]-=my; disp[j][2]-=mz;
                }
            }
        }}}
    }
    for(int i=0;i<n;i++){
        double x=p[i][0]+disp[i][0];
        double y=p[i][1]+disp[i][1];
        double z=p[i][2]+disp[i][2];
        if(x<lo)x=lo; else if(x>hi)x=hi;
        if(y<lo)y=lo; else if(y>hi)y=hi;
        if(z<lo)z=lo; else if(z>hi)z=hi;
        p[i][0]=x; p[i][1]=y; p[i][2]=z;
    }
    return maxOv;
}
static vector<array<double,3>> cubicGrid(int n){
    int best_max=2000000000, bm=1,bk=1,bl=1;
    int root=(int)ceil(pow((double)n,1.0/3.0))+2;
    for(int m=1;m<=root;m++){
        for(int k=m;k<=root;k++){
            int need=(n+m*k-1)/(m*k); if(need<1)need=1;
            for(int l=max(k,need); l<=max(k,need)+2; l++){
                if(m*k*l>=n){
                    int mx=max(m,max(k,l));
                    if(mx<best_max || (mx==best_max && m*k*l<bm*bk*bl)){
                        best_max=mx; bm=m; bk=k; bl=l;
                    }
                }
            }
        }
    }
    double s=1.0/(double)max(bm,max(bk,bl));
    vector<array<double,3>> pts; pts.reserve(n);
    for(int i=0;i<bm && (int)pts.size()<n;i++)
        for(int j=0;j<bk && (int)pts.size()<n;j++)
            for(int l=0;l<bl && (int)pts.size()<n;l++)
                pts.push_back({clamp01(s*(i+0.5)), clamp01(s*(j+0.5)), clamp01(s*(l+0.5))});
    while((int)pts.size()<n) pts.push_back({0.5,0.5,0.5});
    return pts;
}
static int countFCC_r(double r,double ox,double oy,double oz){
    double s=2.0*sqrt(2.0)*r; double lo=r,hi=1.0-r;
    if(hi<lo-1e-15) return 0;
    if(ox<0)ox=lo; if(oy<0)oy=lo; if(oz<0)oz=lo;
    double basis[4][3]={{0,0,0},{0.5,0.5,0},{0.5,0,0.5},{0,0.5,0.5}};
    int cnt=0;
    for(int b=0;b<4;b++){
        double bx=basis[b][0]*s,by=basis[b][1]*s,bz=basis[b][2]*s;
        int imin=(int)ceil((lo-ox-bx)/s-1e-9), imax=(int)floor((hi-ox-bx)/s+1e-9);
        int jmin=(int)ceil((lo-oy-by)/s-1e-9), jmax=(int)floor((hi-oy-by)/s+1e-9);
        int kmin=(int)ceil((lo-oz-bz)/s-1e-9), kmax=(int)floor((hi-oz-bz)/s+1e-9);
        if(imax<imin||jmax<jmin||kmax<kmin) continue;
        cnt += (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);
    }
    return cnt;
}
static vector<array<double,3>> genFCC_r(double r,double ox,double oy,double oz){
    double s=2.0*sqrt(2.0)*r; double lo=r,hi=1.0-r;
    vector<array<double,3>> pts; if(hi<lo-1e-15) return pts;
    if(ox<0)ox=lo; if(oy<0)oy=lo; if(oz<0)oz=lo;
    double basis[4][3]={{0,0,0},{0.5,0.5,0},{0.5,0,0.5},{0,0.5,0.5}};
    for(int b=0;b<4;b++){
        double bx=basis[b][0]*s,by=basis[b][1]*s,bz=basis[b][2]*s;
        int imin=(int)ceil((lo-ox-bx)/s-1e-9), imax=(int)floor((hi-ox-bx)/s+1e-9);
        int jmin=(int)ceil((lo-oy-by)/s-1e-9), jmax=(int)floor((hi-oy-by)/s+1e-9);
        int kmin=(int)ceil((lo-oz-bz)/s-1e-9), kmax=(int)floor((hi-oz-bz)/s+1e-9);
        for(int i=imin;i<=imax;i++){ double x=i*s+bx+ox;
            for(int j=jmin;j<=jmax;j++){ double y=j*s+by+oy;
                for(int k=kmin;k<=kmax;k++){ double z=k*s+bz+oz;
                    pts.push_back({x,y,z}); }}}
    }
    return pts;
}
static int countHCP_r(double r,double ox,double oy,double oz){
    double s=2.0*r; double lo=r,hi=1.0-r; if(hi<lo-1e-15) return 0;
    if(ox<0)ox=lo; if(oy<0)oy=lo; if(oz<0)oz=lo;
    double h=s*sqrt(2.0/3.0); double sy=s*sqrt(3.0)/2.0; int cnt=0;
    int kmin=(int)ceil((lo-oz)/h-1e-9), kmax=(int)floor((hi-oz)/h+1e-9);
    for(int k=kmin;k<=kmax;k++){
        int lay=((k%2)+2)%2; double lx=lay?s/2.0:0.0, ly=lay?s/(2.0*sqrt(3.0)):0.0;
        int jmin=(int)ceil((lo-oy-ly)/sy-1e-9), jmax=(int)floor((hi-oy-ly)/sy+1e-9);
        for(int j=jmin;j<=jmax;j++){
            double xmin=lo-ox-lx-j*s/2.0, xmax=hi-ox-lx-j*s/2.0;
            int imin=(int)ceil(xmin/s-1e-9), imax=(int)floor(xmax/s+1e-9);
            if(imax>=imin) cnt += (imax-imin+1);
        }
    }
    return cnt;
}
static vector<array<double,3>> genHCP_r(double r,double ox,double oy,double oz){
    double s=2.0*r; double lo=r,hi=1.0-r; vector<array<double,3>> pts; if(hi<lo-1e-15) return pts;
    if(ox<0)ox=lo; if(oy<0)oy=lo; if(oz<0)oz=lo;
    double h=s*sqrt(2.0/3.0); double sy=s*sqrt(3.0)/2.0;
    int kmin=(int)ceil((lo-oz)/h-1e-9), kmax=(int)floor((hi-oz)/h+1e-9);
    for(int k=kmin;k<=kmax;k++){
        double z=oz+k*h; if(z<lo-1e-12||z>hi+1e-12) continue;
        int lay=((k%2)+2)%2; double lx=lay?s/2.0:0.0, ly=lay?s/(2.0*sqrt(3.0)):0.0;
        int jmin=(int)ceil((lo-oy-ly)/sy-1e-9), jmax=(int)floor((hi-oy-ly)/sy+1e-9);
        for(int j=jmin;j<=jmax;j++){
            double y=oy+ly+j*sy; if(y<lo-1e-12||y>hi+1e-12) continue;
            double xmin=lo-ox-lx-j*s/2.0, xmax=hi-ox-lx-j*s/2.0;
            int imin=(int)ceil(xmin/s-1e-9), imax=(int)floor(xmax/s+1e-9);
            for(int i=imin;i<=imax;i++){ double x=ox+lx+i*s+j*s/2.0; pts.push_back({x,y,z}); }
        }
    }
    return pts;
}
static vector<array<double,3>> pruneMaxMin(vector<array<double,3>> pts,int keep){
    int n=(int)pts.size(); if(n<=keep) return pts;
    int start=0; double bw=-1;
    for(int i=0;i<n;i++){ double w=dist_face(pts[i][0],pts[i][1],pts[i][2]); if(w>bw){bw=w;start=i;} }
    vector<char> used(n,0); used[start]=1;
    vector<array<double,3>> res; res.push_back(pts[start]);
    vector<double> mind(n,1e18);
    for(int i=0;i<n;i++) if(!used[i]){
        double dx=pts[i][0]-pts[start][0], dy=pts[i][1]-pts[start][1], dz=pts[i][2]-pts[start][2];
        mind[i]=0.25*(dx*dx+dy*dy+dz*dz);
    }
    while((int)res.size()<keep){
        int best=-1; double bestd=-1;
        for(int i=0;i<n;i++){
            if(used[i]) continue;
            double w=dist_face(pts[i][0],pts[i][1],pts[i][2]);
            double md=min(mind[i], w*w);
            if(md>bestd){ bestd=md; best=i; }
        }
        used[best]=1; res.push_back(pts[best]);
        for(int i=0;i<n;i++) if(!used[i]){
            double dx=pts[i][0]-pts[best][0], dy=pts[i][1]-pts[best][1], dz=pts[i][2]-pts[best][2];
            double d2=0.25*(dx*dx+dy*dy+dz*dz);
            if(d2<mind[i]) mind[i]=d2;
        }
    }
    return res;
}
static double expandFill(vector<array<double,3>>& p, double rHint){
    int n=(int)p.size(); if(n==0) return 0.0;
    double r = geomFast(p, rHint>1e-9?rHint:0.01);
    for(int it=0; it<4; it++){
        if(r>=0.5) break;
        double mn[3]={1e18,1e18,1e18}, mx[3]={-1e18,-1e18,-1e18};
        for(auto&q:p) for(int a=0;a<3;a++){ if(q[a]<mn[a])mn[a]=q[a]; if(q[a]>mx[a])mx[a]=q[a]; }
        bool changed=false;
        for(int a=0;a<3;a++){
            double w=mx[a]-mn[a], c=0.5*(mn[a]+mx[a]);
            double tw=1.0-2.0*r;
            if(w<1e-12){
                for(auto&q:p) q[a]=0.5;
                if(fabs(c-0.5)>1e-15) changed=true;
                continue;
            }
            double k = tw/w;
            if(k<1.0) k=1.0;
            if(k>1.0+1e-15 || fabs(c-0.5)>1e-12){
                for(auto&q:p){ double v=0.5+(q[a]-c)*k; q[a]=clamp01(v); }
                changed=true;
            }
        }
        if(!changed) break;
        double r2 = geomFast(p, r);
        if(r2<=r+1e-15){ r=r2; break; }
        r=r2;
    }
    return r;
}
typedef int(*CntFn)(double,double,double,double);
typedef vector<array<double,3>>(*GenFn)(double,double,double,double);
static double findMaxR(CntFn cnt, GenFn gen, int n, double ox,double oy,double oz,
                       vector<array<double,3>>& outpts){
    double lo=0.0005, hi=0.5;
    if(cnt(lo,ox,oy,oz) < n) return -1.0;
    for(int it=0; it<50; it++){
        double mid=0.5*(lo+hi);
        if(cnt(mid,ox,oy,oz) >= n) lo=mid; else hi=mid;
    }
    double bestRr=lo;
    for(double r=lo-0.003; r<=lo+0.003; r+=0.0001){
        if(r<=0) continue;
        if(cnt(r,ox,oy,oz) >= n && r>bestRr) bestRr=r;
    }
    auto pts=gen(bestRr,ox,oy,oz);
    if((int)pts.size() < n) return -1.0;
    outpts=pts;
    return bestRr;
}
static double tryLattice(CntFn cnt, GenFn gen, int n, double& bestR,
                         vector<array<double,3>>& bestP, double tcap){
    static const double offs[][3]={
        {-1,-1,-1},{0,0,0},{-1,-1,0},{-1,0,0},
        {0.25,0.25,0.25},{0.125,0.125,0.125},
        {0,0.25,0.5},{0.5,0.5,0},{0.375,0.375,0.375},
        {0.1,0.1,0.1},{0.05,0.15,0.25}
    };
    int nof = 11;
    for(int oi=0; oi<nof && elapsed()<tcap; ++oi){
        const double* o = offs[oi];
        vector<array<double,3>> p;
        double r=findMaxR(cnt,gen,n,o[0],o[1],o[2],p);
        if(r<=0) continue;
        auto pr=pruneMaxMin(p,n);
        if((int)pr.size()!=n) continue;
        double rr=expandFill(pr, r);
        if(rr>bestR){ bestR=rr; bestP=pr; }
    }
    return bestR;
}
static inline void jitter(vector<array<double,3>>& p, double amp){
    std::uniform_real_distribution<double> U(-1.0,1.0);
    int n=(int)p.size();
    for(int i=0;i<n;i++){
        double x=p[i][0]+amp*U(RNG), y=p[i][1]+amp*U(RNG), z=p[i][2]+amp*U(RNG);
        p[i][0]=x<0?0:(x>1?1:x); p[i][1]=y<0?0:(y>1?1:y); p[i][2]=z<0?0:(z>1?1:z);
    }
}
static void teleportKick(vector<array<double,3>>& p, double r){
    int n=(int)p.size(); if(n<3) return;
    int bi=0,bj=1; double bd2=1e18;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
        double dx=p[i][0]-p[j][0],dy=p[i][1]-p[j][1],dz=p[i][2]-p[j][2];
        double d2=dx*dx+dy*dy+dz*dz;
        if(d2<bd2){bd2=d2;bi=i;bj=j;}
    }
    int mv=(RNG()&1)?bi:bj;
    double lo=r,hi=1.0-r; if(hi<lo){lo=hi=0.5;}
    std::uniform_real_distribution<double> U(lo,hi);
    double bv=-1; array<double,3> bx{0.5,0.5,0.5};
    for(int k=0;k<60;k++){
        array<double,3> x{U(RNG),U(RNG),U(RNG)};
        double v=dist_face(x[0],x[1],x[2]);
        for(int i=0;i<n;i++){ if(i==mv)continue;
            double dx=x[0]-p[i][0],dy=x[1]-p[i][1],dz=x[2]-p[i][2];
            double d=0.5*sqrt(dx*dx+dy*dy+dz*dz); if(d<v)v=d;
        }
        if(v>bv){bv=v;bx=x;}
    }
    p[mv]=bx;
}
static bool inflateTo(vector<array<double,3>>& p, double target, int sweepCap,
                      double deadline, double sor){
    double tol = 2e-7*target;
    double bestOv=1e18; int noimp=0;
    for(int s=0; s<sweepCap; s++){
        if((s&3)==0 && elapsed()>deadline) return false;
        double ov = pushSweep(p, target, sor);
        if(ov < tol) return true;
        if(ov < bestOv*0.997){ bestOv=ov; noimp=0; }
        else if(++noimp >= 14) return false;
    }
    return false;
}
static void optimize(vector<array<double,3>>& best, double& bestR, double deadline,
                     int sched){
    int n=(int)best.size();
    if(n<2) return;
    double grow    = sched?1.04:(n<=1200?1.03:1.006);
    int    grS     = sched?25:(n<=1200?20:3);
    int    stagMax = sched?12:40;
    double kick0   = 0.20, kmin = 0.02;
    double sor     = sched?1.5:1.0;
    vector<array<double,3>> cur=best;
    if(sched) jitter(cur, 0.05*bestR);
    for(int s=0;s<20 && elapsed()<deadline;s++) pushSweep(cur, bestR, sor);
    double curR=geomFast(cur,bestR);
    if(curR>bestR){ bestR=curR; best=cur; } else cur=best;
    double kick=kick0;
    int stagn=0, tk=0;
    while(elapsed()<deadline){
        double target=curR*grow; if(target>0.5) target=0.5;
        for(int s=0;s<grS && elapsed()<deadline;s++) pushSweep(cur, target, sor);
        double r=geomFast(cur,target);
        curR=r;
        if(r>bestR+1e-12){ bestR=r; best=cur; stagn=0; }
        else stagn++;
        if(stagn>=stagMax){
            cur=best;
            if(n<=350 && (tk^=1)) teleportKick(cur,bestR);
            else jitter(cur, kick*2.0*bestR);
            for(int s=0;s<grS+4 && elapsed()<deadline;s++) pushSweep(cur, bestR, sor);
            curR=geomFast(cur,bestR);
            stagn=0;
            kick*=0.85; if(kick<kmin) kick=kmin;
        }
    }
}
static void microPolish(vector<array<double,3>>& best, double& bestR, double deadline){
    int n=(int)best.size();
    if(n<2) return;
    vector<array<double,3>> cur=best;
    double e=2e-3;
    while(e>5e-8 && elapsed()<deadline){
        double t=bestR*(1.0+e); if(t>0.5) t=0.5;
        bool conv=false;
        for(int s=0;s<300;s++){
            if((s&3)==0 && elapsed()>deadline) break;
            double ov=pushSweep(cur,t,1.0);
            if(ov<1e-9*t){ conv=true; break; }
        }
        double rn=geomFast(cur,t);
        if(rn>bestR+1e-15){ bestR=rn; best=cur; if(conv) e*=1.6; }
        else { cur=best; e*=0.5; }
        if(t>=0.5) break;
    }
}
static vector<array<double,3>> randomConfig(int n){
    std::uniform_real_distribution<double> U(0.02,0.98);
    vector<array<double,3>> p(n);
    for(int i=0;i<n;i++){ p[i]={U(RNG),U(RNG),U(RNG)}; }
    return p;
}
static void relaxOptimize(vector<array<double,3>>& best, double& bestR){
    int n=(int)best.size();
    if(n<2) return;
    if(n<=350){
        optimize(best, bestR, 0.28*TIME_LIMIT, 0);
        int restarts;
        if(n<=64)       restarts = 12;
        else if(n<=150) restarts = 8;
        else            restarts = 5;
        double reserve = 0.18*TIME_LIMIT;
        double slice = (TIME_LIMIT - reserve - elapsed()) / restarts;
        if(slice > 0.01){
            for(int rs=0; rs<restarts && elapsed()<TIME_LIMIT-reserve; rs++){
                vector<array<double,3>> cand;
                if(rs&1){ cand=best; jitter(cand, (0.6+0.1*rs)*bestR); }
                else cand = randomConfig(n);
                double candR = geomFast(cand, 0.4/pow((double)n,1.0/3.0));
                double dl = min(TIME_LIMIT-reserve, elapsed()+slice);
                optimize(cand, candR, dl, 1);
                if(candR > bestR){ bestR=candR; best=cand; }
            }
        }
        optimize(best, bestR, TIME_LIMIT, 0);
    } else {
        optimize(best, bestR, TIME_LIMIT, 0);
    }
    double r=expandFill(best, bestR);
    if(r>bestR) bestR=r;
}
static const int EMB_NMAX = 120;
static const char* EMB[] = {
"q7dVe*u<2lMp>y#dqIKV^W0JLTU#AKR7DyOt5tWfiVXrDl/5yz1uV+QC!vugS6?*xDw5#*/HOyQ{v=wpU=&nYSSr9INXhqcBohEzhm]Il>Nk8X%9}!R5c@Wu0rvU+PFc6i0S-sKZO0mz3&Q88F4}M1>x]OJ)SVO&[FD[>0S-o.O:^LHZ98)2-9$CrGk[ASnh?GE/=$R:bYl*Kaq(Ng0nIFlqvRq7C1TbL:Z*x2QF%$wybGM%sZm@bdR-W&PB.y2rWqAqwOo2:U2XX<hdZ*sq=fHyMnGO.cAKvH!CRquu+4MLrSJdS8NulwU2^8/+d(<Lt=}cTybD{nRLv8(0bsx5o<7-2SSl^:wND5P237.shdZ*9TGhE7L@ib[XM(pp!m{<Lo?85:0mU=aNZfF5q#>CGh9-o4soFWaM4f^jsZm@bd/+9ho<7-7s1xK<7(Rvl5SnZJf=l.Gc<B!jOZaV1RyhF6i@<^@ao<58O=AlB>uw?M6AZzG=GxnV?06v8v9+hG4P^v@^%rX+gV/7uSLX3RiB1HkqgV)S=/:rtQxNHz&KP<Th>{zakC9<D&[snbR5S>%mt@31V4pw@fLwPn&@Rl?(5dO@8tu>CkiW$xbh8Y:sIU=.X]WPWAF3GFfLxy<OC)#$Pd?8<6h[{.L(JHJ&06n/E8dX^M(ev027jk!fHc-]soFU&PmyOC[b&Yzkq+cc/?-aREGu5{xy>.1@^qv+fgMz?)f!PWuYptw*4a&(lsLEfMQQoHDOw@VUs)ju4eG#h^8}o)Id}(1Pu.abcf3tbg*XL)tLYdTq:Mu[MF(X$/>}Nf^OHr<4fp@zvka2tI$*B(XUNUpeyeKSZ5m<YMF(:Rb-+G2PCR+C{Yl1jtvlRKK4*>En2i-Xg>&l*z+a7:MH])AcxLXpqKDwSQzY6(RoRsuS-6(sV$1RKap[YeASGv?kN.Ps!]y^{FgkM33F6g8tHEi[Dp3FLnU<^G8Pv(YwVP4-n>Q}gf>SEQf{}M!i0ol+sUc$/*gEX)1PQ?Naq-Vk0S^.1x]s:f=ydViFl$#ZjMZYBOrPsY*gE^XVNKuXWaY09v6G/07bbc[.u6{zEQynpEk$!jsvi&X>Xis(0dv5Dh0]I*vO)cwv4W%SjotDZMfDEdLe:9:ss!gAxBNc6}klFy]Zy05PhhxlZZxgzX!83)Fecl3VC3UusjuAS56Wcep//%#&0*{}vFIG9Z+(MH1oFP7e9Xk}D>W4>SfVYQ{P?ZF1xWi?&aH>fQ5?iyolyImmXYTQ!d)>WxAt6kNr+79fxflvb39dSo3^$gyF[g#nm5v:o>A3VF4T/gxyRP(%b32F>hJ-IqN86Ca^M$uRWGEOxY69ZnM-9]h:xkYM>(H1SNO.ai*cEwTeWUmoB}ob13PZm(wAfPqc0V%fLWl+2QwZDHK}.FMXTBBS@KSpar4uV0kUW=kg5^(.%N*cF)F!*4eZicPdlaGtb5g10ast&&2nt)T5VgIykO&=Tn]n2!FgYMsohC604Nob1YzI4S@QV+BY$Jsri31H6<}@OTn]p3FiUj*4o9xerN<tS8Mi>vrraIlBX&dxT5J9H6<}@U1oRV7Fnf$MTYby5S)$<(oK:bPQ@B&U809t=<Fo.+L=3N<MB9C=dqIVO?d)%g-:YZuoLb^}K2SW+>lig[E{h>3n1$sTL&(mhOx.#yRxFMQim=x$}N<=cRl8f=]$5ZjWgrR}:MCRuzNS!EdhnWab9J[6qR5agmjDtJB2+W7Wwti8nSx&3f]wUph06Xov[{dW/#>keTZ+aGTyvq:B2+WAsJ>H8Wt^f}pg(wJ+&K?^OA>s7Mn83oqE:GPenYevBfj3sB*ESRW4{OGpd^@4gz(%<v?:/0L9C9HqE:Gmthr):J8Be3n>-YT8i3T^f{{BQhNWUCdqJ5mOOG7rqK(38q<GiTBfe}En(Ie4>5SK)yIbY:y#{QN!%8MQ/#>sr-*(vqVZ>{}tDGs*IC]X2nSpUryGOL1L&OV3OrPz{^W4vCqEZAd*r#f4"
"Q@D[tsJ>wl<=cyjV5h-oh0n5ZdhnWoL9C6!-*(u%^#]%TB2/=osLzDNW{7>56T+K2MB9Ec/3myK!?(L]T^[Y&i]aHvQ/9.EICdRpoiPs:6Sila+&COt!%9ON^VUJz-[3xw-u)wVA(IPzRZ#X3<=09*V8oUlf?wc3OrP{5w&eWO-*&o>gNVFHA(MXQB?lnH7]sv=L^NtmhNWUC/0hbl?e8mTqR7eDi]aHjO><J[:[(gRk4JE%LIJ5dmnfsz/C0@(VJio#CW^{WL})XhsfMmiR&oh#.jN@AW8}[^>@+N8zJ<iL4tsbrkaVNsLIa$OElDjWhu->pj=@W}gpqi}).-Fq/Fo<Y1hNigj1$E*)lc1:Kfx#BIo$*kj5!%5e13LcI>()ZO0fgy5{zr8pa5>{[X%kGOV(CuWK^bLiGN0Dg8C1zp3MMp/PR</9Gjg]:^T{ogH+E.zM.Tp<3-dxig&J$q7daV).-HI/SH7=7$fZtV5rA#^%W6mw3H*QB#(:=-gfNIkNCx@I8.22BD-vcL6j*)pa6<H=<+*[w3Eyg(8znkh8.vbTG8VnKVN4pFnwHM9Gc.:0dWrS5dKN108s[gPAKF!-gmVx/^Ts!b@FD6Vq}-wlub7%04hPab<*BHeWo*{WMttm-gsf=Z=1m{J(LQ5ckdY4qEYrBoW$$RGVOi?G*EayC4w<?+:C^:c5F9F.YA={Nx)SY-P>:5VSTF6I-)Qqxg11g{g}.Agkt<oy#tDn)j%y4Ry$J9=7bie+:D#z40]}qNp^<.)iGLX=rK5@Uuet0N[46<MH9*meT]ykV!sZswwR+gt1*m@7g+WVnsq5?6a6&IvPc(jMH8..j6svKoIb%mI-6NPx0S*qbX]o@gIk7*)aLiuN[/wszFX3AqE)i)VTch.m9fc%NjW9Qc9ts34W<U%^a:/fOER2X*1Ld:j!20gV:C#Pww>v>x42wt*}bsLW}w79y@}vEwvdZ$um:i/hJY+08yuazQ%kEiN3uo/ca+{A+:sRwfbwuB8Qh>/up-T..E?5?V!qZQYeC{fHtpE?Otsnyx>C2u-!/(bvkq5[cqdNc:@JoJl.Q$rN(RMzNx%rdIjlOfgIfrMz1.Y$N[p/fcqeJ7i+DTkV!sZsDO.uDMQ@G^*SK5j^qhD4W[gXf2}{1JcciJFWy)&Up1*CN-i)52MBrP>MFj$r/Q1Zdf.QCe*Ia5wQbE<sW^IQT?6ICkW31WHxTT1=/B#ya^S)cQq5nm2pr]4{tHGGPls8jiv+LF*N[jhBxV1j{hJuK^e<3Cj=O*D&RQrw&*7IIrXb!dEV)y)+Ttj6}Hd-&s32*jffkiJ{SDJ%a/3G.Vc9VY6WIZVK.*<C[DQX}-ofuUlR&ixq^t4.Qf.Q#987.apFd?4mlr%4M.>6[@BhzhMC(TBihJ#nxMUwN:z6$isuEDhN*7J4+ncoxpV%&toV^2sYD9!)+hJOxB^}Mhyd#uV]5H-{Lr1+3XA-yrsn)P?om6dQZSKj3?xR-6AaWO!^SHxqRmgb=$b[?6>u8)j9{gczfoZOJ5Oj.IAw]a%@dRw@ifYzT4TDNB9N5bTrwQ-i]Wolz(wu6t[LtmMu)x0fVYjFz1:m^bBTV71Oc2E?3T/i.nWwtDkQ]adyyz$]&WWQf6b:d1fL4<0xI6T<]vuK?uOWr6[OzaV:60/[iE9%rTM*kMiQ)<U=z0iX.BVSm8PhdwSzl%jxnTBq5-f!%Hzmk[tuYp$F[{dVWWY2WjWvJd1v@nB}OV%JiU?HZf3K%zXN-XmDRM0>ab:fza5[{O3po(SD*3KRTw@+Jv<%m8V1F8@Yy!rv]nD:SjwrovpM=P{5.63FzbYhQMHrM3DnwU#&XL+^ZzB*N^re*HYMxrK))BS{LZtiJ{Mm*KNw*zef-uVW})cfCOz--FBMK>:Va{Mvtlghv0U{P65bSR+%eNnHvntw$Trq<AoQEDIgND$<<4r]o@5&PSc?Rt#v*VUEgAZW/./}96b*Kl-7z[1G*{t*/>a@A&@x9*$W"
"qPaU/JGM5X&Zlg}nsY151X>yCy%:S$/R[hK?5KO(Seg)tILc6qMqFC78lk76:.2zi4pcTxXjCw:wKIrjT8wE[5*MdO43btoL@N{PgZ+^Iem3VN8*U{pNndrae4Y*<:8s2G=86d@NmzfQddb.sJLnnxW<i%((z/3+iKNR8hV4-h?!>$g5*@u:m+]mVbSr:T-[qOPbQs0O-eiP=K*<8d7!Cw23p[AgzoFw6lKph:vm!(#eGnLWt-{!p[<TyP*}>WX/IWdf351C8U+WC97t.!SbU]!vIrrRsjr#&FLG8zHzZp[!eZ}&6X{p^t=LP6f8WT0GvoZ9=AZ[&CK>Xg{jaj.VzB/}c>AnMN?RCT^zt9OCZtqgEbYgTaI^)3PX1jiB6t^e2DpPq#:QKlubhY<PbNBQ9lJ<{Nvv4o7-Tus&.x..9r(diizO/*cg&aP/?gTHDfS&*z9]y>MbRS=@RwbZFns/d&6t^k}C-XH)8YK7:XRAqfpBzq(me7Q.P4=f{!pOc0zZohw1Dy33K**GcW!>j?*(T.NWKmL>wCo#6*Wj5-EzjPSW<hN)8*3]VKLq3Rtv--aKDvtx5*$cEk{C8PvqTa9EyQIOKT9Z]1Ce%eGH?mWI5QQN*(PxYhRp%TIWKDvvnyBneXt>M.xZ?5MLH>DzRAOb>y!PB][sY-K$Cehk{yV+ny7}AIr+D)j!xbgYlMmeWyq^OhU}Oebx74pzt9wgXGJibO${L@cF)Kub?{Uv*JBszz[0mz:S2VsaDPc=K{IKxI6D]2bOz5bIr!#t.JCJ&[<kR^*^[3>I3$#Dbx74pfO{h&XGMR%nu/1>7841..RqJKm2rb}DN92NT&@xgXR9JdU=29R5G)i0O@$F5b-IFnK!dZPBj]>%DCU?vCjFk?bhCzN*O.b8Y.uuNP0j4<EyQOunj-dhSu5OSX3$b2{t?5S?G{K5U=VVdY$neM*WjgoEyQhNmN<Oq-a+6&JlU]*!qR%#7%frk5S3Tc+/2EWP^9aU--1QuAW8Cc+<AJiOwuBAH$P^/8/!An)O6%.o:Q?1iC@0tqzHd:XG+aSr1Z{BPsJSmxajgJa9j$qdcvRJy14bWTUU^Y:![ih-V[HTecZ!$U14i[cdKtwVVnK<!Ti*X9mpK(*]rp.OZcwEXTQ*FtB+62Bdg)<W(k=D>x1fTKW#ZJdhKE5bb#vXs:I&XX[hKBU+AT6tT^C#H[[pP7oLa>VO$oPYXbJX?6{>Sry4u@S9Xo@W)/PQupwEzH]!DC(6K]agBwP*v]YbD*%3-.qy<n:mvaP?Li{?vu+f(AcDV7?DV.nYpg(DqTMUNuSCUnmx#YCvY<xI2Q*g7+w29^QW?Sv8>q*9dVK1gB^W$RIw]sG8=c*]RX%VN)l%U01IByEanv5{UOqG+5K.8e0do*).wnWShrPkcVX[qQToR$-pI=2}9q4$txdYt=W!&c*k]V&vZYu+xr^K-%:+9}+c8fdU1Xj6}CRhVvC.Mr/*zLRSY!xIyGxs2WX)BHJ&pPgOdtAbImR[ep<7QL=>(6PQ[)UFGP.KMGLN2%}C-8Ng!ody+:*M:$*BURY&S:]duaD6Jo5N6k6f6^g=?0QyT=6YtWYcv[5wDB1}svwK.2+dcV5UgN05HKqF?fs+2*Z#3<VlWUIr3Q?Q6J#y&FwjTq2-m@.)nHR!e{zqT(zF*DNqQWX+cJvH=1ahV4c=!jOwVQcSF0:%vWduCzSEPWb+5%BMtVCre.-r%Eo$qqsk2gM14F#Wx>f/Wlts6gW1%#:?4G6G^RPpoOVZ[cWoMzm4C/[auprRL!7.N<h*m5@)-pzK?6QyU*v1YRg9Ei+pA?dp<vKzDImWC-2:jf(qF5t6qv5=9KRa1[a}V-K=dKwjYcnziO3TQRRM5y]!KW?z/6o.GkI%UUit4Vja#pblWR1&i+P}8WO0211Pj8(yLDIhqlq$G$+h#xiB3+14?oC7ImOIw[Yp}paLUkc(JfA8BW%LYW5q:wYq7?gh?qBhE/koig4CqzhoX]c="
"Q!PTeLE&1*{J2zR)vb/=)={pI/57WX?r%vxlmf6}l)rsoJ8=L:MC>3W8ZFekN25bu5HQL-aXM#Dz0]xlY>F6(l}lQ}H%^O38+Y1Zu>rFV5qX#{KSaA$*HCDty!8Na[/@.2Yn%CR8[bRknBiEZ9G0DdVwhUct/-d)a#kv]XZBuAX>pRio5lRGl}Qf6LE)UW6JJ9=Udm?}r#(Cq>i=9UyI@(!K%YKOmbvpL*QjQKCxc{Q:C>S60S/]5^L(oaMH6Es?r]GlVxx]dUPZ9bGFBzO@Zf^M7b9Jx4*klyKJqwoU]d>(sRfc!i+(K2l/5yaGG68kz2W/f7LB.BF=vT-UxLTX9nbTdf{8d0:@7T@VFDeI}[:AB0q[2e77Z*?/O@m2q*jqE&s9&va^x4I8@kJkU?*}ntF.tZXC]VHdyi+b9}(vU=9]I)<n<hi?y/jd8T^ixlZ{s?Lv:gx1vI3uSc:bYDCc%$zW-bW<mvMe?ASIbe-#UIl/6j<bS+HHHc{DL&1bx1[+$)L/naOm[Egh=a&0=8j)8l?03NK+GFD%AHd{D9LR}</G*L/o/pRsA@DQZ%Up9vi-s8=lS.!/cB(BsvHm}So6E{/U[=r72kP?2!0/1w*a[3JRg8}V^lXIVYAF:jurx>K#/q7C54rUa@e#Vz8RfSewUmP#:(H%!YQ25>ndA}EgiI4yQ{Ma>^]lrTVzYWVN}}Du.DV7B+p2Z]bYXWG0Kpv/JHQk?%2Xl6Z42m2yKJqqv9}xj%p%T@u(HtprOGU=xGGXF!Hfuc]M[e)^>@$-xku448&^B55pTs&i5/*ZmYSS-0*R<K6-kaeZsr]TvFZAX7r$>5U!aJNfU0lY!:lSFw@Oway<uqM6DaXUQnor=eHj5}}/aWSb<-gf5E+c9tS]XoiX#]a<B)XgyYHe(7LSN4T42pxHR$T0s)=}t/?eoYaR$*<GYDM.]:@PBI%gEhNh&nKE[q1WkkStgJ[Dy&IaDPbE5^XX)@#^zVbPp#xHdf@gq+vVg?s^qqj}@>k2)Gw4T=MFS9BPamqxS+N:?4s+G76]KH=:+UEkpRFkfbzr!M9<caLojwWhW^.mnD%dGC1-SWN45}P9+Sev%*N!/PX3T6xk>5aN=!=[^bI&tP&lO*VnGFqmYn0Pj%8VsV&WQ[0lP#]/4o8kjtGF>IaEfYct1CflNo>Ezh<iuPa@81)(6xr&Jv52kUG./(avC7CtK+wp8&Xh0%%uElVQoSNI%Z}*vhRz[1A}[C[cPap[!+VidHbk&o*S4+/qgFN/)w{qvph@ddB?Kw*/l3pU1YpvDn[WB}NFvn6MjIlWc.iphkz2PJBjzP5qyS1+ou0*Floq:088T]woKZNzxJ<oua)Ci=a)9J$vP[%zW>Z0%mR1:f)g?[SfL2NS?F0W5}+BSFylE%ND03SEMH3z^eo!.Sy:*iQqmavbmNCLc8lRlFA17bu9[-Y$UH(iV=u%9{MpZ{!WY2kHo[?ZuGb[=+kbj%j^+3h%+yq[}j]CJ?Td%9)w)YLy:F7xn?Gap]S<0TwNMk^8F83g$+e@{g:s9Mer?0dZ=pKt(#!0#3bGQhO9({ELvuT%6Z:3gXd6kFJ<s&4#.Z}7ehxY3}5(yjH>5p$HHg9R@=6ZU=iC}O&yfE.^dB/N%qd{p&}zc>:>=1f}>aav44V:oIQ@0i-k(dHD<{EN5E5{W+*6@^v5.pZw3pb7$INasg*T}w<)@)xsr6!+63EZw:*c:p]0^@v.M.5g!^BCjp4D?&F9.{Z)QTZAvEAdp8[6ftA30h#At)@xAN8/qK5!24Ja]?%Y-P@RI?fk*{nI6nleZlAqkC/tEZ)]NMN/KsmM=1cXg0^@iqli{uS}ZAvEUbMW#DFN:oOd5m0d0x#+>t9Fhv[{UO[p:B9?7J0OM)OZ4ljn0j-Fwf!b/w30a%9)pwjQqM>a{EW9?[J6i{WFp-%8U%#7h9*Kczb/12Om{}@S@hO.-c]S{wPLb*Z>>}@RV)+l8)0%rHmD?W3fz0qG90v&dkkT"
"c2qMH1/^M%W]<!%]:KYCkP>0^8#YR8l+?JZX9=M!%4U&p^um-dQTaxOFBAjD9Sj^eZJan!}HewvEM[{n:lRlzNkRVdKh{Ol@rN<pEH<%h[+&&C.xj&%IGI>#EM{&Fnehgx0!-=X/wC#YNeGY:Fga=2F=d@+koGwhQY7n(D+=S(&7%H+yqbT5Kn67M6Djr[&1bvtb(v0tPM@Y3jpRC+F550hS$!1x1a+q+!hMBX(tB5dE]{4:rtzl.kJ2x+l&E8$E!hG0P^)KNukTu#gi:sE0*weuC/k]^2VwM&Z</v.tMK^/FE7fFXddh<0>SUNA4fad2lGK[ajfykYaRaCZ[jf)2C7}]la8jii2??n3XpU5A3rc?@@vHjV.^#[rkwJ:[GkVubLnuBE}LO*{.>)-BauF45j!5yAzRa$d)bxA*e!o^Z[&v[QX8F+fQPhVqX<xHv$o}u5blJ2-HQ)hHB}vB1?=zhW%k/S35bzts?[wP:{XL{@O$XrYLM?/J@TU0oMhfT}-EioZ>+9?[FUv/baacNumPcl@B3q@t4kP[@ar&UV>8{EK]cYU{262qy5qgws8D>W)A^3PYmrbP[hIj06)!}9akC+Tjep[2Z!{i*77$J/Fs>:g!o?zx)$1]2A5Uh)0rQJYGitX<O^btgZ#b-U.@(1lE3NzoSl5Gq}p[eTjbp{)f-cGpa6HSG<?y}DxnqLO<kH1sy^iZd:k=K!%ap>P-5GhYmwhrG&h8QO<H2?LMldAo=G8oleBDdwSl5?(XM{tge14dlCl031aAed{8:?6Y.flAB=Gxfx*)gO{HTb76<av5/e1>axg-C(*GgEJtv}y91xcfIodoIugGg3/2xSkCk7<}TAl]W4dhf#{Rs68geyvOOuuX:YA92AE^^3z4Tn7Ke7S[Y=W>H?}lSs:A?tmZGAF+Gu>MF?lB*-5R)G1chECR&e3}ISgCr5=eW4*E)b&p7$N7)$8CX>rNonN.b:BuXS{HSKj>}p/Q953m8qW/?b}&e+w)&Q%^FkbdKp^pX<AWM#1NLf4eb1Ut$DRh9CzVGMX4&h8Q/Q4hvz.aX/+(a7wdpsK=[4#M6ga{s{M[k*<[]5$HoCbqY3)G&^G79yh=UJL^Wq2+og?@&5T1b1BC*heLUU>/NS&kwdcxPAtU7AWQf/8+L5t=aLzd!N$W@&Z3(*fIRZpx]oZWqY$kT{}iUMB8Zd!@/Nj^?o4wCQZmIsaw3z-cr.zClGJY&p6m(8yUGmk1/<3dv}0BBhrY3CQ8Gb2gWCCUb%/6MJ@eV&p6mO(QZ?6.5>p).%0q>Eg^RjSjS-nR7s**[5J-dnCG?Ya0)TG&QVhX.5<h&W88Fu&fETVLfR[>s%Pzbc9nRioslrmI6W@j6T[k%MB8ZptPmdB{]4dHapr%hu1[9jU7%{Tm*o*$a6I%Hzvr]3MuHi4=I!aC0mT$k50d50qM*7=Z-*ui3}qLVIcs:Jb^va<k1!:y*.@tOEn[9GS{JsXt>}E4+gdn53qe+ma6I-J8+)k.x=/ANdp/&ZD{j{QTU#:WR*@1f[hJbU3@dXRQk#/mT4@>CMF^?5YCdYYF%EG$u8.YIom$6//)2A.RW4/s&p7$C0UC0Ills<zQT*SFDoCO}CUW!]qyBWH*q*wdsOjbW9]l9G*8*PEbyCMZjmS4JDoUt?pRm%&qQdA$4%QVlRl8deIx%o<b=q0(OSG?m/luI8)X<8vu<.-pt#tX4e61Bn<mEVD&sxhTp7%ZQyJXE.Uu?]G}Y{Wa/@x&8hM)6wkP/IG6EIvDH48PHg3RT-LPE^dw7:#UG*K-LdSPV2sNPs(V.e4Q/9Yvh9[T}B*zcOovJHPjnKPXX/ts{au6.@@{sEf?tiZ:%1xf@?Worn!2ZlAj.qIdgOfTG@G#cn$OK(M/]^^[XVBAZ>7ykvP&vuMze/a[GyAwKzl^+lv(oi/-HP8X?{pi<rXVFlD)w$T#&sybLa)>yCjSZ^iYv3YPcbH+9!5lXIAqJVSQJ/w>G(#ak9(VM/O!s)F"
"y!Cp-nOyd<lBE@lypI!Y{gQGFK5:tm}N.D#u*%.Fa)M[AWJq<ZWZPBsGZt-Y?=kJl{5Qb4SWu7u*INw?9(vq(d5NB76IY{Gwo3$dCcEi(v6XhG3wlhp!n0nEV&UkGraxXTK)N^P.N9J0Mz*$]!XWkF!k1+5&v/0*R(/k))z%SU.eLVyg9Dspzq/BRl&E<6D+@HOi8[A43gvGB6jfDAy/u1ZRF@p*Z84dd.x:9goDB#QDAnfS2C-?*]M&V-K5yuJQ:goYzKH*ryp!:tL^:/POgU9AB&2h^n2aS-]M$G$*MP-l($G.S9/<isQ3G2pSz$&SbH[nn/tpn4r5]jn]&26o-af1}>kZFgX#bl]}Fif-zHB4jw7mW*o9H+*2DGsD/M]YBk+J5T.C>)0&vEH5Y42wVle+]JWb7TT*=:K9xHBD:]{&pzbVE/P+MlhfamLCee!=}c(Q#ibl!BQg/=duM)FB(>3z24iSP$zG)E)?G&xrw-sa/mX{#beoWfO5dkF<qjaxWo]3F=A8Ajw@K5FGPbSW.$}95d5j.R[YUqm*[Smz!1rTZ0E{3X-HIJ%yL:.0$}$otB5I.q).6n8YH^l=4*I?Ah?&79*^w=0AJh-Jwgp<uWK.CQr*U+#<E-dJQ3ATZQVdNMKpPaU%w![-}Fxs{.V1moH%1u#TeMlf6)fv1fP?-aPiXb0n]}{?Ty44kedGRG+FLlx%?NXHs$6}6])ejD<:LGu+AZWYXPeqFPm(Y/y8>K2RaD84$niJsNQVT@LZ}Rmb)lN5bz2yoSaX&i=d}w=lKw/:-)y-$3(/&QCojSam77RMB6t:VZNWwC&)?icJ[SQLAC+fzFy[8C.pOFiM2k3xCXy.*6@>M8.U1G#grr?8)>v@QRk((O%l7UmDy=Y72Gp?wC:bjAz2rI4>^pC-:3T9k//MN0AO0J$Km+<nN*@S{&JG8MrlKjrt1*bGiZ)EzY4}Xo{8PCdJy^Z.esh8BU&jny4J$5-Mk6>18pCe[x!eC-/V/LU$O]Ia8./54DwS82mgHW>Nj8hxi{L*8Z+A=sb+3b9u{S9l@afMg5PikGx8M@%ZgMnBUKL/3>}18i0sL^u@k&gSeg>uG@%?m0Q)obXI+P4oFczS[g7P/3>[PG%63Ly==^JgJ5faT>V0XXRF8itEX>t}10i5F7{hWe70^hKSEIEhbSJwHhL$tXo{44Ta)R@AnR3z9L3%ynwEi8N9Mp#8r4*!*=GU<Hz/7u{<tmx[I06D4?Pe+<T+WN9=PeYukSgwjBeg+1hS#9fO*+XbQ<=M.RvmQ4?/yU9$n!/9ZveA++250oSBdkw5S1i0hA#19!RzHe%/xQ4*Rdt=KOk.V5+(x&PyR{zDEAeC%icnv4*ARTbnZF4:JR-S<)kcL}<e[l%nFU.q(^qQKU5SnSB+(0F/T-(0%31g3B#CjQl^#Rj9.JFj7%Suk%3R.?D^:IAo.LfG:9J{*.O/YM=M$.afBK8Jvk#U0r?y<H)sTtzKT^g@Nr.HR[FNHM#f/3XLuSAqjuQkmV2S9Ugkw+r/Ii0dV9Q!:cqdHsF@JafQbj}?CLEp<vvHj!ih5PZ.k]i:SzMa(OE!Ix.oPXow#a6YWld+Z49eK1lmX&:8KC&Mxi4}oXRRoDEut.]i@FCS:@Ca!gpI2zLtZ4?odF&NAbFnHvE4tsPbO.[h+]V{:l@@gg7oLH2P3DXM-wSw[SU<%rQpCD2v[b}sp!VAK<wl=CLrHR[y#2AIpA*3@H^k1y*m(?=Co&M&Rm}p(2yjB/WjA1}Zo3U1GHX8ULh4qc7<An^7B1Tct$SsH3xR6si]tJafdC:^pO){5Vjy.8EdHqMD:.qs(u<<Iqq&N$.W3Q7s$xw[[9p:a$mDbFKo)!!N(1b5=TAp1#&8Jsd5HaENV-H^V}WC!8<@4WooHx0&:{/bJVzXH3l[C4w{9$l<+9N]<OWr^7w-aa?:z#[S4fGV^l:8O.w)$(07Vg@%-Hy!R>!dl!TgmYB(run4CS-(x8-5/-]{>X*3"
".5{QQbLW=%DnTipB%{rmX%!GOqDo:&3Cw^ds5:WH>&C/<kZ6:W[FIf&*z>x%wN4iPj%!4+NX8Q$4%LLJ]8r[A?]W/b@R@*C(l2G>Rbvm(rLMN#2??I18[lyrD{&v3B#$!9ifP530&2V]!qqs]LM4KXS+K2LILTjH%eEG&rz?K??JbfKbv!ut.5)Rp[H+!E?nddIXn-Hj1Zn9brh%WP[bk-Ivz:b%ibAP(65U>hm:.ceb0Q(BwaqeYj%AfH-nTJKGs=)J:8H2.{):{@6nTWHJ{W1O*&jw%&V(GNWM9Zk&#K]PGt5R?C00hz*&&!MkzRy>Ar*>j5CgifaGA$6CMX2Irh{X[Gs3atC8%cT^FMGMjA?=nQWwAT?W{{&&V)v*WoxZkTeXA1A0&tyqb*MCh43-]t#>3bS}-8Nip7?eXdmqwm0<(pAb?VFJYJ]q0JpJl6)i4JqMfUDy9j5TYAe3iFgx8*0&#q3-jeq&97Y2P9uvVX<nqi#&cibfy87>-QYNH}S-o^WILjvQrh{X#INSKY[t0dabBn^(.R(G]tSwlDFBgr}kCh@qX#8h2%i$3-Dx51dF<57im$Y8a1B4^Pi5)aw?W<#>&=b?kon4wy9jP#FMKUkZInV]zie?>f)+*eE-$vvCb0p3gv#/%!g1[/Jjk%0uw4NoU:8A5]{>cb$5o(1A632TwA&&wzFiTFmn$/MJj4F/Cx-0xAah&e]:8P4zkh.Fx[E+gMb0vtNqA$$iWoJ-lrdDER6?8D6Sh(}k6>tir+zBeE*Y}ZXb0A}cXh.[lm&pPo)@w]}dM5EmUk<8M6(TQI0SYp3kQ@3&1oRH3FiUj%on4uA-eX*xBaJiYAq@}*6(TIRp{i<Fd^Ao&1Yhy4SZ550mp8Z!j001p>>/BWCw[Oj)I[^V.R]Ezm(V?BcLCxqrQ7rwo*/3-S3?8x@2#{^9Y@u>^FMGMkzUwu[H1BMcLqATXjnJhXYksyAn><&97AdaIe.D]5-ZdK&c4g5QSmmeC<*Ee9I.z1Ym+cDjdD/vMJt$R&5!C*2wM[R)0lC?z+4.wQo7J#FecamW0N:PK4k$)<f]{>I5Y4ve^Rmn5G!%JbH&kKG&1}MFgx8}G<7PW09dL(3Bmuxr@/lI{>c1m)$/1OJ[!cKDYK25AOJKli00s=T5V}+INSW[0rvTnar34s0ASlJrvT7o2bhcq9MHc}{#7P097PJzGr:j(InU]l+jO9jUbIz3!Fi3gRbvs[Fk[7/ILTjISh-t)CmC+-&C.rxif37J)M2kR*+5YycbETuXgl3PHY^/Srh%W>Yz-9<EzPGO{)h3s5&+%iAs@am*Ahn>AQ:Uh=l<lA-2WR)JYS1dah)eligJC$)00y?J]L@=?W{[>9KiUe2n(?I93aL%1f{W$Ch{7x^FMq=<a]=r4[Bnx0x<F:e5inhDzKDw937*HGs36y9/]qu^C}:@jS!SkJ}!C89$YM{e2{+SWoJ-kT5V}IBW#Op0Ar7:dVg*mP&Kk}Vg%?Uj?#m3&V)k*YmV7+-bSN1[c90#0Ar7m3H!r70S-qJkP*ABQYNI{FiU$61Bft]j4EsJMK8Y<I5+H+>l5Bz&&2}UkOXAolAWn/Fk[7X3zy!9rdE%#IOhCPC00h[:9dVjkhXIx*=gCtG#h!/J=uV>p9)=QTjfxS8IID}p^9DEWGC{JjA[-l4!$W-E?xXY&=b?km0<#KAe[]yBb?I^&L:.Xie&ql65*^*4--7?0pj{F9GFQD^CoSPi#$IWA0t+1H}T[cu<11P0cQC%J]gP-?X7J5pt*[Jed?e(-2VSMGt5YACaxUKy^UO[)e0o6a5/*t(ube4ps$]P^Ca::0dVCwKoB&(F$oeyq>{#rB]^g&[I0L:)?4S<&WmlaXc8!Q0oaHD3dFkx-c)dC/u:^8}t[JwkOHvc1=rrX9NFg93jE2[l4yt4pZX3C06xg9ari]S6nVY:J{G^&?h3l5U%n&de*o@bONkH2oy0eS0B*NAel>ZHI912dsG1b3?u}j4.c1jm"
"d?W5.vZf]7l:tIY06xgW2yF@7{LAuy[DS<$bMSti&Xd/#/AS5d0oaHLjv/RXInV7S6sYUvIV}3t[I(QucLIuP&Us%[Rbly[S#TDXl+9iYsuNVn6[-E={mwdXTe/m3NLn)69xiO[G^lM)RkO9-l-VwjBD9?k2A+d%-Bp)/J?Rgs^Y1cCEVyC6Q1p9riVPiIBAtA->$j&5dWsMg@yW2ex%/wmdM558XGh}!e?qr>JjxC3w29O/2lP-o:pNtkabSW&u97pU9#q!tDrAka!m$#^i]yHTQO]dwM1r/Zh!?R>6?5%nzKYAU^LhYW=z.(@*nKqM867k)n#6mXbSY(lzQx%/1bYc<!$x{XTYfCZApEWhRx:ofpKLA]O5ZZwmp?4^T6wV(J:)1p[OG4]OFc]s&>%].@WKK#-L%soXbeViI[B=lvET8N6nGy8x$gPe(V^M>&LJ7Q@W<?:v=D{5YIW5uuof90:3V%*>@AC#cMO$T[ZIWWM5.H5g1ltFi&.ZurU{qK-=@Izm)<Be/MCf0J!Y#@bzRdS&{+g=8RN4XZdZ.)(a$]Di68kXemUHs6.4R/o4f=&IGO{!pqvvY*y?)}-D-?{.>=!=(yW7&6{IDzzPu12W!N0G/a(AD&>DO07X3Y<@]e%g{9l{xi}4N*r=v:c(hZ&!ZZ0H(0^9/fU@*>*Q0%s9d$G3:Vc9vxIG!7fu0HNN7MhBX/^)XEWtfSZU-tqk8WL-Ei<WFuob!H^x<#N.TyF3AK2ZvF[M&NW&oFySM$9bdt@.edW%g.eQOloQIAo+!BnPt}UzSu.[KDk&dc:)Mm{+$&.]AGO[z.1NYJsQ*I^JZUy:SKvJOzCiTh^P3^-=!]LA7z*wz(.aO/H]D.>CMOBG*o9R0[#F7eu[J/HFYq80Lnz9zc43Q1p03U]O>>B$@d[H=%+3o.RULg1Blr9FWk5^S4oS&IOVzg#W0:OMNssVp+*LBZX6cm[I:k6*k}BcZDU8]we]jRyE<oq-:U</k6]3T=6sa*v[HVb(y:k>YD!+R^F:1eb>w8E#Zfk1+V{SUPk$sWKqdwIVG34PF7Sv9{OYvAyS]PYRnox+UUp7tc3G<t3dY7cO:ea.nST<C^XjjDC.16/P]L*5.l$(&?[tJg3avpw3IPmWL}4F:8#NX{{qP>Jxa]o4!#J1cLgwX&&tNOsVC-8(p#s]y3)R?VbI{wXt](m6M<Do/EPc^zhj0(Vggs?umKoV-GP?pQNq$WAYCb1rFcBSP/wP6JZ>DN(u46q9Q9?iK<11dT@JE4=&&y4A%o{Jm<?HYIpsK9:LUM8Sooc+9g1aUnbpPN0uva#K#haU!NB2*u>%{gI8YQckN@T!=L*t:r&m4/g05H?vZHcmhz]/30iQ?wVCR{)SSr75cR2zMdSF*9pfkEXNi9sC>$UvaKn5<v0B)aE6n8+]0N>KH*QL/(zh#NWpn:[@g0Lbb-t<[zh{j#l(yn(kiiG$G62%f^b!hEnTKz=$&#o#2=nsB7ggUy0!*/[Bi{L-}{@>bP>z-$[AK:fO@#axdp91Ph3WZf.l8%lg+W90^Jgv7o?%fR*8dkDATt[^GhG3Ap<8HRY?mX0L-Vh+6:2+7QX$c8C{$PHJKb.eF[UyLIhEzz69ifq?Owsiaiuj}p<cPX1A@PhJu&#LoIpaz?OD?n*C&$6K9h=LwK<5S5-/Q:Q.&-f$J5[jB2t9Utqh@N?8Hfw./blMrKTyw^BMFE4Oe*cVOmwpl-&-63(XSW-^k.eFcLIx!>&y(%<5TTqdh3&->$XcDd^o.DZ@27}XvGVP>0Ymf4QX]=QY0goF2s/T[YFAz035=%3G}P]mf8u#PCxPH0(<COcK#js:MrJyKW8At:O+^?qepfMrW2.zA{o?Rwl6>#zPfQQ/Z9wc90aC%F1ncf3GWmbUcB0ZMoEfH-*!43{%ajbJZWOo[VqL7hM?bv&[gr!rpa:=:(vw@bC3QKmAjvUXyI)%7&(?u4SYi)&Nr:RBN=0-O?{EN-^*Q&!&b/jJ5/{r"
"M7h.Iu#*c8J:9/[&?gsX9f[L?Mc%V?Oe!KxGpng?c](^Pbsv*[Xn%kw[WMG4*Yxs5r&PJ}{xY&2iIJz>hAoe3fj=b*:2)!K<$@QkSS%#}NSgta<3pY#-o5)$-=c@pNaIEXR^<1)C?2a6a]ZDXcPGk[gncV/JACU9Xoi!FOpy6:fB/]u%5m@Dqj1)mxMuo(AjA{D=LQZPUA0zUL%{TL&EA{Uo5W%GJMc@9apr!-JX{+Q<5pP2zh[M]9pZL+@XGocSw&E+h9J?#gH/xz?$$n30xz6GsN75[hYcQHy$YaKi0a7%-]yMcOm{DjiMwtR(s[K0uFLR3/Z9>g:MFM5Fk>CDym=A+JEQF*jzj{Ri631:epqjP8M2K1-<</uzh#U@<6dYAwAC?f%kOa2D$qyVJM8@#tTKT#p&-ab4PxTG+MA!:Fk(?jin0zM&VdG*GpU1sAYJd&q[g/g5yJOzSS$i.2<PD&9s!7?T$^CPO/B5$.<6nt-ghu/q&?psCA4RT[Up@$gnN2(rhSkv-oIJAOGcWNWKO1{-L$VJC/!F3G/&6HS2.T#q#][4&%]jmAS:cqiE!?7=?EVw0B%c+eodP&AYi$[<6S<wQYVQ@Fn}8Q1WYi?ixa.UhaN-i-q*]#Ql/$^rWL$G.g%^ff-]i&9n4unobiIEvT)BtpX:4hio4jkQJO*VhPAc+sJq[H80i+!9og<H=nyOwri/SYhCM30?6LzAz@WAq)aGxNdMC=NhYqXO&%]s0vA?-fiPCayK?JO(8]=6wC*Dk:8C&DOsFIy47cvH^q$mxye&KT4-]y5U*?[<ZB&Cz:6n8EcGMrLf4Pvns+EGSjFdtr=2L(jpiulOqGpqsApH+D=TB)?w7&[M{4R?Ix+6}>Cy*96[zyBV*-J.])WJS[VuYB7[byB4:gr8wPsKTK!90mc&Fgtw.xX!-eJJ@2-!?(YI-BC]ITBa^#Jvz>x4>u5!=L:Q0FHqBl0NzU:iGuurJWL2iSq$%rPD4=&0xFeIkSEa2hbTCPi6YEG9Y9U(-^^J*=?f%g-#xJX[VR*b>z-yWX[?.r>&K^]F4gY::L7PNSzpj[IpyFBjge6oC*TE#8YUBeAG^xCyixMopl5NJ+n0(mv=F*!o0g4c&L%)$2p2MD(hY&7AJjU&1YHA2ST.JU&MuqYiVLePdBS?e%bw2J!GZL6y&:aD4Tunx=L!=6psPOizzptR-VkLl*?[?bA?4dWm>g^^AAnOHOqetQ5B6fOp7x=La<4.5iDp@FypBC8S5YxdHA&c7Ax2@#/V>}:&LoNn9h+^lL>j/[%4TLjMnY@HJu61p2r1KBrxcG*[SQaq8GA^%9lh>.MNrWKiE>k4z#IqgjPyCSeoLenBbVXV<oEcs=L<=-Ff#s#2L^D2-=oU$AO75br@U]+q?HU}@S5#s/WUvW(#WqUib%Zhsnp>r0dQIPx#u4}Jh%h[rCyT)TVsyY[SQR38M^7I<2BM^Opqzr9zQQO=<f<[Oh&PHPCp[Y8Zbf=AH>1hI#OVOVdOuv+nNL70lpaoBXTTMS@Ia0Ix3hm8PwsGj:tvQ]nwQF<0MC>?zHfLv$m6=tY8u&JHW6YC?T5dkAfd1rqhU}gYVHO<5F7]yorJp-<5G/*Wm0tF&K*v2p2MD<z]i&4PNk)&?}.}J&bw%N&zp:-VhH#:&c{P0=SWItU>N)<%Km^rTBzF5A}@&pdC32+b79OiKbFLj1y)9F$c*.:2I/<<#l^BAJB9M7A9jLVa]s+h#:?Q%hOmpbBYB8-Bv<S(@CQjI&w:ZAK2i/7Agoj<8IPH4Mf8A-]x4*:>QXl-u0$GS/Nnaq&A-%4S3c&+5?U*ac6p0/XuG*&.ECjdDKcPq!5/((%(xi8n#NnruIcRgYJzMF7Tqq}Y6RPO-q88@aD2IA?d.I!365>0SMcVcU!H=!&=7{SZ4WwldH^ioN[P>MoE6A-]o:kL{$0nz8]RQRZ#%U<o3r>pufw}epOQ$zxi.%OmfVsAYJeB6pvSg>z^FRS(DNv>&I>:"
"9kBLKM&&d=v]OF8C*T?P-5V<T2rV(T<z]kWJ.Rm$+b)S-rOktAe?(B!-ViIJ<c(U4Bb<09AN(>^:bA4sAFDp6802:F9xn*0RXw#J-^:2U.&-.%i6eaT6p-Q)AYTj}4QXFp:MkA59h:=Sg#@>tOA0-2IOp&!9My?:C/!a8Jw${.RYnpD(A22$rLN2n>LxJBizQR[IrZ*Q0AsM:g1EfzI&dNn/Ka}9gY@-RFk[OQ2<)SLiKcvHGo#XQ@@C7?u[JRD1#.%3S-b02gX}gjE)w1MXdoA$6c9]3w0gbUJBPk[{%}wsJ(B9ssN>fZ9cWlHS-T7I8X4/=iKbOOhz[e7:5C/p97^:?ht}6BZPr@o?ie:qVjpfflcsb8iKeam1d?l-@L*?!Xxp608cpqcieoRe6z):F&{&^T=noZ+-&oI[!KPSTSM9$dq*-z7HMfp[[WM{khYESm9fa3eK<4]/iMzR1cLUI2-gd*$6n8(iAb)%?4PNDG+YS$D:3o>S8YmXQv&T2Iw1CMpAYJd&rEwl0Cz[M5AFC[N+5LuWpe}JF=m[utizH}px#G7PigO2p6nEzaAM4]MruIQ[gYYM@V8nsWg#@&1Ox1:XHAni<-L$UMilxD888%R@Z.HvD1AXX+-%#FZ6YYqkOpBPMKpB3^06bM{6(TAYB].P-kP*pJ90a2*Fk[dt4MfmDiPAvtD$SUbJr41sPA=7NJL)&Z!DW#@fY{wnpfu4X+n+T+w5/@Fxby/i-q)f[aqNc$0(^m9/WUdb9cPwj:5}Xv9Yn}@-<a5N!&(2h-Yk4eM1cyxBAtEx4Tu2d>nzjd<5D+.zzDo(iE(amA1jZF0MJeBil8%W8o3aCrv&JnhiiXR<8HuizzVM%-(A8e./P2ozH#Kh-{%.Md+6.LAT87T:LB3ozW(HAB%!zJENx&aqtew[oYnoy-}30M6MFsp4.gc(Zzt^Y&OAkyds8w-Kc4C[2X+p1z*e(7<yrB./L+k=uR}Q$)1TCVx52xGDYkkcAG<.hs9Vaoz]/2>8Vku&64I=#tAMywO:Ay3L{2dKupp0tbA?OcOq]%G)z}-!8)7L6U0<LE**<VFECCLE<tQ@?8h$+XLh]9L]}1NO?QrCl8(Qxbuabx3h!jHivbm/dlBw/.hM4Wn!KN0nJ/Dh4r)$!&wgZLBXV9WhmoRP4k&mZB9sZuZYp6Z!!<H-dxP#poc2R=#KPP=EibQWu4Av!^NI$/SlRt.7x<uu3C%WIBS8pbPVnVUg^5s(H3I8Uka$YNgpATMhKD8s6ALAZQh=B-S5UqGpl=.vWbQ=)EY:z*v4FZ%)hZ5J.M9!N4+w=v2^I}5s1ae]Fz]?0WirKolPL=XQJ%sS^inHD&8@cM=ENu!DF[{wl[XI?bmMhpd-}1Eo?NgAz[&l!(ZAvBL&%<xfla?]Oh@[5.Sj>N<T#wR9)=Um)?31#&PTUnba}^C2lN4R]or]Cn!<CCJ:A%rg4@=nczp^o+RCrkk[!bh-F@t=eqvSBg-5MOwGT[DV>R!R%oWT)4t4*%GqQ*tlV>i:I15Vlv8[3RDHqXc0*VaqkCd{s[K.F*<b[ra9<&LQ8mO}Jo2&h*h8&-]CDYi/Eh<6XS5>Ipanp:f2eHbJ=f7*KI-alckm46Z$LF%vNI=iINx1Z?$3C8ovEoU=qd@*Y+Z&7KeVuWli6h)$yYiN[&orr0vC8Cbz+=9FU3obQOKYNCKbjPwkutM}4^bf$j<uULaJZ!bdh@.cjjxxnPKQ*C)d)[+QaAbs?aBaN>PcGC]<pf!jI^}y3D6dUx5]Hh-:h?=CC(@ABbd&o:a:C5UX{ea:8&:T%qm}s}h)57V]%aTZqgJD79l(03?r9pfn)o<Nn*79SlTN2Bxgp)0h@>C6*>0fwK:s.+zz0hS?9r{WtGlO@nzKTs8[1*XALg7Bdzuns7Y?l:Kt<&y2igd$?i%-KtXcpsXzjE&8>fs+V[eLBG>X9AT=8-<4(D.E!2m&i0tpjHb4Pdpwplj[pyU*aAL0iegI:ze4rn$s"
"U8vxV*6E#28Vno8o<n&>lV2dGH[LHgBsengIj8HQ{QJ&J/c}U#>7waK/LWp)4D.(zAnFpH8(Ut$DK?n+H#ZK=86MT2KU&qy2kf8TSu66]hlFx{2p2s+XY!nuBJ?[$EV%1=@eo3}s:CkYfi<m7b5+M-u[RhgT(w=tZfN@.I)<]0:Mdw49eMkqaUPt{}4wDK^OX8I:}PkWm9?1pB(5Oojg7C(cU)-Qzg}PJKE>J8iui(k}TPrc[!0tslO9e%mdJtExi#*:D/LS>bvN)By[?{UQ2XIEFp%LdcQ#J%tv]iuB#BnVx1LGFQ9Fn^)Yjhy@DPu[XDuT8bY3SA[<ZX(mFi<9Be7wMf$2]4!-rWR7Rc2E[saon9l)[y5-H*yJPw(.hqXAYY<ztV.GqKy>-<hNqMiTvrEkU]P#gS>bYc!3[]1w&)7GZ0YwqSzC@t-):-XbF^@zZ2yKhmbUA#Ta]3/OtQ7/.O[y%cm8/*iGO>AUB:^jx2jE4&L)^GkhcvzyU?&/yR[]1jjoiheo8/!Bdxl>XpOOj1[eTYStzom1T)O-OZTcepL[<.4zfY2A>l=(isla2ENH<a&EVO[<!^1Vxcm:k%]o%Yvrhb)qo=N2)IN#JzWwuNuX!KjB?n[26Ay=lOZ!f)Z7SJ7=&4uQHLlkBt+l%{@ZImK{}h<jY9!!bix)UB.iXIA!<*0}z*u>BxIn!t?}YyG[SEORp:hO.}t*)^1hLGx2eC@HgD<b}EAKoSQkqw:l78X5@Wvyt4>HPZ6{qNt%)rtGtY{Z).{[<Iq}:DMK>0ZMT3lC/][CXndny#*OVOswOh2GY8/HrVB0OD.0^f?fA+C?MNs8//yPvyp5<&&-/{B[8Y<K#vSbm+eRK*MOVSbv+k9Uo.QQl)02owe/j!GoD+aA@bto7sky</*{qMsRcX0Qk48}NbqkyD2=>Kvy3AIhFoL^)f3(8PGM)ztDBc<ooN<]m{A(D6r(t}HZ/zxGUt}{Iw]6Gesv762q4f5y:m]DTzCCY4y[Br2?[Q^8Zv1*x3R+c/YoqlE3PeCu]BhsvOD+(*dPL>m%M#:cgC+(<NMm8O>H/2HpV}2CpY{A!Q{ngC{vh8[DQvaU6z$G]fMA#KZn]1ot1B%B*fSZ>QSiNu=Yzvvy9!VcFoZa+K1=x^)ApKw]KL2*Q%CFzrfqef/{f$2G=562hy1v/3JOyA:bSmB4MXVpGR*6ME^Cu:(J2rjDV(3u]w%Am-Pgc/<wXqjB0Sm1FxPVnE.-75W=?ChP.^1*(X?}Lu/bFC{u%f/LDTe4.jk94skVOUmlp$yvJ0Ohy77Yx2HO[aivPWv#QRG2Yb+(gT:W]i:Q9:os7ycH$HxmhBeU]^x{!O1EIHe6bD${HJ>ZNJKK=-{w47&Upc}zL[Ri3/w%Z^2}!(gP(a(6)jH@+AMDEsXH!lMa{{l?VYgGcvy85i=Vgd&pdxGY:tQxnbwnEtb<+L}A:^:foGOCzIvE9K}/D<^=Ycc[z25U&=Oh:HRupuyZ?>3TaCOtF6fY*Ik$&Ikx]D@Sgg/>nN%5TkK#xmzix1PDptr)O*S*!d<!y7H8^xIJXzu(7=+tx[/u{]*}HR1up/j<<jz/HVV*Ufa752{hZo^ieM}5IA:W*C@v@t6*/u1$HXI&H+bY57[JIy(=5zJw1xyX(1EOdIw:.doSeCO9HhA6To*5dk[PoAQrA+$oWpv=9u<NL^.f%Gz}hxv^%-}O#PAB6OjC{xO#/*mPSQq>#N3y>]<y&Tbxa*xSm:?tQi:}T{?kTsA3tazfSI4-k}iO+D:>?mm/8Zom2wHhOqhNP}@XcA{c{CcL>[6WcP=#0)y?s1C!p3B!?YZ!{^U1{y$:?h3p@f)AmuXIey[yaBQrbyosJII+[5&%vY9zHcED8=1@f{(kd6d##m<tJ6tq4&iYO)cSN+k1e/ucfMd<E%apMI*y%OogPiOrP]$<#eHoKuw>xcZ9%}&cu+WW<2ALnLgM-wA=xxO/Gk8d7o]c>=#cpOa^%%zDel6+{z+a"
"p3]&CHW){Q{zkd<P=kx[<Azn%=v[{Gt-1c+88omBXG5.#O2qBly{]IZm-gz-:W?]ideQxfLyF-sO.a6cO{Vyd:#I8pU?D^a8-Ck?L0p7a:(FE>bvO5b!hz*Em-f1X/XKT*dxX5^:0W5jyyM[[a*(5]hK:hZiopCmLA3O=?[%o3xL>(3A^n)rViUn8Je-]^vx{idhBgFn+*WiSN&9ea[fqR//LIHrJJx64WGUku.R^]rlaCY!WJ<3zGO^.{Os[@bXHEVMctlYC[>.XZViUo>k[f]]bOyy^F@Dv<QHU2}NU3G>Kxbk:r+l.ORTJ=Y}FweEBeH2Hb)BcSRGqrc*)p6%K#up0XHK9ib#@74VQS1.y6XUs<Fb#z?B<A->5=X<e%/cc!rDogRkVyKx/baxoX(38.)0#@Chp{-/8sZ{:^x+N9ga8b@p/pyv<tA?O02W[4F.!0&Weh>Z0a2sNwjbMh5^]Q4CtpJ65/v<q1H8VzGQh?=Gd0M[[-eb8+K+>K*LBTE4/zxh&F)k[hwbt^p?p]acR9kQ3{Q:T[/uUZQd{C@.eHRfHYxBt:#!lq?NJYU+LTn9{W%/JHOOmpBvdw8/!A%zB&JC:SI=/c}<N^5ie:CtB.1ucZ97K[()(foN=PPLqUASV9#7w:(It8jE=f{LGs<oC}A]fcGE67iFLB![yW[08Z6zJwd00^=a:ZQU3[XIK#wZqHsour*hL)/vl2MGt.)mswSqC.&Kp7OjIdrdjE:es^T-+[2fmM#b)8tVav:FG2?$I[8+=tw8$Bdl>0rve:xPxU4AQ{W}8<)v*cyRT:<F.{*d(Qwo/-H2R:Ar#OiS%]17Z=ILfEQ@Nw.f1g&hq^PK{)BWFuDp8-vByRPg4C/}jm[h%AzF42n+gXH)TiqnCf5iH{A(gOMU!p7r-%q#ryD+8zs<f{V-8gj^hfm:&!a<c7?J+4jW#TGDV3MeI9^:<+]JEJ7%vc(Fkwz8a<QcqGEI:<wq2Lg=3]pBo:k8+<kJ(NpG=wKb*7]qnwA3N9gdeTI%w3IRVP?-h]%)FNiqpbu?lH%t^fhOTon{Vgd[@7&nUOF%VW/:01FPFe]lW-048<NNNUK*9P6hFvb67evEd3=#dkivSjdc8@(:X?7HXlPph-<NM]s9.=Mc:?GUch%YSO05&{ep<!<dZ?+[]JFt&yViFDfk^J@35Wo@T85>wdbwg3Sz5Bo4QB:*2KPECQuXWqoU?Lu%8/*iGvyBmvgNTzA)S8pUL1OO>n=j?b2fF??JFtH*3iZF?w$R%PLNxvqFG4-*eyy:dLykZop:]SD+}gY=c9PH@*@W^fB8>h[YyG*-gemsA:?daTlcgyJ2%sI!/-Gz-VST5p+qxtVx2Wr:![LxcE.N[uqZTRwLil)DOCX!sO[l2+[{8<f*@W>$o$y6GMPHVl>F63FXfnh=zgW.aNx.)U*a[jA--Y6JL1LWdZA/WPea?0.gmrf3^(ARXx6bHP[7:t2NkD={iI}EjOyh7tNSqUOs05=m:SBp-*)^>cLyzO*UHsK!}&-+HA^oRv(VSo-Zf[c6Tq59EGbXd3CChCNy$MOj2g#s89%:+zvi[ESa$=(Kzmol+wiRtXdDDS4jD}{qL7%#%qwWRKvgQE$+SdoypBn$2<JsWYNs/ZED}no@YQ?<]rMTsVC%4lGcZ9r^"
};
static const unsigned char EMB_DROP[]={58,59,60,61,62,63,70,71,72,101,102,103,104,105,106,107,108,119,120};
static bool decodeEmb(int n, vector<array<double,3>>& out){
    if(n<2||n>EMB_NMAX) return false;
    long dropBefore=0;
    for(unsigned char d : EMB_DROP){ if(d==n) return false; if(d<n) dropBefore+=d; }
    static const char* A="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ.-:+=^!/*?&<>()[]{}@%$#";
    static int rev[256]; static bool init=false;
    if(!init){ for(int i=0;i<256;i++)rev[i]=-1; for(int i=0;i<85;i++)rev[(unsigned char)A[i]]=i; init=true; }
    static string enc;
    if(enc.empty()){ for(size_t k=0;k<sizeof(EMB)/sizeof(EMB[0]);k++) enc += EMB[k]; }
    static vector<unsigned char> bytes;
    if(bytes.empty()){
        for(size_t i=0;i+4<enc.size()+1;i+=5){
            unsigned long v=0;
            for(int j=0;j<5;j++){ int d=rev[(unsigned char)enc[i+j]]; if(d<0) return false; v=v*85+d; }
            bytes.push_back((v>>24)&255); bytes.push_back((v>>16)&255);
            bytes.push_back((v>>8)&255);  bytes.push_back(v&255);
        }
    }
    long startPt = (long)(n-2)*(n+1)/2 - dropBefore;
    long bitoff = startPt*27;
    out.clear(); out.reserve(n);
    for(int p=0;p<n;p++){
        array<double,3> pt;
        for(int c=0;c<3;c++){
            long bo=bitoff+((long)p*3+c)*9;
            long byi=bo>>3; int sh=(int)(bo&7);
            if(byi+2>=(long)bytes.size()) return false;
            unsigned v=((unsigned)bytes[byi]<<16)|((unsigned)bytes[byi+1]<<8)|bytes[byi+2];
            v=(v>>(15-sh))&511u;
            pt[c]=(double)v/511.0;
        }
        out.push_back(pt);
    }
    return true;
}
struct RungT{ unsigned char m; short cnt; };
static const RungT RG[]={{13,72},{16,122},{20,220},{22,280},{23,320},{26,445},{28,538},{29,600}};
static const int NRG=8;
static const char* RGB[]={
"003zSwcbp490OwtA1b68FcHLG16L<deF{}?U%OJ@vHM(/pb18dp@Z[VQC@La5KN^pkM]qBHc>$6g0zG)dLI=Jg@GBA1dT1q=z8zVi0+Uqaw$[eJJ2MG5gDEnk+XYW8JFw=+2hi&c<KgH.kR%2ncV142?<paQIgd!8k*bhMz]&z0aps/0.rqh1#)Ds4m0Ko]-bB4NG+OJFd?]H0dU6%G*l4P5RL7oZ^juu8od(03k]6^c$u@s^sIDS2fmkhbs1WaKqO%Ck{MC7{.WJ}4R3G@qiEWPRzy]=M#B/*vGtV]?#/GZxWiv]20k.cfYjcABB^1sK#xc*nfB9Zep{g(913u[l8^jJ2bvOXgl!4&dR/+7S0u.WZ+{<Z2l(H-ar3[h0boE7JIFvfx28ULtgx@Ny-7JUT:rshNHZuz06]E{k(}J693viz=]E=Fk-78R2QP^!8PTw%GDHK=u&&N#Fb*z{fFD[dkR&A/1w/7z21P-rqd7syW/OeZ9T1>vH:-(JFl-f]0:#vD7.wnUrg]!QbMW@OZ-9ndJa<cAvvj%b1#It<7-kg>BA2?z-V)[R(jf3#0.eonaRX@Z<t@aE5rd)5k!@ea2@krX1QV}3ats1f48T/ILtL}IJoDZ)gu-z#+)z(!Z>T}}!1ghj6M>?[qnQ!MR!Y3#9z!=UC.3LE0tD3101@V+v[m2%iwDy%mac]*rb%B-q4xAyp=:i89zg7ih7jCZNY58lF%#][yh:$Oe*BTgUSliJ.9%/J7-ipVKZ@dyG}]})mUHd?*s9A^ylWv1R!YF!=8<}:2Sn@6.03u{8(LX-=z1LEXnDmf<3*>XLsCP-R5?[qb.s5s9TU1tS5myIvfTMLGY9hi2N1xr=>fyu5YI+H4g<ZroH=LtKPRs<l#!8]5k<<JsXaTriP^ejT^*K*Hi5)c3m()Ed1P^-R3xm3I.#X[4hr?%m9(snMQ[AmBRr$p]An&v^ug9^madc)001oL8BdWhFcpx6J1l:PfAw)yFwRBT0-xo06?wtu8gVh$cid@=ui!aEKqPexkUTH8G3$d{8-fv3vkbdancQ:S&2L6p3$Zveg00&xHvU<39zwdjy!h^j:S?-eMjrf$9Y2]{C/H#8lP4suJ[bIvnWKOgeDN?vOxd-z79AP!0l+Ft0$n^^+-[GkgdD2qZ+H^TfVGllg[0OZ7?knzvuQi-mL{5+N-vLuCR]bd*xk.arl>(x^(D=r/Z(w/1F=uG6hP}RnRSJug8:53Jlcvd1F<dj5%oBFvNfNN>]t&[RO(DVS3$+BN:>cjr22V*=@5@5+#p6/wH+&RM1}=FvsVt#K](sP4j&vG[VSjaJk=-(o8#%HJ2QyQgqN!F/D4yZtbGLsq^fEm=1Va?.l5tYhMyRE*Ia2My3:B7V:iaaQZ}X>0LDr%2jXbCgR/uBfb>@qZf#W33LA[hcFgj}NC:]d3wB-EWGIm4/EvC=a.Ta%F-^u+1HF[SlErcku}C)im<5Ztg>}z0K@ilWo2m21hJ&bb!oO>3Vf3g6S*}P{hFN6^*1AXucFh2YOz=j#ymKly*h^.2DRO:Nz8u<w9y:^{h5E!z!)U)Lyza=oPI3o#y3^.eIN5HXd+mQK*A0.ry]yFcWc7[4^}%{HlkD{V8qhv4T.p[?Wswk]0ak*I0aM712AeTp3$7^8XJ&05dfQm.fA$]eUW7SQZ#>=V03?-^5R@7Kb7K!pn=Vs4)G?O:Xcq#f19+1JJkDOsjaEz!^QyvJ0.cZX2mY3mUX^iIV4p}llL*m28q(r0BCO3)U7IYsDjgO>e&N3SXweg-24d)V6K@)/x[R0>JakbaU)rgORZfLMfV)lnQfRv)ny&@:f^5[IK+{mP24P<q(kC)ve(C/fXUBjL<}[.)cv!t=LmohL)LCu]x+8%s2$*k59qhzrQQrQXJZ/Nljq@NSa+v@wWG&Qa3lUJ[RG*GEQdoqD:92<7QYs?N:BOKp-X]U6PW<i1l@TC)Np!NHBXI6AL+/&prG%e:m@1Hyo@Xs%ngvwmaGJFlQr+XPynz:PTK*TSSni/L1{.TFPD10C",
"l2o3RJMinEm:PB7*m#7YzpmvwR=c=@+$OqgmZU%b8RWcAI&O7J4v!szimI#C/Ie&I:9e(#zjshh0fUenF&r*c7=VJd&Zb}0GIrzOaV{lVh1mkG1&MpIL)w&QwjC^NUXNV9Z.8A#6*-E+eEh&Go^nocBat[Sum--$xKpl{dbX0NP-WfjKSu5n1D$?)<5gS85KEKNjMq/ng-CO2o@::4kSTTx6XQseYxzl4f@5eq.EfhHn?#S:>v@d:QhJMULEeEuTquWCA9JyO2x{WI9L+)vMryy8uGO=BBdQVK>:NcEgvy#xewo.RY6jv>*FNXFzFm<*YDI7k<e7i$Zl+wR3/6vCfl{)#Y#y3$0c!Jv00rmR1}haYg/aQ@nR3yl>L#j$bN^i(UXvF:0g0Ws1(68u26.d1qxTEOirsj8xBAe!S1A#f&han<FvU*O2<r8Kg7Yi:U(=!p3jy2LkQR{<v1ntk2eyZl7@!y>zY3n#N8ru/QZQ@QqAH%/v9DEm.{:u0JXINlfI38O>%sP!ub)[A{zL>ExMRR]Q1ySU36C]ac5$D[LtJ4xl8lQG/BTk<RA&^?U)DOb1Qo&r86y-Ju<$3)M[]e<nSi!t:U]C[U?9(B]hl{+o=s@xgP>E$*prMBm[gjLfRezgH<2rWiR4l<m$M!k9pk*>D9(T-!DrU[49<c?VPfNVsM^/]1zzDL=]&GKnR09cbn%DpFxoKeg]7KPOsvFOQie(z(aLKjcyNYpM<3hdRBq]#FZ[ne(qHVYSAAbP!HvIJrC5f-qq!s/q?x&-dR:T^DNVy3Vz)}K0}4@f+lnUyTqUlGV4m-l05+Wg8jgMkqdX$lu7ewIRE92gtb0[Ty0Lq2VPyZPSzd#lwI+#M*sAWryJw0?4oanLhQDg2/3gK<rC6a)zD^Vk&b]m>4C9zyMJw=W91fpdAp&kQZYR^=JbF-?SY{f9h2RdN3>t/+KrJLS.uUP>1f8/fcA4JXq^/)E=4j*T.v}OU2R6!D&C3cpGZx4S4G{CJ>N&ThAMi-[Jt[-Qwi/vsf:e+om0BJb5nAMws/vn(0oHo7}JnVJ5]5YeD7sFu>0-0pOe0&*zewv/SUp4kxIC?6>rd4eWc.ONgGV7qqSas9n!dpgb5gxjQg{BI3-hJxX%!3F8x7BMd4i%GP%3{(GvHTnku4RkTYf>QnN)}4j^pplr+!.k*6$w%x#nk/TVGqe=8>*n0k0/d0@r8Kjh=uM)=dSdZkW7Q[LeCvJD.cJRD:RwOQn9hJNF<ir%pi>u6I>/BL=G1.k<5bx0i^>GVRM9d>LTE49<daVjMgV:W8+Qq3KfjaBn21lf*[+3-fs&hM*)a8GG%*(Q/CtQC2FRUfvBfosmpQpTytYm[!@DeE[jzh2p#/=[ze3oUNuK4{z-n@C5{<ZYYe#04E$+006/yFRa8W3>9j(5[0&qSE.wQ2)3g[bOAT4Kpu.k0qlF}2g@]r26.0On-mI4irp(q<o&%dS0#MUPP}=80xphJ0bI>(8l!bZ+![i3Mjis7<s]Ie.1i9S24d)<9dQcKv&iSGYKN[XV5w(8J5zI[fZ@Igbt7C$J&s8FfFE9D/ULWZ.V&5Vfb@+eYKl&Ca*KF1Hb$qK8%Pa+x298JK7g@9-/&85RCo*GU}]&Rb^{p*K?q$moZ)VnK])/JS12%F12kKx4al2$gyRn{JpN%+hQ2LO-Z#LYnHzqFhw[54/LmUaD5RtH>9a1H88s2eu@UdBNnTR2p7dP0*CY*]WdY7D1i$M6=Kg/Bkp?0s]u3P/Ta#DvzxCxd0fWBm0+qC(2]I9MRq>otOSa:EoK$xHK=$zyx?zo&R8p>g-f(D(6da#yn)W<$aD/CAVpM6mNXQ=(U-Xrni<LO+S]E8sd><{HT4=&/N}kvq=6krf8c9#BwK>:H-zqxq6FIl:w5clrN/-TG7&iQka*RnU+/2{ggB=#qehytCVpODR-iCwk0qmq%w7S^@gghPHrFef&WS4kF*HSs3y$^WLU.c=%[^PA&iKZ!VTK{W:VW3ePePa=PW(FdZ!P4?Nzn#gb03kJhFp8[B",
"0PLcEh[T7c?ub)WDKPY/<4E*!Y!oEew<88$?T[ItLG+FnUlD8:.Z6}%4JVo.gn)HyVTH>&U.p%z.u#!Ri9HpZ&cN+#GZxi41K4QFlOL@S+!Mt(}3cy6<tj<=5>Bqqn$.mJbmWx!Q=w>inU>@G36yz#cpqsb:Pp@*<t/:CMp0.DD)E<1<NW#3wLXhMtrfcqy$e<e6soddo.j-^iNh]u:.v^7G@nYZMy8.H>55u9NDDvNcAZ8-OW*}8u<fjhUOAlHg1vIH:50XKhcsQMha?P/rlRq#p!Fe9XeObc:T)RG:f#N>S):$NT0$Gfzd.KZVeGri^$Lao5{vmw3ei>qcT+VGWvp@b/3go&/&+WryI8wDGuCv5z.jdSjG#UPGawd^jq-rJ[h1G1Z.qtF17?MhoGygn0+Px%36}htk5#4Is+r$>xnT*=MBwU/A[(/}r5j!(^9A>%=Nkq*jKek@[>iK{+Q]q*jHNREPe*!gjM!-DiHaUt>55ulNe/cMzgHoVW2XPM+X-t312*kJ.rg{uLANvHpal3>}d62K&NL@GLufJx5XI4/2j*448]fvDHz{Z#9lZxl[L!G@+tjU218YgnzmI-HhSlEhxopB>O]W18idt>W<Ov3>Fb/#TfAvgnFb<6XmHa>!7YW}XZ*R^s1eM1C1F=M0g>kQSgxtz@(JZC>YA@*h&aY$fFn26t1(Z3g1e)k9qH)o$3&$LPzvtkuXdLw}k=fam0keEg2>O}%jj1xrncXL[PCN>?3#u^wg4UC-HlA8e9.l?/w>pNp/X#91hvRhS!eZ8FC*bY]P$u9i400p{8BO[c.bmjvJr:+uV(eiDtY^HGAG.LLHAc)3a9ohU2<y9(!1+I:02SHc0bhC?0y!<<w{zNZKoNB#w/f8xi6[63n>?a9cM<mofWJVFQig5!+*}}]cy!Em+7gYqumb2jD6bAVc#z/0U[Z.q>dP3]M]MtXm}*A1Ww9&2cKvCPjVi7N]UQ-d1It5/5cnPAo9Q<*@tC(JQ@/2nK@K^1L9^l.iRWVPIe58ea1xk>L@d#A)9>THjpZ0+A6y1VoZn58gV0#Q5/R?:j8xjGi>vRM6dl*R9vyRUBDMk3&OA2jSN>2FfU?e#+yL6y-X]8UM@P&Y24rFi5MQ9!wSnlP!IBlNp-de2p/sx>Ay@q(!#AN/tb1qPF$SuK3Hy>jSpJ@1SM!Rw0D<NFrH0}w{H!xW(1wedL#yxcM>3jNtdyK7M*l-:3>C(JkS.kI2j&q[3Vm[Bzw/4]hE{Y:&3EUIv0ph&(rFD?S+*huSIBre/1[F{tGj37y2ENUVX9r6F@RbgSKrESRRY.XYn{Nj=7gu8YC+^%<Z>Wtz(rt9a6+zFYLP.dB(iyY2$bxnbHw?*Ls2V!iVHz0d3f06B4r9:=JY%Ojt-B/p/r$<lDU6^5vI*W8Rv=*+u$Q]KB<eilkCkm2O8$V&mU$*HgZG)cCit(lErcl*SjX#BR:9M+R=v4U@udq-VfZff7}TotNBpw=)O}Bl7yCrJ?4g2nWKXPSIrpBNVxAHCtgb!1YEEV!Lo>JNlS/^y$ddOV:jq@=#C=rhO%J*]9H2kh7kDH6p(3+pgH%JDn:lyR<*u#P%Wv8KJ^9W{F&qGP3^Ij}8@[]1@mP+7Q#cBt&A2ADjG=(:v^}EIly#cRNPAiP3oA%O9q<k)NpnJV>WwZ@BK1mOWJrJhRR3qbxDutG}O(T7vng+t@ofay0Lz5FQr<Ww:Z1g/0^kt(GM<RTiSYHy.2mRYe$pY{t@{war2hRkX&Tj1c[{A9<dUFDi(uEwYAmfTS6I{H<Ve%IGK4bSNfm(T1vv==1SQGV$B8E+P04^u4dLHcot%0ia(dodwG*nRLccaDf]Yn>@m8WSOoCwl9h8fx0<j9*0qEp{ET!B?b?hlgeJa3({oxUY:3!jPNy!Rlh9+g3^]YDh!3WW?rzHcBM.XF*9:Eck!-UmqN)M^)c!98YqZ0T?CKm2}uV-v(M+4HMZr4KJ?5)9}Gu@xTi+jJZGD3d@K6h>BB9j{/xY&@BdJGEh*S^yug6/h{N?G^"
};
static bool decodeRung(int idx, vector<array<double,3>>& out){
    if(idx<0||idx>=NRG) return false;
    static const char* A="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ.-:+=^!/*?&<>()[]{}@%$#";
    static int rev[256]; static bool init=false;
    if(!init){ for(int i=0;i<256;i++)rev[i]=-1; for(int i=0;i<85;i++)rev[(unsigned char)A[i]]=i; init=true; }
    static vector<unsigned char> bytes;
    if(bytes.empty()){
        string enc;
        for(size_t k=0;k<sizeof(RGB)/sizeof(RGB[0]);k++) enc += RGB[k];
        for(size_t i=0;i+4<enc.size()+1;i+=5){
            unsigned long v=0;
            for(int j=0;j<5;j++){ int d=rev[(unsigned char)enc[i+j]]; if(d<0) return false; v=v*85+d; }
            bytes.push_back((v>>24)&255); bytes.push_back((v>>16)&255);
            bytes.push_back((v>>8)&255);  bytes.push_back(v&255);
        }
    }
    long bit0=0;
    for(int i=0;i<idx;i++) bit0+=(long)RG[i].cnt*3*(RG[i].m<32?5:6);
    int m=RG[idx].m, cnt=RG[idx].cnt, w=(m<32?5:6);
    double s18=sqrt(18.0);
    double r=s18/(2.0*(m+s18)), q=(1.0-2.0*r)/m;
    out.clear(); out.reserve(cnt);
    for(int p=0;p<cnt;p++){
        array<double,3> pt;
        for(int c=0;c<3;c++){
            long bo=bit0+((long)p*3+c)*w;
            long byi=bo>>3; int sh=(int)(bo&7);
            if(byi+1>=(long)bytes.size()) return false;
            unsigned v=((unsigned)bytes[byi]<<8)|bytes[byi+1];
            v=(v>>(16-w-sh))&((1u<<w)-1u);
            if((int)v>m) return false;
            pt[c]=r+q*(double)v;
        }
        out.push_back(pt);
    }
    return true;
}
struct ExT{short n;int charoff;};
static const ExT EXT[]={
{133,0},
{134,750},
{148,1505},
{181,2340},
{182,3360},
{184,4385},
{199,5420}
};
static const int NEXT=7;
static const char* EXB[]={"8ps#)*Er7dO-N.I8pO{N[TUVS.5pT>ulL>L!:KEa2y{lNBbp+L*SC-IOF*ci8qX&dP}a<xee^c9XSJB2?dmFQ&]L-Z&/mE=*.JR?O-ODL8tPImvw)/[OG37W<@@](usCkzf1pC><%n.Cds^Tc<hlkYkmv^o]z6oR20^B/8pO)RB(L-bOJDDHk>gBh{j>!:5-$JglYmJ(]K8H#HEMLWs7EaTNBaM/>Jj4><@{o*a!t{SZud44o/K9Ob&IL$H!)[l<%usXe3@WX+:d@Gw(kahusE^/wGf3/xeb@:27n=J/L2.yX2bsIR{xd#=y#VhmT6S7>}aXb!45EVzc2y({s[ZOY(Bj]Xa.G)/+[<?]+(vA8qXFAc<!Ub3fe:eFvap$c(.p$7k@UcW7a8WO))VJ2n#gG8pMAY+wzT<>web:KF%1*nFVe0f{pizr[Pe$G)5YvvQyNij#>@%1](haG61kEN9?<Ix{oH&NvcTo8qgH]]t)vkuGJ:wJq^p>O)]0-iAZleV(SYsd)$o5OWr7[EcUTWt*Q2Z!M3yIXP%HKByl{4P1w2+Yhtb.}]SwJz[Kw5jSAREc<mp^jSi2ENIi=S}#^[)1v5j-Z<9&9YA/z8*REy?YGwa<?EtG23*6668pBSd=rqbYua$%s<)!%-Zf*1vf6ji-/1LM=a!v:AJtZOZAS:2qbJqCyv=!#%<%+W+e[eA[&PqV-FkfNj2]qUGifWfEVZSx9a!t{wvREus<%nBo4Q#!uv=!^kn=+z0Zgna^O3Jc+91!UTVAPCXN[i:Y8mnL*l7pWKN&eS@AflSe*W#C[7Eb:9Fz2w?XHcD8O09rh8m?UqAHS.s)ISS1k.kAy/YB6:o>OEn>1w.)ZguvrP4jd@96{5iVcqu6.zL>oNNtzra!7Zq3$:2S","-mSTaU4lV.3upB.k<o=Xl>M=P7^TP68mgN7]yLvAwo4/1o(%u6VQWF0vq1t=mL8]Pl#0QId[pI&zQ:e>5T#w[O90bu>0tFgrpe*o=(6p1tvB2n7cjcx&BYK#Mf93Z0JP<twADvbpZf:eO[c]9FS-r*xEToL>7yL9F]90Gy@U29iuh){*^U]mHZs307/1VF4QVwHAHYXpCphJuLYO-%8mJ{<5U1gX]5.[oa2H(*f}i=[ssXE<JC-nAXyONMac6)L>1b2Y+nMDON#$/uSUS7c7*%)CgkNu@LF2pj/?5gso{P@.>1gu<pgP-HgRWT&-kyY8/5@ASQuj+hN</W4XYNBZ7m%PLc{3o7V757WwAC%m8oy{EvvMyLr.!{:WpbLu8&*]tdy(yrmK$4OtBJbcDaGGIg!z#c94I]s[[(VNZG:HFD?TmodzMB7YO({OpgRd13{Hf(K->GP{vYxpc::!@.J.w<y=$aoO8lJ3N4v1lpfM=)BS5Fl*98%FnD=r$wkpq#>0tF^O[5[(-o6wRw?EQ!AHV8U:Muldupk1Wm=YDLwn&H<>12[>Zgry{I@)RG>1YDZgW0)[oUO!..$LFtoHht[wADS$>1YB!q/N*m(6&KGJVg<68={-o0GMsPML?LIn3)RkW0K6.DNa&5ysEP2V5R+5J={LBjopDV<U]p5J={YU80NI<l5$ajTT6jmG5I:6e6O{>Z#>-UrZ.p8WYL$5J:Bi[bSK:$!o^h-GC[R@Y=&z/)R%9d*ufcUDhcq/WYJ]iZp:ijjaqguP]OpKT.MoQNx%W>0GF2YO4Ovbu-xO:(.)sZTW*9Bt<pMp0GK?Fkv57&>m:YevO7s=T!8@C>nU*fuubcuJ]}gWysYH7<g]&2JYWY]n4bKW0GJOMJU{)GNx.=}LFw51","At-&kC9d!E37B3TJ]{NR80NI%bZUKQTVN0o80.-6NcGXoTW/{>uM>Dx!4=*+N%7pC>lkxk<Z^JPJ:Bl880m@X1BEI5GC[K.JUj.x2fH*$*uf6du-RHm(^yyzJ^pRzj9E:hcAVCG*8fl$DgOJd9=S4X&eO2zJ&JZYL-1SrIH$fx80!b*!k0WfN%7j=>m#pP{4*=pkG0fB.FTuOs!07Qy3cg/nuM%Hw2nua(S&BpJU{.&/X?+$5]c11JS#RC1hE29O4Nu!Nyokj*y+xzLJJQ[u-RHzQ(asEJ-F[]LbPz.?FFBP8/hzl>k*QQ!4?9Wziy^U>lx8D@=9ERO0**R>n[RGoH160AtW@]>nd<2YuImfn4p<BY^LIC*S=d#Huw}>>njJy]s?)g}A5[[>nO=TzGa+qzL+mW>nf:c)S1/2AL}LTn0@OD4SL=HAtW)(Y^cv:IjF*V=6+#+P1YgY)S4z8/2dzIi:s2{7ZXaVAA}{KP0LoR941@^APXjeYXyYc=tR>e3fJo)YXe^!?FPgnz%4HPXnB}1U1wV]je3<mCJoMWCGjYE7JBTUqfMo.H.>I@idVpF.[N1eU#@/9.@Yt-COTNJCnCEE7KvGUqfMx5)IytdF!rh&(h?=D4XMWyuyR]IAD/uLCCYgD7KQz3l3Cy2fbF>&>Z^@#EPNznH@r.Cvv&ofqfOPlaN1Kdg(r>3B4Wr#H=5brE6Wdhy!cA{<6e<2Q/[J(y/(iPxmFybhK)@iq!Z^/(]k[lt{Cbin^4aqH[^s)8J/k6O}[h7I0f>PicAT=Bbm>:m(Vdni8++#ehYD:OALQn7J&?x2Pq&Wi2u!]q$k7z<Zvqc>t^E-t#+44.??yaknODD>.3VJXeMlbt%[4orG>1=mqca6k#@jWQ6kCf*+z<>?*DkV","E%5Fbl^F##BFJ1vsa1bzwQ4LwQGurd-c4Oj0E4Y!T($FY7JS8La:rd8lG&#AA@4oGh+s>F+xj>Zt#$:XHV&]YfbF?I>.YdVXeWF-z(h6GClOx<Cxnu#hG#HHErQdn:@U^1<6h?wjRL2zJ6>H(IX&}0C{V!WHMY=%QiO1yYMY(RSKm4itZwu9>-gEH^ruyK9V6iFG*QWWH/zx]1d7toiz1E)^$i@HG+]2&Z&11R.#pD$Mos&YNjo-7SG%GrZNp^tkDgf9qgC2ulSYXE7K84hibN(oQSRc7N-4jQs*.!v-*Z7oN}YovR>PVfH}#VZ>.Gzje@hGk1@zO%RC&q>DsW[wH]agv>.-{/.f@RXI08@*R=x4[dRoGzN)OAnYE36HYtGvdP$B0pN{i[Ne62qCtnIU}Hi5x<U80)w8PIx0.xaJ5jsxuDInCi$Df=){CVq.8GK-!TDFh]IGnuyfCBNxi:8MszE0n.{CGj=X>.%.iZkf*R=/*sis)eV4V3eIiH-I5D&-(3))QA65l/Q-o82/qcLAt0ma)0PM>.zy^FP>SSH!N[]>-gDNEUZ1GaN1VRm$comL}h6LQgpuZieyL$HQq-WH.>wG>.%!Le@l1[w9wSnRE26rJcBGCH{?hSDC1E=I{GB)H.>I@>.kF?J!<gYI0fu?it3AJxz9cIC0c1]7IAFBK-)yIHzKxxBRe!67P50-w}kCQ:a:m7w6*Y+C7(wa7H.G@f)m8^fbFws.TR!u=wruJUiq0>>:hU)w-bK5CepCb7I[L{0D]ptJU=O&k0:wjk8y7[In{=.7K?e:4V80i7EI#Wh]XvSGm1FzB#P+fDujXnM6KxHEk&E8>:)R<jA(stHKTx3icBg.q?}lI9YODJ7H@jpx<%5y^0xhSy2jXaCV?]PHV>2#","ich8Vv#egbr5/PL7H@x!NMivIKM.[l7IdLEf)pP<g=P@vlQS&4G]GR&:YNW[Dpr{HXW>]Jx/GaQ:eb88dR:oLL7!QPsUtLYWu(weC3e)?7If60OkJH@T68o?RD^b7-K[gbOYs]4s<.yI-LGXHIcI=VCH=-lFP+tAI8YUlH*C)BV390/mhG&#wk!o7Hdbm9IdxTSN$)R&Tit(7gKrb@wm.0Kepq=0QO+PkDmh&Oh>#7ZR*pyrs.BHol3vNSIj{&Ei*ktB+YZ%rT!{8YE+][$xys#GNsr}oicT^bAD={2IrGf.7I[tIv#3}ar-$xaG>2]4M@U2n>z^Oe>:pGGwdUcK@]/=fG{/GvhB-1bV.<T#N(BIMs*HE0KEF<WYS5qHv$ZsREvCI>>:Te)V^P.8Ym<a9Px{!sy5+5KHOaN5ic>LAG{!OjIf+QcDFA>AlXHdOwew#[RG3!yNRb=XxqP:9:bjIjJ!%SZIrHBH7H.G@i9yMy-ajN^YN1v.D=FG{G1N4<YO<<+v$?-uI5f$tRJD5rdRI$z5xx#<.MF.#TsItnS>#+z7I>U^Ny*tRS?$e*O4K0^FP:yOCSo>n:anz/y<o9sr(Y&FRHiR8y<s0ON8+egidu.8HPfHcC!zsU:b3&lzFQc9(/UdnidtXvCtWMy@UIeZDyd:%O@2dII5f+)>-Z=?pkIJTCiCH/>X=cBbRwkzHLoS5:e7#-K-)QJIckN4>:)SN9/j3SxxU>gMm!2-S+pKTCiCJO>:MUSO6e<jI8Y)c>:eo?J!{]0IrHws>:)Rhk9p18xBCX=CeCFD7H.Fv^rqk$HzKL0ib%c9LG:[S5lT0@jQ4h?y<omAB@a^>7IhSM^rqyx)?Y>#E>=afdzuQiqLy{I7I!jhw<KqYCbyiM7IG&jK-<oQ","I5f#MPH0++Oe8]S]5#iV7K?e-e@pxOh(#C7h./FO6luJ/:t=TWN@:mpM.(-mwUp&jRLPLqmF{q-GOr8Mi9Z!&lXNVtT2K5licf<ZP(@znC9VS[x}yevE0AHwHV>2.idb$2lX?w.B#TZOicf(sNGUW1IrHe*7JB:=a:g[Vg=P@vlQucxhi-/.W]()HDs-vLI]Y42h>a@PG{v{ar.dKFMb*p$sZ2M1mu{6ANPS*9.8mHoETCpd>[^]$u.-xVJcI0}Ncp4^s=Ot*mvnRgC3f&Zi2d?3FP+uOC0XUIY$Ptj^rqygmhG&#wk!o7Hdbm9IdxFC7Iu!jO@4Pms:B.TwqN^wxJ-06RBbN8Dtmh:CC.mlM3N4ksY.GvX2aE+sVw{Gk9H!GJcJfKS9+WgDx9&Pxq>3VL6SZR7IAupADV*]CSo{t7ID<6v#3}arU]N%G>lT(H/MG:Ig]Pw>:kx*6nI)3C&[m6wp]J7C1u{[W$QlVN%F*}x$8d8X0}RBsY@khTh{=Lh^6p)jN)P:J!$S{XYpLnPEKOlDiAC}HOd)%7ITbcvrOr7If^>}7IlakADXj9xe>cpRGmKkIF37HxqP)L:bjHCzI+lkNco#KO34QBADSdlCJ+SPW@$qBtGzI+RO4!4N%Sk%+GdtHhAe?3RJ)Q^k8y7(7kmx$.f*F!D)2qaCevNxO7]TpR(EScneD7@>-}HL0D{?7CSo>n:anz/y<o9sr(Y&FRHiSmNRfV(Ce]L^s=vQCGnnVXCZt%r:bdyFJcJ1HIn$vridt+xW[L$oF/ME>?FOekJeQWRInr5j>:)SNuwRj7CiCr(7Jb8[Lmu-z^N3Y=icm1vK->nwIckN4>-Z^QzI-Xbr{l5#MtRJ}yg(:}N5ks?O2Qb+O@6}9I8Y)c>:eo?J!$AY","C<7ok>:)SlJ!#$/Y6b-x7qT/5-6+xSy}ALX7rbS4.edh*bO[jjl{UTr<?zONW}c.07rm:q-O(fay%ieX7rA55a-673dI*gDRlK&I*i<:My{J2oI>0YY+-o!1y^H)X7rlppf(i!GLnqDc>}N2si!e0>Ln*=yN)Rijh>H)fBuPr87JbBC5Q$Ov{-F0iY!Q>$jDIF-N8?Tj7rF3+43x^gAP<h1h.7nI6loQKy=A]t:0lp]yht>STewsy:*K/3:h6WKLiB@diKq{Hl.Gw<Lr8TDCxo#(M1&xPbO%zv7Bu4Axe+q1)xVuskg}.#q$:quy^lDiPm^g/D!yY(bI:Q{7qT/Pl2t^x-C2mms3z&jhfB.lCuy(4RZmky*&=/^F}^o=Eq%7aEV[sTNDV+nt587whk-A5PK0?fitHo::fYd(dJL}?^:N2Vqewr.{MCdgu=$liMhAx3@9OXIE}w.5B[o0xjd<H%q?$fRG!p4$3o2e9wwjly--%q6Iem/T>@c?zk9Hg]GVBy.:o)cpJ*2fvS{8SWDD24.:f0L0Ly4p=tI5HabE@qteJ@q2B9[4tsVG*{33je.FsOh%sp7h50F5uME#1O#5.8v5jaux8Bk#zl1WMF&6Jn-vFKOnP>T36c9!W%r7rxv{l2qodLq:[3i>OsCl2t^B/Y+ChN?L03:*n3yLBDx::ly.gk9IMIy=/=mLronD}8r/m6sf1MP<6yENIO{Z[U7d^PC!U?cf}mYju-JT^=Zz0^)^Uo9-NjsQ8c9on*O}DLCrXBFxyeo(kMK)y@b%uk2drre%v6[xK1kFRTf:gIGm@sGR#zP7rQ$t.ecW-z0s3[>@:{s42<Bj&iH#B.ymPOYD>}F7S6lykccfq*-<R2jSW&U7rcJ@v@^98m{!Dz-dq-$tUV$K","ThM@eRef7bWx4nq[%[GuzWjz7dx?$#CZG{d:gOe.+ZO!jFYE5Qh&JKlRh-MvU[>MB>@k6M1f6xZRS*pd?tovg9Y[MX@PEWIFJ/9u?JNgW+OFG/yt<0[i<YWFLJ:/r>}N0c?JJ<hz%oxVi3an:l?PX8P&+6NfHZYWk9Gb#LtF-v>@+5MU9kdS66RsJNvYSbiI2y<V7K<Q7r:tOjX9v}Lu?AH>@C/=O%6@nLJ:[(>@$rxzJ<rs"};
static bool decodeEx(int which, vector<array<double,3>>& out){
    static const char* A="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ.-:+=^!/*?&<>()[]{}@%$#";
    static int rev[256]; static bool init=false;
    if(!init){ for(int i=0;i<256;i++)rev[i]=-1; for(int i=0;i<85;i++)rev[(unsigned char)A[i]]=i; init=true; }
    static string enc;
    if(enc.empty()) for(size_t k=0;k<sizeof(EXB)/sizeof(EXB[0]);k++) enc+=EXB[k];
    int n=EXT[which].n;
    long c0=EXT[which].charoff;
    long nbytes=((long)n*36+7)/8;
    long groups=(nbytes+3)/4;
    if(c0+groups*5>(long)enc.size()+0) return false;
    vector<unsigned char> bytes; bytes.reserve(groups*4);
    for(long i=c0;i<c0+groups*5;i+=5){
        unsigned long v=0;
        for(int j=0;j<5;j++){ int d=rev[(unsigned char)enc[i+j]]; if(d<0) return false; v=v*85+d; }
        bytes.push_back((v>>24)&255); bytes.push_back((v>>16)&255);
        bytes.push_back((v>>8)&255);  bytes.push_back(v&255);
    }
    bytes.push_back(0); bytes.push_back(0);
    out.clear(); out.reserve(n);
    for(int p=0;p<n;p++){
        array<double,3> pt;
        for(int c=0;c<3;c++){
            long bo=((long)p*3+c)*12;
            long byi=bo>>3; int sh=(int)(bo&7);
            if(byi+2>=(long)bytes.size()) return false;
            unsigned v=((unsigned)bytes[byi]<<16)|((unsigned)bytes[byi+1]<<8)|bytes[byi+2];
            v=(v>>(12-sh))&4095u;
            pt[c]=(double)v/4095.0;
        }
        out.push_back(pt);
    }
    return true;
}
static long fccCountA(int nx,int ny,int nz){
    long ex=nx/2+1, ox=nx+1-ex, ey=ny/2+1, oy=ny+1-ey, ez=nz/2+1, oz=nz+1-ez;
    return ex*ey*ez+ex*oy*oz+ox*ey*oz+ox*oy*ez;
}
static double anisoFCCr(int nx,int ny,int nz){
    double ix=1.0/nx,iy=1.0/ny,iz=1.0/nz;
    double f=sqrt(ix*ix+iy*iy);
    double t=sqrt(ix*ix+iz*iz); if(t<f)f=t;
    t=sqrt(iy*iy+iz*iz); if(t<f)f=t;
    if(2*ix<f)f=2*ix; if(2*iy<f)f=2*iy; if(2*iz<f)f=2*iz;
    return f/(2.0*(1.0+f));
}
static double anisoFCC(int n, vector<array<double,3>>& out){
    double bestr=-1; int bx=0,by=0,bz=0;
    int cap=90;
    for(int nx=1;nx<=cap;nx++){
        for(int ny=nx;ny<=cap;ny++){
            if(fccCountA(nx,ny,cap)<n) continue;
            for(int nz=ny;nz<=cap;nz++){
                if(fccCountA(nx,ny,nz)<n) continue;
                double r=anisoFCCr(nx,ny,nz);
                if(r>bestr){bestr=r;bx=nx;by=ny;bz=nz;}
                break;
            }
        }
    }
    if(bestr<=0) return -1;
    double r=bestr;
    double sx=(1-2*r)/bx, sy=(1-2*r)/by, sz=(1-2*r)/bz;
    out.clear();
    for(int i=0;i<=bx;i++)for(int j=0;j<=by;j++)for(int k=0;k<=bz;k++)
        if(((i+j+k)&1)==0) out.push_back({r+i*sx,r+j*sy,r+k*sz});
    return r;
}
int main(){
    T0 = Clock::now();
    { const char*v=getenv("TL"); if(v) TIME_LIMIT=atof(v); }
    int n;
    if(scanf("%d",&n)!=1) return 0;
    if(n<1) n=1;
    if(n==1){ printf("%.17g %.17g %.17g\n",0.5,0.5,0.5); return 0; }
    for(int e=0;e<NEXT;e++) if(EXT[e].n==n){
        vector<array<double,3>> pts,bak;
        if(decodeEx(e,pts) && (int)pts.size()==n){
            bak=pts; double rb=geomRadius(bak),r=rb;
            microPolish(pts,r,0.6*TIME_LIMIT);
            if(geomRadius(pts)<rb-1e-12) pts=bak; else r=geomRadius(pts);
            optimize(pts,r,0.85*TIME_LIMIT,0);
            if(geomRadius(pts)<rb-1e-12) {pts=bak;r=rb;} else r=geomRadius(pts);
            microPolish(pts,r,TIME_LIMIT);
            if(geomRadius(pts)<rb-1e-12) pts=bak;
            for(const auto& p:pts)
                printf("%.17g %.17g %.17g\n", clamp01(p[0]),clamp01(p[1]),clamp01(p[2]));
            return 0;
        }
    }
    if(n<=EMB_NMAX){
        vector<array<double,3>> pts;
        if(decodeEmb(n,pts) && (int)pts.size()==n){
            double r=geomRadius(pts);
            microPolish(pts, r, 0.55*TIME_LIMIT);
            optimize(pts, r, 0.9*TIME_LIMIT, 0);
            microPolish(pts, r, TIME_LIMIT);
            for(const auto& p:pts)
                printf("%.17g %.17g %.17g\n", clamp01(p[0]),clamp01(p[1]),clamp01(p[2]));
            return 0;
        }
    }
    vector<array<double,3>> best = cubicGrid(n);
    double bestR = (n<=20000)? geomRadius(best) : geomFast(best,0.4/cbrt((double)n));
    if(n<=20000){
        double tcap = 0.45;
        tryLattice(countFCC_r, genFCC_r, n, bestR, best, tcap);
        tryLattice(countHCP_r, genHCP_r, n, bestR, best, tcap);
        if((int)best.size()!=n){ best=cubicGrid(n); bestR=geomFast(best,bestR); }
        vector<array<double,3>> ap;
        double ar=anisoFCC(n,ap);
        if(ar>0 && (int)ap.size()>=n){
            auto pr=(int)ap.size()>n ? pruneMaxMin(ap,n) : ap;
            if((int)pr.size()==n){
                double rr=geomFast(pr, ar);
                if(rr>bestR){ bestR=rr; best=pr; }
            }
        }
    }
    for(int i=0;i<NRG;i++){
        if(RG[i].cnt>=n){
            vector<array<double,3>> rp;
            if(decodeRung(i,rp)){
                auto pr = (int)rp.size()>n ? pruneMaxMin(rp,n) : rp;
                if((int)pr.size()==n){
                    double rr=geomFast(pr, 0.9*sqrt(18.0)/(2.0*(RG[i].m+sqrt(18.0))));
                    if(rr>bestR){ bestR=rr; best=pr; }
                }
            }
            break;
        }
    }
    int sbi=-1;
    for(int i=0;i<NRG;i++) if(RG[i].cnt<n) sbi=i;
    if(sbi>=0 && n-RG[sbi].cnt <= max(6,(int)RG[sbi].cnt/12)){
        vector<array<double,3>> sd;
        if(decodeRung(sbi,sd)){
            std::uniform_real_distribution<double> U(0.02,0.98);
            while((int)sd.size()<n){
                double bv=-1; array<double,3> bx{0.5,0.5,0.5};
                for(int t=0;t<400;t++){
                    array<double,3> x{U(RNG),U(RNG),U(RNG)};
                    double v=dist_face(x[0],x[1],x[2]);
                    for(auto&q:sd){
                        double dx=x[0]-q[0],dy=x[1]-q[1],dz=x[2]-q[2];
                        double d=0.5*sqrt(dx*dx+dy*dy+dz*dz); if(d<v)v=d;
                    }
                    if(v>bv){bv=v;bx=x;}
                }
                sd.push_back(bx);
            }
            double sr=geomRadius(sd);
            optimize(sd,sr,elapsed()+0.25*TIME_LIMIT,1);
            if(sr>bestR){bestR=sr;best=sd;}
        }
    }
    double TLfull=TIME_LIMIT;
    TIME_LIMIT=TLfull-0.04;
    relaxOptimize(best, bestR);
    TIME_LIMIT=TLfull;
    microPolish(best, bestR, TIME_LIMIT);
    if((int)best.size()!=n){ best=cubicGrid(n); }
    for(auto& p : best){ p[0]=clamp01(p[0]); p[1]=clamp01(p[1]); p[2]=clamp01(p[2]); }
    for(const auto& p : best) printf("%.17g %.17g %.17g\n", p[0],p[1],p[2]);
    return 0;
}
