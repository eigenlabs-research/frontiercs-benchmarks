// Sphere Packing in a Cube (C++17). Embedded best-known configs (n<=120,
// Specht/Gensane tables, 12-bit quantized+Z85) + lattice warm starts +
// push-apart inflation with kicks + anisotropic expand-to-fill + fine polish.
#include <cstdio>
#include <cstring>
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

// Spatial hash grid (head/next linked list) for O(n) neighbor queries.
struct Grid {
    int G;            // cells per axis
    double h;         // cell size
    vector<int> head; // head[cell] = first point index or -1
    vector<int> nxt;  // nxt[i] = next point in same cell
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
    if(mind2>1e17) return dface; // no neighbor within 2 cells
    return min(dface, 0.5*sqrt(mind2));
}

// Jacobi push sweep; returns max pair overlap seen (pre-move).
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

// ---------------- FCC lattice (count + generate) ----------------
static int countFCC_r(double r,double ox,double oy,double oz){
    double s=2.0*sqrt(2.0)*r; double lo=r,hi=1.0-r;
    if(hi<lo-1e-15) return 0;
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

// ---------------- HCP lattice (count + generate) ----------------
static int countHCP_r(double r,double ox,double oy,double oz){
    double s=2.0*r; double lo=r,hi=1.0-r; if(hi<lo-1e-15) return 0;
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

// ---------------- greedy max-min prune (farthest-point) ----------------
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

// Recenter + anisotropic expand-to-fill (factors>=1 cannot reduce r).
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
        {0,0,0},{0.25,0.25,0.25},{0.125,0.125,0.125},
        {0,0.25,0.5},{0.5,0.5,0},{0.375,0.375,0.375},
        {0.1,0.1,0.1},{0.05,0.15,0.25},{0.2,0.3,0.4},
        {0.45,0.05,0.3},{0.075,0.275,0.45},{0.325,0.125,0.05}
    };
    int nof = 12;
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

// Inflate-bisect optimizer with basin-hopping kicks.
static inline void jitter(vector<array<double,3>>& p, double amp){
    std::uniform_real_distribution<double> U(-1.0,1.0);
    int n=(int)p.size();
    for(int i=0;i<n;i++){
        double x=p[i][0]+amp*U(RNG), y=p[i][1]+amp*U(RNG), z=p[i][2]+amp*U(RNG);
        p[i][0]=x<0?0:(x>1?1:x); p[i][1]=y<0?0:(y>1?1:y); p[i][2]=z<0?0:(z>1?1:z);
    }
}

// Move the bottleneck sphere into the biggest hole (targeted kick).
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

// Try to reach zero overlap at radius 'target'. Returns true if the max pair
// overlap fell below tol (config then supports r ~= target).
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

// Continuous Lubachevsky-Stillinger-style evolution with basin-hopping kicks.
// sched 0 = gentle (fine polish), 1 = aggressive (fast growth, frequent kicks)
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
    if(sched) jitter(cur, 0.05*bestR);   // break lattice symmetry ties
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

// Fine radius ladder (undo quantization / squeeze last digits).
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

// Random configuration inside the cube.
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
        // gentle warm polish, aggressive multi-starts, gentle final polish
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
    // final expand-to-fill (never hurts)
    double r=expandFill(best, bestR);
    if(r>bestR) bestR=r;
}

// Embedded best-known packings for n=2..120 (unit cube), 12-bit coords,
// Z85-encoded. Sources: E. Specht's sphere-in-cube tables (2013), incl.
// configurations by Gensane, Schaer, Goldberg et al. A short relaxation
// polish removes the quantization error.
static const int EMB_NMAX = 120;
static const char* EMB[] = {
"q6=&Ld^:y))mN2Oo8zC)1H2ioUEaIpo7Yr{3S6E)p!HQYo8zC)1uw&liE?HSlZ#7.M&GHv0k99!YLCJV1rG>&9?JMRZ8M3s]D3r.[v/UGZ7}l[8@WlI06*U!NV)/:YsQrh06*<&Yr?NIV&E.IYtamyZW}H>1ouje06*U!ZW}Jv}#n[:06*}(ZW}Jv]z(V6%g!EhZW[(}w46stWMw0c-fmTFOjL!Q0dUsPj7C.qOjL!UnJQ-W-gfAqOjLT0zz{=JhSFWk)2<o(nxA8]twa1!43DbD05&>P:TP?t)2<o[Lf2gA:T&Pjbt{xEOP.tqh7U<!euN7zqP^D^x#xy&2MP.pT+/zyh7U(4B5HfoT^h8Kh7U(L*]-Hpp=exVh3lW.15NebUqDNHhH1zmOs1+}UyD(Th3lW$<50[@Ufcuj+j>{^46N?KUGDDkh3lW$]h4pI0dUeFh4lVC46N?Bp>eTAFdeSi2MPw-UyDZKh4lXh]h9Q80dUGYh4lXh@iaovUqDzUh4lVC15KHN05T<P+jwDL@i4%-05U2Zh3:hj46Lf!0dUGYFceR946N?Bp%eLkFdeRs46N?KUGD^Df/6E$ixw(J6vaV0f/6GbQz5d^([HS?=C)/F.XmL56vaV0Flh#9qXH(Ci!dqC=DLS>oZBF[0m5Q5=DLU2+<rS)1U-<U=.yT$s1nx*1*>[&fJee3QxOk:qe>YkfJP-y(.[}t1-*IQ=-D#3tHCifolz3XX}pf[Yv5xi1<:pDI^lapt(8^AnK<A$fJ++qtAjRK1:3phfKjDdrb<P>ADEKsI!tQnt(5MiJ)BQ{I^({Gz2emvWLg*(XAuDfTbYW-1<:pDfKjEjxgZBJ}NC[*fCsEv8lbw70CbvkfDy7x0%77}@Xbj&=/k4qqVRyS@Xbe?fC<!u395HY0a/BmFdhC6>M0pI0NH2gfDy76socXl@^UBRfDy8r(8J7@@X0i#f6IvG&{Lhg3OHDoJh2-Yt(Dwy5Z.m5FjGor67jc.[xv7nE$^bt>LANz[{<?k^ht[o>FJRS4:fn3=v2/x67vUt5SEB>^h-*qEgZVKQ0TMOeWBb+>hXb1*LaUw^N<S]Xb:(+@M{f3fK&!+mF?!#*IP.h^NQ1qHRFiuE19KmeWBb+>B{r}8}47cPHe}f?2h:F1KIhzu+E-kn1[v3P])2lX#^malE).$/3Y-jeN)qVLgF6Ldj]ikqLL4U3E{>L7K#Vs=whD#R7aiA0iGfoFcVd92MN=7>ZT2<=whDPCMs<J>ZT2<FcvBFHYBLZ>Q*n{^W.*sz7d5NdBK+o^W.?562(Fi-+Mjb!0S@Zx?)cu-)e6if{rH2v&o!506@@>em#foOJgE506%cQk$z0a3*SrG0dUEyFdjH&x?*i@iH5.cFd7V$LnRUZ>fv/9L&xZiOJgEaiv[(t!0vn/OJek5Z:]:ekR*c&:l1]Q1a=SJERX<j4l@tjZ%?>&lO=Ea]Fq.LZ<(KIk6N7R{pbROZ)66z.hGjY2dSdS0dT+zMhJ]#3^pdhkFj>GedB*{LMPzi@c&GEkSxUg{VxlnkHYM}!a-[[QpI>JbRb:TFc9qjTj]ea5Q)S>!d1X]>Jpnc)9[}Alpenl29}E3Zj6j0e2SD$r*Fmw0dUsPe3hl<4mm(V8H@n=MyQC{Tj]eS8O)Siea2PMQWuLKmsU?oFc/@tr3Y$0mzlx6!kAF1[#)Nhr4vK-DgFe?2L6wnUGlpr^>axC{YFm8^:XA9kUhDiH4r+3+Y.BrF5!7s60D^LjjKPiIlfQ7LxekG0m.mkI5oVQ8-Vu(p@pENd/yuq19h][UpxGEh-H//4n.P@p+8Au!D=Sm4nXKWT17!ZdZD9eso6?X0mdyud.O!hY&rJGTagnLfKmUd0%KnY?:rR*F^pg+3=WR)}yn#PdZGqe4ogVL05dXTFcT$0M(%4O?*#9<F^YNFXb:<N0m:^tdZDb9o&JQQrdC5O!Lb*&4oa*<rmKYx!Lf0)VA2E+04LpP!K3q^VA5GSTagDQd.O/i"
"0&)/:r4txT!K3s!o&JQFT17QVFc8^44oa*!04LFUFcVd92MSyC0dUcK!KPW<4od&!rmKItdZDb9nb$cMrmKso!K3r!4od&-0m:!yFdkEcY&rJM0m:!y!Lf2[mLug(Kspw/dg9!1&/I:!=Nx5wO-k:3@CKa]zZ+:HdgHhOOP09MfXcBZv/?-v1#6Ivh!dcOOPV@Z0*3{M@mpSvdg2T84r0W4=^}{rOCyloDdMHrKE+TL/7[eMzcUpR=NxrwdfB>D@CMfPz:Ohp/7*wS@CQo=z}sGvdftQ3[@Lp8fK>fC/7H^D0*5@yKEVc)df.GS[@Nt.KyF-rdf.EP1#6J[=TNekv*QE>4r4^gfEWDHv*QGj{aN]%KEVc)dgasg1#8OCz^%a1dgas<}oLPb=TNA5OC^xS}oJKLKK&#:v?0vD{aU7>fR58W/7)Uy}oNT<zZ+:HdgarH0*3)@fEWDHv/?-v3c$4b=Z:M/v?0u94r6?*KyFuNv?0u91#4FgfK>Bnv*jR%3d08<z>cIPdgHi93d4h@fEW6+OCyFW3d08=fR58Wv*QEg3d2dq=Z:*Sv?0vD@CQo.=TN3r/7a]j4r4^gfXk}KOC^xS}oNT/=^@5kdftOW3c$44KyFFGdgHi90*1&q=NxgDOBR+K4r0W8z}sGvv/?:44r0W1fEWhWv*QE>1#6J[=TNA5OC1S=4r2..=^}!zdftQ3}oHF@=TNekOC^wo3c$4b=Z:*SOCyGv4r4^nz}s:gv/?+y@CKa]zZ+(A/7)T44r2..=^@5kOC^w@@CQoTKspl)OBR:<1#6J*KK<aVdftQ3[@Pyq=Nx5L/7a{N@CKa?fEWhW/7a]j0*3)@fR4?<v*jSS{aU7[=TN3r/7H!/@CMfEKspHZOBR=i{aP$.=^}!zOBR=i}oLPb=Z:X./7)Uy@COk1fEWDH/7a{N[@Nt<z>cIP/7H!/@CQo.=^}VG/8oJS@CQoTC{!=2CcVxvORxk%D6]DocZ5:#VuY$gN<TT&/KxcUy1+lxDdVSwvZW^@ymuPplE(Q2wCL>?ctWzfHg+{?xlG:(NuqyDC{!}]zSk@f69cW0*gFF0O?kX>n:q(]y[Bq:/LhX:vTj:JD1ny]/K^}g3m8jBHs*u$/K^}g4tD$nDdV&+fB-+HCamb[D87[7cH#74{&[dRz}erncIqE>0=M-%od6C=OapG6eS.sTJ[AuhcI!r4/Y3h!?6S8LuG*1]20=lRz!5]/wemEpSKMF=^WL$Y/-D8lY]Z[MW7N!RN$YNbDaki%Zmjka/:kAMA=mDrzZTlJJW=lEb5Z/PG4P9</YyKd@G91Zk?MQdI+PpmnX-{Wz}LOYHMdzG4ozc>Ks5y)cFeqB^RWUYz6)XFBAhkG2s*b%L5$.YcF:?i/z%OGzp{I7cE>2K!y}0)zb3&+o@j2OL?Gm:zp{pbFdmV-EuJnhm<=<*FcVd92MSKa0dUb/cFk?U0=r<PK%X?>q!)XR*r3T7zbU2eVrCvMd<VK<zp{rqcF:<rb5M0SzhTypFc-WreS#3KK%Y1F/^D<NeT1^sLg.3&cphme[#Tlc)oW8e!H3w8ls=A(wt8QdcpS<BVsIU4N<H>0cpCi/o}a*CjE!n!Nz@-qqE<a:O1RB6s*^Mxmv&by6J%cgcp^<:WC(W/wtdSXum?T%A+5H/vXDI=Rx$QRn[uR^5#FcdNAy$*hAX7iN<H=&cqx#pVsCpRO0NO{/%QcL?j99aIp]hpMI)qDdIQNAKKU}L/$A)dVskd8vWrf^cmd!S.EziHMl[4]JBDuzdqpaORp<paMyTM[eUfqzw7/ptn(xb8wT)7va6DNMcmd/p]{/t$v+Gs#*1&ikxdIV^N5iE]zvJdJHBTY2xPlO[utgHvvQ+TPvKrKr*0Tt}MF!Sj8M3/hP$hdW3s}&z8(uIy*1PXAzjPwhMl[hOWw4f:^QC*OO.qOW*18KhY}3G0y1.hVRD]-tb4q(Fy9>clus/PMeUcPNvWrzicl=]+Qgh*gsTCz>*0Tt}T^{^<"
"OHbRVcmd!SWDrRAN{04autgG:WDn!/OOqx2*18BeT^%fvwtbsE*18EfQgfBQ!MGrxK[^DNdqpaMvWsc!cmJXE3s]<0wzbVOP{BSlo}v#Dwf@[0Wu^78eUfqjv+Ggsn]pcL!%j}tMu69#cmJY.*8l+ZOOqx2*1&jL?jlN=tQvJ^cbJ]9rKj0@QGYGUcc2i@+1N[>t?RqGccl[Qgf!4{pBCdl*9mnHqE0A?QB0)=Q@xVAj![anVN9-Ru}80!huWru=$%7xQis)f.DZ#XQUmS>cc-iBS.s.X3FEwzKFk#/n-#.VQN(-:z3?T{VvnRFO22)RtBFA=9IT}8vj2E$*ariReUXw9t?R7Mu})AIhTo0Kt-L}YPsKdsgf&J@QJ6Om*a(?TQfr6Gs%ZqEc8aNv!/Dcq^Qzcsc8]&HR#Sz4th0crc8AiEqD*gOR6R@&c8I-=b3FBmt4cLmFcE66)wjmN9-yXIvD#pj79j^7Rd5M/*eYoGu8p#?s%ZfbFc5]C9I>6Y&zOkFO/l3=(er}/&HQ>^O/cIvs(C}!e.Tv^TU3xv6!H-lRl7}B*eq9.)wksq9[<[@vEo)?*2pT*^W?FQc9s2$sod5XRp]jK*eSkv5(Cuiq6^l=q$r^i6$LXeNrFZfu79h*5(E8uql25Eb]%aZQe>VRl@o4zMgW4Io%dCl7PvOAb[:=[T/*yt!7S#?weCm(7iFnZrCctH*ud:&WPF9Tp)}UFb]dCP}s@Xda>(kINdNWg8miG3Ub6)RsSgXB0-t&OqkAG2NRwz5{etJ0A[NN@!FO8n?k</5ytd?:b]%ex&.!B?q6*7ZR8s[0&.:S+Ug&9!*tU#-j*3a!n*Hqdb+^]f0.MOyn>YBnvs?9/F>B^1^iLKB*FQT]30]gxn$Qo*b}VJi!RZnso^TP.yB#0FX=FZjA)(5pPlYRNY%VP@%cqOoMO2S1qC8]Z{1J=WN-oe)r-bn$&HUF9y>*rheW:tqkvN+wv>Ca7QdTfrn#A$+*F%km])+Z@JB]R0n#X5zQ]JfEn*Hwx*FQUs!QeTAwSk5f*1<4YT*I&!WjeNub=ezSgd{q9aLP^^*FQS?lpr7*+P]Iy*P{P1n+Dtel[^c8bU+[k^Nm]Jh6+HR.Et{oYSiqEGwiX1Q@2=5r^v]*dX]6!qTI87Wvb^.Yio!OMq#^do$sZ)/w0.ZbT=SZ>5M<vO7>HgFo*Y<i1$5Lm1IjXb?Le?Y%KmMm.EB1bU6..r/$}P4<TFfQJ#u%DvBRS6@YnibUEVv<{[PW8c}W5y63/TY$k(vYd!vi*R1C/)yKSnH1<5SO]8874zgx9gUYQP*R7kNVoY-Gl54cXbPL[(5<YIzJWbl6JJ$<)NE?:uI>/n$OEGP+5#jdZzKTW!cYL&HQ&CeRZde@o*Vod}6}wy/ll)qoMms4.1pZ*Jpi1oCw?y.x7oVAR]LLpc*VF*bgc{fXq}{IPbOAz7JLn/h9]37@L[$Bfs!d?eqKXLkufTTNeXX8ql8DX.wFdqiT?D51J[-*GXAQz/9FV7}IRDn1dbzQJWzyR:SBiANbO}-X<X?#%ZrQ2!*U$JN)y[rrZt4TjbMZl55<NYKZEV^rbNh{E3e]-:k?(8YbMZlub0..sAWhae*WWQr6PBZ(5T3R.vi99neGD8v5Itt*vhRp2!!crX/<LJML4Cq&}6GI<9K^vce?K?D89=)H5N=xkLLbPKx*l2DZMXjNbNJC35<S3eZEW5pbNh}[8r2C2&OjePO!sq$dDc<qA?jq^*XG/p5<T<}5Pa6Vvh<mm?m>iAZA-@w*XG/0?m}.r*1V6.L5m$987abpZQC/fbN}Jy}qH9bZMtfNbKJq(&Y:!K5f{OUO!sjTdD8Joa!Fl*bKJq(>eMf?kdxv:bJKO7&YXQH4>Yi+vBBfX&]jwd:s!qtbK4o#3UpNB5zO3jvfzhfRvX[h.c)L%LQmfP8r+HS.q2!EbJYG1W:OV0lp+MHO=cwH(er(O.jzu)O/cH6znz!.ABBdRK.K3JI6jbd[7X/H"
"P8i9O8sH9y5yELbbKJq(&Y=ck-V%(!K}[R:9L(gK.6wj6*.N]39L]f![7XRe*-sj)T?&^k.q36EbHO.1?nhZ*jTnLhbIzy1L0dOHcl}7QL43OCYl2pCj?0:vvnjP-?nm}>jZ*:pbKaIs^MO]pI^X}nL4CyM?niO{.xWRCvoh8<s)OY072gr.voI]VH85F<*Hja$bI0PsoorQO@m&TD!7})jloPgN(Tr{AP0yI*n:vRXj!xv7nv!5TVos!/6?DEkP7.6k@LB:&]kbI+PulfHeJ)Jej>}TabI*k(cFve#jXio6P0TBhckyyh.xW.X*:R&Uo#psyChaW^XC:n??nm}>jTnxrP0yH7!!yVXj>}rvbHO.(L0bzB/^25ZL43Pfdjcp/j?0:vP1x7?!!D)-j>}TabI0LN!!D)-jZ*Vunvxk7W2nuZ6M6cI*:pY.loQGz6Rq&Zvmy%umIKSdByvH0*:i-bY#1hAjTnUCP0yIJWIhoy(!A+zP0yH@ND{DM/}bR^L4<mqdjgT(j?12GP1x9oWImHC(!B8eP0/sT?nm}>j>}kAnwh)iW2rY^6M6zT*:p.fb0z9a.K!p+*:i-Ao#tWoBLF66*+3ymY#5LGfcY%umY2L$[?ly4b{J53v<UqQmSh]LfetcabjWxi&Xf9]IZ4<KP+In^}S8iSv{Y4#bj^36{hO)<7IuhxbkFX(@od=WfwMW.NeLw!IpL7L7K#hWizE>k10]e-yOq=[bkO*Ws}z[H7h56x//6Ciw{8]ljkyv)bj#}AWQ<5(Q1{3LT+H)wnOat@^3UU7bkO$d1]&l%KGjJ$?3KjIdV[p/xwARB?3QV3(:9p}Rd!Yi?33j51XmwT^4vC!O2IQ<uc{)IIROY[?3d&O3ZEy{+?S-Pbgn4!/I:0xeY#G[bgO}T0XG+8uN1sOw+4z%/I:0BwTsO)*@^bJ^K&R&eVT(PbgZbSa$YpvWOJEx?6C6X0XEvMeMq)suyNV}D.$bC(x^U=w=p43!)xSS^Q?#3MQUwKeZ!lyS+mSz?07.IK$v)gq!w]VMM!1[33s#O<cMs+x9>A7BO.Fa]4++ax^}h{S{ua:^<0IgSzgOX*piUU^Yr69uvey<^K&R$+B4VwxG:/0INT>^+*8^5?3djg?o$a9iO+(hw=p5Cct86V[%==fNr(W8a$.Aq^JT&6q/TrhNtAbMhZjzlN35[+4B#.^Uh9Lrwdv9!n?:M%kQ8=Ixjrr3r-S^jru^(*NVY9!}wrKZrJp6*bauQezpU[8NXGS#YB0]orP$Y}SHxbMNYUe@oUI-m![%$bwdvaQw}#+0Ho1zVFgYun4C6uHOMCJvxngWQMuI83![cyZSfQGuzopHq5ujHpNy2g+DkLvw![c*urGkMNOZ%(@hZjzlN35)xNDE0m6wUAFdKb3t@MOn?/4>4%?eg-x.yo7$kQ8V.wcah]Nq:vv/4K53?3}}iHgP6udvG2DbsiUQ&^J7rhDzNUN5jyrp0WOlbXGb+b0)MzVlXnq{98$x*^io!p1{<-N<7g{bw(f>(mFZ<c$3SYnyw+cmdn)P*F>(vm?%h5XvC{*n.{*yNh$3)qE(*$nG{SwwvvR<w]Kt@C4RNob1u:MAWkH<-hdcTb2c*P.xZr%h@g++MtLyjXSU^SA*X-JuB?E?tKOrrQK.5QP{gLk[*6:pG{>-<SmEgSo1rY*c+B$D^P@:@.C5W=-uF3VtP5:Xr{nwM*lZ!A?l>)hzqi)PyU=K(z%=p6m&Zu@pp2x6*tgJtXT<#?a%FjPa%+]E&V:Z8wfJRYxU>@qdr)3-RKtpBVNOA*e-29yb4-B.a%Fk=L4is@oEAo=a%XOh9O<5Rb5b1??pVaml%=UIT?9S9tpox0hpr.Pb4.Zxa%{&k@Oq6W/4o4{bI?f<ucTSvX2kY9B!%7MCu#EyHM]HdsgIF4mLHffF1oeO?pYGfg9Qm=TD1sSQ9?0adKhwmbbomK?p:<A?qb7^a%FIjQIq9!Ro&Ia-V?%:/knllo[RoIbINlG?qc^W"
"VlFbTe*L(3a%*EH/(?4*?dATs?qc[!j(6%b.qi&W?qc[J5/&dKz6vsFa<#h<!?75b9I@DYa<}xE&Vl(8f^k8uOGY8BdQKvELhmx6z0*j9?qW25a6H#KmXeKxl$)*#tsL5=?xV/=2MN-c9u<1G?xA$Q[/uNlZUQN4yXW(L/8u&A0Dmw+KllBX9Pz20&./CTT{9*d1S%C&v4(rjB{</-*e5%xw2fT]XYN-rj(SpwUozR3?x/w4AVH$-seuw=RfF>mo5@&ts/$x.IFzi{N^maJ&/):cyI&7/dm>Ys&V!HB?yW.i4DpQ0zR=[atFo2t/1fY^0dUKrm^ry#dNZ//&VWd7?yW7ge-FiF4P8nua^FcAqw*Tv8gRh{sv-J@/w[9r8nQtha^9Y^M&AKge):/%sPPEylk=W8G@jzU?Emveg8UQ3@*L>DucRIrg8UjC0E%{K?4/hhvI$(YVTe!Nf[Ht2zrvLe8cx[2?EmvL@Pkn]>70Y#a!v:$nct:+VM?HJsxYX2xxm]E88tHWa^d$>rslC8<a3-P?ER!+EDzYtuALazswv}Q3a3txjh6.}FbcCa:@GnK09=O/?EUZ@mgMHw>h6?(N(yzxiKr4]:Iz)m?F6r&=bq4H]@?V-uoNBap3gqB<>:>p?E([+oC=M5>h6?(?EN9(F!$QV2{Bs9F6]c-2Mt({>ZROdp}DlwDH[r}7A[j<a:r$yF/0p6>PLfJa-.5GF!]:T7A)*d?IQ]3+jBFz3@BCip%eAW<QuUn]h9puUqD.+Z35ql08QTnp{H)qDH@w6s0qX{rGlis!MyW/>ZTPzUq=Mk4D@rS3:hI>a:@N8a@qzEoX$vqGW{cmH@C#@2BfPmB9Tx:.wO#.7Hph[?I6HOT(#^x]mdgCUqc]]IU12C>=XLa?IR57Wb]RW]g&M+a-XyQHYnJC>?Z[c?I{-Ce:dYUV!*oIDP*O@n#QQ2462rJ?I[8.HYDb:>PLwB?IhjSAU)3H7A[80p{H(XClx7q>ZSdtp}uuADI1yc7A[g[a:r$yF/0p6>PK#)eoo6AH6PH%!>.LaT%-TkM*R14}a)NjpNQ#3BQGvr{!XV/Use&5p3tAZ2R!jyf/xqBF!]OO1AWZDBxWW1M*SWx)TBl)T%-Uud=L{u)+IKz?HVL[[!Zer%k5+0Fcur7Dh5L93.JjxpNZosM*R1a&^QgLE^k{}e.m6O&&T(:T$65!DI36F3:h(kUs!B-Z38T^7A[g[pN7Xhxyplo!#*HZ?eg<3T>wi{7V7ena:@NTmgm$77K#zE?Iqf*G:U7U>PLxr?Iqf*JP.Ka7A[26a-XyQF[c3E3SabBUshOqAU<lc7M?&na:3}>lkH$y>PLxr?HVN0d2/BrJOCZbXbgUC3ICTd34J:1p)zGu6:V4L{e41Za:1eoBQI3V>V9(QrH@tIZ37OV(xbxozMQUc>u%11sB8}ca^-g8!Ev2[{e6xua^-pW<QFza@[-kyFasAjHYnJr7K#$UUs+W}p3wP*4i.hha^AIuc&)e#[>^%tUpxl@8Uey%{ob$ra:S!qG:U7J7Q3E6p}=7@p3yC[N]4SA?E(}u9P<Lu{tfl}a=UYn[!Zex>PLxrUse&QF!@}n>?.9}?H#uvp3r:G>?.9}?I{.xEDR&00n-qzp}DoxEDYKt7Q3E6a-XxG[!SHk7A)]xp{H[Q0V9285*qoi?IqaOCWkTj7V6?(oU4=e!xH?0%4{0qp{*WY3Imfa]XLofp}=6?VkoYk3/l^kp{?t/BQJ]3RfVlHa:p3x2Gq6h>PLfza:r%o]:sQ9@#(z0FaThMnb]++]NDSo?Iqc/JP+(]7V6@F?H#uuS}mDX]SH^+UsFT&Y7w?d%f1<0a:S=7BQLxt7Q3E6Fc$}Lmgm$f]XLF$Us!DzIUb+S7V6@Fa:@L$AU)O^3:ip*a:@OhHYz@6>=Wz.?I{-H4D[%W0n-Ii?I{:5o7Z#1>?.r.p{*X/F!]OJ7K#M6a:@L$AU?MW3N6&*p{H(XCMn[i0dUa!?IqdPBQLxt"
"7Q3E6a-XyQF!@}D>PLxr?HVLuVkluG]Dwt+a:1e?XbYvk3/lNBp{*WY1Q?I808P%qUse?H3Imfd>=Wz.p}=6?Y7B<k3N6KbFdpY=CMpNF3XePRp}cGGCMrm6>ZS{sp}DpH3IpI)03M^OFdp-0HYynS7Q3E6a:S=7BQLxz08Qo0?H#uvmgjRM]SH^=Ur<5eo7WSv]SI0MUs!B-Z38T(]XLF$a-XyQHYwPD>?Z[c?H#tlZ38T*3/lEKFdpZ)2MS.v>ZTcap}DqQIUa8D>PLxr?Iqf*JP/kU>?Z[ca-XxH0V92d>PK%.p{*VODH$4t7A[26a:s0q0VfZH03ME)Ur<5eo7ToE7V6-WFdpY=EDYKt7F}yNp}DoxBQGvn3Satkp{*X/G:Q.+3Xex*UsFT&VkoYk3XePRp}=7@mgh@9>UO>.?H#sXAU<li03MN!?HVL]1Q!ep0iX*0Fc3Kq2MS.q0dUjXUsFU%o7WSn7V6@FFdpZ)4D#V808P%qUs!BgDH@w63:i83p}DqQG:U7M3XePRUsFU%lkM*$]NEqvUtak:o7U@[>=WRJ?IQ}GAU)O*03M^O?HVN2G:Q.<>=Wz.?HVL]4D#V80n-zqUtak:o7Z$(7V7enFdp-0JP/kU>PK%.a:r#gAU)O:7A)]ep{H[vmggnV7K#M6FcVb#EDVg.7F}Zna:S^:lkH![>PLfJ?HVN2JP.Kd3N6T2Ur<44VkluB3Xex*p}=6?Wf]*v08Qo0Fcur7DH$4z0iX*0?IQ@4Wf$&B3:ip*p{H(XCMrl}7K#$Ua:r#-Wf}G(]IA3$a:s0q4D@rF08Qo0Fc$}Lnb%/Y>ZS{s?Iqf*JP+(@3N6&*Ur<5elkJFf]NE8MUtajSY7w?e0iY2!Fc$](DI1yi0iY2!?IQ%emgm$f]SI0MUr<3GEDVg.7V7ena:@MJWf$&G]IAl+a:s1zHYz@m>?.9}?H#uvp3w=M]XLfo?I{-H4D#Vd>?.r.a-XwxCMn[c7V6-Wa:1e?VkoYk3N6KbFc3Kq2MM2@>ZSZJp}DoxDI1yc7A[j<a:s1zF/0p6>PL6Ra:1f@lkJF77F}HEFcurSZ35q23Satkp}DpH1Q=G1]IA3$Use&QIU6-M3/lNBUtajSVkjWe0dUa!FdpY=CMpNF3XePRp}cGGCMrm10dUsP?IqeZ1Q?I80iY2!Fc3LzHYz@m>ZTca?IqedVkoYp]DwLMa:1f@nb#F$]XLofUse?H3IpI)0iXZ9Us!Cq3IpI@>=WRJp}=7@o7Z$$]Dwk>?IqdPEDYKt7V6&Np}=7@mgm$77K#$UFcVeiJP+)4>UO>.?IQ%ep3r:G>?.9}?I{:QJP/kL]=OweE)Kw.ATc<a2oX3(E$v[G0TZ5@DqdoQb7nU$BKK6r2yAx7o?Iyjb{Da^2oXkJ?yvw8BNpi0xcqr*SD5k]!Bum<%4QunCRPJYF^pg/<wi.rEcp9HViL>E!VdnO?*#9qQ6EKl7KGY>be&@KZ2*CT2Ek&qaBrmJc=w?F2DRkaaA@ZFb)!(^6mpZ.?*#9<JRF0$2oXtbU5T={4FX?m08e8k?dgJexwVwB@I}jkreHFeG{2OQ{[>g*s&oT:IH19m?*fd{?u2fRa{SFB}})zd??{w//uQ^Z9dS*MEc{Zhyhe./0]F)Hm[]xdJRvX(c/oZ:VuK0ZOdzI({]wobPn&Run{E21{$}8azxnJfJRH&M[GpZk?>yn%quY^?1:mp7ayj&jQ6k..!f4S)Vt{caT*ocM]wunCQ$tR?GYKsO3p}VcaJiB1Q9XTS/eitlU+YmlrnQA!1T3!aCe%UUvG>U+}aUSwqpJ9Lb{dr&1YBs-Ff^0$EFM@X}ROdU?(iU:)Fm@a4=sbnB^!Cy2x1]s*<TMYaxN(Q-vuw=1*j^7n]sAPL&Epc3p[4#FDu0N2P/}D-!5TnY&/Nmbgi5n5}HySCYC*C2AdGt/?TV*T$*fi+iabr]rat2qRIUT8WATx1(Dv1aw@xHX}S#CjM0pq=mji)j>[>*1(Aqu"
"?(VYIQ72}N1*Q(0PaRs+Q6jt{}GS/E?(V-.)FoA(4=sz+P]F[[*sUq.7#IE(atK3wh%?.y1k*@$auGf/F=!/4Y=YolB$%dxh*lxa2CuUJp]UJUj]AIF17AqGauy[c.unYop9qqf/od1EF=!/2154es?]Vj2HZx+j0E[T^p-xsQG-L0^zbp&EQ16HNquG:$@7qUI?{878H5GO5@HXJov1FHr-VVDW1aF#6P/qWw6TJo#19MT6Fc2?62MP3c0n>?bFc#oGi0&z16m=UKaug1>=hsc{2&A+#PE<-Fa{l-BcS3JtUtLD{IUk(H9J@VCCg#}CL^vNP0D<q0E}/U2&SC<i2J9FkUZlP$i%TRVDer0jBNw#Gg6uBP5*fvn?[/GM-r]j%@igpR?]C!vj]BiTCI*vlPPe(Ms@(Iq1i-v{u5!^t:>>}1@bNw/amh0e&R]uT@VykcrL!W&74rs7@^KzpDt3Oa(mYQo@=hrCamhbE5=LwA)&QJ<?#S(o3j[rI0W5V*Vks&Ls$Z#ULon(.n$(l(G.F7j]?SYC&1Vk-F=nGPzld>cCnhJ8BQWgyvbI7=uJ4Z#74yb%3swk.POrU%0SExMbTuOfamoeU<ErvU]Din@F2tqS2MRh40!pAzpV})5&R[w[aIU#ADx52Xj]]2o{tu?=PPc(i(jq!.0hfM1aBo&bysYFMx&wi=nUb)+6QTzy@HsG0UI]EM*w-nX1gqq*&2u^zh}#9P391EjaqaXi^Gns[@Z*Vd?#iQZ@S3$U@[[!Xamof8b-*>1QiS)xalTl[9MjSs0H}*v&1)V%9S:+i0PjN=acDvvK[&t4NnWUEu2[7fK[&t1}jPLsad7}dL%MNP}jP$Rvg&qSM)CXYq?d:vac*8YVi3%q58.goacD&iR#yPxPc6sOf312?/v24j}kttkacDxSw:yIg1%ob%&a/2eM>+cI2^60+mG&GdO^5n<5q8EQ&be-vTy^2e{Ke#6pZrhhDD>w=242ZM&a^*u3C@/:-]7.sW/(Pbp5Mf0(}KIkq}HMz/tZB<6%{1yad-XhsYcuk1%I$3D7Ckzc}WHWTe>C(H]!0{EG>d922WR?sTPk^rK9.0s6b}e&aI+:iUQBW}ksh@GGYO*)GNNLQ#l1OYr@sgyiU5x2h<k7&9(CP(R1T$5MG8y&a.!BZ6lZR1%o9g*sMtZzu/1=:dV/zP/@XsxMthC{#2PQab-lP5+]K].5[WxC.bTBh&nYgU8rtQV*7yfJT9$p{!Cc2aa^#BK{Ad/C}5<&uEDf8yD%ILqq=GZ&cFl]6@-pd.g%>-oE&nJL!%=y{!ipBoA-lPJT5.<.l-mkabC5zHL$8kXtqnItJZ)0rL.uzT6w}A&c1GgF}mrOB03f}&cg35+rnJnJnOZJC*%N/C.3@{j+$Od&cp34+rH/KrejOkoGmiK:bj6gJur>eI0+IP/TsB){(9@nVPyY1JT8}+{>qvtDD(=6g$$/:tc1G3WW).NN%6uW{{Vm*ab+N9EH0yA6T5Exaa:i#2Ku@0Qi$F.&b]?Vzu)zKy2BBBpTfg!G=Ye]+ZE+wU^^9VEE17O2oSm^&cm8+zu&(62Jf}3&c>]n0R!nYZ/]Z(HHG%Y.t2Sq{-C7ZV?05%G%i6MSnIW<ti9G9[+shWU1POktey0(BIS-O.5[FSaa*y-xO28Y{/+@Fn}Kv86T/Grj%bP6oDchrG@=wLCc*-l&bOKECzQNiCoq-fR6Sp1!x0A}Dq1c1R7j+LI:6(cDlm]XWpe<p]UzSf2oO.U&bOK*e^C5Tr^a3doDcivRlb0fDwwzPtg9bK}CINf2y8RWR5jA3T3on%ky}U@elT:=B/}L>lG0gEHw?B5wVBWmfF3A$C>=ZoJR:<$2JfheHF^&M>YGMe2JfZtHGJu&yD/ao:fYATCY0=ACPBcY0]HWw&bOK#{?^dU2oO%hC+#hkui=nqKf&Dk&cFqVBYEvDiiA6qoDBBj9Tw&=2JfF6T=61kEG@KE{-C7ZoAB]Z2BV-={[/!v!t:@t"
"HL.(.{<1ZiabToEF+R/+qY*#+n*i#6?cq4VtowB>?<%%2=iFIXT970soa1K)6U>/CIa!ehC:c1v]Usc9I64z^aa:hX1I^XlZ!t8KoB.pHGU[)(2tD/QV?pGc@S=S]nk2>9orym1i-m5PH-.HLC-2vyMX/%QqfXtvC:BniL!.Ka{!g{HoAB{lT{r)fCc?6gHI^uLxM/eG2zWdxabSGeJt!SfIfLE^R4*.w6U{K-q6gNHR6kQQIEQ)[{[/)-oBtQx@S/!e2EBp(V?OZI{&:1vku^EWV&ghk(PVFGqphhJCS!8VARO6B9uQnCe)X!vEG}Gi2zWg4&cer%4Hj8e^+]Eptj8SJxM&r.2zWg4V?(#BB/@7I2Jg4o&bOJlF+TFM<Sxe>th/mM6%nTf{>JqMC:JT@Jt^NX{!h5WabS}4[+sh^2oO.T&bOLIEG{%MkoCHs&b(:!rjjovWae}-of@tR<Il29U1O#lHH!h%7<nNN{{oRdC:*]>1I&@qZ[&7SoBb!-G%6uN{/=#FoA?PSCzPIZ{:#/8^F:w}m8i0?I64K%C:BnxxnGp&kx@R]V?Gp#JT6/}IqSPmth/miyD&[%C86RzHI^t<2ZiMGIhb%NR5jv$G%4l^Z>VV4ab-vW7M2G{qa}D8R5Aj?DqMp({{o:sab-wwARM)kU1Pi/V?h5t{&Z[KkJ3:xV?($-8:Ane8NF#eHIW@OO^iVsWy05nac3P:B/(*D2EA%}&bW{M4Hm=1Z#fUDoB.qh>zwbzkx@+8HIxYAEG#X<2EBvE&bOugg4PP^]b//@of@fABHibA]lj76Cun]eg{<2$]wE!+a3i0J.sH(>3)ebbHTC+uVgGbA1Zs<%HPrGUK[v0P}RYCDo2#Lv-icE.nj$0f&kz#}lgW-whY<(-a2Ht$m%dD)V}A]iof@fABHnXXAAs1WGeETVd!D^-0V%qeapXGXEHJOxzl$t)CldpSh*6sw]ljfvogKaDh*6sA.R&FKKtR3ghBr3QgRu*CGe<4qnhAd@@!5@*a4QXn-h[+E4b/nQa3Hac<G{0&hY<(-QJ?EI+7Bg1]SH.Jr#QZ)i5v:9wp@RnHR-R&e^#{Tl:QyeW2LLpj4(%TqNXyMa6!XNOW)F/9:PrCtt+7KN[1L:42t4?&lzQ3p7jGCAnXZLD[C&=:%Vb-]p>TYH!rrT.Fw8P45egV&lz+?=j2LX[GeH@o7zCVg{!tp[M=:%&obn96T9t][PKj]9$!4I[:PlZ4J=-4&nKj:&QF0uCCSZVQrTk<nv0Kh{:@iAadeo=1G@)*4OD@dWdHNA7)4w]!/LQh&n?jqj5WjnWS*NroFc*:<Lb#TW]o-fFaV5<h9>gIlv&cnW9xd^.snPZ=(Ti6Cj-2%>w]It}GE+7ad5mUR$LQ{QeR%ey8T*o]w+t}0o0nPH}uD&>-nvq[RD$6&oljz6T9u64Tz:1a07]!d#tt^yNl[nT25!p(B}0X%eUDZI0E=V0R9LN4Z+K%Wh]*+<<SeNpd>.[WdiVd8=4K(zNxh4WeiUS<K}2FTL+Yr&mK1D<GJ?b4OFLD&o?6%9Ua{X]X!Z0o8YrBQ%v]:pH7QzIaAV+&QJVg4Z^lJ&n@E@9Ua344=D:+&o?6%5+imf[K&=<9#ii@[:RCLzxxZT9$!4q?vpRDp!H&V&o*Q/7JcRI[H(S*9#7@o5+c}T^JoIu9$!4L3s+e04Tp%Rt{n-g{IO{IUzA%%o7Yr{1H2irzC6D?I4)cn9U7*M4ODw$We]!>@TMw9k%uOXr2hcS1WOBIZfiQ6&n@xe5+f%b4U+^q9#!41<GEI8zNxRJCiC5g7KUf1[QH?pr2tQ0c&!=4PSvnRMg$ul)HFgsJ!]v#I5Ps0{(Ohs[RE>Fo8paY9Ua$NruAcQ.!nFt>aEkiZcj3XTncd5d<K0kKf[7?9#!0Yj}s%%p#HUo9#w?51G$C{p-*qtWfR]4}+=qqpt/Qp&jYAN>X]ziO}Xo@Qs5I?<<QQE4J=y*&nKj2?vrfLzS6PVa07]n<GEI9"
"K(kTqWehDk(T0:E4=DsC&o?9k>Z<%3((K1R9)mXUK)DJX7AQAVuA@>Qh*#*^Q(3[[9)m+e:6u3S(}g1^9)qo%F:x8A6ctmmE6Ja#e=cRh6ct0.q/Rd/EIkyyQo*MWn?]X>g]9-pts}iFIzRbK=j-/B6s9><KCbI!e5lRL]/E<}r3(vg6?juoXM0vr&wD38.r+F#)&GGv&w6XuGRWM=tPR6)B<4J$>C)G2Ymcd]F9lYv<WBWsyo9Z2j17bV>-:%+osm#mH0yv<ig7=9RIRaB9(%?I[:8Lh)bo#t&wHUm>@qTmVo9xBTep8>i>[%b6r>jB9)kdL+wapzpM1DN9(moLg[)QyrjinIWpeABoiMEn)6TLvWqffIg40[jRqRsZWhRc:)G[T&5f[[2&uUjxuk2]D)bpbwmiia+m5w&%(}kLpTFO>Ni2W]26ctpp&wC8SO!H6D6x7-G9!AK6AP!E1kyfP#//ceSGPy#a>=OVk9!#RXF-[ekuML6%&C[8pp8s$fcqjE79!&I[JP}S-[KUtHah%NxvM&W@&xC#89!Xg+BOefy(uoL99/dr#h+@sNQE&7ho7zc7BHu+ZYIJ:R9!+66*?*1b(Zz%39!AJQ]!yxYA]zD-u.e6=w%BA>FOKOY9!}DKI^HS0QV3AaT7h6QteuU](WP%?F5K^AH.+AU?H.8R9/adg=khuF7Q.[*pe)EqK#jqqjkHZ6L(Z92Ie18(7OEuo9/t(R7y)Gp1dRA.r4I%?dO}gA0]XOAW2nL8yF{5T7QWM=I^6B-s8kkX-7K(P9!:l@-wnM[MMz4y.a>!DjSCTD>?.kp&Eg]Xg4%NNW&QLCo%?:}DtVIN7KLty&C@^Z4I#O>PY72iRD$6$N{s{37Q.K8bMnIrB}iFW7AMN=9-DC/qv/VG<SmpXLf@R/zw(aMpAAPSzJ^EftdyTb[K-}SW%Hq@9VA#b)+Jsg9-zkzqrqn(<SmSR9:rXYrg1.e!Mg?H9-q2<H$9LFvyy>Zrw}%3&Pi$kK-^/Bo]d4<TIo)G<Y4Yb9:1C:ggx3m2{QmqY42P}R2d4LtP:z]EkMen.r1pX8D8}JNOdikrG:fJUPuW6L.Y%6/yiOVS2/m-&IzIPsk6=VVl0fU&J085t5Jl)Qph(w?.BC<5!r(]<.5dXIdzSzoeoh[0cB!Ja6]<SRio/nV3Z+0oi@K}Q2*{f8N.STWCa#@O-Tn78MbhaqPgX#jz#oB(A4?-&IT1Q+C]/1VeIR3&I{@%X}CUMU1$ScRQ%)W3#p.g8N&ngWK@GST61tPak-T=S!B3Z3xSdG<ZaI%&J3CDSm!UnzCZrJ&IFxrJVc(V8Sv-Detq8Bh9aAd<5[a%osLDMg2<rmD(+46&LBN8qr5YR?GKLh9X(XS&P0/rOl6Z<Sf>66.q=N@HGZJ/COE}D.q<8&&#WC:9YkM[HoI.<9hZ13o&o@7(4TVS<9G@u9Y(+pTdDWpVpP5@u>ad$jei&RxHiTAE?-/#lm)Y<-yL1er]a?P-XOV1982Hd&IsJMAPs#qd]o]{XM6B7-RCX?ZtV]*LtiaeAPu2${o1Nq&Mfwq<!><zDRM8J&LBNk:AMQWaHzI!JwoF3{Tn3xN$A1oW/dl/jxLYOd4nXvEl1@(IL=]3+ypEYo$IMe<T^F%9jE9Q&M-Amx3VilmS}=/TS85+ws8n(g*8i/WYp+l+whl=k&G)ys[R^MC!h>LO5TUgFevn8XpHPY9j8^/C(9imo4H7v}@7&$W@9gcJVv31%9TRP&MU#[=9R!*9o250Y3qqO+534M&Voo5B8*gXb-Cpg&*{<y9WfJ:APh4-hY:F<F05.Y:(ao*}l@9rS$2Z.e!Lb+&(Fy:9WfK:VeO:Z9uKR?rpPzn?w[KE[2(pi&N!u9:oHRMJ7y!AnxG*ib-A=&J7y#Fnw{DT:2uM9&*{<yY>4:L2}QC15k-?^q33az!LCeK)O%faW(bqpj@[7}mxVvk9WfJ:Cq*M<h+p(CB8*hXCq&w9epjBs9W*+8d%3ii9z1lQ"
"nx@layvu10mD8[uJeI+^ilm5s9Gl10&N(j#R4w*dJ2<Cb9WDYo5-R:/U&chNS%XqZN^U8MmxVniUkmRC!/IpZBbV2#JeI+^j@}9]B7aFSW(XV&*Jd>b&<QV4F17HR1RwvQ&VoUe9XazM8N7nn5-IhX&Nq]c?w)W0pVo&19{$50wq]H4B7aFSnxq1]/Vo2i9PtU(&NP6Vj@}9]Bk+=)&OCQl?w$O5&7K:v9S+g!rBW6Stojc89U6nk2RC39&4u+xDnF5:5-sYk.aULi9S+hX]P4]6FWVlHsp%E+7olu)a4fj19Tzz$[.?L14O/5$S@.eq[.^CetyV>#vyi<%g2L-stbzTW9S+gv^CS*<ag6}:Vbb$cB[3H.a07/fsp%GI@gJ3vJ-W25p[gK*6=]COag7af&Q{19}:tYFvth#Gl%z<e.qDvumsO>&FcD%D6VG5a0dUFx&RG7vVgk9wxl(zJ9TwCY1vWRAihM@YFcZXiR(]1N0iN7RH08U31S</WWOdp>ynDpo[.&y0DMAjC&Rj]A1[yNA9-qcKT@4dEXn+fwvxFnrH&?Xn(sS^)qq}G@D-cMU3B5H*6x5murq)##{Py<9ag7wS9S+iT]5E85Y)=sGFcil-ywbRJkh(jnowq^ZM+8eI{tL^Kl%p8V@VwGjk15Z0PZA]3=l5)%al02a&T$2/<QS${?MwVD9PwetAO=1MNJp/ts*u?^h*uw7XQQsf&U+&s>Ly.I?w$((9Qt.$a(l99oc90F&T$44e*hn65AjqU9Pt<mng>}yRxa<znC7he-hH3:aXMb+YtRAUSV$oqewECmU/Seah{lLU]X{K$qmdLa6^aK-ohVlyCCMT*Wu1js}Yj^:UU<N(leC+NM/o5fDc$E?>CI9CWToOlZ^(G!(cT$i{tM[Z&U*{e{Y3/F0du3EFjG(ES0De0}gu^6D(a-#&OGoqaZ8r+9QU89i)nZcmP1fZSs!0hS[gcGNxsS!&U{.l980B7?AfDLpNNiDxud!(3:vh@U^=/xN&9m$RN2&%&U%&o7ZEyVoJ=c{&T$4-i1Wg$^&E<1?]9mv?wrIy0o66$9QU8T<2b*h8-#[1H9.-p5-cW8a(Tlo&Ur?+T%1?IsMgd%O{?<wDsV1N?]9kLDeo.79WfI:aVuQX9M<G+.q5bg]C#j4w8Vi*/A02f7k$06XescA:L:Mu*#pzD9Myq5g2fnWbosRY!m4O2gmoOWiAI@P9MolM!HEPEVt&:%grwzkg>8}NFK/bMJS:*<j%JOUxMUO(C4zdW^CgZ$7jOxfrFs=ec0]7Eurco:FaWcyLb.PVE&VIX9Nk)t)JS:(j@}NXwM%-BhY}3vFPptfPbf.tCFnB.*<DQ4C(=YImiboI@d1ViFbIY-{Y4liP)#{jrP<x5+yLVRE=ruBXdNIR4KjHc&2PU?KuK]GiLXU81e#K%rF)1/<VYV0%eE0hSVsBa9yu+AP@F}eSUDN7i&Tjw?):wHXeQ0V0O>g7FK/3S9Nk[t9WqYGFPpVlXe((m}*?=euEQD#SVsAQ{Y5.{(!4^G9NSZ57cS+:bE!uv&Yb<^)JZZFE=r9qw9$OK7A.S[EVs).&Yj*lg2fnWbxrQi&XVbzj%MNcbJpQC9L{Kv5.[zF&78+Y&XV98.q5L<*{!dx9Miz!0O(V{Z-ynDXdNIo}*/^aE=r9qOdl+:=jF@O(MbJMaEcnk&>o<F04V=%SUi190O*ikbosZKXd<xZ)JV-zEJuI/n92ce6N?7Fbs*]s&XVam1ZUi<uvSh=SU+.o>M4Q$uEQD#rQcpcc0<QhY&B79rPNMDj%KOc*-sy0XdNIo@V:}@a05KQARa.thY]4muvSFKFcVcP(y&.nP@F}fFdj?}2MDsb0dSS0Ff!A1Xb[d?@d3%QXesel9WuWC*^+Ggnmyqm)JS:]Z-yYBOdd!+(y>i{0ic>oSV4Nh)lW/=Y@A0ASUi1B6<DgwQ<L3rU1*j43Xt/bFW=MP&Yz.Jj%NMH*&mW#kv+hqhd0p6lvX+f!gdd("
")lY.=%kqq^&xv![^Cr!e2P.hnnm}$.iLVVXbos[w&YXO%@V!]&&nM-OXe()OjaVxHkJj<KXe()OhA@j#*{^%K9MWee&OrAHEJunWwawu:7ATWoEVsD:&Yj^{.q8Kj*&rjx9MiA^g2eoebot5A/(ydn0Pbn8urc0#FcVc9*m+IVY]AcGna4Px-d2igEN&/JARaY2:LYOHEZ*8v9Nk)1>o5XuEOa7vqyzJB<Z7%J09erTFcxnf3zyF!lAg=bSVsAa&OEHAKM5CGwby/$-#<p$ESs]kXeckS)JV-w*-sBZA-7RtL]rKp/dkFNFc}#I>M69X0dUkVSU+-o8LIV!0mSy%FcVd94KnF2&j6*MJT^oO:?O}xbB/2srQcnT>M4R3P!1&<!fOC38LIV*P!1&<SU+-o8LHf#lAgTM9M{le-c/4Xa05KQ9MGni9WuWMFW=UIXese-jaVxPbE!U=rQAaM(y(YQ0dUIBcW(#^4J3u}av>zC&XxjU)JV-zE=rhkASdhM0O(V#bv*Y$Xese-g2eoebArd.&YXP%9Wtg]bJp!o&YXOS&OkD}Z-ynD9MiyE-#^?g*#pjR9MizD[.yR>&2PU?Ob}A93XotzbE!uvXe(>n-u@8b*-vHO9StvEg5d:-umTUdFc9z8.q4Mza05pGn92a?&OnCBE=ruBn9+-o6N(5Lbs?eD&XV9Zzx?WpP!1V1FcxmV)lT?murcg?rQcod1ZS.kP&HpRSUF>1g2cpe*#pjRn8-ly[.uT$al2wkAR>sDhY)LKY){n)rP<B.6<=y(Y){v*FcVd94mkNJ09erTFc@0I7YY6}lvYC{9MWf#g>7C=*{!iVna4QZ0O!ZqZ#vK5OcrjH*mZ4EuAa?xrQcnT>M4R3P[0GzFc@0I9yu+AP@Gb1SU<Um+yOUc*{!qOn9+-=hY]Kw*(q+[n9+-o5.#wcFPpVlXd@sI)lVse0mSG(SVsBW?t(GJwGGXAwawvljaOA$bAq[Q&YXOS)JS:(kJj+QASdgM>?$C:bB/k*9M%2[hY]KJbs?eD&Yj*lj%KOc*-sN?9L{J>1B.xpbs*.F9Myq:+yKW6*^*O^9Nk[10O!ZsFGq=MXdNH^)JS:]Z#v979R5nkJS!orC2^Tge7:&L0>>PLP!1^N9MiA^g2bp+FK!^I&X}@H&OoB@EN&ZPn8-m.g>6Dzbs*<4e7E@R!M#nt09eb!rP<A.(y&.iuvSxRrQAbM6<^)Xlri119Myt6j%KOmFGq=M&XVa7[.yR<a4L-zARGI:/z[}kuvSxRrP<Ak/z{Bc0dUsP!f>rE1ZXhC0ic%hFcVe8c&==@*(qA<n9F)<jaSz1F-nVE&XxkU0O(V{Z]w3XSU+.o(y(YNuAa-DFc@0I7YZMGY@A0AFc@0I6</xtP=H#8Xeckc+yJgAEZ*0BAR>sDin?p=FW=UIJT^qz9WxVaFVoC>!f>rE4mj7WQ0$mG/IdZ+g4Kf7EJuYUOb}B8g2bp!bxrv8&YXO%@V:}@al2wkJT^oO:?V)}&j6*M&Yz.JhB3g+bE!^t&XxabK>Cl=/=uUh9G/2I.@.#e6fY<T9HljpM&HFF*2uPm9Hb]v0OK.q03diPPUYri3jE@DcznidprTN@N?JS2/=uUh9HD&^vCfTTclnuq9Hb]}g1/sr//&%kFcj}e=l+.=(}UK>prKqTh0=lI08yPnFctk#xAg/C)4zq-U@Sp*.p=U?)0]GZpr>i5L<ze<cqI.$&:?x@K>Cl<03dBb&:HI<1NFr=0i[uxFc1tgh#[Ej]mdq$uPn3X:n*.#6mZt*9HWy!2MQKs0od-4FctkPM&KxN8gR+v&+esx4KMYocGnXspqLR*i$&T76ciq0U@z.EyzkYq>cl=fuPO}ahs[7c/)<CzFcVdWj%>^/6hDWUU@-R(+m+[^8l(d3.j8a$y-fGm)bz^]U%6Itg1^Oscv=cI&+*7-O*EjI>hHd&k5IX/MgK0G0ad]p&:}-vwBa*RcB2Jg&:HJk=l!5G4b#dKPW1(V+O/zQ"
"6tZm!9IaRqxAm4fcGn[&&:O#fqpnC/^ND8{9vITya>R(N^Pc$y9u<GA5Z(+b)5cWnSG?UeIR2zmwGrxP9vyNc7SGu=@Ejz(S7/5*a>0]=Eh%OFXTJ<B{@<:8HGHm#&)Dhi5!b<cq:$izEU2}EQ0:@le-Lga9v3^cdNc3Ccveo+qE}hP8N@L]@.?wp9w6h0!mmH/w>1J9XS:{/b[P7!{-gR&XSLsJQ0=zgVvPIaHJU:3{&W+MiPwE-SxmL]:D{]meWlLoDmHD7!nmX9r!jb@&9>6[^B9PQ9a.hZq258!R1zVQeM26WpFGIXpaWFiv[7t/=Lm^OslMVT^Y/a:9vk=brAqjd4#Bfiz[W*nxQu3I1(^}TFX6u-ANH(-Un]kEAq?C*geIBYLIQ53Eoj%duflz<q=9uvFW^Uvzy]g:9gju!X9+&Pogg^$db@{V&[eHAN0vSvio5^D-=S0P*Y1MPe-J[]G]N:u{<XXgbJoR{9w68b*z<K=]=kUp&S@AsD&A3Vl#?J9&[=oD?yP-G=@YhM9sl3)0NPr==}T%z9r+qz&M$u>RYu3jMeQ8t]{lW84&>X$S&^#}3LmL>^6rhT9s:.y0NP<0UPji9MHSuy0NJ1:4=W5!&@f#PIWDPoxg8T?oC=9=Wbh(ipFH@qCi4O^[Z72vffxx}p6>l4BM=&l>ZOwlpkxcph0zh!QFv%BC)9Oh[Z75j=]MpB&}(hCy9.*j/S*[dz?^3a{y[RTflflf&}IFZnR}]$^3r!-&{>xBIX!Dl{!DP>FbR97CMif.^0#x3s:f]!1^<*Z>q#itV3gIA[Z5oZ590y>FMf>4XAp2P<d.:E&}!:Ve?TjTfxh&NL<oZ37+*qL590hUU(Si@IX9b.^6j+vV*eBjB=9CC>VQQ}E$L8j[Z6xaai[SGVhMz1*Afe>^6wumM}RU1O?FNK^bL2uCXNDx[Z5o^fxhJG&@gQS?y)):ue>RUs:&1A1*E1R^3s2Lm*[hk]X.?wfhzm%&}CE)@X5:4^bK]#U*a]a]qyQBa9MG]&}m:u&M-C8Oe*IkLABcg<f.RK=uu2K9p0/7&M-Ca[eQLxFc0@VHZ$G[f{wrRiM&0O>twI==q}859o/Gz5Zzx:f=YOx9G/3W}YD}zcGnE$&$Z&n+DD:9wow}Fp4L+j<LH%r?+A$&&L:Tl8Y.V#SEGP:sP)[.C{/>c>]x13puMvi<Nb4Nf&43s9ouJ&Y9}E]Yl[-HRwi9hunjmP(F6[2yO0ov7r9(ZN%aLrqISqRJXM1:=k=3pi*w/SG[Fz4%dUqGFap*m2MN%^0dUl99oWw-N&3V}=F]-LLy.2f{Q^LLij!-)I4hlW>@NOLmWd^B+SlGo8MG*nYmV0)sHvX7.<?.j5WldcU[}fn8XB7J)Y^SSpmA243Mcd3SEPsj&$QHT*AdE#aroX&VidI{8Z6nn=AW?YCkE=?8wQETa5h.O&:HK1zlz>[clh0&&$VqG)Lf]/=F]Zrp5V>c-o#Ubf(geJLzwU/6[.2nf=Yx/&#nic9X?&s5-TLjVjnPe)Li-pOz$UBLB:I?5ZuD%+]CcQ9mr2ivA/!zlqw/(iG.Ma!0S$r(M1Bu9Gqr!5ZqrF34IR>9mr2jxAg^PllaSh^dbv2!2xRb6M>Oz9n8x.yAkD#=eQk]9mAr@clN/AP%o%>9mr2a86hJ@=2>JJFap^l3LTPn6I-eiFjG/76Y{y(gmo#Yr]&y)x=Xnw=93vFV}-wkyz*j6>70RUFbRjDFZ3%0qH/2HUqCwY^AQX%>2ZpPpaDbX^AN*MgjM=9php.d7Ye.=cXbzw9mS}xa<EpA1AZWZnv@%6d&zyN8l(:OppfOP(M3#0}+wg/9mr2MO9mkZgcgVV9uBw6{G(>obE9OcFjbNXe8U-wRR}Xd<0e6]c1=f^Z$$Gzwl}3HobaLl+$+-C<0:@/ck7/p=eQq{Bs]uec3s.pglH=bp9V:#d$>yvNwv$jS?0!.cm[g2=eQ8iUq=ud](/vFgfvz#9mY=6d(36n?lt5R"
"<1r9fb<Tk7)xM7m9m%]AbBL/igoO1IU%T@te&4J09PF+cN/k78)Lr/NYw>6U9k#Fe.oya?+PaEy9k{&T.ovy%6DCf[r@{J}vAY.m7iCbKVd2G*=neLr+&p*/9k}84FY]aQ6DCyOKWmjKg0ySVVeWN-V0ccc1M()yQ*(OAoy>Tkg@mBygFPjPKWOgG:2)9E9ur2S9lNRye&bKR>cpg{pg%HcwAFy#O>-}@Br{0zwOK#5+Sgvv<2V9[N??bO(eKY[pl:#.-D4jhgJjVM&!N3.*#A)D+-LA><1X=c+nAt6@7Z(>F5HGLxA5%<Adwh]9lnbg1-O$$0fG!L<27!ai#Lhm1f[gVFhDyG>MS3y+Z:N!KW{06eAq/A)=6a-yAtppyzE]k797v=V1i9YobnI&&iveQiG1=WM1)N>19u4hF5v>.M&<Hb8bs:^Vc*=Roc94N)6k719lPw.h:1w19G-#UgoJBu47UyYtB-${&}<Q%j#b)k&]q{(<32ms.oy{]+-LTR<2E?li#f&pgAswY<1X^2?zm[KQ#Kq@BDeebK<8&B09oT^Fn2.Ti#@y/gVHjl<32n1unHR7-2Maa9i)XaQ03hhkod5U9j0i%VccA-<-43S9jn4bqoL%e+sVi^FaR&ZqoKEXBSv@d9kjn&lcALU0f6$tr[f/(R#H.f+n^!gQlK+US#u%keq}a(9j0f@Zbc?m+z/G.V86Y9Wb[pV<Xw}d<47dm>?]ivy8<@AAg$DItXX=(</ea}u2!6Fn#rPpmT$dt<3-1?mciC$ktC80JAWsng*bPr8GqIlFbRnAXb+LW0dT<q<4cmvowu.I!F$<=gcMJmlV<x.g-T5xUC6dyQ03yi8xui[FcNLz](q4Y3GB(rJA21&1WC*!8Jv?#9k32YqoM$W8H}QQ9kjn&nb}Qt3v?.vh(^-mX$w1N1=PF@<4t}[@Pkl92#jQMFd0o{}X[G&g-Ta@y[0<kX$sGQ<X#dn<4Kh.unNTC8CTP}V8X!x]Yta3Y%?8z-&y2.Tze-W+I#tOFc}{04M7f:g{Kb%FbRozT#7UqyedTH<46fPT#9z^g]Z2*V8DjkpbG$>8NlwU<4BI/?zv#jyu4?:<4.CoT#5)Wg#?s/9hgI][Yz&KQ8Y#(rn>${C&oTv1MJw4Fj=KQS#WgZ+hVLl9iIf0a<e?Kh5@sj<59@]E03BcaVdP?Ve[v[Xcl<0=-#nw9p.Jq.o?jvBTA{R9hp==lcq}M7&YcbFdcl?c&/l>>I]/?pjt6X6ZaHRo?]YJroiiZSnQG[>C//<9isoTDL&+d>DQ48pitWVro^kl4nMFf9hJQT&Mtda1xr{h<59%gIN%#&GICy!F4bTwsnNIuQjrQ(-(Fg=:K+}/<=d)]ponz+pbP411MJ/CFkd$v{YKLxpK>:Nh](awS?L2f=t-KY9r4*MtF9J4UzJJXJO36cSug7$1*Fb]9if{Y+a=e09ogluFbBz[HYrwN98m<(VcCj^S#Qe8zl.S[ixFl3Rad-q+8I(h9ir$Y!z/H)9zO8UVo6aD(M+B&jgu7:k?>Os&MyfI*bnpPp9&x&snR2:hgL]n<59%r@XUW}Qt]a[<6BPrdb0SW&(on3th^S:snP0IhiXY}<69vX9YpiYB(MWJq)51rr$Bs$h8Ytf<65hqunZ+wEJp3#9fJY+g07?F2]]h:cEG}f>t-H*BXSpsFdY)Z{Yk(dE$ecM9gBiB[Yv$$=R%cK9&NR{5.O^ihxn4*<6*6fj#KyB{Q9KPyZTowun**A:&xz%V3tA>d&fg{dC8jH<7dbdT#o2{xTa?>SYIm%>?rr92fbysFmRoR)LZ6p9]?6NVf<KgyAuml5WuW7SYH{1)L.sy9.G-%pgrg-n+0^(hGkG/FKFNs7Zp$KhH]kTJTb.Zm?Hn+1far:<2mN2un/?L15qPHFlh{lxA^BhQKW1Th#0EmhXa>X{O{znF1wNUxzCJ8@&$bhF2:Hc{XIex:(vca9g6y7BXLSE:=kxFF9dzs<M8r)<S6KVppHk5"
"h0vX2:XF$xrv3!*xYI{GhPG}6<8c+IGYHe/(WOnvpo?]X6YDrOQFq:B<7dbNG&v2E7DqsK<70MUK&.OP^st<w<8c=Pg17O!QVlkISXIQnOoM3r>Y4iY9f<CxL&$6@::.JK9gF0GFYUgX*8MO6<6*9<xaVb^hJzj/pkU!^T#trG^su9)9tcANzz{J-:ScM6AXUVQM0be(&8MPh9fVK9K&.Pe^0$Pcp65DRN?^0sI*@1NAIRZ3Rcny16SLP-Fc@X7K&.R7<4.tJ9gF3f-n:PGGb>i%9f/ywvAzps@t.g5pd6?jtn?Jp:(v739f/w3N&^egUUE4n9g/2.MM0b[<-eur9gl8$vAAKNhH>?IFgzyP:nB@((}WrNJ(RpbX#-f3<(o[xmqNU$L?^sJhNdTTFmNRVc<Lmr:*3pMAQaQH:N1Aa<YekSpdT.D0M<wg}*X&N9gF3IT#e$w((qGpi4^9eNY.dB:VfRiF5vbFL&:$H:(u@oV0>:OR0&/-t22L=J)oHSO&i>nsxhyr<7x7^vADNThCP7{p7xg.nc:Vehyy/jiw9ju3BH>BhPNoCV5Q.^5Y)2g@nODIVcab!Wb[M?tcW&7:kLdB}LAjX:/r$&<8c.59Y3*0hNdVgrvs)dxSDQ/72%f&<7deVxZB+{7MMVXFgf=?EM3k0jMynxTp+)2zh-%XF8:-H^mP>:DBq^D060e2<7+HLzz[K69i)^EV5jO4ztAjfhSkm&9fU%]K&.Pe:ScH@FkKA*L&>^BoYGN6J?>?TL$!KD}8x#+9g4x1wz{Iu27T4eFu^joncR<!:+!fWg5w)pa<3#^^L>F<p6/.prn%]Z9cYlvVj#KyOTD:l:(9@)Az4w[On+:(+8LZdp98r#{Xn-6hH&biVb=*6d&:50hSEsa:4Ere5Y}qo7i8zK<7hw#DL-X^7VPh{V3[QK+nmBjhNdY?<8c^nv$J-R7GE4+<7x7^wAH]r2c&Bs<6*unsoY+lRq(1a<6*6m>?WqP>nYSXaec!KBM[!W08xB9pn<W5L>w1J3@<nErtULr5Y]mxhGqmZ<7+AuJX]^e{#t>oF14s@nbYZl&^YFC<7dOc@X/o8Nh?EQFc+4cM&Rg$0E{06<7dtMroCn0DRO!^S-PpLynz{&6=-xCFmRCI8Z19J&J-a>p98q)-nGVT?CaLs-}nh]+b4&Y2fc(epjdyd}XNq#hSFp^9pt.+un-^(R5USt<7FS(3#+fm:(vm%pixXVL&[72::^rw9gc@<vAB4j6DAO-rme*)(#i@ohGq48<7+zgO&kRq:Sb}/9d=f<xA5%9swzf6rF}!oHzF]&:F%XF9eZUz[Ym.I:q*d:pc($wN&ru1Hu0o{<9&{bI0qgH:q*9S9vzQxFZXQK@#B&[Faq7HAMxe5>UsB@9e9xFpb+{o<C]Pbpc(%J]Yd0y:wbp<9e9x)yA32t(eM)D<9I.umoFO(<J0&Jpc($YBMxDG8+@+PV9-UkBMuiUsr91.rFQPH@cr!:[w*WjpEx>imd)N4<BvPnF5?SgN?g8P7uiynrF}As1bt&X0onIUFc@jmCMr+[8?m$gVa6!lvAr%jGNZGri6>SoRXHhc3YTl7Fve9C>NlA{]NKy^E[f%!]XLbzhYNTz9e9BUP#XKH8@(n:V9^Y]DMc%K8l{D/&^hS+O&r{:>&eSZS^t8TVXCY@aDs^&Vir@xIY(6K]S/yip3S%xrnK]2z*$z6JPY@t9YD*Vbf5P.Vt)%ej/fiL/Yc$)Va-dRlce8cGNV-QibQXnBXLQfh=U+F<98gSJYg@V{2yli&[=RN!zGqh(9m@YSZj.^DzVh>h(.r=<9U-s?z-GA5959$9f8#rCK%<zh@ZKL9d=ghFYF*!>]tn}:c1B.Gn9b!:z>+K9eZT?lce8d<q--(Fc0.>HYwAJ(3@/Q:cV&DGncUo:z>#m9eZV0Vb>%L8YTBI<8JduJYd@l(eM^8:dtjfGn7w0<z!RHpd?A7GYIfk:wb4xpdjeyGYC)0>$TUc:ctTWlXZ9]"
"<z!&kpd?Bm1MXP/08t&3V9-N(IYkzZ(9mRR:d12ylX:Sy<z/5@pd?CABMzFhh+(VW9eBO{FYJuE>]tGV:c1C[1boMw0c8MgFc@hHnc4W(<BvoLFcVcPHYz@h(3%3t:cV<S1bs7<0c8=[Fc@iVXb+Mo8?m$g<9gM7JYhE](eN0*:dtku1bm!N8/YrcVau?MIYpX@:G-1>pd(KaGYGzF>$T><:ctU&VXA#r8/YJ>Vau&-3MFbJ0i@/IVay@PIYn}h(9m&t:d13MVXEH!8/Y:OVau<[DMh0]h)GTe9f8%TFYM>@>]tZx:c1E6Bb0B+h/Sx?<9EQ#pb?iP<L$m3FdsIrHYDEW(3%m5:cV>!Bb3%lh/SQL<9ESdZbL808@>}V<9&{!JYl0y(eNjK:dtlIEMUl&.@FOJ96MG30MhLajrcgAYmE>o2%64/Y*K3Y96Mq3jh3<F.Za[q95+Q).nzG5Y=n]Wq$m:BijG}*ko8z}K+PhxhawMu-2HIHa^=I3whDWrc{#$hY0K!/@YAO0Y$q=X<h<INe<egjlC7*:l[QO&bTK)KjvyZk95H0C}9dm]ZQ}AKsCu@S755X9iVH1Cq$m^ChkcevO[Vqd95+FIisGZdZ#JCILL-4f3+.+$Z}nC=yZ@Ot/1-m?dqPFRyZdX!Q^w(cjrcn{l<]-97h]ir4mkY.CTokQH]>yc.&YI/962WosU#sli7b+*yZW[t<u?[GkFbsXLM0Y8:V)z*Z)1VeYv/Aj@YDzzjD]O9<hsoF[X.>[k#DjLN9]#hyDlb$uB&4)L6:2Bg^FUFiVHk(l[#P>XJT*aYP.V{YvK^T+E{p]Y$r4jl[t%./*0a9jIfvY96MG%?Ah=tZ2NbFYw73V+E[(NZ#AhF96Mt#0Mi%olx?yhLK}i4:V[?hjIfgW<hsnLe<9}TZ2NqHYwtU^2eGon.[j+/Yv1KR-6SIfjD]O9<ix=J@YGk^jMBZHl[QPr4M>Z)Y*K3Yl[t}<ix$d6.!YDx95G#4f#E/lko8+#<h5RMk0f/GjIf1ULLh^ik0iSZ.[jO^963k6f#B$>ko8HDyYRA[g^SfCjrc8[l[QO&e<egjlx?yhl)oD:1vs7..:C+pl)L9n3+X}PZ-kf$<h<FMdl<3SlC80=l[QO&bTPb]Z^GCpLF/)fyQ#]m+4>JYL:1c}+G7UWXG0^}yZdYn->YjyjzUc!95+OLg^V0*jIfgXLL-4f6H/Th.:ELsy+-}[-(2D6gI>7WlGV+S-&K2D<b}vwMl@eupcoHPly4U&BO<AM7Vln{Y{57P<h<JH@YDzzjMBKF95ku{}[oR{ksx3DLQ5XY:Wc$<Y{6>^Ye-Ey:T-SpvuF9fz*3=1:9PCY4C$:oQPh6Gi<[1+orm)$q(W!Rix$d9N&BJ6l[t%./1*s?Y=oi#Yw73V+E#@VY$r4CF=0eHx1+?T.:BspL:/1(CghR$^rA1?^loV(4yyeH-n!yoIW?(fFY3J[jD55US$Anf.nA[wjv#[D<i2boT#)(P.@FW4LM0V:hO!v)Z}nC=YO/nt=ogdjjvy!:<h<KC=o973.YgOI95kvy^zD}Xlpr0mLK}fYk0hjR.YgOI95+PG0MhL2Z-kf$LKVLO}9b&.jD]rM96p?NhO=@L.[jWp95H1x.nzGojrcKE<h5S3a&uX%ko8HD<h5S3a&uX[ltM=>l[t@sbTK)CZ^GRoyY)6Aix][UZ)1ySYv/zpbTPb)lx?qWl[QQl!iUcsjvyKi<hso@?AjgH.YgOIyYRB&3+WJrZ^GCmyZAr]cCXVWltMXuyZdW>2eC3Tlx?bULLEBV4M]%GZ-kCHyZdX!/1!]LY>^*r<hO(#@YC0bjMBC%l)oEi}9b&VkJxO$LKVKUdl?RKlx?qWl[7t*cCU&vkwQS<<hO(#[XZN^kwQS<LL-4S}9dm$Y{5mRYv/y*k0k55Z#JK3YwtVl[XVtdY>^.!Yw71sjh13clC7]nl[t$h-6R8-.!YZ]yZdX!*R4zhY*KqlYv/Aj}[q4ikFbsX"
"Yw72n3+:gbY$r4jl[QQl/*1JuZ2Ny2Yvod!k0k50.@FOJl[QPr1vzd(Y*KFnl)/.f:V<OKjMBC%<h5S%@YC06kJx@2LM0Wje4%#+kFb=m96MHB->ZSWjIfDj<h5RMf#E/d.Yg+K95ku0bTJFhksuWj<hsn7f#E/d.!YK)96p&I0MdqOko8WFLKVMJ:V*+1.[jslyZW]@jh5nN.:Dat963mz.nCrUjrcgA95H0C[XW:EjrcgA95G#IcCS1{lC7-ll)L9n1vwsCZ^Gu-LK}gT2%Q*dY{57Pl[t@sdl>C]lGun5l)oFd-6OnIjvyZk<hsn7f#E/iZ-knFYv1Jj1vtH6Z)1ySyZ@M}cCQPVZ?:j<yZdWydl?RPkA>/R<h<J40MhL2Z?:y(yZW}*/1:WajzUr*<hO){=ob[Cll50ILKVLO@YzeXko8HDLLh^Wdl>C]lx?4dyZAs&{p#)>jzUr*963lE{p#)/kFbPkLK}hN//$X$Y{507<h<J44M]%LZ2Nj095kvy/1-m!Y$qXfYvoene4}eslC88pl[t}<hO=@QZ?:y(YvK^f!iYw@Y{5mRYwtVl}[i$hY$q>hYw73h*R7kGlGun5LM0XR?Ao&wY*K3Y<h<I9k0m(w.@FOJl[QPr1vzd(Y*KFnl)/.f:V<OKjMBC%<h5S%@YC06kJx@2LM0Wje4%#+kFb=m96MHB->ZSWjIfDj<h5RMf#E/d.Yg+K95ku0bTJFhksuWj<hsn7f#E/d.!YK)96p&I0MdqOko8WFLKVMJ:V*+1.[jslyZW]@jh5nN.:Dat963mz.nCrUjrcgA95H0C[XW:EjrcgA95G#IcCS1{lC7-ll)L9n1vwsCZ^G!qyY)6Aix][UZ)1ySYv/zpbTPb)lx?qWl[QQl!iUcsjvyKi<hso@?AjgH.YgOIyYRB&3+WJrZ^GCmyZAr]cCQPVZ?:j<yZdWydl?RPkA>/R<h<J40MhL2Z?:y(yZW}*/1:WajzUr*<hO){=ob[Cll50ILKVLO@YzeXko8HDLLh^Wdl>C]lx?4dyZAs&{p#)>jzUr*963lE{p#)/kFbPkLK}hN//$X$Y{5BTYv/y*k0k55Z#JK3YwtVl[XVtdY>^.!Yw71sjh13clC7]nl[t$h-6R8-.!YZ]yZdX!*R4zhY*KqlYv/Aj}[q4ikFbsXYw72n3+:gbY$r4jl[QQl/*1JuZ2Ny2Yvod!k0k50.@G3896MF8g^ZlxZ^G!qYvogc-6OnIjD]GO<ix=64M>Z>kJx=0LM0XR/*1JnlC88pLL-5N.nCrUjzUG&<ibc@=o973.YgOI95kvy^zD}Xlpr0mLK}fYk0hjR.YgOI95+PG0MhL2Z-kf$LKVLO}9b&.jD]rM96p?NhO=@L.[jWp95H1x.nzGojrcKE<h5S3a&uX%ko8HD<h5S3a&uX[ltM=>l[t@sbTK)CZ^GRoyY)6Aix][UZ)1ySYv/zpbTPb)lx?qWl[QQl!iUcsjvyKi<hso@?AjgH.YgOIyYRB&3+WJrZ^GCmyZAr]cCXVWltMXuyZdW>2eC3.Y>^[*l[7up2%TTEZ#JK3yYRC=/1=8fY>!7&YvK+=e<9}RkJx=0<ix=61vs7*kA>SPLM0V:ix[Gulx?bUl[7vj!iYw>ltM@)LLh!Q@YC06ko8WFLLh!Q}9hHAkA>.aYv/Aj}[ni)jIfoh96MGG{q2.gkFb=mLKVKgjh2Ch.[jWp96p&4bTPb@ksv3n<hso@^zG=dY>^[*Yw72n3+X}SkFbPkLL-5N:V)z+.[jWpyZW}**R7kNY*KFnYv/Be+E#@YjIf9f96MGG[X:zjjrcKE95H1x->T${jIf9f<hso24M]%GZ#JZ5yZ@MEix#MzZ)1:WYv/Be-6R8)jD]VQ<h5T[=ob[MjMBR#<ix+bf#B$?Z-k8DYv1Iog^TON.&%ZV96MFMa&toTZ-knFyYRC=/1-mZlx*@SLLEA.k0k50.Yg+K95+QA^zJPJY=n<}<h5S%@Yze:jMBn}l)oDo"
"dl!!elGt]1l)oDobTMqIlpruql)L8!cCS1#kwQD?963kKcCU&qlC7]nl)Lah!iYw@Y*KFnYvod!k0hjWZ#Jv1YwtUra&toTZ?:4?yZW]@hO.ZaltM@)l)/.f-6Onv.!YK)yZdW>3+ZuXZ^GRoyZAs&{q2.gkwQwrYvK=l3+Zu:Y$q>hl[QPr2%WE&Z#JZ5yY)6Aix@.$.&$6Z963k6g^WA1Z^GRoYvofh]G)v2kA>/RLM0WX2%Py&kA>/RLLECP//$X)lC7]nLL-4S@YE*Cko8<HLLh/L:V)z]jD]rM96p&I0MhL7Y=obEl)Lah//}>PZ2N3$Yv1Jj3+ZuXZ}nvmyZW{Adl>C$kA>#T<h<J$^zJPJY>!7&Yw73h*R30cjIfoh<ibc@=oe.*ll5fKLKVMJ=o97gjrcgA<hO(5e<c+}lGu0JyZ@N>{q2.ljzUG&963mz->W/qjIfoh<hso@?Am20Z2Ny2YwtT&k0m(BZ#JZ5YwtWg.nvlE.!Yv>96p?Nf#Dy8ll50Il)oFd-6LB#.!Yv>yZdV{jh2Ch.:C}r963lE[X:zeko8z}Yv1Jj3+WJwY$qXfl[QOwix#Mu.@G3895G#4g^Wz@.:Dat95G#4g^TOSZ^GCmYvoenbTMqIlx?bUl[QPr1vs7*ksu<lLK}hN//}>IlC7-lLL-3Ye<c+}ll5fKl)/.f-6R8)jvyvg95+PG0MeZYY=n@Cl)L9n2%Q*8Z#Jv1yYRB&2eC3VZ?:y(yZdWydl?RPkA>/R<h<J40MhL2Z?:y(yZW}*/1:WajzUr*<hO){=ob[Cll50ILKVLO@YzeXko8HDLLh^Wdl>C]lx?4dyZAs&{p#)>jzUr*963lE{p#)/kFbPkLK}hN//$X$Y{5BTYv/y*k0k55Z#JK3YwtVl[XVtdY>^.!Yw71sjh13clC7]nl[t$h-6R8-.!YZ]yZdX!*R4zhY*KqlYv/Aj}[q4ikFbsXYw72n3+:gbY$r4jl[QQl/*1JuZ2Ny2Yvod!k0k50.@G3896MF8g^ZlxZ^G!qYvogc-6OnIjD]GO<ix=64M>Z>kJx=0LM0XR/*1JnlC88pLL-5N=oe.}jrcKE<hO){=oe.}jMBn?ms&2?Q.EGdX/mQT8(cAum5-P+mh@iv9c5M9H}^<]Y5V-1vQqFxTmFYOni(zrz9+:%m*IyxYbN.Gxs6cc}[s3d^=ln/8)g^5OgiBKZh1=7M:vK?rc%Z{6xAaBK^RsJTgC#wX@<neHPv{7S#4ab(Qi?dzPJi{rd5q8H2=axXCoMETbi6RK*#ByK$i1vU0]C$l#gTFz%!%Wjv4$iCsF7qX/ovrC]@qJ9fJ)@<xi>>sS-*vfP(]m<w%EoO<=7%mq*bhKX5]xnFO[YIJJ]/X2WefSd0-c(BqVUN3?7qrfVDC*z^ft<xk<0{mie9bw}QcXC0BH!eh@])Vu1{MlA?OovDkER@Q)Iyf<sFRL=}cM1&/R<hOl!!%Qn6=xbGd<xD>ze>a&6Sa+Gm<wAIHupt>?X<9TUX=+8o/04L0w%726m1dM@j9/)/gEJ$ex9rk]P$lIeX/mKgx%YUea?y2EsWFuaMdVEg}]}DUX]8Qw8(cnCK?bBcX<A=il)JIW>9FA}l:(J@b^6MPRT6tWZ2RF>8[[+J5Xp)tmDvN}<x+op@Zw]KmuA-=Y/QNc!&&7>Y5We-lQ0}9Q+[v/N&wtpW$hcKr*jM7!r1aIZ79NJU0}zTmDwD#8&l{3FW!6!?5y{PlXh:YBxxGuMpn6*<Au[fGCniwXvXZ38?%umFW=UQiu9/oV&NXsJZ>4Gn91W18&^AcHx3p+Xi]2D8&k2t^ywm-}ZQ92<A6T5^yqDDXvZ/h8?CH+dhLv865W<Wlg=v:l)$y]p^RQ^YehlSHo%rRXn9]WXZv]ON$>.)osC(bm#o{6W11b7m[B:(H&<A)EJhl]XA4*ulY[(<V)icjIw}i+x?1kG29*c:eno++yU2k/bPgvVHc#/XKoRB{H$AV3.doYsy=HE!r[MER"
"Y>#Ii<AyyXDRnwkGoN-9<AnH-k1BWU={FW3KPz$dNgQx?Xrf7!XxbB@P$3ll+&B3)J(iYogZVS1c:.phL7)hOsUeN@I&#f1<ASQ@hLVF6W3l$EWml{RtqJ@:vDqMqSlN)Pi[&UwH]])^8&Gv*I7b%OoxpC^k<1]N]z(H)S#DZnY7>c13+81KY$o21LuU{>K*@4ZLslCfU!ASuI+=S2^6AEYY+d^nD!^XbbGTmZ<B3tlFW<Z<q/53Xy^CO(CYq.)n3duL<BfPrJZ(Fe4DI^fnUeDkBTMl&XmMjvWR2*9BiG=(m#Orz<BfR7@ZQF0ISw=!Y)-FP8%#Zddcb?E8<G$lwB{*HW=zFb8&jDj5W[?OWI@#tasea]tN1#ya-*YvH7Q$w>c/s>nO.09L^AU8jhsj0]sQqBBzG2rP$xw.WSs5)8/b/vO7UIjcB?8dm)1vWvyCl^DMc&]BKvJ-MS>zvXs9Db8!)5ax?u3UWNsv68!NQXK*.{pnGh9D<Dm2=<JFO.I<m-i8!3fa>JT0@y-v1&lv0n9LQQtfjWMP@ZF)cVD96DL98q.^8/bNzGHs*+ooQlLmD<.AV#O@dWO!oS8!IbIx}e%[Wb^K.8!x>4Q/?!KYjZ+!y/ApbhP]]JUX.Z}<D!X1vyz0wbtP:HllnZ.Q!c#@n/IN5L^dVtXL0(IWWHTl<ErN%wBNiJWWC[bCoKz#O>f(j9<9?#U!AX9wjbBA5Pxp+8!Y#ENv]^[nGhBM&KI#Vpd-fyWZ:&gA2}QiLR/@?Jbm]CoV-BfI$Jpof=L14Y*>kK8}Yq8nXG1kx)@Ty@Z(gv1Ip.LN+L]}Wctp^X4<Q/Yg$@f=pc/=n-VEmLdrRA8gXkOW=A&ABzE+rEO2^mmnY7pLZe!WO>j7p![XAUYj<3&2@#^cnXKo9pB>kYL>s3JaMPa=Ls$IWQ-(e.nKOm3<D-{/zB)l/n-USkB?QB.owX-(Wh=DV8=XZlrui2QWtjq$8^ry/6-qSz(kJ>.B+X7@.l[CA+dXftyO${irt3Lg)l/G9YD/hCov1LUq@<TWZ8/(=qP[d7n/H.JLZz.iI>a{Nv.rA@8^y5qRMn^O9yGeINypzNg$MN%m@{3E<F0(z4O182XSY=XLuU/qsVi%K>r:J9cRlUr7ETIvg/lnlM1=+NyVcjfWta^<TJ<ASth<qho1S}Lpn@LQR4QAaZ)<8[LwyAn4Oy[gWhSnX8=<Ot(^p#0r>Z9yp}L+}a*]CMlR#62pxv(:S976Co5#X?<FHayup{Ma)C$Axq+/@%>xM6-lyEml<F/=LSwHazo1VbC<FF#[EO7L=n/HTnCz5m3rpavIWDaj9L$RA4up]B1ZpwJm8^ce3Q^qZon/H[yB0n5R!yN7ff}2x>8^f{(V{w[!Wq8b1yl(eEl[Sq8n/H*AcQ:FvOcqp+SzLN%y{3xDW-+fgWkFS8o5rgQCb[1-n%I?38=@7tqmgZX98K0-ldYdLQ=#NGWl@!yl/Y>*>%N]p1-fz*O76W-SLsS.WhSrX8=WrU-6ElHhrp-ZmUd%&sU}7514hGG<Fi?os!rmZo29K8<GqasFWJpsr2YdoZ7.THU1B[ro604?<G3eUU1E=H0VJq*mi<G8cEA+QVBPIh8.n=&1u&{$&86kOkk7+pnv96uC&hXVzLMfMr{tZ57eHI@UW&nF(WIlCrLr%(HhGN4bv](.VC7=ZZiC&cJ.ISAVkUdx8.0^kFWi1keGl4=xySo{7)t)>o]gQg8.g7>cL&).VC3XE8ZhJF^xX7oo=CO&<K7o]?cR9Dx]Skao3Qc*f%Tg(<WL@Hl3yQvGCY[)Vo:s}UvFRwHjR:>tXdCezS[z{rdzdQqx45onX}*!6:e!mo&ku88ZVRy5WD6j8NHofmL>V&lf-5eHI(lGyp6h=/0!4H)@6}^V4(%he>}cRVoqsaon*6n1xcl/]gjQfKRa)I(2m>u4Zs{:8.kqnccA@wVGf^6Akw&>=p]l6v7{S^y(RjQHq@r9eKty>xzB9C5WEO/o]g=S"
"<K7n@)Of(!(eZO1lhP<+e0YhLo&27&<K^R@!gOt?m2?Uv+fA.?jq{Qup2%Zc?}NrWI7.eG+h?sqUv-xKmOP0:.EWtOYqS5#8)cvMFs0ROXYT*49.JPbo[}uW8.K<[iv%mAFv3fPY&G-@D3A{N6{cJ9YqXuZ*IC%8(eZG:V5RKAJ.OP:PbPx<AR)2FCrOe@o(x8!W%0*-)OcE@4=1=)zgki6=q0xRp2$e98XGtSql.5Yp8:>HGAntV+yHPGJPe{4Fn4E(89ZXUV98Yp8YFb8&JZJdPbO?D8XGvo(G]JipupxsB}I51m20VA^D6FIaBrQZa@EHrU{1X6UqdEAf%&h0y?owGn=zEP)Ol<4I&H?)nyOv[6S%.>CUnt3V.Ugs-28efCA7a!JLR>)[V+wFZz!3BoziP@m8IrGmnc]YJUI#VAJ@uRVefQjUpYJC?t1f166eEQq74p.>dfNp25BBb<MvB{yJKNZi<!S6JP/wxcv6l<bd5bfAX[F&:0m}(V2RO78Y6.r7i<<R6HtklxmWsdf%O7<bAoRQD&u=3>LpuWpkpCO<MHZD(G%DfV/G!+W[bRFw{Gx?Hk<=fM#<4H35t]56QS)Bk$M?J7#gH26Dm{Bk$l9^C6DqAHMfcPN26O98iXTlZ)XMAx3bCX&JYsyF-Mc8<M@LFZd%)?DPj+[Zi-qP+JZ2=ALw%rJMZT{26^>!U]sUF8Y0XjV>L/0(P[7%Zp5n98%.c>C!QYyk$*KJ0KowEC^OI*GBQ+$EOii}owONk<I1of=p]lL/w#mwl0^q96CPS!V5N47k%.^p<p(l!pc?Jr<Mh9R)OmzF6Y!z3JM$pU.lJ7+*?hTxJMTPF(HS>Ypp7qZ<N8clk2ghbp8-^v8U72)0K0JqX/gf+E/V.)o7T4{3/}Au8U71cP%6GLtPiPrp@Lvi[VM:Yp)*Eh8U711k2w%gEJ4>!p^MB7a{3^np?<rk<P=I*NJ0gtj9jW+J?JUjS}mk4otXSOqx2iD[VPv.S+M[+AG>4FXS{#n?m)U[ZknZ9}>RfNY$Fj?<QoNk}9R*]p^%Z3<QL7{EO^Xo8GT+iUAu2POgda5k/nDeW.e*.{sjYLV65?tAp!4n2qm92U9*/NNZzx$XSraJPQmf(<Q5OB5WdlVp^^({ZuEXB3*$M[uf8L)M7a7RXPW6qfxABjz*y!wQ#DHCU5lB1n]vEJ.loc4p^yQGB@<Gw2zb^+59i@48VAi]K%0r>l:TkcplHc6b>h(83/}H{AxeY-hiO%GUwY[(<P:z6x>7}?Q)HkT<Qo%+(/9Y5p-w$MAzxt:6:U/K08IPpAC3qkBsZ8d)q+51o!ng@e(ogObY-LknQ0$(1FY!:!>lG1y/.)}U24WrAodjCYib*e3&VrhfN3g)l*lR6bUkSQk&&kA8U72mj6-=#^b+:lbw{9)9-7>#-FB5/8VAL<0K40#p-w(cW&}}l4P18pgp^xXOk8]x2q*idUyntA8V56M.lorgp#@aep{FWm40iA+tZ)u?:)dH%@.Q9Sp#}.X8SRwD&JwdfU0y%G8SRxu&JtDifb&$<k&F>T>*7ggtFM@nJpB4?>*9(1UcNWf8TE/NP@#nl3(.vc8/-Ifj*MxPQ):RqOe$JP1P^Y$APE]ndB7!VlG^CKq9:9y8SRwD>59DyUcNWfk&-zu<pmdv4RBX5k&#@@<pmdv4Nxtk8SW&+tLH{S/$T]EZOGKSNvebX&=xyjpZ=]-!I873qbYUQBhZPf5W29bUh]*WZhu!Fe(xy6zw}NCx5<Bu>q.Mqyicg/<RH42s21/W2xVq8BFnON]Q.x2qi)DkBU}Mwh=Mf3O{2Ec94bwq(Wr@*Z&e%:z%L!%k2EY/SHcBTDKG8Mxt9yFqh>4F8Tf6t5W4sqXAaP*Zi0brg=JbPL&SXbNUFPs8l+4u)3{MgMu*IsxayC#2$6M>Mg<798]X8g?&$:%PcS2b(gt[Uo}Mb{8Tm^)+I]l7=Gc1..A>uql9w6z?Bkw2Z9C=38%38{j?&JoyodbD"
"2hWb*Jj}yO<RYwWD5a)u7lGUn.o2JrYD/>zA(b8:8TGR#*PwDhU5vMK<SmptCmUXXq5XY7<Q[VS9-j3ukiiF1ZZ:jff#Se/x(laixZO:V7qrp/qez7T<RW{d)OPf>.SqV(8TZFkpeQ9(&/H(KZVm#F9-i@(qrnO4<Smpfl9h6Sg-B=zZAGZA/])u74DO+nk?ap4![&8AT.bQg8RAl{a*4dQTS2()U%CHZm!EXPcGMzM-Mn(l?CN}3T@u>L8VAIL^CzHhjS{s]YPMMxGA@}[<R(RTx3Kpt![*:h4i0:Ck*=Z}bN>WBLr^RbZA:mN/e2Y%{zTk^.g[<Pe(C27jMOg[mld(tra+kvksVWTZ0@EuBtV=Kmh<X9BxiiLrn{fY/Ed]BppWy^ctWqq[[b/H<Ta^Z/ed9T4vGa<8R:*}J-3!DX8!%=y*H2I>f6ykT.=+/LznCwbS]#+r]>X6OTb%*dyV]pgJUCZB5$X!R.}j7T-L:okC}dFWO&6?T(Q$Ox3KqkV<<+#+]kn#K+}0qdkQ6mo{F%zyuX])g!0XN8[^a$OMMh+jm:xzGYKD*Ni$xv/]]ep:EbK-Vz<7g7hxnitXi-d<T%D1AJHejT4J0PVAZKDd>r)wqOfG!yd(N->[7H}-i@sZYz-XZ9-xVV<k?C^ZpTkcoy+*ZxUL]WJh#Jrf%ehyqFiSkp3Tdvlvs.N{Jfc4ZMCS:Tnl[A1I[4?8R{t<)OS5Vp2#<.8QUXGa*2A}qyy=Q<S>Ha?CN8y:k[P9/LTP4^xhHTlB9gXBs7!+duiBwqKPAK<TH7<k2PcxD/V^LZj5$}a*6&UqBePa<TC^1^xbWITS2<n8Q]h3P@*P1y-p!j8QUW0Xxt/wyE3pHBO(vZ/VCbWT=fb[8RVYa^xgeB[/3tv8QUXGduURaqSP798Vt<BD7yw.nX-qfb}k}Sb3W]YqywRBo({/HP@*O@TS2W$k/s6@^xd2/4m51Sk?oZ.V<^@ddNqN38Q]kTeaC)S+rT1BZSA0s8i}@[Vb17Ko)J1Bwy>)KqCzn+ZElZ}x{zz:qGDy5ZElTG)OUsG%0k(/x2gPm-+W]WyI7z+JhyzR.l4s<+nO@Fx0.$^3?Ia(h<)*lB9b:&cPDfWhGK)RNo8wouq+<8qSPn=o>/&kNuQAUi]}ThBeH!dxfB1+f)nG?8RnTv*.kTKT.b8^8QUWCV<Yjv+vXO)Nm9mNiHmA9L:M/XZD4W@du{lXqOLkVNmxc5MO+a4NzgvLNyt4Fx}d[0g[4(5NHb<W5?6ieqGo*TALttxe(A=u&#Xg$ZCIeoea?{GMa}3><TZ8!dt$@&qSPCbZnyni-]$db6$>Nn.i]kc@2={j/1}*D8R=FO/2tGw-tWIRlZfQFqlp+(d3Wha.JXt>tTrkT*#PYYl+VCr:XH}wfxB/:lW3GhbS=}3qB?hS<Td#mEO$79FAncH<T%C}[Y*eAqnj.YM)RjZe8%e}qOIM@<TY>^peZfY)=u60bvjpfcIwweG=nCtl6no<x{fE/R:4C*B8gwj-0%2p/1Tkex2@bNMO}bn:A8Cmw@h>k-+gZ1qCnP7K#9PZg.h<v{t(s@<8!/I)OPjhT(R88ZA:m.I$m>w]g.JMNmx7L?COzHT-F}?8R:>Kf%c<F[/3oP<SQ+/eaKRbC&YKu<L>4}[Vr.S^{$+Db9(UH5<--lTTaU&Jb&[R*A[o!]k^qd<S@tv=qH5!4vHv.k/6M?g-wMuTW{DkkG5$Gl9i&:3tbKI8RAhbhF!su!R{}Eo<o8jrmZdOTY2ZcB4f%:^xdxhMM1Eik^Q@vctsyIgAa?YB2$^ZR-5eHxld+gk.?qT!d3V*E)AmYk^vCViu/x{LP53XNkO}FIi^lDT!9e18R:V#^xgew.Y4q3x69*xAJLhwL8N15Nyl}qR-L*KeosW</YpZWJ-5c.hSVo4No&=R?CSCX4i183M.T$](6rmx<?nOv/YLka9LN%c0OV%twkgUpdv4Uay}&76Jj@i8ct.=]]rM11<Ta^Z?yN+k"
"(9&&OwQw48ik&c4qSNmx/+{4od9HZ8qJe/<<TR&DXxR=IMYfi/<TDyLV8%6JqSPyO<TR*^Obc-:+&H8G<TYG0H@(5=zJ70J8Re-td$qf=[-5hvZA:h0=qG/LyVM$*JeNB/EP0%ETW7ffy+#NZe2bxKqSPkt<SQ->.l1Td+nO^i*r)8?Zevk63(WJ}8QUVL-0!]H3{.*F8ReW?-+QKe3YH@#x0.%+![WK]dl$G]k/NOKYyQ8F]w}fu<T%CLuq=-X:=99jJhywQ+s9t%T*z4o8R}(lAJLhdy8UH4Jhc(&Hhf5<T.a&uk/NR8J-7^zT*kETw$7ldQZJcG+vWDIJj9FkN9d+43Y!.#w%**7ea.<Hskd.%8R{)/r+^O]du7csVxLZ4N<<qUT*jU5JeGhK![.V.yE36x8Q[C*r2lFc-x.Jbk/U:%jm!>iL!P6a<TY@dl9n!QT*73W8Re3BcspvpZL?dP*WRx:?wEUDqrnF^I/8B8FVX2V2%=^lZCKK6du[)WqGD09<TY}<peT/4]w}A=ZDLLc9L(t8:@nWfx0FE6-+HEqT.bgEk/U+zdu{kZ[%II5<TY}zuq/(n-tWkxx2C8@/Tx(^{(7aY8RAirMt8q/:W0wV8ReUZWRouOrx>G9Jk/MI/.dOq3Z<:/B6xyjV8$s0+4vX6w#J>TN<*R6dJl%YVBE%]*F?^(4yFOjx1NVHjm*iEL:Lho<TDCejm^I$]w}P5ZDLK&zhYy/b?Z]=8R}YgFVWY2vRbIMkSDHpK/YLUqrnB4<TvM4J-3!DxYQGFk+p>FyXcTnhjmiq"
"Nmxc@du}OhqOL6t<TY}<peZg3hruC*ZC2S/"
};
static bool decodeEmb(int n, vector<array<double,3>>& out){
    if(n<2||n>EMB_NMAX) return false;
    static const char* A="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ.-:+=^!/*?&<>()[]{}@%$#";
    static int rev[256]; static bool init=false;
    if(!init){ for(int i=0;i<256;i++)rev[i]=-1; for(int i=0;i<85;i++)rev[(unsigned char)A[i]]=i; init=true; }
    // gather full text
    static string enc; 
    if(enc.empty()){ for(size_t k=0;k<sizeof(EMB)/sizeof(EMB[0]);k++) enc += EMB[k]; }
    // decode all bytes lazily once
    static vector<unsigned char> bytes;
    if(bytes.empty()){
        for(size_t i=0;i+4<enc.size()+1;i+=5){
            unsigned long v=0;
            for(int j=0;j<5;j++){ int d=rev[(unsigned char)enc[i+j]]; if(d<0) return false; v=v*85+d; }
            bytes.push_back((v>>24)&255); bytes.push_back((v>>16)&255);
            bytes.push_back((v>>8)&255);  bytes.push_back(v&255);
        }
    }
    long startPt = (long)(n-2)*(n+1)/2;   // sum_{k=2}^{n-1} k
    long bitoff = startPt*36;
    out.clear(); out.reserve(n);
    for(int p=0;p<n;p++){
        array<double,3> pt;
        for(int c=0;c<3;c++){
            long bo=bitoff+((long)p*3+c)*12;
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

int main(){
    T0 = Clock::now();
    { const char*v=getenv("TL"); if(v) TIME_LIMIT=atof(v); }
    int n;
    if(scanf("%d",&n)!=1) return 0;
    if(n<1) n=1;

    if(n==1){ printf("%.17g %.17g %.17g\n",0.5,0.5,0.5); return 0; }
    if(n<=EMB_NMAX){
        vector<array<double,3>> pts;
        if(decodeEmb(n,pts) && (int)pts.size()==n){
            double r=geomRadius(pts);
            microPolish(pts, r, 0.55*TIME_LIMIT);  // undo 12-bit quantization
            optimize(pts, r, 0.9*TIME_LIMIT, 0);   // opportunistic improvement
            microPolish(pts, r, TIME_LIMIT);
            for(const auto& p:pts)
                printf("%.17g %.17g %.17g\n", clamp01(p[0]),clamp01(p[1]),clamp01(p[2]));
            return 0;
        }
        // fall through to solver if table malformed
    }

    // ---- warm starts: grid + FCC + HCP ----
    vector<array<double,3>> best = cubicGrid(n);
    double bestR = (n<=20000)? geomRadius(best) : geomFast(best,0.4/cbrt((double)n));
    if(n<=20000){
        double tcap = 0.45;
        tryLattice(countFCC_r, genFCC_r, n, bestR, best, tcap);
        tryLattice(countHCP_r, genHCP_r, n, bestR, best, tcap);
        if((int)best.size()!=n){ best=cubicGrid(n); bestR=geomFast(best,bestR); }
    }

    // ---- seed from embedded 120-pack + maximin insertions (n slightly >120) ----
    if(n>EMB_NMAX && n<=EMB_NMAX+48){
        vector<array<double,3>> sd;
        if(decodeEmb(EMB_NMAX,sd)){
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

    // ---- local relaxation within time budget ----
    double TLfull=TIME_LIMIT;
    TIME_LIMIT=TLfull-0.04;
    relaxOptimize(best, bestR);
    TIME_LIMIT=TLfull;
    microPolish(best, bestR, TIME_LIMIT);

    // safety: ensure exactly n valid points
    if((int)best.size()!=n){ best=cubicGrid(n); }
    for(auto& p : best){ p[0]=clamp01(p[0]); p[1]=clamp01(p[1]); p[2]=clamp01(p[2]); }

    // output
    for(const auto& p : best) printf("%.17g %.17g %.17g\n", p[0],p[1],p[2]);
    return 0;
}
