// Sphere Packing in a Cube — solver (C++17)
// Approach:
//   1. Warm-start constructions: balanced cubic grid (always valid) + FCC + HCP
//      lattices clipped to the cube (binary-search spacing so >= n points, then
//      greedy max-min prune to exactly n). Keep the best geometric radius.
//   2. Local relaxation (Lubachevsky-Stillinger-style growth): repeatedly try to
//      grow the target radius and re-jam by pushing overlapping pairs apart and
//      clamping points into the shrunk feasible box, using a spatial hash grid so
//      each sweep is O(n). Keep the best *actual* geometric radius seen, within a
//      wall-clock budget (config allows 1s).
//   3. Hardcoded near-optimal configs for n=2,5,6.
// The geometric radius r = min(0.5*min_pair_dist, min_wall_dist) is what the
// checker uses.
//
// Build: g++ -O2 -std=gnu++17 solution.cpp

#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <random>
#include <cstdlib>

using namespace std;
static double ENVD(const char*k,double d){const char*v=getenv(k);return v?atof(v):d;}
static int    ENVI(const char*k,int d){const char*v=getenv(k);return v?atoi(v):d;}

using Clock = chrono::steady_clock;
static Clock::time_point T0;
static inline double elapsed(){
    return chrono::duration<double>(Clock::now()-T0).count();
}
static const double TIME_LIMIT = 0.86;

static inline double clamp01(double x){ return x<0.0?0.0:(x>1.0?1.0:x); }
static inline double dist_face(double x,double y,double z){
    double df=x; if(y<df)df=y; if(z<df)df=z;
    double t=1.0-x; if(t<df)df=t; t=1.0-y; if(t<df)df=t; t=1.0-z; if(t<df)df=t;
    return df;
}
// Exact geometric radius, O(n^2). Only used for small candidate sets.
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

// ------------------------------------------------------------------
// Spatial hash grid (head/next linked list) for O(n) neighbor queries.
// ------------------------------------------------------------------
struct Grid {
    int G;            // cells per axis
    double h;         // cell size
    vector<int> head; // head[cell] = first point index or -1
    vector<int> nxt;  // nxt[i] = next point in same cell
    void build(const vector<array<double,3>>& p, double cell){
        int n=(int)p.size();
        h = cell>1e-9 ? cell : 1e-9;
        G = (int)floor(1.0/h)+1; if(G<1) G=1;
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

// O(n) approximate min-pair + wall min. Searches +/-2 cells (safe for packings).
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
    if(mind2>1e17) return dface; // no neighbor found (shouldn't happen)
    return min(dface, 0.5*sqrt(mind2));
}

// ------------------------------------------------------------------
// One relaxation sweep (Jacobi): push pairs closer than D=2r apart, clamp
// into the feasible box [r,1-r]. Grid cell = D so +/-1 neighbor block suffices.
// ------------------------------------------------------------------
static std::mt19937 RNG(12345);
static void pushSweep(vector<array<double,3>>& p, double r, double relax){
    int n=(int)p.size();
    double D=2.0*r, lo=r, hi=1.0-r;
    if(hi<lo){ for(auto& q:p){ q[0]=0.5;q[1]=0.5;q[2]=0.5; } return; }
    Grid g; g.build(p, D);
    int G=g.G;
    static vector<array<double,3>> disp;
    disp.assign(n, {0.0,0.0,0.0});
    double D2=D*D;
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
                    double f=0.5*overlap/d;
                    double mx=f*ddx, my=f*ddy, mz=f*ddz;
                    disp[i][0]+=mx; disp[i][1]+=my; disp[i][2]+=mz;
                    disp[j][0]-=mx; disp[j][1]-=my; disp[j][2]-=mz;
                }
            }
        }}}
    }
    for(int i=0;i<n;i++){
        double x=p[i][0]+relax*disp[i][0];
        double y=p[i][1]+relax*disp[i][1];
        double z=p[i][2]+relax*disp[i][2];
        if(x<lo)x=lo; else if(x>hi)x=hi;
        if(y<lo)y=lo; else if(y>hi)y=hi;
        if(z<lo)z=lo; else if(z>hi)z=hi;
        p[i][0]=x; p[i][1]=y; p[i][2]=z;
    }
}

// ------------------------------------------------------------------
// cubic grid fallback (always valid; achieves the baseline radius)
// ------------------------------------------------------------------
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
                         vector<array<double,3>>& bestP){
    static const double offs[][3]={
        {0,0,0},{0.25,0.25,0.25},{0.125,0.125,0.125},
        {0,0.25,0.5},{0.5,0.5,0},{0.375,0.375,0.375},
        {0.1,0.1,0.1},{0.05,0.15,0.25}
    };
    int nof = 8;
    for(int oi=0; oi<nof && elapsed()<0.45; ++oi){
        const double* o = offs[oi];
        vector<array<double,3>> p;
        double r=findMaxR(cnt,gen,n,o[0],o[1],o[2],p);
        if(r<=0) continue;
        auto pr=pruneMaxMin(p,n);
        if((int)pr.size()!=n) continue;
        double rr=geomRadius(pr);
        if(rr>bestR){ bestR=rr; bestP=pr; }
    }
    return bestR;
}

// ------------------------------------------------------------------
// Relaxation driver: grow the target radius and re-jam, keep best actual.
// ------------------------------------------------------------------
static inline void jitter(vector<array<double,3>>& p, double amp){
    std::uniform_real_distribution<double> U(-1.0,1.0);
    int n=(int)p.size();
    for(int i=0;i<n;i++){
        double x=p[i][0]+amp*U(RNG), y=p[i][1]+amp*U(RNG), z=p[i][2]+amp*U(RNG);
        p[i][0]=x<0?0:(x>1?1:x); p[i][1]=y<0?0:(y>1?1:y); p[i][2]=z<0?0:(z>1?1:z);
    }
}

// Continuous Lubachevsky-Stillinger-style evolution with basin-hopping kicks.
// Evolves 'best'/'bestR' in place until wall-clock 'deadline'. Starts from the
// given 'best' (already a candidate). Returns iteration count.
static int lsRun(vector<array<double,3>>& best, double& bestR, double deadline){
    int n=(int)best.size();
    if(n<2) return 0;
    double grow    = ENVD("GR", 1.006);   // per-attempt target growth
    int    grS     = ENVI("SW", n>2000?3:6);   // sweeps per growth attempt
    int    stagMax = ENVI("KS", 40);      // stagnation before a kick
    double kick0   = ENVD("TP", 0.20);    // kick amplitude (fraction of 2r)
    double kmin    = ENVD("TF", 0.02);

    vector<array<double,3>> cur=best;
    for(int s=0;s<20 && elapsed()<deadline;s++) pushSweep(cur, bestR, 1.0);
    double curR=geomFast(cur,bestR);
    if(curR>bestR){ bestR=curR; best=cur; }

    double kick=kick0;
    int stagn=0, iters=0;
    while(elapsed()<deadline){
        iters++;
        double target=curR*grow; if(target>0.5) target=0.5;
        for(int s=0;s<grS && elapsed()<deadline;s++) pushSweep(cur, target, 1.0);
        double r=geomFast(cur,target);
        curR=r;
        if(r>bestR+1e-12){ bestR=r; best=cur; stagn=0; }
        else stagn++;
        if(stagn>=stagMax){
            cur=best;
            jitter(cur, kick*2.0*bestR);
            for(int s=0;s<grS+4 && elapsed()<deadline;s++) pushSweep(cur, bestR, 1.0);
            curR=geomFast(cur,bestR);
            stagn=0;
            kick*=0.85; if(kick<kmin) kick=kmin;
        }
    }
    return iters;
}

// Random configuration inside the cube.
static vector<array<double,3>> randomConfig(int n){
    std::uniform_real_distribution<double> U(0.05,0.95);
    vector<array<double,3>> p(n);
    for(int i=0;i<n;i++){ p[i]={U(RNG),U(RNG),U(RNG)}; }
    return p;
}

static void relaxOptimize(vector<array<double,3>>& best, double& bestR){
    int n=(int)best.size();
    if(n<2) return;
    // For small/mid n, spend part of the budget on random multi-starts to escape
    // the lattice basin (helps irregular / near-cube / exact-cube counts), then
    // polish the winner. Larger n rely on the strong lattice warm start.
    if(n<=500){
        double warmDead = 0.40*TIME_LIMIT;
        lsRun(best, bestR, warmDead);   // polish the lattice warm start a bit
        int restarts;
        if(n<=30)       restarts = 24;
        else if(n<=80)  restarts = 12;
        else if(n<=200) restarts = 6;
        else            restarts = 4;
        double slice = (TIME_LIMIT - elapsed()) / (restarts+1);
        for(int rs=0; rs<restarts && elapsed()<TIME_LIMIT; rs++){
            vector<array<double,3>> cand = randomConfig(n);
            double candR = geomFast(cand, 0.5/pow((double)n,1.0/3.0));
            double dl = min(TIME_LIMIT, elapsed()+slice);
            lsRun(cand, candR, dl);
            if(candR > bestR){ bestR=candR; best=cand; }
        }
        // final polish of the global best with remaining time
        lsRun(best, bestR, TIME_LIMIT);
    } else {
        lsRun(best, bestR, TIME_LIMIT);
    }
    if(getenv("DBG")){
        double ex = (n<=8000)? geomRadius(best) : -1;
        fprintf(stderr,"[dbg] n=%d believedR=%.9f exactR=%.9f\n", n, bestR, ex);
    }
}

int main(){
    T0 = Clock::now();
    int n;
    if(scanf("%d",&n)!=1) return 0;
    if(n<1) n=1;

    // ---- hardcoded / near-optimal special cases ----
    if(n==1){ printf("%.17g %.17g %.17g\n",0.5,0.5,0.5); return 0; }
    if(n==2){
        double a = sqrt(3.0)/(2.0*(1.0+sqrt(3.0)));
        printf("%.17g %.17g %.17g\n", a,a,a);
        printf("%.17g %.17g %.17g\n", 1.0-a,1.0-a,1.0-a);
        return 0;
    }
    if(n==5){
        double a = 0.2639320225, rad = 0.33380;
        double u1[3]={1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0};
        double u2[3]={1.0/sqrt(6.0), 1.0/sqrt(6.0), -2.0/sqrt(6.0)};
        printf("%.17g %.17g %.17g\n", a,a,a);
        printf("%.17g %.17g %.17g\n", 1.0-a,1.0-a,1.0-a);
        for(int k=0;k<3;k++){
            double t = k*2.0943951023931953;
            double x=0.5+rad*cos(t)*u1[0]+rad*sin(t)*u2[0];
            double y=0.5+rad*cos(t)*u1[1]+rad*sin(t)*u2[1];
            double z=0.5+rad*cos(t)*u1[2]+rad*sin(t)*u2[2];
            printf("%.17g %.17g %.17g\n", x,y,z);
        }
        return 0;
    }
    if(n==6){
        double pts[6][3]={
            {0.257641909958, 0.742661633859, 0.378489876371},
            {0.742358090042, 0.257338366141, 0.621510123629},
            {0.621507841464, 0.742696876521, 0.742323942460},
            {0.378492158536, 0.257303123479, 0.257676057540},
            {0.742662776614, 0.620832032508, 0.257304266073},
            {0.257337223386, 0.379167967492, 0.742695733927}
        };
        for(int i=0;i<6;i++) printf("%.17g %.17g %.17g\n", pts[i][0],pts[i][1],pts[i][2]);
        return 0;
    }

    // ---- warm starts: grid + FCC + HCP ----
    vector<array<double,3>> best = cubicGrid(n);
    double bestR = geomRadius(best);
    if(n <= 6000){
        tryLattice(countFCC_r, genFCC_r, n, bestR, best);
        tryLattice(countHCP_r, genHCP_r, n, bestR, best);
        if((int)best.size()!=n){ best=cubicGrid(n); bestR=geomRadius(best); }
    }

    // ---- local relaxation within time budget ----
    relaxOptimize(best, bestR);

    // safety: ensure exactly n valid points
    if((int)best.size()!=n){ best=cubicGrid(n); }
    for(auto& p : best){ p[0]=clamp01(p[0]); p[1]=clamp01(p[1]); p[2]=clamp01(p[2]); }

    // output
    for(const auto& p : best) printf("%.17g %.17g %.17g\n", p[0],p[1],p[2]);
    return 0;
}
