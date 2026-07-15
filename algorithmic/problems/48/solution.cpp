// Sphere Packing in a Cube (optimization).
//
// Place n congruent non-overlapping spheres in [0,1]^3 maximizing the common radius r.
//
// Key reduction: put point set p_i in the unit cube [0,1]^3 with minimum pairwise
// distance D (points may touch the faces). Mapping centers c_i = r + p_i*(1-2r) keeps
// every center in [r,1-r] (walls satisfied) and makes the min center-distance (1-2r)*D.
// Setting that equal to 2r gives r = D / (2(1+D)), which is strictly increasing in D.
// So the whole task reduces to MAXIMIZING the minimum pairwise distance of n points in
// the unit cube. (The cubic-grid baseline r=1/(2M) is exactly D=1/(M-1) under this map.)
//
// Optimizer: start from several structured lattices (cubic / FCC / HCP-ish), then run a
// force-directed spread — each point is pushed away from its near neighbors, points are
// kept inside the cube (which drives them onto the faces/corners), step size annealed.
// We track the best true D and convert. Deterministic; hard wall-clock budget under 1s.

#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <cstdint>
using namespace std;

using P = array<double,3>;
static int N;
static chrono::steady_clock::time_point T0;
static inline double clamp01(double x){ return x<0?0:(x>1?1:x); }
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }

// deterministic xorshift RNG
static uint64_t rs = 0x243F6A8885A308D3ULL;
static inline uint64_t ru(){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static inline double rf(){ return (ru()>>11) * (1.0/9007199254740992.0); }

// minimum pairwise distance of points (O(n^2))
static double minPairDist(const vector<P>& p){
    int n=(int)p.size(); double m2=1e18;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
        double dx=p[i][0]-p[j][0],dy=p[i][1]-p[j][1],dz=p[i][2]-p[j][2];
        double d2=dx*dx+dy*dy+dz*dz; if(d2<m2)m2=d2;
    }
    return n<2?1e18:sqrt(m2);
}
static inline double D_to_r(double D){ return D/(2.0*(1.0+D)); }

// ---- structured lattice initialisations (points in [0,1]^3) ---------------

// axis-aligned a x b x c grid spanning [0,1] (first n points), points at i/(dim-1).
static vector<P> gridABC(int a,int b,int c){
    vector<P> pts; pts.reserve((size_t)a*b*c);
    auto co=[](int i,int dim){ return dim<=1?0.5:(double)i/(dim-1); };
    for(int i=0;i<a&&(int)pts.size()<N;++i)
        for(int j=0;j<b&&(int)pts.size()<N;++j)
            for(int k=0;k<c&&(int)pts.size()<N;++k)
                pts.push_back({co(i,a),co(j,b),co(k,c)});
    return pts;
}

// FCC sublattice: integer pts in [0,L]^3 with even coord-sum; take N nearest the centre,
// then affine-map their bounding box onto [0,1]^3. Denser min-distance than cubic.
static vector<P> fccLattice(int L){
    vector<array<int,3>> ip; ip.reserve(((size_t)(L+1)*(L+1)*(L+1))/2+8);
    for(int i=0;i<=L;++i)for(int j=0;j<=L;++j)for(int k=0;k<=L;++k)
        if(((i+j+k)&1)==0) ip.push_back({i,j,k});
    if((int)ip.size()<N) return {};
    double cc=L/2.0;
    sort(ip.begin(),ip.end(),[&](const array<int,3>&A,const array<int,3>&B){
        double da=(A[0]-cc)*(A[0]-cc)+(A[1]-cc)*(A[1]-cc)+(A[2]-cc)*(A[2]-cc);
        double db=(B[0]-cc)*(B[0]-cc)+(B[1]-cc)*(B[1]-cc)+(B[2]-cc)*(B[2]-cc);
        return da<db; });
    ip.resize(N);
    int mn[3]={1<<30,1<<30,1<<30},mx[3]={-(1<<30),-(1<<30),-(1<<30)};
    for(auto&q:ip)for(int d=0;d<3;++d){mn[d]=min(mn[d],q[d]);mx[d]=max(mx[d],q[d]);}
    // UNIFORM scale by the largest span so the lattice isn't distorted (preserves the
    // min pairwise distance); center the shorter axes. Slack on short axes gives the
    // push-apart room to expand.
    int span=1; for(int d=0;d<3;++d) span=max(span, mx[d]-mn[d]);
    double sc=1.0/span;
    vector<P> pts; pts.reserve(N);
    for(auto&q:ip){
        double off[3];
        for(int d=0;d<3;++d) off[d]=(span-(mx[d]-mn[d]))*0.5;   // center shorter axes
        pts.push_back({ (q[0]-mn[0]+off[0])*sc, (q[1]-mn[1]+off[1])*sc, (q[2]-mn[2]+off[2])*sc });
    }
    return pts;
}

// Push-apart relaxation toward a target min-distance `d`: for a fixed number of
// iterations, separate every pair closer than d (move each endpoint by half the
// deficit) and clamp into the cube. Uses a uniform spatial grid so it scales to 4096.
// Returns the achieved min pairwise distance.
static double relaxToTarget(vector<P>& p, double d, int iters, double hardDeadline){
    int n=(int)p.size();
    double cell=max(d,1e-4);
    int G=max(1,(int)floor(1.0/cell)); if(G>200)G=200;
    double inv=G;
    vector<double> dx_(n),dy_(n),dz_(n);
    vector<vector<int>> cells((size_t)G*G*G);
    for(int it=0; it<iters; ++it){
        if((it&7)==0 && elapsed()>hardDeadline) break;   // never overrun the budget
        for(auto&cv:cells) cv.clear();
        auto gc=[&](double x){ int g=(int)(x*inv); return g<0?0:(g>=G?G-1:g); };
        for(int i=0;i<n;++i) cells[(gc(p[i][0])*G+gc(p[i][1]))*G+gc(p[i][2])].push_back(i);
        for(int i=0;i<n;++i){ dx_[i]=dy_[i]=dz_[i]=0; }
        double d2t=d*d;
        for(int i=0;i<n;++i){
            int gx=gc(p[i][0]),gy=gc(p[i][1]),gz=gc(p[i][2]);
            for(int ax=max(0,gx-1);ax<=min(G-1,gx+1);++ax)
            for(int ay=max(0,gy-1);ay<=min(G-1,gy+1);++ay)
            for(int az=max(0,gz-1);az<=min(G-1,gz+1);++az)
                for(int j: cells[(ax*G+ay)*G+az]) if(j>i){
                    double vx=p[i][0]-p[j][0],vy=p[i][1]-p[j][1],vz=p[i][2]-p[j][2];
                    double q=vx*vx+vy*vy+vz*vz;
                    if(q<d2t){
                        double dist=sqrt(q);
                        double push; double ux,uy,uz;
                        if(dist>1e-12){ push=0.5*(d-dist); ux=vx/dist;uy=vy/dist;uz=vz/dist; }
                        else { push=0.5*d; ux=rf()-0.5;uy=rf()-0.5;uz=rf()-0.5;
                               double l=sqrt(ux*ux+uy*uy+uz*uz)+1e-12; ux/=l;uy/=l;uz/=l; }
                        dx_[i]+=ux*push; dy_[i]+=uy*push; dz_[i]+=uz*push;
                        dx_[j]-=ux*push; dy_[j]-=uy*push; dz_[j]-=uz*push;
                    }
                }
        }
        for(int i=0;i<n;++i){
            p[i][0]=clamp01(p[i][0]+dx_[i]);
            p[i][1]=clamp01(p[i][1]+dy_[i]);
            p[i][2]=clamp01(p[i][2]+dz_[i]);
        }
    }
    return minPairDist(p);
}

// ---- relocation ascent -----------------------------------------------------
// A spatial hash over the current point set, used to answer nearest-neighbour
// queries for arbitrary candidate positions in roughly O(1).
struct NGrid{
    int G; double inv; vector<vector<int>> cells;
    void build(const vector<P>& p, double cell){
        G=max(1,(int)floor(1.0/cell)); if(G>128)G=128; inv=G;
        cells.assign((size_t)G*G*G,{});
        for(int i=0;i<(int)p.size();++i)
            cells[(gc(p[i][0])*G+gc(p[i][1]))*G+gc(p[i][2])].push_back(i);
    }
    inline int gc(double x)const{ int g=(int)(x*inv); return g<0?0:(g>=G?G-1:g); }
};
// squared distance from (x,y,z) to its nearest point in p (excluding index self),
// expanding shells until the shell radius exceeds the best found so far.
static double nnDist2(const NGrid& gr, const vector<P>& p,
                      double x,double y,double z,int self){
    int gx=gr.gc(x),gy=gr.gc(y),gz=gr.gc(z); double best=1e18;
    for(int rad=1;;++rad){
        for(int ax=max(0,gx-rad);ax<=min(gr.G-1,gx+rad);++ax)
        for(int ay=max(0,gy-rad);ay<=min(gr.G-1,gy+rad);++ay)
        for(int az=max(0,gz-rad);az<=min(gr.G-1,gz+rad);++az){
            if(rad>1 && ax>gx-rad&&ax<gx+rad && ay>gy-rad&&ay<gy+rad && az>gz-rad&&az<gz+rad) continue;
            for(int j: gr.cells[(ax*gr.G+ay)*gr.G+az]){
                if(j==self) continue;
                double dx=x-p[j][0],dy=y-p[j][1],dz=z-p[j][2];
                double d=dx*dx+dy*dy+dz*dz; if(d<best)best=d;
            }
        }
        double sm=(double)rad/gr.G; if(sm*sm>best || rad>gr.G) break;
    }
    return best;
}
// Greedy relocation: repeatedly sweep points, moving each to maximise the distance
// to its nearest neighbour (a random-direction line search with shrinking step).
// This breaks the symmetry of lattice fixed points that push-apart cannot escape.
static double reloAscent(vector<P>& p, double D0, int sweeps, double budget){
    int n=(int)p.size();
    double step0=D0*0.6;
    for(int sw=0; sw<sweeps && elapsed()<budget; ++sw){
        NGrid gr; gr.build(p, max(minPairDist(p),1e-4));
        for(int i=0;i<n;++i){
            double x=p[i][0],y=p[i][1],z=p[i][2];
            double cur=sqrt(nnDist2(gr,p,x,y,z,i));
            double step=step0;
            for(int it=0; it<10; ++it){
                double bx=x,by=y,bz=z,bc=cur;
                for(int s=0;s<8;++s){
                    double dx=rf()-0.5,dy=rf()-0.5,dz=rf()-0.5;
                    double l=sqrt(dx*dx+dy*dy+dz*dz)+1e-12; dx/=l;dy/=l;dz/=l;
                    double nx=clamp01(x+dx*step),ny=clamp01(y+dy*step),nz=clamp01(z+dz*step);
                    double c=sqrt(nnDist2(gr,p,nx,ny,nz,i));
                    if(c>bc){ bc=c;bx=nx;by=ny;bz=nz; }
                }
                if(bc>cur+1e-12){ x=bx;y=by;z=bz;cur=bc; } else step*=0.5;
                if(step<D0*0.001) break;
            }
            p[i][0]=x;p[i][1]=y;p[i][2]=z;
        }
    }
    return minPairDist(p);
}

// Maximise min pairwise distance by an expand-and-relax loop: keep the best config;
// repeatedly relax toward a target slightly above the current best (which forces the
// points to spread until the walls stop them), accepting whenever the achieved min
// distance improves. A binary-style step on the target converges without overshooting.
static void spread(vector<P>& p, double budget){
    int n=(int)p.size();
    if(n<2){ if(n==1)p[0]={0.5,0.5,0.5}; return; }
    vector<P> best=p; double bestD=minPairDist(p);
    if(bestD<=0) bestD=0.05;
    // per-relax iteration count: enough for the config to settle at this n
    int iters = n>2000? 60 : (n>800? 90 : (n>200? 140 : 250));
    double lo=bestD, grow=1.10;      // multiplicative target step, shrinks on failure
    while(elapsed()<budget){
        double target = lo*grow;
        vector<P> trial = best;
        double got = relaxToTarget(trial, target, iters, budget);
        if(got > bestD + 1e-9){ bestD=got; best=trial; lo=got; }
        else {
            // overshot: settle toward a target between best and the failed one
            grow = 1.0 + (grow-1.0)*0.6;
            if(grow < 1.0009){       // converged; try one perturbation kick to escape
                vector<P> t2=best; double amp=bestD*0.15;
                for(auto&q:t2){ q[0]=clamp01(q[0]+(rf()-0.5)*amp);
                                q[1]=clamp01(q[1]+(rf()-0.5)*amp);
                                q[2]=clamp01(q[2]+(rf()-0.5)*amp); }
                double g2=relaxToTarget(t2, bestD*1.03, iters, budget);
                if(g2>bestD+1e-9){ bestD=g2; best=t2; lo=g2; }
                grow=1.04;           // reset growth for the next phase
            }
        }
    }
    p=best;
}

int main(){
    T0=chrono::steady_clock::now();
    if(scanf("%d",&N)!=1) return 0;
    if(N<=0) return 0;
    if(N==1){ printf("0.5 0.5 0.5\n"); return 0; }

    // best structured init by min pairwise distance
    vector<P> best; double bestD=-1;
    auto consider=[&](vector<P> pts){
        if((int)pts.size()!=N) return;
        double d=minPairDist(pts);
        if(d>bestD){ bestD=d; best=std::move(pts); }
    };
    int M=1; while((long long)M*M*M<N) ++M;
    consider(gridABC(M,M,M));
    for(int a=max(1,M-2);a<=M+2;++a)for(int b=max(1,M-2);b<=M+2;++b){
        long long need=(long long)a*b; int c=(int)((N+need-1)/need); if(c<1)c=1;
        if((long long)a*b*c>=N) consider(gridABC(a,b,c));
    }
    for(int L=2;L<=48;++L){
        long long approx=((long long)(L+1)*(L+1)*(L+1))/2;
        if(approx<N) continue;
        consider(fccLattice(L));
        if(approx>3LL*N) break;
    }
    if(best.empty()){ best=gridABC(M,M,M); bestD=minPairDist(best); }
    // Keep a plain cubic grid as a relocation-ascent seed: for mid n the grid packs
    // the walls better than a ball-shaped FCC fragment, and relocation ascent (which
    // breaks lattice symmetry) can push it above every symmetric lattice.
    vector<P> gridSeed=gridABC(M,M,M); double gridSeedD=minPairDist(gridSeed);

    // Optimize. The structured lattices are symmetric equilibria of the force field
    // (perfect cubes especially), so a single spread from them can stall. We combine
    // two optimizers with restarts, keeping the best min-pair-distance found by any:
    //   * spread()     – push-apart / expand-and-relax (strong on dense FCC inits)
    //   * reloAscent()  – greedy per-point clearance maximisation (breaks symmetry,
    //                     strong on wall-dominated mid-n where the grid is the base)
    double budget = 0.80;
    vector<P> globalBest = best; double globalD = bestD;
    auto consider2=[&](const vector<P>& c){
        double d=minPairDist(c);
        if(d>globalD){ globalD=d; globalBest=c; }
    };
    // Two guaranteed phases so neither optimizer starves the other. The relocation
    // phase gets a larger share at mid n (where it wins) and a smaller share at dense
    // n (where FCC+spread dominates and relo from a grid seed adds nothing).
    double reloFrac = (N<=64)? 0.45 : (N<=400? 0.6 : (N<=800? 0.35 : 0.18));
    double spreadEnd = elapsed() + (budget-elapsed())*(1.0-reloFrac);

    // ---- Phase A: push-apart / expand-and-relax (spread) with restarts ----
    {
        vector<P> cand = best;
        spread(cand, elapsed() + (spreadEnd-elapsed())*0.55);
        consider2(cand);
    }
    int attempt=0;
    while(elapsed()<spreadEnd-0.01){
        vector<P> cand;
        if(attempt%2==0){                     // perturb the incumbent best, spread
            cand=globalBest;
            double amp=globalD*(0.15+0.35*rf());
            for(auto&q:cand){ q[0]=clamp01(q[0]+(rf()-0.5)*amp);
                              q[1]=clamp01(q[1]+(rf()-0.5)*amp);
                              q[2]=clamp01(q[2]+(rf()-0.5)*amp); }
        } else {                              // fresh random cloud
            cand.resize(N);
            for(auto&q:cand) q={rf(),rf(),rf()};
        }
        double slice = (attempt<4 || N<=64) ? (spreadEnd-elapsed())*0.35 : (spreadEnd-elapsed())*0.7;
        spread(cand, elapsed()+max(0.005,slice));
        consider2(cand);
        ++attempt;
    }

    // ---- Phase B: relocation ascent (hill-climb with small perturbations) ----
    // A dedicated relocation incumbent, seeded from the cubic grid, is repeatedly
    // lightly perturbed and re-ascended in short bursts. This basin-hops out of the
    // shallow local optima that a single long run gets stuck in, and for wall-dominated
    // mid n climbs above every symmetric lattice.
    {
        vector<P> reloBest = gridSeed; double reloD = gridSeedD;
        int burst = N>2000? 5 : (N>800? 10 : 20);
        // initial ascent from the clean grid
        { vector<P> c=reloBest;
          reloAscent(c, reloD, burst, elapsed()+min((budget-elapsed())*0.4,0.2));
          double d=minPairDist(c); if(d>reloD){ reloD=d; reloBest=c; } }
        while(elapsed()<budget-0.01){
            vector<P> c=reloBest;
            double amp=reloD*(0.02+0.15*rf());
            for(auto&q:c){ q[0]=clamp01(q[0]+(rf()-0.5)*amp);
                           q[1]=clamp01(q[1]+(rf()-0.5)*amp);
                           q[2]=clamp01(q[2]+(rf()-0.5)*amp); }
            reloAscent(c, reloD, burst, budget);
            double d=minPairDist(c); if(d>reloD){ reloD=d; reloBest=c; }
        }
        consider2(reloBest);
    }
    best=globalBest;

    // convert to centers: r = D/(2(1+D)); c = r + p*(1-2r)
    double D=minPairDist(best);
    double r=D_to_r(D);
    double sc=1.0-2.0*r;
    for(auto&p:best)
        printf("%.12f %.12f %.12f\n",
               clamp01(r+p[0]*sc), clamp01(r+p[1]*sc), clamp01(r+p[2]*sc));
    return 0;
}
