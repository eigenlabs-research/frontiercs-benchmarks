#include <bits/stdc++.h>
using namespace std;

struct P { double x,y,z; };

static double radiusOf(const vector<P>& a) {
    double r=0.5;
    int n=(int)a.size();
    for (const P& p:a) r=min(r,min({p.x,1-p.x,p.y,1-p.y,p.z,1-p.z}));
    for(int i=0;i<n;i++) for(int j=0;j<i;j++) {
        double x=a[i].x-a[j].x, y=a[i].y-a[j].y, z=a[i].z-a[j].z;
        r=min(r,0.5*sqrt(x*x+y*y+z*z));
    }
    return r;
}

static vector<P> cubic(int n) {
    int m=1; while(1LL*m*m*m<n) ++m;
    int ba=m,bb=m,bc=m, best=INT_MAX;
    // A rectangular grid is useful for counts just below a cube as well.
    for(int a=1;a<=m;a++) for(int b=1;b<=m;b++) {
        int c=(n+a*b-1)/(a*b);
        if(c<=m && a*b*c<best) ba=a,bb=b,bc=c,best=a*b*c;
    }
    vector<P> v; v.reserve(n);
    for(int i=0;i<ba && (int)v.size()<n;i++)
      for(int j=0;j<bb && (int)v.size()<n;j++)
        for(int k=0;k<bc && (int)v.size()<n;k++)
          v.push_back({(i+.5)/ba,(j+.5)/bb,(k+.5)/bc});
    return v;
}

// D3 is the face-centred cubic lattice: integer triples of even coordinate sum.
// It has nearest-neighbour distance sqrt(2), and is asymptotically the densest
// possible equal-sphere lattice packing.
static vector<P> fcc(int n) {
    int w=0;
    for(;;w++) {
        long long q=1LL*(w+1)*(w+1)*(w+1);
        long long cnt=(q+1)/2; // one of the two parity classes
        if(cnt>=n) break;
    }
    vector<array<int,3>> q; q.reserve(n);
    for(int x=0;x<=w && (int)q.size()<n;x++)
      for(int y=0;y<=w && (int)q.size()<n;y++)
        for(int z=0;z<=w && (int)q.size()<n;z++)
          if(((x+y+z)&1)==0) q.push_back({x,y,z});
    int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for(auto p:q) for(int d=0;d<3;d++) lo[d]=min(lo[d],p[d]),hi[d]=max(hi[d],p[d]);
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double s=1.0/(span+sqrt(2.0));
    vector<P> v; v.reserve(n);
    for(auto p:q) v.push_back({.5+(p[0]-(lo[0]+hi[0])*.5)*s,
                                .5+(p[1]-(lo[1]+hi[1])*.5)*s,
                                .5+(p[2]-(lo[2]+hi[2])*.5)*s});
    return v;
}

// Alternating projection onto the non-overlap constraints.  This is only used
// for small instances, where boundary effects dominate lattice asymptotics.
static bool relax(vector<P>& p, double r) {
    int n=p.size();
    double side=1-2*r;
    if(side<0) return false;
    for(P& a:p) {
        a.x=min(1-r,max(r,a.x)); a.y=min(1-r,max(r,a.y)); a.z=min(1-r,max(r,a.z));
    }
    const double need=2*r;
    for(int it=0;it<900;it++) {
        double worst=0;
        // Reverse/cyclic order changes the otherwise rather strong lattice bias.
        for(int ii=0;ii<n;ii++) for(int jj=0;jj<ii;jj++) {
            int i=(ii+it)%n, j=(jj+it)%n;
            double dx=p[i].x-p[j].x,dy=p[i].y-p[j].y,dz=p[i].z-p[j].z;
            double d=sqrt(dx*dx+dy*dy+dz*dz);
            if(d<need) {
                double ux,uy,uz;
                if(d<1e-14) { // deterministic nonzero direction
                    unsigned h=(unsigned)(i*1103515245u+j*12345u);
                    ux=(h&1)?1:-1; uy=(h&2)?1:-1; uz=(h&4)?1:-1;
                    double t=1/sqrt(3.0); ux*=t;uy*=t;uz*=t; d=0;
                } else ux=dx/d,uy=dy/d,uz=dz/d;
                double move=(need-d)*.505;
                p[i].x+=ux*move; p[i].y+=uy*move; p[i].z+=uz*move;
                p[j].x-=ux*move; p[j].y-=uy*move; p[j].z-=uz*move;
                worst=max(worst,need-d);
            }
        }
        for(P& a:p) {
            a.x=min(1-r,max(r,a.x)); a.y=min(1-r,max(r,a.y)); a.z=min(1-r,max(r,a.z));
        }
        if(worst<2e-10) return radiusOf(p)>=r-2e-9;
    }
    return radiusOf(p)>=r-2e-8;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<P> best=cubic(n), t=fcc(n);
    if(radiusOf(t)>radiusOf(best)) best=t;

    // The two-ball optimum is the body diagonal, not a lattice cell.
    if(n==2) {
        double r=sqrt(3.0)/(2*(1+sqrt(3.0)));
        best={{r,r,r},{1-r,1-r,1-r}};
    }

    // Try to remove finite-box boundary loss.  A candidate is retained only
    // after its actual geometric radius has been recomputed.
    if(n>=3 && n<=28) {
        double cur=radiusOf(best);
        const double pi=acos(-1.0), delta=pi/sqrt(18.0);
        double hi=min(.499999, pow(delta*3.0/(4*pi*n),1.0/3.0));
        for(int step=0;step<13;step++) {
            double goal=(cur+hi)*.5;
            vector<P> q=best;
            // Contract about the cube centre to make room for the new margin.
            double old=max(1e-12,cur);
            double factor=(1-2*goal)/(1-2*old);
            for(P& a:q) { a.x=.5+(a.x-.5)*factor; a.y=.5+(a.y-.5)*factor; a.z=.5+(a.z-.5)*factor; }
            if(relax(q,goal)) { double rr=radiusOf(q); if(rr>radiusOf(best)) best=q; cur=radiusOf(best); }
            else hi=goal;
        }
    }
    cout<<setprecision(17);
    for(const P& p:best) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
}
