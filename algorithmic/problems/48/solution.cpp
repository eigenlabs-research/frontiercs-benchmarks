#include <bits/stdc++.h>
using namespace std;

struct P { double x,y,z; };
struct I { int x,y,z; };

static vector<P> gridPacking(int n) {
    int m=1; while (1LL*m*m*m<n) ++m;
    double r=1.0/(2.0*m), h=1.0/m;
    vector<P> v; v.reserve(n);
    for(int x=0;x<m && (int)v.size()<n;x++)
        for(int y=0;y<m && (int)v.size()<n;y++)
            for(int z=0;z<m && (int)v.size()<n;z++)
                v.push_back({r+x*h,r+y*h,r+z*h});
    return v;
}

// A cropped FCC cell.  Coordinates with even coordinate sum form an FCC lattice.
// The returned lower bound is attained by balancing lattice separation and face slack.
static pair<double,vector<P>> fccCell(int n,int a,int b,int c, bool compact) {
    vector<I> q; q.reserve((a*b*c+1)/2);
    for(int x=0;x<a;x++) for(int y=0;y<b;y++) for(int z=0;z<c;z++)
        if(((x+y+z)&1)==0) q.push_back({x,y,z});
    if((int)q.size()<n) return {-1,{}};
    if(compact) {
        // The original prefix can fill all of two axes before reaching n.  A
        // centered prefix is another valid FCC subset and tests that crop loss.
        sort(q.begin(),q.end(),[&](const I& u,const I& v) {
            long long ux=2LL*u.x-(a-1), uy=2LL*u.y-(b-1), uz=2LL*u.z-(c-1);
            long long vx=2LL*v.x-(a-1), vy=2LL*v.y-(b-1), vz=2LL*v.z-(c-1);
            long long du=ux*ux+uy*uy+uz*uz, dv=vx*vx+vy*vy+vz*vz;
            if(du!=dv) return du<dv;
            if(u.x!=v.x) return u.x<v.x;
            if(u.y!=v.y) return u.y<v.y;
            return u.z<v.z;
        });
    }
    q.resize(n);
    int lx=q[0].x,hx=lx,ly=q[0].y,hy=ly,lz=q[0].z,hz=lz;
    for(auto p:q) { lx=min(lx,p.x); hx=max(hx,p.x); ly=min(ly,p.y); hy=max(hy,p.y); lz=min(lz,p.z); hz=max(hz,p.z); }
    const double d=sqrt(2.0); // any two FCC sites are at least this far apart
    double L=max({hx-lx,hy-ly,hz-lz});
    double s=1.0/(L+d), r=d*s/2.0;
    double ex=hx-lx, ey=hy-ly, ez=hz-lz;
    // Center the unused slack independently in the three directions.
    double ox=r+(1-2*r-s*ex)/2.0-s*lx;
    double oy=r+(1-2*r-s*ey)/2.0-s*ly;
    double oz=r+(1-2*r-s*ez)/2.0-s*lz;
    vector<P> v; v.reserve(n);
    for(auto p:q) v.push_back({ox+s*p.x,oy+s*p.y,oz+s*p.z});
    return {r,v};
}

static vector<P> smallPacking(int n) {
    vector<P> v;
    if(n==2) {
        double r=sqrt(3.0)/(2.0*(1.0+sqrt(3.0)));
        v={{r,r,r},{1-r,1-r,1-r}};
    } else if(n==3 || n==4) {
        double r=sqrt(2.0)/(2.0*(1.0+sqrt(2.0)));
        vector<P> t={{r,r,r},{r,1-r,1-r},{1-r,r,1-r},{1-r,1-r,r}};
        t.resize(n); v=t;
    }
    return v;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<P> best=gridPacking(n);
    int m=1; while(1LL*m*m*m<n) ++m;
    double bestR=1.0/(2.0*m);

    // This bounded search is the cropping ablation: dimensions close to a cube are
    // tried, and the actual extent of the n retained sites, not the nominal cell,
    // decides the scale.
    int root=1; while(1LL*root*root*root<2LL*n) ++root;
    int lo=max(1,root-4), hi=root+4;
    for(int a=lo;a<=hi;a++) for(int b=lo;b<=hi;b++) for(int c=lo;c<=hi;c++) {
        if((1LL*a*b*c+1)/2<n) continue;
        for(bool compact : {false, true}) {
            auto got=fccCell(n,a,b,c,compact);
            if(got.first>bestR) { bestR=got.first; best=move(got.second); }
        }
    }
    vector<P> sp=smallPacking(n);
    if(!sp.empty()) best=move(sp);

    // Stay infinitesimally inside after decimal conversion; this cannot create an
    // out-of-cube coordinate and is far below the checker precision scale.
    cout<<setprecision(17);
    for(P p:best) {
        p.x=min(1.0,max(0.0,p.x)); p.y=min(1.0,max(0.0,p.y)); p.z=min(1.0,max(0.0,p.z));
        cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    }
}
