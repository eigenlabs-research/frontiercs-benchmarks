#include <bits/stdc++.h>
using namespace std;

struct P { int x, y, z; };
struct Candidate { vector<array<double,3>> p; double r = -1; };

/* Center an integer lattice fragment in the cube.  Its shortest lattice
   vector has length d in the supplied integer coordinates. */
static Candidate make_candidate(vector<P> q, double d) {
    int lx=q[0].x, ly=q[0].y, lz=q[0].z, hx=lx, hy=ly, hz=lz;
    for (auto a:q) {
        lx=min(lx,a.x); ly=min(ly,a.y); lz=min(lz,a.z);
        hx=max(hx,a.x); hy=max(hy,a.y); hz=max(hz,a.z);
    }
    int sx=hx-lx, sy=hy-ly, sz=hz-lz;
    int widest=max(sx,max(sy,sz));
    double h=1.0/(widest+d);
    Candidate c;
    c.r=h*d*.5;
    c.p.reserve(q.size());
    double ox=(1.0-h*sx)*.5, oy=(1.0-h*sy)*.5, oz=(1.0-h*sz)*.5;
    for (auto a:q) c.p.push_back({ox+h*(a.x-lx),oy+h*(a.y-ly),oz+h*(a.z-lz)});
    return c;
}

static void central_order(vector<P>& v, int cx, int cy, int cz) {
    // Twice coordinates are used, so this also handles half-integral centers.
    sort(v.begin(),v.end(),[&](const P&a,const P&b) {
        long long ax=2LL*a.x-cx, ay=2LL*a.y-cy, az=2LL*a.z-cz;
        long long bx=2LL*b.x-cx, by=2LL*b.y-cy, bz=2LL*b.z-cz;
        long long da=ax*ax+ay*ay+az*az, db=bx*bx+by*by+bz*bz;
        if(da!=db) return da<db;
        if(a.x!=b.x) return a.x<b.x;
        if(a.y!=b.y) return a.y<b.y;
        return a.z<b.z;
    });
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    Candidate best;
    auto consider=[&](vector<P> v, double d) {
        v.resize(n);
        Candidate c=make_candidate(move(v),d);
        if(c.r>best.r) best=move(c);
    };

    // Cubic grid: this retains the elementary construction as a safe fallback.
    int g=1; while(1LL*g*g*g<n) ++g;
    vector<P> v;
    for(int x=0;x<g;x++) for(int y=0;y<g;y++) for(int z=0;z<g;z++) v.push_back({x,y,z});
    central_order(v,g-1,g-1,g-1);
    consider(move(v),1.0);

    // FCC: integer triples of one parity have nearest-vector length sqrt(2).
    int m=0;
    while(true) {
        long long cnt=0;
        for(int x=0;x<=m;x++) for(int y=0;y<=m;y++) for(int z=0;z<=m;z++) if(((x+y+z)&1)==0) ++cnt;
        if(cnt>=n) break;
        ++m;
    }
    v.clear();
    for(int x=0;x<=m;x++) for(int y=0;y<=m;y++) for(int z=0;z<=m;z++) if(((x+y+z)&1)==0) v.push_back({x,y,z});
    central_order(v,m,m,m);
    consider(move(v),sqrt(2.0));

    // BCC fills the small-n and boundary transition cases well.  In units of
    // half a cubic cell it is the even sublattice plus the odd sublattice.
    int t=0; while(1LL*(t+1)*(t+1)*(t+1)+1LL*t*t*t<n) ++t;
    v.clear();
    for(int x=0;x<=2*t;x+=2) for(int y=0;y<=2*t;y+=2) for(int z=0;z<=2*t;z+=2) v.push_back({x,y,z});
    for(int x=1;x<2*t;x+=2) for(int y=1;y<2*t;y+=2) for(int z=1;z<2*t;z+=2) v.push_back({x,y,z});
    central_order(v,2*t,2*t,2*t);
    consider(move(v),sqrt(3.0));

    cout<<setprecision(17);
    for(auto a:best.p) cout<<a[0]<<' '<<a[1]<<' '<<a[2]<<'\n';
}
