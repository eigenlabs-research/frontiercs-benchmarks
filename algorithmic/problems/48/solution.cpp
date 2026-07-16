#include <bits/stdc++.h>
using namespace std;

struct P { double x,y,z; };
struct Candidate { vector<P> p; double r=-1; };

static Candidate gridCandidate(int n) {
    int m=1; while (1LL*m*m*m<n) ++m;
    Candidate c; c.p.reserve(n); c.r=0.5/m;
    for(int i=0;i<m && (int)c.p.size()<n;i++)
      for(int j=0;j<m && (int)c.p.size()<n;j++)
        for(int k=0;k<m && (int)c.p.size()<n;k++)
          c.p.push_back({(i+.5)/m,(j+.5)/m,(k+.5)/m});
    return c;
}

/* Points are integer coordinates.  FCC is i+j+k even; BCC is the
   union of the even and odd sublattices.  A compact radial prefix is
   subsequently fitted, after rotation, into the cube. */
static vector<array<int,3>> compactSites(int n, bool bcc) {
    int L=(int)ceil(cbrt((double)n))+4;
    vector<array<int,3>> a;
    for(int x=-L;x<=L;x++) for(int y=-L;y<=L;y++) for(int z=-L;z<=L;z++) {
        int parity=(x+y+z)&1; if(parity<0) parity+=2;
        bool take = bcc ? ((x&1)==(y&1) && (y&1)==(z&1)) : parity==0;
        if(take) a.push_back({x,y,z});
    }
    sort(a.begin(),a.end(),[](const auto&A,const auto&B){
        int da=A[0]*A[0]+A[1]*A[1]+A[2]*A[2];
        int db=B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
        if(da!=db) return da<db;
        return A<B;
    });
    a.resize(n);
    return a;
}

static P rotateP(array<int,3> q, double az, double ay) {
    double cz=cos(az), sz=sin(az), cy=cos(ay), sy=sin(ay);
    double x=cz*q[0]-sz*q[1], y=sz*q[0]+cz*q[1], z=q[2];
    return {cy*x+sy*z,y,-sy*x+cy*z};
}

static Candidate latticeCandidate(const vector<array<int,3>>& sites, double near, double az, double ay) {
    vector<P> q; q.reserve(sites.size());
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(auto v:sites) {
        P p=rotateP(v,az,ay); q.push_back(p);
        lo[0]=min(lo[0],p.x); lo[1]=min(lo[1],p.y); lo[2]=min(lo[2],p.z);
        hi[0]=max(hi[0],p.x); hi[1]=max(hi[1],p.y); hi[2]=max(hi[2],p.z);
    }
    double span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double scale=1.0/(span+near), rad=near*scale*.5;
    Candidate c; c.r=rad; c.p.reserve(q.size());
    for(P p:q) c.p.push_back({rad+(p.x-lo[0])*scale,rad+(p.y-lo[1])*scale,rad+(p.z-lo[2])*scale});
    return c;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    Candidate best=gridCandidate(n);
    vector<pair<double,double>> angles={{0,0},{M_PI/4,0},{0,M_PI/4},{M_PI/4,M_PI/4},{M_PI/6,M_PI/6},{M_PI/4,M_PI/6}};
    // FCC preserves the incumbent's high-density asymptotic packing.
    auto fcc=compactSites(n,false);
    for(auto [a,b]:angles) { Candidate c=latticeCandidate(fcc,sqrt(2.0),a,b); if(c.r>best.r) best=move(c); }
    // BCC is deliberately the additional finite-cluster representation.
    auto bcc=compactSites(n,true);
    for(auto [a,b]:angles) { Candidate c=latticeCandidate(bcc,sqrt(3.0),a,b); if(c.r>best.r) best=move(c); }
    cout<<setprecision(17);
    for(auto p:best.p) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
}
