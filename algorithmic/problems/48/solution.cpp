#include <bits/stdc++.h>
using namespace std;
using P = array<double,3>;

// The ordinary grid is retained as a certified floor.  The other candidate is
// an FCC lattice: integer triples with even coordinate sum have nearest
// neighbour distance sqrt(2).
static vector<P> grid(int n, double &rad) {
    int m=1; while(1LL*m*m*m<n) ++m;
    int a=(n+m*m-1)/(m*m), b=m, c=m;
    vector<P> q; q.reserve(n);
    for(int i=0;i<a && (int)q.size()<n;i++)
      for(int j=0;j<b && (int)q.size()<n;j++)
        for(int k=0;k<c && (int)q.size()<n;k++)
          q.push_back({(i+.5)/a,(j+.5)/b,(k+.5)/c});
    rad=0.5/m;
    return q;
}

static vector<P> fccBox(int n, double &rad) {
    int w=0;
    for(;;++w) {
        long long cnt=0;
        for(int x=0;x<=w;x++) for(int y=0;y<=w;y++) for(int z=0;z<=w;z++)
          if(((x+y+z)&1)==0) ++cnt;
        if(cnt>=n) break;
    }
    vector<array<int,3>> v; v.reserve(n);
    // Interleave slabs rather than take a face-first prefix.  This makes a
    // partly occupied box have a substantially smaller bounding box on many n.
    vector<array<int,3>> all;
    all.reserve((w+1)*(w+1)*(w+1)/2+2);
    for(int x=0;x<=w;x++) for(int y=0;y<=w;y++) for(int z=0;z<=w;z++)
        if(((x+y+z)&1)==0) all.push_back({x,y,z});
    double mid=w*.5;
    stable_sort(all.begin(),all.end(),[&](const auto&A,const auto&B){
        double da=(A[0]-mid)*(A[0]-mid)+(A[1]-mid)*(A[1]-mid)+(A[2]-mid)*(A[2]-mid);
        double db=(B[0]-mid)*(B[0]-mid)+(B[1]-mid)*(B[1]-mid)+(B[2]-mid)*(B[2]-mid);
        if(da!=db) return da<db;
        // deterministic shell order, deliberately not a lexicographic slab
        return A<B;
    });
    all.resize(n);
    int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for(auto a:all) for(int d=0;d<3;d++) lo[d]=min(lo[d],a[d]),hi[d]=max(hi[d],a[d]);
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    const double d=sqrt(2.0);
    double scale=1.0/(span+d);
    rad=d*scale*.5;
    vector<P> q; q.reserve(n);
    for(auto a:all) q.push_back({(a[0]-lo[0])*scale+rad,(a[1]-lo[1])*scale+rad,(a[2]-lo[2])*scale+rad});
    return q;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    double rg, rf; vector<P> ans=grid(n,rg);
    // The regular tetrahedron is a known useful exception to a grid.
    if(n==4){
        double r=sqrt(2.0)/(2.0*(1.0+sqrt(2.0)));
        ans={{P{r,r,r},P{r,1-r,1-r},P{1-r,r,1-r},P{1-r,1-r,r}}}; rg=r;
    }
    vector<P> f=fccBox(n,rf);
    if(rf>rg) ans.swap(f);
    cout<<setprecision(17);
    for(auto p:ans) cout<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
}
