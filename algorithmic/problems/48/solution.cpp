#include <bits/stdc++.h>
using namespace std;
using P = array<double,3>;

// Put a point set in its largest possible uniformly scaled copy of the cube.
// Keeping half a nearest-neighbour distance at every face makes the two
// possible limiting constraints equal whenever the lattice has a contact.
static vector<P> normalize(vector<P> a) {
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for (auto &p:a) for(int d=0;d<3;d++) lo[d]=min(lo[d],p[d]),hi[d]=max(hi[d],p[d]);
    double md=1e100;
    for (int i=0;i<(int)a.size();i++) for(int j=0;j<i;j++) {
        double x=a[i][0]-a[j][0], y=a[i][1]-a[j][1], z=a[i][2]-a[j][2];
        md=min(md,sqrt(x*x+y*y+z*z));
    }
    double w=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double s=1.0/(w+md), pad=md*s*.5;
    for(auto &p:a) for(int d=0;d<3;d++) p[d]=(p[d]-lo[d])*s+pad;
    return a;
}
static double value(const vector<P>& a) {
    double r=1;
    for(auto p:a) for(int d=0;d<3;d++) r=min(r,min(p[d],1-p[d]));
    for(int i=0;i<(int)a.size();i++) for(int j=0;j<i;j++) {
        double x=a[i][0]-a[j][0],y=a[i][1]-a[j][1],z=a[i][2]-a[j][2];
        r=min(r,.5*sqrt(x*x+y*y+z*z));
    }
    return r;
}

static vector<P> bcc(int n) {
    int q=(int)ceil(cbrt((double)n))+3;
    struct Q { int x,y,z; long long key; };
    vector<Q> v; v.reserve(2*(2*q+1)*(2*q+1)*(2*q+1));
    // Coordinates with all parities equal are a BCC lattice (in half units).
    for(int x=-q;x<=q;x++) for(int y=-q;y<=q;y++) for(int z=-q;z<=q;z++)
        if ((x&1)==(y&1) && (y&1)==(z&1)) {
            long long k=1LL*x*x+y*y+z*z;
            v.push_back({x,y,z,k});
        }
    sort(v.begin(),v.end(),[](const Q&a,const Q&b) {
        if(a.key!=b.key) return a.key<b.key;
        if(a.x!=b.x) return a.x<b.x; if(a.y!=b.y) return a.y<b.y; return a.z<b.z;
    });
    vector<P> a; a.reserve(n);
    for(int i=0;i<n;i++) a.push_back({(double)v[i].x,(double)v[i].y,(double)v[i].z});
    return normalize(a);
}
static vector<P> fcc(int n) {
    int q=(int)ceil(cbrt((double)n))+4;
    struct Q { int x,y,z; long long key; };
    vector<Q> v;
    for(int x=-q;x<=q;x++) for(int y=-q;y<=q;y++) for(int z=-q;z<=q;z++)
        if(((x+y+z)&1)==0) v.push_back({x,y,z,1LL*x*x+y*y+z*z});
    sort(v.begin(),v.end(),[](const Q&a,const Q&b) {
        if(a.key!=b.key) return a.key<b.key;
        if(a.x!=b.x) return a.x<b.x; if(a.y!=b.y) return a.y<b.y; return a.z<b.z;
    });
    vector<P> a; a.reserve(n);
    for(int i=0;i<n;i++) a.push_back({(double)v[i].x,(double)v[i].y,(double)v[i].z});
    return normalize(a);
}
static vector<P> grid(int n) {
    int bestA=1,bestB=1,bestC=n, bestMax=n;
    // Enumerating the two smaller dimensions finds the balanced rectangular grid.
    for(int a=1;a*a*a<=n*2;a++) for(int b=a;b<=n/a;b++) {
        int c=(n+a*b-1)/(a*b); if(c<b) c=b;
        if(a*b*c>=n && c<bestMax) bestA=a,bestB=b,bestC=c,bestMax=c;
    }
    vector<P> p; p.reserve(n);
    for(int x=0;x<bestA && (int)p.size()<n;x++) for(int y=0;y<bestB && (int)p.size()<n;y++) for(int z=0;z<bestC && (int)p.size()<n;z++) p.push_back({(double)x,(double)y,(double)z});
    return normalize(p);
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<P> ans=grid(n), t=bcc(n);
    double best=value(ans), v=value(t);
    if(v>best){best=v; ans.swap(t);}
    t=fcc(n); v=value(t);
    if(v>best) ans.swap(t);
    cout<<setprecision(17);
    for(auto p:ans) cout<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
}
