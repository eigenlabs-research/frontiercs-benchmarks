#include <bits/stdc++.h>
using namespace std;

struct Point { long long x,y; };
static vector<Point> a;
static vector<char> primeId;

static uint64_t hilbert(uint32_t x, uint32_t y) {
    uint64_t d=0;
    for (uint32_t s=1u<<30; s; s>>=1) {
        uint32_t rx=(x&s)!=0, ry=(y&s)!=0;
        d += (uint64_t)s*s*((3*rx)^ry);
        if (!ry) {
            if (rx) { x=s-1-x; y=s-1-y; }
            swap(x,y);
        }
    }
    return d;
}
static inline double edge(int u,int v,int t) {
    double dx=(double)a[u].x-a[v].x, dy=(double)a[u].y-a[v].y;
    double z=sqrt(dx*dx+dy*dy);
    return (t%10==0 && !primeId[u]) ? z*1.1 : z;
}
static double cost(const vector<int>& p) {
    double z=0; for(int i=1;i<(int)p.size();++i) z+=edge(p[i-1],p[i],i); return z;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    a.resize(n); for(auto &q:a) cin>>q.x>>q.y;
    primeId.assign(n,true); if(n>0) primeId[0]=false; if(n>1) primeId[1]=false;
    for(int i=2;i*(long long)i<n;i++) if(primeId[i]) for(int j=i*i;j<n;j+=i) primeId[j]=false;

    vector<vector<int>> candidates;
    vector<int> ids(n-1); iota(ids.begin(),ids.end(),1);
    // Both monotone directions are cheap safeguards for instances where x order is already ideal.
    for(int rev=0;rev<2;rev++) { auto v=ids; if(rev) reverse(v.begin(),v.end()); candidates.push_back(v); }
    const uint32_t LIM=(1u<<31)-1;
    for(int mode=0;mode<4;mode++) {
        vector<pair<uint64_t,int>> z; z.reserve(n-1);
        for(int id=1;id<n;id++) {
            uint32_t x=(uint32_t)(a[id].x+1000000000LL), y=(uint32_t)(a[id].y+1000000000LL);
            if(mode&1) x=LIM-x;
            if(mode&2) y=LIM-y;
            z.push_back({hilbert(x,y),id});
        }
        sort(z.begin(),z.end());
        vector<int> v; v.reserve(n-1); for(auto e:z) v.push_back(e.second);
        candidates.push_back(v); reverse(v.begin(),v.end()); candidates.push_back(v);
    }
    vector<int> p; double best=1e300;
    for(auto &v:candidates) { vector<int> q; q.reserve(n+1); q.push_back(0); q.insert(q.end(),v.begin(),v.end()); q.push_back(0); double c=cost(q); if(c<best) best=c,p.swap(q); }

    // Small exact 2-opt moves retain the spatial ordering but repair short crossings and
    // account for every position-dependent carrot multiplier in the tested interval.
    for(int pass=0;pass<2;pass++) {
        for(int i=1;i<n-1;i++) for(int len=1;len<=8 && i+len<n;len++) {
            int j=i+len;
            double oldc=0,newc=0;
            for(int k=i;k<=j+1;k++) {
                oldc+=edge(p[k-1],p[k],k);
                int u = (k-1<i || k-1>j) ? p[k-1] : p[i+j-(k-1)];
                int v = (k<i || k>j) ? p[k] : p[i+j-k];
                newc+=edge(u,v,k);
            }
            if(newc+1e-7<oldc) reverse(p.begin()+i,p.begin()+j+1);
        }
    }
    cout<<n+1<<'\n'; for(int x:p) cout<<x<<'\n';
}
