#include <bits/stdc++.h>
using namespace std;
struct Point { long long x,y; };

static uint64_t hilbert(uint32_t x, uint32_t y) {
    uint64_t d=0;
    for (uint32_t s=1u<<19; s; s>>=1) {
        uint32_t rx=(x&s)!=0, ry=(y&s)!=0;
        d += (uint64_t)s*s*((3*rx)^ry);
        if (!ry) { if (rx) { x=s-1-x; y=s-1-y; } swap(x,y); }
    }
    return d;
}
static uint64_t morton(uint32_t x, uint32_t y) {
    uint64_t z=0;
    for(int b=0;b<20;b++) { z|=(uint64_t)((x>>b)&1)<<(2*b); z|=(uint64_t)((y>>b)&1)<<(2*b+1); }
    return z;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<Point> p(n);
    long long minx=LLONG_MAX,maxx=LLONG_MIN,miny=LLONG_MAX,maxy=LLONG_MIN;
    for(auto &a:p){ cin>>a.x>>a.y; minx=min(minx,a.x); maxx=max(maxx,a.x); miny=min(miny,a.y); maxy=max(maxy,a.y); }
    vector<char> prime(n,true); if(n>0) prime[0]=false; if(n>1) prime[1]=false;
    for(int i=2;i*(long long)i<n;i++) if(prime[i]) for(int j=i*i;j<n;j+=i) prime[j]=false;
    auto dist=[&](int a,int b)->double { return hypot((double)p[a].x-p[b].x,(double)p[a].y-p[b].y); };
    auto weight=[&](int source,int step)->double { return (step%10==0 && !prime[source])?1.1:1.0; };
    auto cost=[&](const vector<int>& q)->double { double s=0; for(int k=0;k<n;k++) s+=weight(q[k],k+1)*dist(q[k],q[k+1]); return s; };
    auto routeFromCycle=[&](vector<int> a, bool rev)->vector<int> {
        int at=find(a.begin(),a.end(),0)-a.begin(); vector<int> q; q.reserve(n+1); q.push_back(0);
        if(!rev) for(int z=1;z<n;z++) q.push_back(a[(at+z)%n]);
        else for(int z=1;z<n;z++) q.push_back(a[(at-z+n)%n]);
        q.push_back(0); return q;
    };
    vector<vector<int>> candidates;
    vector<int> xord(n); iota(xord.begin(),xord.end(),0);
    candidates.push_back(routeFromCycle(xord,false)); candidates.push_back(routeFromCycle(xord,true));
    const uint64_t S=(1u<<20)-1;
    long double dx=maxx-minx, dy=maxy-miny;
    vector<pair<uint64_t,int>> keys; keys.reserve(n);
    for(int mode=0;mode<2;mode++) {
        keys.clear();
        for(int i=0;i<n;i++) {
            uint32_t X= dx==0?0:(uint32_t)((long double)(p[i].x-minx)*S/dx);
            uint32_t Y= dy==0?0:(uint32_t)((long double)(p[i].y-miny)*S/dy);
            keys.push_back({mode?morton(X,Y):hilbert(X,Y),i});
        }
        sort(keys.begin(),keys.end()); vector<int> a; a.reserve(n); for(auto z:keys)a.push_back(z.second);
        candidates.push_back(routeFromCycle(a,false)); candidates.push_back(routeFromCycle(a,true));
    }
    // A polar cycle supplies a different basin for rings and spoke-like point sets.
    long double cx=((long double)minx+maxx)/2, cy=((long double)miny+maxy)/2;
    vector<pair<pair<long double,long double>,int>> polar; polar.reserve(n);
    for(int i=0;i<n;i++) {
        long double ux=p[i].x-cx, uy=p[i].y-cy;
        polar.push_back({{atan2l(uy,ux),ux*ux+uy*uy},i});
    }
    sort(polar.begin(),polar.end()); vector<int> pa; pa.reserve(n);
    for(auto z:polar) pa.push_back(z.second);
    candidates.push_back(routeFromCycle(pa,false)); candidates.push_back(routeFromCycle(pa,true));
    // A geometry-scaled vertical serpentine order is useful for layered point sets.
    int bands=max(1,(int)sqrt((long double)n * max((long double)1.0,dx) / max((long double)1.0,dy)));
    bands=min(bands,max(1,n));
    vector<pair<pair<int,long long>,int>> stripe; stripe.reserve(n);
    for(int i=0;i<n;i++) {
        int b=dx==0?0:(int)(((long double)(p[i].x-minx)*bands)/(dx+1)); if(b>=bands)b=bands-1;
        long long yy=(b&1)?-p[i].y:p[i].y; stripe.push_back({{b,yy},i});
    }
    sort(stripe.begin(),stripe.end()); vector<int> sa; for(auto z:stripe)sa.push_back(z.second);
    candidates.push_back(routeFromCycle(sa,false)); candidates.push_back(routeFromCycle(sa,true));
    vector<int> q=candidates[0]; double best=cost(q);
    for(size_t z=1;z<candidates.size();z++){ double v=cost(candidates[z]); if(v<best){best=v;q.swap(candidates[z]);} }
    auto edge=[&](const vector<int>& a,int k)->double { return weight(a[k],k+1)*dist(a[k],a[k+1]); };
    // Exact deltas: all edges whose endpoint or weighted source can change are recomputed.
    for(int pass=0;pass<2;pass++) {
        for(int len=3;len<=8;len++) for(int l=1;l+len-1<n;l++) {
            int r=l+len-1; double before=0,after=0;
            for(int k=l-1;k<=r;k++) before+=edge(q,k);
            reverse(q.begin()+l,q.begin()+r+1);
            for(int k=l-1;k<=r;k++) after+=edge(q,k);
            if(after+1e-7<before) best+=after-before;
            else reverse(q.begin()+l,q.begin()+r+1);
        }
        for(int i=1;i+1<n;i++) {
            double before=0,after=0; for(int k=i-1;k<=i+1;k++) before+=edge(q,k);
            swap(q[i],q[i+1]); for(int k=i-1;k<=i+1;k++) after+=edge(q,k);
            if(after+1e-7<before) best+=after-before; else swap(q[i],q[i+1]);
        }
    }
    cout<<n+1<<'\n'; for(int v:q) cout<<v<<'\n';
}
