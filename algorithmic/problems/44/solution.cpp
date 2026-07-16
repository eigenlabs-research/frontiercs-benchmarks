#include <bits/stdc++.h>
using namespace std;

struct Point { long long x,y; };
static vector<Point> p;
static vector<char> primeCity;
static int n;

static inline double edgeCost(const vector<int>& r, int k) {
    int a=r[k], b=r[(k+1)%n];
    double d=hypot((double)p[a].x-p[b].x, (double)p[a].y-p[b].y);
    return ((k+1)%10==0 && !primeCity[a]) ? d*1.1 : d;
}
static double score(const vector<int>& r) {
    double z=0;
    for(int k=0;k<n;k++) z+=edgeCost(r,k);
    return z;
}
static vector<int> cutAtZero(const vector<int>& cyc) {
    int q=0; while(cyc[q]!=0) ++q;
    vector<int> r; r.reserve(n);
    for(int i=0;i<n;i++) r.push_back(cyc[(q+i)%n]);
    return r;
}
// Exact local cleanup.  Reversals are deliberately short: their complete,
// position-dependent carrot cost can then be checked rather than approximated.
static void improve(vector<int>& r) {
    for(int pass=0;pass<2;pass++) {
        bool any=false;
        // Adjacent exchanges are the length-two case, handled at the same time.
        for(int i=1;i<n-1;i++) {
            for(int len=2;len<=8 && i+len-1<n;len++) {
                int j=i+len-1;
                // Only edges i-1 through j are changed by reversing [i,j].
                double old=0;
                for(int k=i-1;k<=j;k++) old+=edgeCost(r,k);
                reverse(r.begin()+i,r.begin()+j+1);
                double now=0;
                for(int k=i-1;k<=j;k++) now+=edgeCost(r,k);
                if(now + 1e-9 < old) any=true;
                else reverse(r.begin()+i,r.begin()+j+1);
            }
        }
        if(!any) break;
    }
}
static uint64_t hilbert(uint32_t x, uint32_t y) {
    uint64_t d=0;
    for(uint32_t s=1u<<29;s;s>>=1) {
        uint32_t rx=(x&s)!=0, ry=(y&s)!=0;
        d += (uint64_t)s*(uint64_t)s*((3*rx)^ry);
        if(!ry) { if(rx) { x=(1u<<30)-1-x; y=(1u<<30)-1-y; } swap(x,y); }
    }
    return d;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>n)) return 0;
    p.resize(n);
    long long minx=LLONG_MAX,maxx=LLONG_MIN,miny=LLONG_MAX,maxy=LLONG_MIN;
    for(auto &a:p) { cin>>a.x>>a.y; minx=min(minx,a.x); maxx=max(maxx,a.x); miny=min(miny,a.y); maxy=max(maxy,a.y); }
    primeCity.assign(n,true);
    if(n>0) primeCity[0]=false;
    if(n>1) primeCity[1]=false;
    for(int i=2;(long long)i*i<n;i++) if(primeCity[i]) for(int j=i*i;j<n;j+=i) primeCity[j]=false;

    vector<int> best(n); iota(best.begin(),best.end(),0);
    double bestCost=score(best);
    // Reversing a cycle preserves its geometric edges but not which city is
    // the source at the globally numbered carrot steps.  Test both directions.
    auto consider = [&](const vector<int>& r) {
        double v=score(r);
        if(v<bestCost) { bestCost=v; best=r; }
        vector<int> back(n); back[0]=0;
        for(int i=1;i<n;i++) back[i]=r[n-i];
        v=score(back);
        if(v<bestCost) { bestCost=v; best.swap(back); }
    };
    // The input-order baseline is also a spatial cycle direction; its reverse
    // can differ under the position-dependent penalty just like every other tour.
    consider(best);
    // Alternating y bands retain the input's x order inside every band, so
    // building many resolutions is linear rather than a succession of sorts.
    long long yr=maxy-miny;
    for(int bands: {2,4,8,16,32,64,128}) for(int flip=0;flip<2;flip++) {
        vector<vector<int>> box(bands);
        for(int id=0;id<n;id++) {
            int b = yr==0 ? 0 : (int)((__int128)(p[id].y-miny)*bands/(yr+1));
            if(b>=bands) b=bands-1;
            box[b].push_back(id);
        }
        vector<int> cyc; cyc.reserve(n);
        for(int b=0;b<bands;b++) {
            bool rev=((b&1)^flip)!=0;
            if(!rev) cyc.insert(cyc.end(),box[b].begin(),box[b].end());
            else cyc.insert(cyc.end(),box[b].rbegin(),box[b].rend());
        }
        consider(cutAtZero(cyc));
    }
    // Two Hilbert orientations supply a locality-preserving alternative when
    // horizontal band sweeps cut across compact two-dimensional clusters.
    for(int trans=0;trans<2;trans++) {
        vector<pair<uint64_t,int>> a; a.reserve(n);
        long long lo1=trans?miny:minx, hi1=trans?maxy:maxx;
        long long lo2=trans?minx:miny, hi2=trans?maxx:maxy;
        long long r1=hi1-lo1, r2=hi2-lo2;
        for(int id=0;id<n;id++) {
            long long u=trans?p[id].y:p[id].x, v=trans?p[id].x:p[id].y;
            uint32_t X=r1? (uint32_t)((__int128)(u-lo1)*((1u<<30)-1)/r1) : 0;
            uint32_t Y=r2? (uint32_t)((__int128)(v-lo2)*((1u<<30)-1)/r2) : 0;
            a.push_back({hilbert(X,Y),id});
        }
        sort(a.begin(),a.end()); vector<int> cyc; cyc.reserve(n);
        for(auto z:a) cyc.push_back(z.second);
        consider(cutAtZero(cyc));
    }
    improve(best);
    cout<<n+1<<'\n';
    for(int x:best) cout<<x<<'\n';
    cout<<0<<'\n';
}
