#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };
static vector<Point> pt;
static vector<char> primeId;
static int n;

struct KDNode { int id, l=-1, r=-1; long long lx, rx, ly, ry; };
static vector<KDNode> kd;
static vector<char> usedKD;
static int buildKD(vector<int>& a, int lo, int hi, int dep) {
    if (lo >= hi) return -1;
    int m=(lo+hi)/2, axis=dep&1;
    nth_element(a.begin()+lo,a.begin()+m,a.begin()+hi,[axis](int u,int v) {
        return axis ? pt[u].y<pt[v].y : pt[u].x<pt[v].x;
    });
    int z=(int)kd.size(); kd.push_back({a[m],-1,-1,pt[a[m]].x,pt[a[m]].x,pt[a[m]].y,pt[a[m]].y});
    int l=buildKD(a,lo,m,dep+1), r=buildKD(a,m+1,hi,dep+1);
    kd[z].l=l; kd[z].r=r;
    for(int q:{l,r}) if(q>=0) { kd[z].lx=min(kd[z].lx,kd[q].lx); kd[z].rx=max(kd[z].rx,kd[q].rx); kd[z].ly=min(kd[z].ly,kd[q].ly); kd[z].ry=max(kd[z].ry,kd[q].ry); }
    return z;
}
static inline long double boxDist2(const KDNode& z, long long x, long long y) {
    long double dx=x<z.lx?z.lx-x:(x>z.rx?x-z.rx:0), dy=y<z.ly?z.ly-y:(y>z.ry?y-z.ry:0);
    return dx*dx+dy*dy;
}
static void nearestKD(int z, int from, int& ans, long double& best) {
    if(z<0 || boxDist2(kd[z],pt[from].x,pt[from].y)>=best) return;
    int v=kd[z].id;
    if(!usedKD[v]) { long double dx=(long double)pt[v].x-pt[from].x, dy=(long double)pt[v].y-pt[from].y, d=dx*dx+dy*dy; if(d<best) best=d,ans=v; }
    int a=kd[z].l,b=kd[z].r;
    if(a>=0 && b>=0 && boxDist2(kd[b],pt[from].x,pt[from].y)<boxDist2(kd[a],pt[from].x,pt[from].y)) swap(a,b);
    nearestKD(a,from,ans,best); nearestKD(b,from,ans,best);
}
static vector<int> makeNearest() {
    usedKD.assign(n,false); usedKD[0]=true;
    vector<int> alive; for(int i=1;i<n;i++) alive.push_back(i);
    int root=-1, remain=n-1, rebuild=remain, cur=0;
    vector<int> r; r.reserve(n+1); r.push_back(0);
    while(remain) {
        if(root<0 || remain*2<=rebuild) { vector<int> fresh; fresh.reserve(remain); for(int v:alive) if(!usedKD[v]) fresh.push_back(v); alive.swap(fresh); kd.clear(); kd.reserve(remain); root=buildKD(alive,0,(int)alive.size(),0); rebuild=remain; }
        int v=-1; long double best=numeric_limits<long double>::infinity(); nearestKD(root,cur,v,best);
        usedKD[v]=true; r.push_back(v); cur=v; --remain;
    }
    r.push_back(0); return r;
}

static inline double edge(int a, int b, int step) {
    double d = hypot((double)pt[a].x - pt[b].x, (double)pt[a].y - pt[b].y);
    return (step % 10 == 0 && !primeId[a]) ? d * 1.1 : d;
}

static double cost(const vector<int>& r) {
    double ans = 0;
    for (int t = 1; t <= n; ++t) ans += edge(r[t-1], r[t], t);
    return ans;
}

// One exact 3-edge descent pass.  This also accounts for the fact that the
// carrot multiplier belongs to an edge position, rather than to a city.
static void adjacentDescent(vector<int>& r) {
    for (int p = 1; p + 1 < n; ++p) {
        int a = r[p-1], b = r[p], c = r[p+1], d = r[p+2];
        double before = edge(a,b,p) + edge(b,c,p+1) + edge(c,d,p+2);
        double after  = edge(a,c,p) + edge(c,b,p+1) + edge(b,d,p+2);
        if (after + 1e-7 < before) swap(r[p], r[p+1]);
    }
}

// Adjacent swaps cannot repair a crossing whose endpoints are separated by a
// few cities.  Unlike ordinary 2-opt, reversing a segment also moves city IDs
// across carrot positions, so evaluate every affected edge at its true step.
// The small window keeps this a linear-time cleanup even at N=200000.
static void windowedTwoOpt(vector<int>& r) {
    const int W = 8;
    for (int p = 1; p + 1 < n; ++p) {
        int bestq = -1;
        double bestGain = 0.0, before = 0.0;
        int covered = p - 1;
        for (int q = p + 1; q <= min(n - 1, p + W); ++q) {
            // Extend the unchanged segment cost by one edge as q grows.
            while (covered < q + 1) {
                ++covered;
                before += edge(r[covered-1], r[covered], covered);
            }
            double after = 0.0;
            for (int t = p; t <= q + 1; ++t) {
                int u = (t-1 >= p && t-1 <= q) ? r[p+q-(t-1)] : r[t-1];
                int v = (t   >= p && t   <= q) ? r[p+q-t]     : r[t];
                after += edge(u, v, t);
            }
            double gain = before - after;
            if (gain > bestGain + 1e-7) bestGain = gain, bestq = q;
        }
        if (bestq >= 0) reverse(r.begin() + p, r.begin() + bestq + 1);
    }
}

static unsigned long long hilbertKey(unsigned int x, unsigned int y) {
    // Coordinates are 30-bit.  The classic Hilbert rotation formulation
    // gives a locality-preserving total order without floating point bins.
    unsigned long long d = 0;
    for (unsigned int s = 1u << 29; s; s >>= 1) {
        unsigned int rx = (x & s) != 0, ry = (y & s) != 0;
        d += (unsigned long long)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = ((1u<<30) - 1) - x; y = ((1u<<30) - 1) - y; }
            swap(x, y);
        }
    }
    return d;
}

static vector<int> makeHilbert() {
    vector<pair<unsigned long long,int>> v;
    v.reserve(n-1);
    long long minx=pt[0].x, maxx=pt[0].x, miny=pt[0].y, maxy=pt[0].y;
    for (auto &q: pt) { minx=min(minx,q.x); maxx=max(maxx,q.x); miny=min(miny,q.y); maxy=max(maxy,q.y); }
    long double sx = maxx == minx ? 0 : ((long double)((1u<<30)-1) / (maxx-minx));
    long double sy = maxy == miny ? 0 : ((long double)((1u<<30)-1) / (maxy-miny));
    for (int i=1;i<n;i++) {
        unsigned int x=(unsigned int)((pt[i].x-minx)*sx), y=(unsigned int)((pt[i].y-miny)*sy);
        v.push_back({hilbertKey(x,y),i});
    }
    sort(v.begin(),v.end());
    vector<int> r; r.reserve(n+1); r.push_back(0);
    for(auto z:v) r.push_back(z.second);
    r.push_back(0); return r;
}

// IDs are already x-sorted.  Within short consecutive x slabs, visit one
// horizontal direction at a time; reversing direction in each slab makes a
// robust lawnmower route for alternating-row and banded data.
static vector<int> makeSlabs(int width, bool firstUp) {
    vector<int> r; r.reserve(n+1); r.push_back(0);
    int slab=0;
    for (int l=1; l<n; l+=width, ++slab) {
        int u=min(n,l+width);
        vector<int> q; q.reserve(u-l);
        for(int i=l;i<u;i++) q.push_back(i);
        bool up = firstUp ^ (slab&1);
        sort(q.begin(),q.end(),[&](int a,int b) {
            if (pt[a].y != pt[b].y) return up ? pt[a].y < pt[b].y : pt[a].y > pt[b].y;
            return a < b;
        });
        r.insert(r.end(),q.begin(),q.end());
    }
    r.push_back(0); return r;
}

int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if (!(cin >> n)) return 0;
    pt.resize(n);
    for (auto &p:pt) cin >> p.x >> p.y;
    primeId.assign(n,true);
    if(n>0) primeId[0]=false;
    if(n>1) primeId[1]=false;
    for(int i=2;i*(long long)i<n;i++) if(primeId[i])
        for(int j=i*i;j<n;j+=i) primeId[j]=false;

    vector<int> best; best.reserve(n+1);
    for(int i=0;i<n;i++) best.push_back(i);
    best.push_back(0);
    adjacentDescent(best);
    double bestCost=cost(best);
    auto consider = [&](vector<int> r) {
        adjacentDescent(r);
        double c=cost(r);
        if(c < bestCost) { bestCost=c; best.swap(r); }
    };
    consider(makeHilbert());
    consider(makeNearest());
    // Geometric progression gives both very narrow row-following strips and
    // wider strips for noisy coordinates, without depending on test values.
    for (int w: {4,16,64,256}) {
        consider(makeSlabs(w,false));
        consider(makeSlabs(w,true));
    }
    // This is intentionally after selection: it is an exact-only refinement
    // of whichever incumbent constructor won, rather than another heuristic
    // whose quality could displace the preserved pool.
    windowedTwoOpt(best);
    cout << n+1 << '\n';
    for(int x:best) cout << x << '\n';
}
