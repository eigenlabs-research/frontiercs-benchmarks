#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; int id; };
struct Node {
    int p, l = -1, r = -1, par = -1, cnt = 1;
    long long lx, rx, ly, ry;
};

static vector<Point> pt, city;
static vector<Node> tr;
static vector<int> where;
static vector<char> alive;

int build(int lo, int hi, int dep, int parent) {
    if (lo >= hi) return -1;
    int mid = (lo + hi) / 2;
    nth_element(pt.begin() + lo, pt.begin() + mid, pt.begin() + hi,
        [dep](const Point& a, const Point& b) {
            return (dep & 1) ? (a.y != b.y ? a.y < b.y : a.x < b.x)
                             : (a.x != b.x ? a.x < b.x : a.y < b.y);
        });
    int v = (int)tr.size();
    tr.push_back({mid, -1, -1, parent, 1, pt[mid].x, pt[mid].x, pt[mid].y, pt[mid].y});
    where[pt[mid].id] = v;
    tr[v].l = build(lo, mid, dep + 1, v);
    tr[v].r = build(mid + 1, hi, dep + 1, v);
    auto pull = [&](int u) {
        Node &a = tr[u];
        a.cnt = 1;
        a.lx = a.rx = pt[a.p].x; a.ly = a.ry = pt[a.p].y;
        for (int c : {a.l, a.r}) if (c != -1) {
            Node &b = tr[c]; a.cnt += b.cnt;
            a.lx = min(a.lx, b.lx); a.rx = max(a.rx, b.rx);
            a.ly = min(a.ly, b.ly); a.ry = max(a.ry, b.ry);
        }
    };
    pull(v);
    return v;
}

// Bboxes remain conservative after deletions; cnt lets exhausted subtrees vanish.
static inline long double boxDist(int v, long long x, long long y) {
    const Node &a = tr[v];
    long double dx = x < a.lx ? (long double)a.lx-x : x > a.rx ? (long double)x-a.rx : 0;
    long double dy = y < a.ly ? (long double)a.ly-y : y > a.ry ? (long double)y-a.ry : 0;
    return dx*dx + dy*dy;
}
void nearest(int v, long long x, long long y, int &best, long double &bd) {
    if (v == -1 || tr[v].cnt == 0 || boxDist(v,x,y) >= bd) return;
    Node &a = tr[v];
    if (alive[a.p]) {
        long double dx=(long double)pt[a.p].x-x, dy=(long double)pt[a.p].y-y;
        long double d=dx*dx+dy*dy;
        if (d < bd) bd=d, best=a.p;
    }
    int u=a.l, w=a.r;
    if (u != -1 && w != -1 && boxDist(w,x,y) < boxDist(u,x,y)) swap(u,w);
    nearest(u,x,y,best,bd); nearest(w,x,y,best,bd);
}
void erasePoint(int idx) {
    alive[idx]=0;
    for (int v=where[pt[idx].id]; v!=-1; v=tr[v].par) tr[v].cnt--;
}

static vector<char> primeFlag;
static inline long double edge(const vector<int>& r, int t) {
    const Point &a=city[r[t-1]], &b=city[r[t]];
    long double dx=(long double)a.x-b.x, dy=(long double)a.y-b.y;
    long double z=sqrtl(dx*dx+dy*dy);
    if (t%10==0 && !primeFlag[a.id]) z*=1.1L;
    return z;
}
long double cost(const vector<int>& r) { long double z=0; for(int t=1;t<(int)r.size();++t) z+=edge(r,t); return z; }
void adjacentImprove(vector<int>& r) {
    int n=(int)r.size()-1;
    for (int pass=0; pass<3; ++pass) {
        bool changed=false;
        for (int i=1;i<n-1;++i) {
            long double old=0; for(int t=i;t<=i+2 && t<=n;++t) old+=edge(r,t);
            swap(r[i],r[i+1]);
            long double now=0; for(int t=i;t<=i+2 && t<=n;++t) now+=edge(r,t);
            if(now + 1e-9L < old) changed=true;
            else swap(r[i],r[i+1]);
        }
        if(!changed) break;
    }
}

// A bounded 2-opt pass repairs small crossings left by a constructive ordering.
// Unlike ordinary TSP 2-opt, reversing a segment changes which city starts each
// tenth edge, so recompute every affected edge exactly before accepting it.
void shortTwoOpt(vector<int>& r) {
    int n=(int)r.size()-1;
    const int MAX_SPAN=8;
    for (int i=1; i<n-1; ++i) {
        int last=min(n-1, i+MAX_SPAN);
        for (int j=i+2; j<=last; ++j) {
            long double old=0;
            for (int t=i; t<=j+1; ++t) old+=edge(r,t);
            reverse(r.begin()+i, r.begin()+j+1);
            long double now=0;
            for (int t=i; t<=j+1; ++t) now+=edge(r,t);
            if (now + 1e-9L >= old) reverse(r.begin()+i, r.begin()+j+1);
        }
    }
}

// Map a square-grid coordinate to its Hilbert position.  Unlike a row/column
// sort, consecutive positions remain local even when the data bends or forms clusters.
static uint64_t hilbert(int64_t x, int64_t y) {
    // The standard iterative rotation deliberately permits temporary negative
    // coordinates after a reflected quadrant, so keep these signed.
    uint64_t d=0;
    for (int64_t s=1LL<<30; s; s>>=1) {
        int rx=(x&s)!=0, ry=(y&s)!=0;
        d += (uint64_t)s*(uint64_t)s*((3*rx)^ry);
        if (!ry) {
            if (rx) { x=s-1-x; y=s-1-y; }
            swap(x,y);
        }
    }
    return d;
}

static vector<int> hilbertCycle(bool reverseDirection) {
    int n=(int)city.size();
    long long minx=city[0].x, maxx=city[0].x, miny=city[0].y, maxy=city[0].y;
    for (const Point &p: city) {
        minx=min(minx,p.x); maxx=max(maxx,p.x);
        miny=min(miny,p.y); maxy=max(maxy,p.y);
    }
    unsigned long long sx=(unsigned long long)(maxx-minx), sy=(unsigned long long)(maxy-miny);
    const unsigned long long LIM=(1ULL<<31)-1;
    vector<pair<uint64_t,int>> keyed; keyed.reserve(n);
    for (const Point &p: city) {
        int64_t xx=sx ? (int64_t)((unsigned long long)(p.x-minx)*LIM/sx) : 0;
        int64_t yy=sy ? (int64_t)((unsigned long long)(p.y-miny)*LIM/sy) : 0;
        keyed.push_back({hilbert(xx,yy),p.id});
    }
    sort(keyed.begin(),keyed.end());
    int at=0; while(keyed[at].second!=0) ++at;
    vector<int> r; r.reserve(n+1); r.push_back(0);
    for(int k=1;k<n;k++) {
        int pos=reverseDirection ? (at-k+n)%n : (at+k)%n;
        r.push_back(keyed[pos].second);
    }
    r.push_back(0);
    return r;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    pt.resize(n); where.assign(n,-1);
    for(int i=0;i<n;i++){cin>>pt[i].x>>pt[i].y; pt[i].id=i;}
    city=pt; // pt is later permuted while building the KD tree; routes use city IDs.
    primeFlag.assign(n,true); if(n>0) primeFlag[0]=false; if(n>1) primeFlag[1]=false;
    for(int i=2;i*(long long)i<n;i++) if(primeFlag[i]) for(int j=i*i;j<n;j+=i) primeFlag[j]=false;

    vector<int> best; best.reserve(n+1);
    best.push_back(0); for(int i=1;i<n;i++) best.push_back(i); best.push_back(0);
    adjacentImprove(best); long double bestCost=cost(best);

    // A cheap stripe ordering is useful on horizontally banded data, while retaining the baseline fallback.
    vector<int> stripe; stripe.reserve(n+1); stripe.push_back(0);
    vector<int> ord; for(int i=1;i<n;i++) ord.push_back(i);
    sort(ord.begin(),ord.end(),[](int a,int b){return pt[a].y!=pt[b].y ? pt[a].y<pt[b].y : pt[a].x<pt[b].x;});
    stripe.insert(stripe.end(),ord.begin(),ord.end()); stripe.push_back(0);
    adjacentImprove(stripe); long double cc=cost(stripe);
    if(cc<bestCost) bestCost=cc,best.swap(stripe);

    // A cyclic Hilbert traversal is a separate spatial-decomposition candidate.
    // Cutting its cycle at city 0 permits both directions while respecting the output contract.
    for (bool backwards : {false, true}) {
        vector<int> curve=hilbertCycle(backwards);
        adjacentImprove(curve); cc=cost(curve);
        if(cc<bestCost) bestCost=cc,best.swap(curve);
    }

    tr.reserve(n); int root=build(0,n,0,-1); alive.assign(n,1);
    vector<int> nn; nn.reserve(n+1); nn.push_back(0);
    // erasePoint takes a point-array index; locate city 0's stored point.
    int zeroIndex=tr[where[0]].p; erasePoint(zeroIndex);
    int cur=zeroIndex;
    for(int step=1;step<n;step++) {
        int q=-1; long double bd=numeric_limits<long double>::infinity();
        nearest(root,pt[cur].x,pt[cur].y,q,bd);
        nn.push_back(pt[q].id); erasePoint(q); cur=q;
    }
    nn.push_back(0); adjacentImprove(nn); cc=cost(nn);
    if(cc<bestCost) best.swap(nn);

    // Keep the portfolio choice intact, then monotonically refine that exact route.
    shortTwoOpt(best);
    cout<<n+1<<'\n'; for(int v:best) cout<<v<<'\n';
}
