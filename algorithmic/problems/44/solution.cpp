#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x, y; };

static uint64_t hilbertOrder(unsigned x, unsigned y, int pow=21) {
    uint64_t d = 0;
    for (int s = pow - 1; s >= 0; --s) {
        unsigned rx = (x >> s) & 1u;
        unsigned ry = (y >> s) & 1u;
        d += (uint64_t)((3u * rx) ^ ry) << (2*s);
        if (ry == 0) {
            if (rx == 1) {
                x = (1u << pow) - 1 - x;
                y = (1u << pow) - 1 - y;
            }
            swap(x, y);
        }
    }
    return d;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    if (!(cin >> N)) return 0;
    vector<Pt> p(N);
    for (int i=0;i<N;i++) cin >> p[i].x >> p[i].y;

    vector<char> isPrime(max(2,N), true);
    isPrime[0]=false; if (N>1) isPrime[1]=false;
    for (long long i=2;i*i<N;i++) if (isPrime[i])
        for (long long j=i*i;j<N;j+=i) isPrime[(int)j]=false;

    auto dist = [&](int a, int b)->double{
        long double dx = (long double)p[a].x - p[b].x;
        long double dy = (long double)p[a].y - p[b].y;
        return (double)sqrt(dx*dx + dy*dy);
    };
    auto stepCost = [&](const vector<int>& r, int t)->double{
        if (t < 1 || t > N) return 0.0;
        double m = (t % 10 == 0 && !isPrime[r[t-1]]) ? 1.1 : 1.0;
        return m * dist(r[t-1], r[t]);
    };
    auto routeCost = [&](const vector<int>& r)->double{
        double s=0;
        for (int t=1;t<=N;t++) s += stepCost(r,t);
        return s;
    };

    vector<vector<int>> candidates;
    auto makeRouteFromMiddle = [&](vector<int> mid){
        vector<int> r; r.reserve(N+1); r.push_back(0);
        for (int v: mid) if (v!=0) r.push_back(v);
        r.push_back(0);
        candidates.push_back(move(r));
    };

    vector<int> mid;
    mid.reserve(max(0,N-1));
    for (int i=1;i<N;i++) mid.push_back(i);
    makeRouteFromMiddle(mid);                         // strengthened baseline fallback

    mid.clear();
    for (int i=1;i<N;i++) mid.push_back(i);
    sort(mid.begin(), mid.end(), [&](int a,int b){
        if (p[a].y != p[b].y) return p[a].y < p[b].y;
        return p[a].x < p[b].x;
    });
    makeRouteFromMiddle(mid);

    long long minx=p[0].x, maxx=p[0].x, miny=p[0].y, maxy=p[0].y;
    for (auto &q:p) { minx=min(minx,q.x); maxx=max(maxx,q.x); miny=min(miny,q.y); maxy=max(maxy,q.y); }
    const unsigned M = (1u<<21) - 1;
    auto norm = [&](long long v, long long lo, long long hi)->unsigned{
        if (hi == lo) return 0;
        long double z = ((long double)v - lo) * (long double)M / ((long double)hi - lo);
        if (z < 0) z = 0; if (z > M) z = M;
        return (unsigned)(z + 0.5L);
    };

    for (int mode=0; mode<4; ++mode) {
        vector<pair<uint64_t,int>> ord; ord.reserve(N);
        for (int i=0;i<N;i++) {
            unsigned X = norm(p[i].x, minx, maxx);
            unsigned Y = norm(p[i].y, miny, maxy);
            if (mode & 1) Y = M - Y;
            if (mode & 2) swap(X,Y);
            ord.push_back({hilbertOrder(X,Y), i});
        }
        sort(ord.begin(), ord.end(), [](auto &a, auto &b){ return a.first < b.first || (a.first==b.first && a.second < b.second); });
        int pos0=0; for (int i=0;i<N;i++) if (ord[i].second==0) { pos0=i; break; }
        vector<int> r; r.reserve(N+1); r.push_back(0);
        for (int k=1;k<N;k++) r.push_back(ord[(pos0+k)%N].second);
        r.push_back(0); candidates.push_back(r);
        r.clear(); r.reserve(N+1); r.push_back(0);
        for (int k=1;k<N;k++) r.push_back(ord[(pos0-k+N)%N].second);
        r.push_back(0); candidates.push_back(move(r));
    }

    int bestIdx = 0; double bestCost = routeCost(candidates[0]);
    for (int i=1;i<(int)candidates.size();++i) {
        double c = routeCost(candidates[i]);
        if (c < bestCost) bestCost = c, bestIdx = i;
    }
    vector<int> r = move(candidates[bestIdx]);

    auto swapDelta = [&](int a, int b)->double{
        if (a<=0 || b<=0 || a>=N || b>=N || a==b) return 0.0;
        if (a>b) swap(a,b);
        int ts[4] = {a, a+1, b, b+1};
        double old=0, neu=0;
        vector<int> used;
        for (int t: ts) if (t>=1 && t<=N && find(used.begin(), used.end(), t)==used.end()) {
            used.push_back(t); old += stepCost(r,t);
        }
        swap(r[a], r[b]);
        for (int t: used) neu += stepCost(r,t);
        swap(r[a], r[b]);
        return neu - old;
    };

    // Put prime-numbered cities on as many expensive 10th-step source positions as is locally profitable.
    for (int pos=9; pos<N; pos+=10) {
        if (isPrime[r[pos]]) continue;
        int best=-1; double bd=0.0;
        int L=max(1,pos-60), R=min(N-1,pos+60);
        for (int j=L;j<=R;j++) if (isPrime[r[j]]) {
            double d = swapDelta(pos,j);
            if (d < bd) bd=d, best=j;
        }
        if (best!=-1) swap(r[pos], r[best]);
    }

    // Cheap local swap descent.  It preserves validity and improves both geometry and carrot placement.
    for (int pass=0; pass<2; ++pass) {
        bool any=false;
        for (int i=1;i<N;i++) {
            int bj=-1; double bd=0.0;
            int R=min(N-1,i+10);
            for (int j=i+1;j<=R;j++) {
                double d = swapDelta(i,j);
                if (d < bd) bd=d, bj=j;
            }
            if (bj!=-1) { swap(r[i], r[bj]); any=true; }
        }
        if (!any) break;
    }

    // Bounded 2-opt style reversal descent.  Unlike arbitrary swaps, a short reversal can
    // remove local crossings created by space-filling orders while keeping the search linear.
    auto revDelta = [&](int l, int rr)->double{
        if (l < 1 || rr >= N || l >= rr) return 0.0;
        int a = l, b = min(N, rr+1);
        double old=0.0, neu=0.0;
        auto mappedCity = [&](int pos)->int{
            if (pos >= l && pos <= rr) pos = l + rr - pos;
            return r[pos];
        };
        for (int t=a; t<=b; ++t) {
            old += stepCost(r,t);
            int u = mappedCity(t-1), v = mappedCity(t);
            double m = (t % 10 == 0 && !isPrime[u]) ? 1.1 : 1.0;
            neu += m * dist(u,v);
        }
        return neu - old;
    };

    for (int pass=0; pass<2; ++pass) {
        bool any=false;
        for (int i=1; i<N-1; ++i) {
            int bestR=-1; double bd=0.0;
            int lim = min(N-1, i+8);
            for (int j=i+1; j<=lim; ++j) {
                double d = revDelta(i,j);
                if (d < bd) bd=d, bestR=j;
            }
            if (bestR!=-1) { reverse(r.begin()+i, r.begin()+bestR+1); any=true; }
        }
        if (!any) break;
    }

    cout << N+1 << '\n';
    for (int v: r) cout << v << '\n';
    return 0;
}
