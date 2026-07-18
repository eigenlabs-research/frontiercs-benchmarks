#include <bits/stdc++.h>
using namespace std;

struct Pt { double x,y,z; };
struct Cand { double r; vector<Pt> pts; };

static double coord(int i, int d, double r) {
    if (d <= 1) return 0.5;
    return r + (1.0 - 2.0*r) * (double)i / (double)(d - 1);
}

static long long parityCount(int a,int b,int c,int par) {
    long long e1=(a+1)/2,o1=a/2, e2=(b+1)/2,o2=b/2, e3=(c+1)/2,o3=c/2;
    long long even = e1*e2*e3 + e1*o2*o3 + o1*e2*o3 + o1*o2*e3;
    long long total = 1LL*a*b*c;
    return par==0 ? even : total-even;
}

static double parityRadius(int a,int b,int c) {
    double bestA = 1e100;
    int d[3] = {a,b,c};
    for (int dx=-2; dx<=2; ++dx) for (int dy=-2; dy<=2; ++dy) for (int dz=-2; dz<=2; ++dz) {
        if (dx==0 && dy==0 && dz==0) continue;
        if (((dx+dy+dz)&1) != 0) continue; // same checkerboard parity
        int del[3] = {abs(dx),abs(dy),abs(dz)};
        bool ok = true;
        double s2 = 0.0;
        for (int t=0;t<3;++t) {
            if (del[t] >= d[t]) { ok=false; break; }
            if (del[t] && d[t] <= 1) { ok=false; break; }
            if (del[t]) {
                double v = (double)del[t] / (double)(d[t]-1);
                s2 += v*v;
            }
        }
        if (ok && s2 > 0) bestA = min(bestA, sqrt(s2));
    }
    if (bestA > 1e50) return 0.5;
    return bestA / (2.0 + 2.0*bestA);
}

static Cand smallCorners(int n) {
    vector<array<int,3>> v;
    for (int x=0;x<2;++x) for (int y=0;y<2;++y) for (int z=0;z<2;++z) v.push_back({x,y,z});
    int bestMask = 0, bestH = -1;
    for (int mask=0; mask<256; ++mask) if (__builtin_popcount((unsigned)mask)==n) {
        int mh = 10;
        for (int i=0;i<8;++i) if (mask>>i&1) for (int j=i+1;j<8;++j) if (mask>>j&1) {
            int h = (v[i][0]!=v[j][0]) + (v[i][1]!=v[j][1]) + (v[i][2]!=v[j][2]);
            mh = min(mh, h);
        }
        if (mh > bestH) { bestH = mh; bestMask = mask; }
    }
    double A = sqrt((double)bestH);
    double r = A / (2.0 + 2.0*A);
    Cand c; c.r = r;
    for (int i=0;i<8 && (int)c.pts.size()<n;++i) if (bestMask>>i&1) {
        c.pts.push_back({v[i][0] ? 1.0-r : r, v[i][1] ? 1.0-r : r, v[i][2] ? 1.0-r : r});
    }
    return c;
}

static Cand cubicGrid(int n) {
    int q = 1;
    while (1LL*q*q*q < n) ++q;
    double r = 1.0 / (2.0*q);
    Cand c; c.r = r; c.pts.reserve(n);
    for (int i=0;i<q && (int)c.pts.size()<n;++i)
        for (int j=0;j<q && (int)c.pts.size()<n;++j)
            for (int k=0;k<q && (int)c.pts.size()<n;++k)
                c.pts.push_back({(i+0.5)/q, (j+0.5)/q, (k+0.5)/q});
    return c;
}

static Cand parityGrid(int n) {
    int root = (int)ceil(cbrt((double)max(2,2*n)));
    int maxD = max(4, 3*root + 12);
    double bestR = -1.0; int ba=1,bb=1,bc=1,bp=0;
    for (int a=1; a<=maxD; ++a) for (int b=1; b<=maxD; ++b) for (int c=1; c<=maxD; ++c) {
        long long cnt0 = parityCount(a,b,c,0), cnt1 = 1LL*a*b*c - cnt0;
        int par = (cnt1 > cnt0) ? 1 : 0;
        if (max(cnt0,cnt1) < n) continue;
        double r = parityRadius(a,b,c);
        // Prefer slightly squarer boxes on ties for robustness of the truncated prefix.
        double tie = r - 1e-12*(a+b+c);
        double btie = bestR - 1e-12*(ba+bb+bc);
        if (tie > btie) { bestR = r; ba=a; bb=b; bc=c; bp=par; }
    }
    Cand res; res.r = bestR; res.pts.reserve(n);
    for (int i=0;i<ba && (int)res.pts.size()<n;++i)
        for (int j=0;j<bb && (int)res.pts.size()<n;++j)
            for (int k=0;k<bc && (int)res.pts.size()<n;++k)
                if (((i+j+k)&1) == bp)
                    res.pts.push_back({coord(i,ba,bestR), coord(j,bb,bestR), coord(k,bc,bestR)});
    return res;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n; if (!(cin >> n)) return 0;
    vector<Cand> cs;
    cs.push_back(cubicGrid(n));
    cs.push_back(parityGrid(n));
    if (n <= 8) cs.push_back(smallCorners(n));
    int bi = 0;
    for (int i=1;i<(int)cs.size();++i) if (cs[i].r > cs[bi].r) bi = i;
    cout.setf(ios::fixed); cout << setprecision(12);
    for (int i=0;i<n;++i) {
        const Pt &p = cs[bi].pts[i];
        cout << min(1.0,max(0.0,p.x)) << ' ' << min(1.0,max(0.0,p.y)) << ' ' << min(1.0,max(0.0,p.z)) << '\n';
    }
    return 0;
}
