#include <bits/stdc++.h>
using namespace std;

struct Candidate {
    vector<array<int,3>> p;
    double unit;                 // minimum lattice-point separation in integer coordinates
    int span = 0;
    double radius = -1;
};

static long long fccCount(int q) {
    // Points with x+y+z even in a q by q by q integer box.
    return ((long long)q*q*q + (q & 1)) / 2;
}
static long long bccCount(int q) {
    long long e = (q + 1) / 2, o = q / 2;
    return e*e*e + o*o*o;        // x,y,z all have the same parity
}

static void compactOrder(Candidate &c, int q, int n) {
    const double mid = (q - 1) * 0.5;
    sort(c.p.begin(), c.p.end(), [mid](const auto &a, const auto &b) {
        double da = (a[0]-mid)*(a[0]-mid) + (a[1]-mid)*(a[1]-mid) + (a[2]-mid)*(a[2]-mid);
        double db = (b[0]-mid)*(b[0]-mid) + (b[1]-mid)*(b[1]-mid) + (b[2]-mid)*(b[2]-mid);
        if (da != db) return da < db;
        return a < b;
    });
    c.p.resize(n);
    int lo[3] = {INT_MAX, INT_MAX, INT_MAX}, hi[3] = {INT_MIN, INT_MIN, INT_MIN};
    for (auto a : c.p) for (int d=0; d<3; ++d) lo[d]=min(lo[d],a[d]), hi[d]=max(hi[d],a[d]);
    c.span = max({hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2]});
    c.radius = c.unit / (2.0 * (c.span + c.unit));
}

static Candidate makeFCC(int n) {
    int q=1; while (fccCount(q) < n) ++q;
    Candidate c; c.unit=sqrt(2.0); c.p.reserve(fccCount(q));
    for (int x=0;x<q;++x) for (int y=0;y<q;++y) for (int z=0;z<q;++z)
        if (((x+y+z)&1)==0) c.p.push_back({x,y,z});
    compactOrder(c,q,n);
    return c;
}
static Candidate makeBCC(int n) {
    int q=1; while (bccCount(q) < n) ++q;
    Candidate c; c.unit=sqrt(3.0); c.p.reserve(bccCount(q));
    for (int x=0;x<q;++x) for (int y=0;y<q;++y) for (int z=0;z<q;++z)
        if ((x&1)==(y&1) && (y&1)==(z&1)) c.p.push_back({x,y,z});
    compactOrder(c,q,n);
    return c;
}
static Candidate makeGrid(int n) {
    int q=1; while ((long long)q*q*q < n) ++q;
    Candidate c; c.unit=1.0; c.p.reserve(n);
    // Center-first ordering gives the same safe grid separation and can reduce a partial box span.
    for (int x=0;x<q;++x) for (int y=0;y<q;++y) for (int z=0;z<q;++z) c.p.push_back({x,y,z});
    compactOrder(c,q,n);
    return c;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    Candidate best = makeGrid(n), f = makeFCC(n), b = makeBCC(n);
    if (f.radius > best.radius) best = move(f);
    if (b.radius > best.radius) best = move(b);

    int lo[3] = {INT_MAX, INT_MAX, INT_MAX};
    for (auto a : best.p) for (int d=0; d<3; ++d) lo[d]=min(lo[d],a[d]);
    // Margin in lattice coordinates is half the minimum separation.  This makes
    // the face clearance equal to half the nearest-neighbour clearance.
    double margin = best.unit * 0.5;
    double den = best.span + best.unit;
    cout << setprecision(17);
    for (auto a : best.p) {
        cout << (a[0]-lo[0]+margin)/den << ' '
             << (a[1]-lo[1]+margin)/den << ' '
             << (a[2]-lo[2]+margin)/den << '\n';
    }
}
