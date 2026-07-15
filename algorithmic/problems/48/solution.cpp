#include <bits/stdc++.h>
using namespace std;

struct P { double x, y, z; };

// Number of even and odd integers in the inclusive interval [l,h].
static pair<long long,long long> parityCount(int l, int h) {
    if (l > h) return {0, 0};
    long long total = (long long)h - l + 1;
    long long ev = total / 2;
    if (total & 1) ev += ((l & 1) == 0);
    return {ev, total - ev};
}

struct Phase {
    double q[3];
    double d = 0;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    cout << setprecision(17);

    // The very small cases are clipped lattice simplices.  In particular this
    // avoids paying for unused lattice layers at the boundary.
    if (n <= 4) {
        vector<array<int,3>> a;
        double sep;
        if (n == 2) {
            a = {{0,0,0}, {1,1,1}};
            sep = sqrt(3.0);
        } else {
            a = {{0,0,0}, {0,1,1}, {1,0,1}, {1,1,0}};
            sep = sqrt(2.0);
        }
        double step = 1.0 / (1.0 + sep);
        double margin = step * sep * .5;
        for (int i = 0; i < n; ++i)
            cout << margin + step*a[i][0] << ' ' << margin + step*a[i][1]
                 << ' ' << margin + step*a[i][2] << '\n';
        return 0;
    }

    // A guaranteed baseline candidate: a cubic grid with ceil(cuberoot(n))
    // positions on each side.  It is retained if a clipped FCC lattice does
    // not actually beat it (important for partially filled small layers).
    int k = 1;
    while (1LL*k*k*k < n) ++k;
    double gridR = 1.0 / (2.0*k);

    // FCC points are integer triples of even coordinate sum, scaled by
    // d/sqrt(2).  For a proposed nearest-neighbour distance d, retain only
    // points in [d/2,1-d/2]^3.  Trying translations of the lattice cell is a
    // cheap boundary-condition experiment rather than assuming one arbitrary
    // alignment is best.
    vector<Phase> phases;
    const double vals[] = {0.0, 0.25, 0.5, 0.75};
    for (double x : vals) for (double y : vals) for (double z : vals)
        phases.push_back({{x,y,z}, 0.0});

    Phase best;
    for (Phase ph : phases) {
        double lo = 0.0, hi = 1.0;
        for (int it = 0; it < 58; ++it) {
            double d = (lo + hi) * .5;
            // In units of d/sqrt(2), this is the available interval after
            // reserving a face margin d/2 on both sides.
            double B = sqrt(2.0) * (1.0 / d - 1.0);
            int lx[3], hx[3];
            for (int t = 0; t < 3; ++t) {
                lx[t] = (int)ceil(-ph.q[t] - 1e-11);
                hx[t] = (int)floor(B - ph.q[t] + 1e-11);
            }
            auto X = parityCount(lx[0], hx[0]);
            auto Y = parityCount(lx[1], hx[1]);
            auto Z = parityCount(lx[2], hx[2]);
            long long cnt = X.first*Y.first*Z.first + X.first*Y.second*Z.second
                          + X.second*Y.first*Z.second + X.second*Y.second*Z.first;
            if (cnt >= n) lo = d; else hi = d;
        }
        if (lo > best.d) best = ph, best.d = lo;
    }

    if (best.d * .5 <= gridR * (1.0 + 1e-12)) {
        // Exact-n prefix of a full grid; every retained point still has the
        // advertised face margin and grid separation.
        for (int i = 0; i < n; ++i) {
            int x = i % k, y = (i / k) % k, z = i / (k*k);
            cout << (x + .5) / k << ' ' << (y + .5) / k << ' ' << (z + .5) / k << '\n';
        }
        return 0;
    }

    // Move infinitesimally inside the feasible value found by binary search,
    // so floating point formatting cannot put a boundary point beyond a face.
    double d = best.d * (1.0 - 2e-12);
    double s = d / sqrt(2.0), r = d * .5;
    double B = sqrt(2.0) * (1.0 / d - 1.0);
    int lx[3], hx[3];
    for (int t = 0; t < 3; ++t) {
        lx[t] = (int)ceil(-best.q[t] - 1e-10);
        hx[t] = (int)floor(B - best.q[t] + 1e-10);
    }
    int put = 0;
    for (int x = lx[0]; x <= hx[0] && put < n; ++x)
        for (int y = lx[1]; y <= hx[1] && put < n; ++y)
            for (int z = lx[2]; z <= hx[2] && put < n; ++z) {
                if ((x + y + z) & 1) continue;
                double X = r + s*(x + best.q[0]);
                double Y = r + s*(y + best.q[1]);
                double Z = r + s*(z + best.q[2]);
                // These clamps are only a final guard against a one-ulp roundoff.
                X = min(1.0, max(0.0, X));
                Y = min(1.0, max(0.0, Y));
                Z = min(1.0, max(0.0, Z));
                cout << X << ' ' << Y << ' ' << Z << '\n';
                ++put;
            }
    return 0;
}
