#include <bits/stdc++.h>
using namespace std;

// Periodic lattice packings.  We choose the best of cubic, FCC, and BCC
// lattices after cropping it by the cube; every emitted point has the
// corresponding lattice radius as its distance to the boundary.
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    long long n;
    if (!(cin >> n)) return 0;

    enum Type { GRID, FCC, BCC } bestType = GRID;
    double bestR = -1, fp = 0, bp = 0;

    // A balanced ordinary grid is useful at the discontinuities of a lattice.
    long long g = 1;
    while (g * g * g < n) ++g;
    double gr = 1.0 / (2.0 * g);
    bestR = gr;

    // FCC: integer triples of one parity, nearest-neighbour distance sqrt(2)q.
    long long k = 1;
    auto fcount = [](long long t) -> long long {
        __int128 v = (__int128)t * t * t;
        return (long long)((v + (t & 1)) / 2);
    };
    while (fcount(k) < n) ++k;
    fp = 1.0 / ((double)(k - 1) + sqrt(2.0));
    double fr = fp / sqrt(2.0);
    if (fr > bestR) bestR = fr, bestType = FCC;

    // BCC: integer-cell corners plus cell centers.  Its nearest distance is
    // sqrt(3)a/2, so the face margin is sqrt(3)a/4.
    auto bcount = [](double a) -> long long {
        const double h = sqrt(3.0) * 0.5;
        long long kc = max(0LL, (long long)floor(1.0 / a - h) + 1);
        long long kb = max(0LL, (long long)floor(1.0 / a - h - 0.5) + 1);
        __int128 v = (__int128)kc * kc * kc + (__int128)kb * kb * kb;
        return v > LLONG_MAX ? LLONG_MAX : (long long)v;
    };
    double lo = 0, hi = 2.0;
    for (int it = 0; it < 80; ++it) {
        double mid = (lo + hi) * .5;
        if (bcount(mid) >= n) lo = mid; else hi = mid;
    }
    bp = lo * (1.0 - 1e-10); // stay on the safe side of a count discontinuity
    double br = sqrt(3.0) * bp * .25;
    if (br > bestR) bestR = br, bestType = BCC;

    cout << setprecision(17);
    long long out = 0;
    if (bestType == GRID) {
        // g by g by g contains enough sites; omitting sites cannot hurt.
        for (long long i = 0; i < g && out < n; ++i)
            for (long long j = 0; j < g && out < n; ++j)
                for (long long l = 0; l < g && out < n; ++l) {
                    cout << (i + .5) / g << ' ' << (j + .5) / g << ' ' << (l + .5) / g << '\n';
                    ++out;
                }
    } else if (bestType == FCC) {
        // Shrinking q infinitesimally avoids accidental outward rounding.
        double q = fp * (1.0 - 1e-13), r = q / sqrt(2.0);
        // k remains a sufficient bound after the shrink (and can only gain sites).
        for (long long i = 0; i < k && out < n; ++i)
            for (long long j = 0; j < k && out < n; ++j)
                for (long long l = 0; l < k && out < n; ++l)
                    if (((i + j + l) & 1) == 0) {
                        cout << r + q*i << ' ' << r + q*j << ' ' << r + q*l << '\n';
                        ++out;
                    }
    } else {
        double a = bp, r = sqrt(3.0) * a * .25;
        const double h = sqrt(3.0) * .5;
        long long kc = max(0LL, (long long)floor(1.0 / a - h) + 1);
        long long kb = max(0LL, (long long)floor(1.0 / a - h - .5) + 1);
        for (long long i = 0; i < kc && out < n; ++i)
            for (long long j = 0; j < kc && out < n; ++j)
                for (long long l = 0; l < kc && out < n; ++l) {
                    cout << r + a*i << ' ' << r + a*j << ' ' << r + a*l << '\n';
                    ++out;
                }
        for (long long i = 0; i < kb && out < n; ++i)
            for (long long j = 0; j < kb && out < n; ++j)
                for (long long l = 0; l < kb && out < n; ++l) {
                    cout << r + a*(i+.5) << ' ' << r + a*(j+.5) << ' ' << r + a*(l+.5) << '\n';
                    ++out;
                }
    }
    return 0;
}
