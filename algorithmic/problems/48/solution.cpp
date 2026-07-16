#include <bits/stdc++.h>
using namespace std;

struct Box {
    long long a = 0, b = 0, c = 0;
    long double r = -1;
};

static long long ceil_div(long long x, long long y) { return (x + y - 1) / y; }

// FCC points are the integer triples of even coordinate sum.  For a box with
// spacings h_i, its short vectors include both (+-1,+-1,0) and, when a
// dimension has three layers, (+-2,0,0).  The latter is important for thin
// rectangular boxes.
static long double certified_radius(long long a, long long b, long long c) {
    long long d[3] = {a, b, c};
    long double inv[3];
    for (int i = 0; i < 3; ++i) inv[i] = 1.0L / (d[i] - 1);
    long double q = 1e100L;
    for (int i = 0; i < 3; ++i) {
        if (d[i] >= 3) q = min(q, 2.0L * inv[i]);
        for (int j = i + 1; j < 3; ++j)
            q = min(q, sqrtl(inv[i] * inv[i] + inv[j] * inv[j]));
    }
    // Each actual lattice distance is q*(1-2r), and it must be at least 2r.
    return q / (2.0L * (1.0L + q));
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long n;
    if (!(cin >> n)) return 0;
    cout << setprecision(21);

    // The two opposite body-diagonal corners are the exact general optimum
    // for two congruent spheres in a cube, and avoid a small-case regression.
    if (n == 2) {
        long double r = sqrtl(3.0L) / (2.0L * (1.0L + sqrtl(3.0L)));
        cout << r << ' ' << r << ' ' << r << '\n';
        cout << 1.0L-r << ' ' << 1.0L-r << ' ' << 1.0L-r << '\n';
        return 0;
    }

    // Three vertices of this parity tetrahedron already have the same
    // separation as all four, so it is an exact small-n improvement over a
    // cubic grid for both n=3 and n=4.
    if (n == 3 || n == 4) {
        long double q = sqrtl(2.0L);
        long double r = q / (2.0L * (1.0L + q));
        long double lo = r, hi = 1.0L - r;
        cout << lo << ' ' << lo << ' ' << lo << '\n';
        cout << lo << ' ' << hi << ' ' << hi << '\n';
        cout << hi << ' ' << lo << ' ' << hi << '\n';
        if (n == 4) cout << hi << ' ' << hi << ' ' << lo << '\n';
        return 0;
    }

    // Baseline cubic grid.
    long long D = 1;
    while (D * D * D < n) ++D;
    long double best_r = 1.0L / (2.0L * D);
    Box best;

    // A rectangular finite FCC box contains ceil(abc/2) points.  It suffices
    // to examine sorted, non-dominated boxes: for fixed a,b the least useful
    // third dimension is the first one that supplies enough parity points.
    const long long need = 2 * n - 1;
    long long amax = (long long)ceill(cbrtl((long double)need));
    for (long long a = 2; a <= amax; ++a) {
        long long bmax = (long long)ceill(sqrtl((long double)need / a));
        bmax = max(bmax, a);
        for (long long b = a; b <= bmax; ++b) {
            long long c = max(b, ceil_div(need, a * b));
            long double r = certified_radius(a, b, c);
            if (r > best_r) {
                best_r = r;
                best = {a, b, c, r};
            }
        }
    }

    if (best.r < 0) {
        long long made = 0;
        for (long long i = 0; i < D && made < n; ++i)
            for (long long j = 0; j < D && made < n; ++j)
                for (long long k = 0; k < D && made < n; ++k, ++made)
                    cout << (i + .5L) / D << ' ' << (j + .5L) / D << ' '
                         << (k + .5L) / D << '\n';
    } else {
        long double h[3] = {
            (1.0L - 2.0L * best.r) / (best.a - 1),
            (1.0L - 2.0L * best.r) / (best.b - 1),
            (1.0L - 2.0L * best.r) / (best.c - 1)
        };
        long long made = 0;
        for (long long i = 0; i < best.a && made < n; ++i)
            for (long long j = 0; j < best.b && made < n; ++j)
                for (long long k = 0; k < best.c && made < n; ++k) {
                    if ((i + j + k) & 1) continue;
                    cout << best.r + i * h[0] << ' ' << best.r + j * h[1] << ' '
                         << best.r + k * h[2] << '\n';
                    ++made;
                }
    }
    return 0;
}
