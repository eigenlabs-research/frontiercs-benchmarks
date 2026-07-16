#include <bits/stdc++.h>
using namespace std;

enum Kind { CUBIC, FCC, BCC };

static long long fccCount(int L) {
    long long e = L / 2 + 1;       // even integers in [0,L]
    long long o = (L + 1) / 2;
    return e * e * e + 3 * e * o * o;
}
static long long bccCount(int L) {
    long long e = L / 2 + 1;
    long long o = (L + 1) / 2;
    return e * e * e + o * o * o;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // A conventional cubic grid is useful for very small cases.  The two
    // centered cubic lattices have better asymptotic density; all candidates
    // are scaled so their lattice nearest-neighbor distance equals 2r.
    int m = 1;
    while (1LL * m * m * m < n) ++m;
    double bestR = 1.0 / (2.0 * m);
    Kind kind = CUBIC;
    int span = m - 1;

    int lf = 0;
    while (fccCount(lf) < n) ++lf;
    double rf = 1.0 / (sqrt(2.0) * (lf + sqrt(2.0)));
    if (rf > bestR) {
        bestR = rf;
        kind = FCC;
        span = lf;
    }

    int lb = 0;
    while (bccCount(lb) < n) ++lb;
    double rb = sqrt(3.0) / (2.0 * (lb + sqrt(3.0)));
    if (rb > bestR) {
        bestR = rb;
        kind = BCC;
        span = lb;
    }

    cout << setprecision(17);
    int made = 0;
    if (kind == CUBIC) {
        double r = 1.0 / (2.0 * m), step = 1.0 / m;
        for (int x = 0; x < m && made < n; ++x)
            for (int y = 0; y < m && made < n; ++y)
                for (int z = 0; z < m && made < n; ++z, ++made)
                    cout << r + x * step << ' ' << r + y * step << ' ' << r + z * step << '\n';
    } else if (kind == FCC) {
        double step = 1.0 / (span + sqrt(2.0));
        double r = step / sqrt(2.0);
        for (int x = 0; x <= span && made < n; ++x)
            for (int y = 0; y <= span && made < n; ++y)
                for (int z = 0; z <= span && made < n; ++z)
                    if (((x + y + z) & 1) == 0) {
                        cout << r + x * step << ' ' << r + y * step << ' ' << r + z * step << '\n';
                        ++made;
                    }
    } else {
        double step = 1.0 / (span + sqrt(3.0));
        double r = sqrt(3.0) * step / 2.0;
        for (int x = 0; x <= span && made < n; ++x)
            for (int y = 0; y <= span && made < n; ++y)
                for (int z = 0; z <= span && made < n; ++z)
                    if ((x & 1) == (y & 1) && (y & 1) == (z & 1)) {
                        cout << r + x * step << ' ' << r + y * step << ' ' << r + z * step << '\n';
                        ++made;
                    }
    }
    return 0;
}
