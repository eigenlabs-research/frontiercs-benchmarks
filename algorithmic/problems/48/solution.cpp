#include <bits/stdc++.h>
using namespace std;

// A cubic grid is useful for the few small instances where the boundary cost of
// an FCC cell outweighs its better bulk density.  Otherwise use the even
// sublattice of Z^3, whose nearest-neighbor distance is sqrt(2) lattice units.
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    cout << setprecision(17);

    // The two antipodal points are better than either finite lattice choice.
    if (n == 2) {
        const long double r = sqrtl(3.0L) / (2.0L * (1.0L + sqrtl(3.0L)));
        cout << (double)r << ' ' << (double)r << ' ' << (double)r << '\n';
        cout << (double)(1-r) << ' ' << (double)(1-r) << ' ' << (double)(1-r) << '\n';
        return 0;
    }

    int p = 1;
    while (1LL * p * p * p < n) ++p;
    long double gridR = 1.0L / (2.0L * p);

    int L = 0;
    auto fccCount = [](long long a) {
        // Among [0,a]^3, the even-parity class has ceil((a+1)^3/2) points.
        long long q = a + 1;
        return (q*q*q + ((a & 1) ? 0 : 1)) / 2;
    };
    while (fccCount(L) < n) ++L;
    long double fccR = 1.0L / (sqrtl(2.0L) * L + 2.0L);

    if (gridR >= fccR) {
        // Centers of p equal intervals in each coordinate.  Taking any prefix
        // cannot reduce its guaranteed distance to a cube face.
        int made = 0;
        for (int x = 0; x < p && made < n; ++x)
            for (int y = 0; y < p && made < n; ++y)
                for (int z = 0; z < p && made < n; ++z, ++made)
                    cout << (double)((x + .5L) / p) << ' '
                         << (double)((y + .5L) / p) << ' '
                         << (double)((z + .5L) / p) << '\n';
    } else {
        // Spacing s makes the FCC nearest-neighbor half-distance s/sqrt(2).
        // The L lattice-unit span plus two such margins exactly fills a side.
        long double s = 1.0L / (L + sqrtl(2.0L));
        long double margin = s / sqrtl(2.0L);
        int made = 0;
        for (int x = 0; x <= L && made < n; ++x)
            for (int y = 0; y <= L && made < n; ++y)
                for (int z = 0; z <= L && made < n; ++z) {
                    if (((x + y + z) & 1) == 0) {
                        cout << (double)(margin + s*x) << ' '
                             << (double)(margin + s*y) << ' '
                             << (double)(margin + s*z) << '\n';
                        ++made;
                    }
                }
    }
    return 0;
}
