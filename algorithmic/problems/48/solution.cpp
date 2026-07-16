#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long n;
    if (!(cin >> n)) return 0;
    cout << setprecision(17);

    // The two antipodal body-diagonal points are better than the first
    // two points of a lattice cell.
    if (n == 2) {
        double r = sqrt(3.0) / (2.0 * (1.0 + sqrt(3.0)));
        cout << r << ' ' << r << ' ' << r << '\n';
        cout << 1-r << ' ' << 1-r << ' ' << 1-r << '\n';
        return 0;
    }

    // For 5 through 8, a subset of the eight corners of the inset cube
    // has radius 1/4; it bridges the small-size lattice count jump.
    if (n >= 5 && n <= 8) {
        for (long long q = 0; q < n; ++q) {
            cout << ((q & 1) ? .75 : .25) << ' '
                 << ((q & 2) ? .75 : .25) << ' '
                 << ((q & 4) ? .75 : .25) << '\n';
        }
        return 0;
    }

    // FCC lattice: integer triples with even coordinate sum, scaled so
    // nearest neighbours are distance d.  The M+1 by M+1 by M+1 block
    // is exactly contained in the cube after a radius d/2 inset.
    long long M = 0;
    while (true) {
        long long L = M + 1;
        long long e = (L + 1) / 2, o = L / 2;
        __int128 cnt = (__int128)e * e * e + (__int128)3 * e * o * o;
        if (cnt >= n) break;
        ++M;
    }
    double d = 1.0 / (1.0 + (double)M / sqrt(2.0));
    double r = d * 0.5;
    double h = d / sqrt(2.0);
    long long made = 0;
    for (long long i = 0; i <= M && made < n; ++i)
        for (long long j = 0; j <= M && made < n; ++j)
            for (long long k = 0; k <= M && made < n; ++k)
                if (((i + j + k) & 1) == 0) {
                    cout << r + h*i << ' ' << r + h*j << ' ' << r + h*k << '\n';
                    ++made;
                }
    return 0;
}
