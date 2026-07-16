#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    cout << setprecision(17);

    // The two-sphere case is the body diagonal packing.
    if (n == 2) {
        long double r = sqrtl(3.0L) / (2.0L * (1.0L + sqrtl(3.0L)));
        long double t = 1.0L - r;
        cout << (double)r << ' ' << (double)r << ' ' << (double)r << '\n';
        cout << (double)t << ' ' << (double)t << ' ' << (double)t << '\n';
        return 0;
    }

    // D3 (FCC) is the densest lattice in three dimensions.  Integer points
    // with even coordinate sum have nearest-neighbour distance sqrt(2).
    int M = 0;
    while (true) {
        long long side = M + 1;
        long long cap = (side * side * side + 1) / 2;
        if (cap >= n) break;
        ++M;
    }
    long double fccR = 1.0L / (2.0L + sqrtl(2.0L) * M);

    // Keep the universally available balanced cubic grid when the finite FCC
    // box is in one of its capacity gaps.
    int q = 1;
    while (1LL * q * q * q < n) ++q;
    long double gridR = 1.0L / (2.0L * q);

    if (fccR > gridR) {
        long double step = sqrtl(2.0L) * fccR;
        int made = 0;
        for (int x = 0; x <= M && made < n; ++x)
            for (int y = 0; y <= M && made < n; ++y)
                for (int z = 0; z <= M && made < n; ++z)
                    if (((x + y + z) & 1) == 0) {
                        long double X = fccR + step * x;
                        long double Y = fccR + step * y;
                        long double Z = fccR + step * z;
                        cout << (double)X << ' ' << (double)Y << ' ' << (double)Z << '\n';
                        ++made;
                    }
    } else {
        // A q by q by q grid has clearance and half nearest distance 1/(2q).
        int made = 0;
        for (int x = 0; x < q && made < n; ++x)
            for (int y = 0; y < q && made < n; ++y)
                for (int z = 0; z < q && made < n; ++z) {
                    long double X = ((long double)x + 0.5L) / q;
                    long double Y = ((long double)y + 0.5L) / q;
                    long double Z = ((long double)z + 0.5L) / q;
                    cout << (double)X << ' ' << (double)Y << ' ' << (double)Z << '\n';
                    ++made;
                }
    }
    return 0;
}
