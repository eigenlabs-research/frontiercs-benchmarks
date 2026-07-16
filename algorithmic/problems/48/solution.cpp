#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    cout << setprecision(17);

    // The two antipodal corners are better than the smallest FCC cell.
    if (n == 2) {
        long double r = sqrtl(3.0L) / (2.0L * (1.0L + sqrtl(3.0L)));
        cout << (double)r << ' ' << (double)r << ' ' << (double)r << '\n';
        cout << (double)(1-r) << ' ' << (double)(1-r) << ' ' << (double)(1-r) << '\n';
        return 0;
    }

    // A balanced Cartesian grid is a safe fallback (and equals the stated baseline).
    int t = 1;
    while (1LL * t * t * t < n) ++t;
    int ga = t, gb = t, gc = (n + t * t - 1) / (t * t);
    long double gridR = 1.0L / (2.0L * t);

    // FCC points are integer triples of even parity.  Their nearest neighbors
    // differ in two coordinates, giving distance sqrt(2)*h.
    int q = 0;
    auto fccCount = [](long long side) -> long long {
        return (side * side * side + (side & 1LL)) / 2;
    };
    while (fccCount(q + 1) < n) ++q;
    long double fccR = 1.0L / (2.0L + sqrtl(2.0L) * q);

    if (fccR > gridR) {
        long double h = (1.0L - 2.0L * fccR) / q;
        int made = 0;
        for (int i = 0; i <= q && made < n; ++i)
            for (int j = 0; j <= q && made < n; ++j)
                for (int k = 0; k <= q && made < n; ++k)
                    if (((i + j + k) & 1) == 0) {
                        long double x = fccR + i * h;
                        long double y = fccR + j * h;
                        long double z = fccR + k * h;
                        cout << (double)x << ' ' << (double)y << ' ' << (double)z << '\n';
                        ++made;
                    }
    } else {
        int made = 0;
        for (int i = 0; i < ga && made < n; ++i)
            for (int j = 0; j < gb && made < n; ++j)
                for (int k = 0; k < gc && made < n; ++k) {
                    cout << (double)((i + 0.5L) / ga) << ' '
                         << (double)((j + 0.5L) / gb) << ' '
                         << (double)((k + 0.5L) / gc) << '\n';
                    ++made;
                }
    }
    return 0;
}
