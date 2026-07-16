#include <bits/stdc++.h>
using namespace std;

// All constructions below have face margin d/2 and pair distance at least d.
// type: 0 = cubic, 1 = FCC, 2 = BCC
struct Candidate {
    double d;
    int m, type, phase;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    long long n;
    if (!(cin >> n)) return 0;

    int mc = 1;
    while (1LL * mc * mc * mc < n) ++mc;
    Candidate best{1.0 / mc, mc, 0, 0};

    int mf = 1;
    while ((1LL * mf * mf * mf + 1) / 2 < n) ++mf;
    Candidate f{1.0 / (1.0 + (mf - 1) / sqrt(2.0)), mf, 1, 0};
    if (f.d > best.d) best = f;

    int mb = 1, bp = 0;
    for (;; ++mb) {
        long long cap = 0;
        int bestPhase = 0;
        for (int p = 0; p < 8; ++p) {
            long long a[3][2];
            for (int ax = 0; ax < 3; ++ax) {
                a[ax][0] = (mb + 1) / 2; // indices 0,2,...
                a[ax][1] = mb / 2;
            }
            long long c1 = 1, c2 = 1;
            for (int ax = 0; ax < 3; ++ax) {
                int bit = (p >> ax) & 1;
                c1 *= a[ax][bit];
                c2 *= a[ax][bit ^ 1];
            }
            if (c1 + c2 > cap) cap = c1 + c2, bestPhase = p;
        }
        if (cap >= n) { bp = bestPhase; break; }
    }
    Candidate b{1.0 / (1.0 + (mb - 1) / sqrt(3.0)), mb, 2, bp};
    if (b.d > best.d) best = b;

    cout << setprecision(17);
    long long made = 0;
    const double q = best.type == 0 ? best.d :
                     best.d / sqrt(best.type == 1 ? 2.0 : 3.0);
    const double margin = best.d * 0.5;
    for (int i = 0; i < best.m && made < n; ++i)
        for (int j = 0; j < best.m && made < n; ++j)
            for (int k = 0; k < best.m && made < n; ++k) {
                bool use;
                if (best.type == 0) use = true;
                else if (best.type == 1) use = ((i + j + k) & 1) == 0;
                else {
                    int bits = (i & 1) | ((j & 1) << 1) | ((k & 1) << 2);
                    use = (bits == best.phase || bits == (best.phase ^ 7));
                }
                if (use) {
                    cout << margin + q * i << ' ' << margin + q * j << ' '
                         << margin + q * k << '\n';
                    ++made;
                }
            }
    return 0;
}
