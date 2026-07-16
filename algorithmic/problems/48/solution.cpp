#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    cout << fixed << setprecision(15);

    // For a few balls, binary cube codes give much better boundary behavior
    // than starting a bulk lattice.
    if (n <= 8) {
        vector<int> code;
        int minHam;
        if (n == 2) {
            code = {0, 7};
            minHam = 3;
        } else if (n <= 4) {
            // The even-parity tetrahedron has mutual Hamming distance two.
            code = {0, 3, 5, 6};
            minHam = 2;
        } else {
            for (int v = 0; v < 8; ++v) code.push_back(v);
            minHam = 1;
        }
        double d = sqrt((double)minHam);
        double r = d / (2.0 * (1.0 + d));
        for (int i = 0; i < n; ++i) {
            int v = code[i];
            cout << ((v & 1) ? 1.0-r : r) << ' '
                 << ((v & 2) ? 1.0-r : r) << ' '
                 << ((v & 4) ? 1.0-r : r) << '\n';
        }
        return 0;
    }

    // FCC points are integer triples with even coordinate sum.  Their
    // closest spacing is sqrt(2) lattice units.  A q-by-q-by-q bounding box
    // therefore fits with face margin r=1/(2+sqrt(2)q).
    int q = 0;
    while (true) {
        long long s = q + 1;
        long long cap = (s*s*s + ((s & 1) ? 1 : 0)) / 2;
        if (cap >= n) break;
        ++q;
    }
    const double r = 1.0 / (2.0 + sqrt(2.0) * q);
    const double step = sqrt(2.0) * r;
    int printed = 0;
    for (int x = 0; x <= q && printed < n; ++x)
        for (int y = 0; y <= q && printed < n; ++y)
            for (int z = 0; z <= q && printed < n; ++z)
                if (((x + y + z) & 1) == 0) {
                    cout << r + step*x << ' ' << r + step*y << ' ' << r + step*z << '\n';
                    ++printed;
                }
    return 0;
}
