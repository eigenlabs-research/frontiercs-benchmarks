#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // The integer points with i+j+k even are an FCC lattice.  Its shortest
    // vector has length sqrt(2).  Crop a nearly cubic block of that lattice.
    int q = 1;
    while (((1LL * q * q * q + 1) / 2) < n) ++q;
    vector<array<int, 3>> fcc;
    fcc.reserve(n);
    for (int i = 0; i < q && (int)fcc.size() < n; ++i)
        for (int j = 0; j < q && (int)fcc.size() < n; ++j)
            for (int k = 0; k < q && (int)fcc.size() < n; ++k)
                if (((i + j + k) & 1) == 0) fcc.push_back({i, j, k});

    int lo[3] = {q, q, q}, hi[3] = {-1, -1, -1};
    for (auto p : fcc) for (int d = 0; d < 3; ++d) {
        lo[d] = min(lo[d], p[d]); hi[d] = max(hi[d], p[d]);
    }
    int span = max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
    const double rt2 = sqrt(2.0);
    double fccRadius = rt2 / (2.0 * (span + rt2));

    // This is also a certified fallback matching the stated grid construction.
    int m = 1;
    while (1LL * m * m * m < n) ++m;
    double gridRadius = 1.0 / (2.0 * m);

    cout << setprecision(17);
    if (fccRadius > gridRadius) {
        double scale = 1.0 / (span + rt2);
        double margin = 0.5 * rt2 * scale;
        for (auto p : fcc) {
            double x = margin + scale * (p[0] - lo[0]);
            double y = margin + scale * (p[1] - lo[1]);
            double z = margin + scale * (p[2] - lo[2]);
            cout << x << ' ' << y << ' ' << z << '\n';
        }
    } else {
        int made = 0;
        for (int i = 0; i < m && made < n; ++i)
            for (int j = 0; j < m && made < n; ++j)
                for (int k = 0; k < m && made < n; ++k, ++made)
                    cout << (i + .5) / m << ' ' << (j + .5) / m << ' ' << (k + .5) / m << '\n';
    }
    return 0;
}
