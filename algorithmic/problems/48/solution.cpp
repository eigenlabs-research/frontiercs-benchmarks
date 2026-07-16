#include <bits/stdc++.h>
using namespace std;

struct P { int x, y, z; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // The integer points with even x+y+z are an FCC lattice; its nearest
    // neighbor distance is sqrt(2).  Use the smallest cubic block holding n.
    int side = 0;
    while ((long long)(side + 1) * (side + 1) * (side + 1) / 2 < n) ++side;
    // Odd cubes have one more even-parity point than half their volume.
    while (true) {
        long long total = 0;
        for (int x = 0; x <= side; ++x)
            for (int y = 0; y <= side; ++y)
                for (int z = 0; z <= side; ++z)
                    if (((x + y + z) & 1) == 0) ++total;
        if (total >= n) break;
        ++side;
    }

    vector<P> fcc;
    fcc.reserve(n);
    for (int x = 0; x <= side && (int)fcc.size() < n; ++x)
        for (int y = 0; y <= side && (int)fcc.size() < n; ++y)
            for (int z = 0; z <= side && (int)fcc.size() < n; ++z)
                if (((x + y + z) & 1) == 0) fcc.push_back({x, y, z});

    int lo[3] = {fcc[0].x, fcc[0].y, fcc[0].z};
    int hi[3] = {lo[0], lo[1], lo[2]};
    for (const P &p : fcc) {
        lo[0] = min(lo[0], p.x); lo[1] = min(lo[1], p.y); lo[2] = min(lo[2], p.z);
        hi[0] = max(hi[0], p.x); hi[1] = max(hi[1], p.y); hi[2] = max(hi[2], p.z);
    }
    int span = max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
    const double root2 = sqrt(2.0);
    // s*sqrt(2)/2 equals the face margin (1-s*span)/2.
    double fccRadius = root2 / (2.0 * (span + root2));

    int q = 1;
    while ((long long)q * q * q < n) ++q;
    double gridRadius = 1.0 / (2.0 * q);

    cout << setprecision(17);
    if (fccRadius > gridRadius) {
        double scale = 1.0 / (span + root2);
        double off[3];
        for (int d = 0; d < 3; ++d)
            off[d] = (1.0 - scale * (hi[d] - lo[d])) * 0.5 - scale * lo[d];
        for (const P &p : fcc)
            cout << off[0] + scale * p.x << ' ' << off[1] + scale * p.y << ' ' << off[2] + scale * p.z << '\n';
    } else {
        // A centered q-by-q-by-q grid is the stated baseline construction.
        int made = 0;
        for (int x = 0; x < q && made < n; ++x)
            for (int y = 0; y < q && made < n; ++y)
                for (int z = 0; z < q && made < n; ++z, ++made)
                    cout << (x + 0.5) / q << ' ' << (y + 0.5) / q << ' ' << (z + 0.5) / q << '\n';
    }
    return 0;
}
