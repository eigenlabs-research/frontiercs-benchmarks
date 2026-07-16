#include <bits/stdc++.h>
using namespace std;

struct P { double x, y, z; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // The simple grid is retained for sizes where the boundary of a finite
    // FCC block costs more than its better bulk density.
    int g = 1;
    while (1LL * g * g * g < n) ++g;
    int gz = (n + g * g - 1) / (g * g);
    vector<P> grid;
    grid.reserve(n);
    for (int i = 0; i < g && (int)grid.size() < n; ++i)
        for (int j = 0; j < g && (int)grid.size() < n; ++j)
            for (int k = 0; k < gz && (int)grid.size() < n; ++k)
                grid.push_back({(i + .5) / g, (j + .5) / g, (k + .5) / gz});
    double gridR = 1.0 / (2.0 * g);

    // D3/FCC is the integer lattice with even coordinate sum.  Its nearest
    // neighbours have distance sqrt(2).  A uniform scale and a face margin
    // make both the lattice contacts and cube contacts attain the same r.
    int q = 1;
    auto fccCapacity = [](int a) -> long long {
        long long v = 1LL * a * a * a;
        return (v + 1) / 2;
    };
    while (fccCapacity(q) < n) ++q;
    struct I { int x, y, z; };
    vector<I> raw;
    raw.reserve(n);
    // Sorting by coordinate sum starts with a nearest-neighbour pair, while
    // remaining symmetric among the three cube axes.
    vector<I> all;
    all.reserve((size_t)fccCapacity(q));
    for (int x = 0; x < q; ++x)
        for (int y = 0; y < q; ++y)
            for (int z = 0; z < q; ++z)
                if (((x + y + z) & 1) == 0) all.push_back({x, y, z});
    sort(all.begin(), all.end(), [](const I& a, const I& b) {
        int sa = a.x + a.y + a.z, sb = b.x + b.y + b.z;
        if (sa != sb) return sa < sb;
        int qa = a.x*a.x + a.y*a.y + a.z*a.z;
        int qb = b.x*b.x + b.y*b.y + b.z*b.z;
        if (qa != qb) return qa < qb;
        if (a.x != b.x) return a.x < b.x;
        if (a.y != b.y) return a.y < b.y;
        return a.z < b.z;
    });
    raw.assign(all.begin(), all.begin() + n);
    int lx = raw[0].x, hx = lx, ly = raw[0].y, hy = ly, lz = raw[0].z, hz = lz;
    for (const I& p : raw) {
        lx = min(lx, p.x); hx = max(hx, p.x);
        ly = min(ly, p.y); hy = max(hy, p.y);
        lz = min(lz, p.z); hz = max(hz, p.z);
    }
    int span = max({hx - lx, hy - ly, hz - lz});
    const double near = sqrt(2.0);
    double scale = 1.0 / (span + near);
    double fccR = near * scale * .5;
    vector<P> fcc;
    fcc.reserve(n);
    for (const I& p : raw)
        fcc.push_back({fccR + scale * (p.x - lx),
                       fccR + scale * (p.y - ly),
                       fccR + scale * (p.z - lz)});

    const vector<P>& ans = (fccR > gridR) ? fcc : grid;
    cout << setprecision(17);
    for (const P& p : ans) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
    return 0;
}
