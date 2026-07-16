#include <bits/stdc++.h>
using namespace std;

struct Point { double x, y, z; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;

    // An FCC lattice has the best known bulk density.  Make a moderately
    // oversampled finite piece of it, then use farthest-point insertion to
    // choose exactly the requested number of sites (which also handles small
    // n, where the cube corners are particularly useful).
    int g = (int)ceil(cbrt(8.0 * n)) + 1;
    vector<Point> cand;
    cand.reserve(5 * n + 32);
    for (int i = 0; i <= g; ++i)
        for (int j = 0; j <= g; ++j)
            for (int k = 0; k <= g; ++k)
                if (((i + j + k) & 1) == 0)
                    cand.push_back({double(i) / g, double(j) / g, double(k) / g});
    // The other parity's corners complete the useful boundary configuration.
    for (int x : {0, 1}) for (int y : {0, 1}) for (int z : {0, 1})
        cand.push_back({double(x), double(y), double(z)});

    vector<double> nearest(cand.size(), numeric_limits<double>::infinity());
    vector<Point> chosen;
    chosen.reserve(n);
    Point first{0, 0, 0};
    chosen.push_back(first);
    for (size_t i = 0; i < cand.size(); ++i) {
        double dx = cand[i].x, dy = cand[i].y, dz = cand[i].z;
        nearest[i] = dx * dx + dy * dy + dz * dz;
    }

    double minD2 = numeric_limits<double>::infinity();
    for (int take = 1; take < n; ++take) {
        size_t best = 0;
        for (size_t i = 1; i < cand.size(); ++i)
            if (nearest[i] > nearest[best]) best = i;
        Point p = cand[best];
        minD2 = min(minD2, nearest[best]);
        chosen.push_back(p);
        for (size_t i = 0; i < cand.size(); ++i) {
            double dx = cand[i].x - p.x;
            double dy = cand[i].y - p.y;
            double dz = cand[i].z - p.z;
            double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < nearest[i]) nearest[i] = d2;
        }
    }

    double d = sqrt(minD2);
    // Mapping p to r + (1-2r)p puts every raw boundary point at distance r
    // from a face.  Equality of pair and face constraints gives this r.
    double r = d / (2.0 * (1.0 + d));
    r *= (1.0 - 1e-13); // avoid losing a touching constraint to output rounding
    double scale = 1.0 - 2.0 * r;
    cout << setprecision(17);
    for (const Point &p : chosen)
        cout << r + scale * p.x << ' ' << r + scale * p.y << ' ' << r + scale * p.z << '\n';
    return 0;
}
