#include <bits/stdc++.h>
using namespace std;

struct Point { double x, y, z; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    vector<Point> ans;
    if (n == 2) {
        // The two opposite body-diagonal positions are better than a grid.
        const double s3 = sqrt(3.0);
        const double r = s3 / (2.0 * (1.0 + s3));
        ans = {{r,r,r}, {1-r,1-r,1-r}};
    } else {
        // FCC points are integer triples of even coordinate sum, scaled so
        // their nearest-neighbour distance is d and their face margin is d/2.
        int m = 0;
        while (((long long)(m + 1) * (m + 1) * (m + 1) + 1) / 2 < n) ++m;
        const double d = 1.0 / (1.0 + (double)m / sqrt(2.0));
        const double fccRadius = d * 0.5;

        int q = 1;
        while ((long long)q * q * q < n) ++q;
        const double gridRadius = 0.5 / q;

        if (fccRadius >= gridRadius) {
            const double a = d / sqrt(2.0);
            ans.reserve(n);
            for (int i = 0; i <= m && (int)ans.size() < n; ++i)
                for (int j = 0; j <= m && (int)ans.size() < n; ++j)
                    for (int k = 0; k <= m && (int)ans.size() < n; ++k)
                        if (((i + j + k) & 1) == 0) {
                            // Spell endpoints using the margin, avoiding a
                            // possible last-bit excursion beyond the cube.
                            double x = (i == m ? 1.0 - fccRadius : fccRadius + a*i);
                            double y = (j == m ? 1.0 - fccRadius : fccRadius + a*j);
                            double z = (k == m ? 1.0 - fccRadius : fccRadius + a*k);
                            ans.push_back({x, y, z});
                        }
        } else {
            ans.reserve(n);
            for (int i = 0; i < q && (int)ans.size() < n; ++i)
                for (int j = 0; j < q && (int)ans.size() < n; ++j)
                    for (int k = 0; k < q && (int)ans.size() < n; ++k)
                        ans.push_back({(i + 0.5) / q, (j + 0.5) / q, (k + 0.5) / q});
        }
    }

    cout << setprecision(17);
    for (const auto &p : ans) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
    return 0;
}
