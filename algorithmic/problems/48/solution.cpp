#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    cout << setprecision(17);

    // The best two-point member of this family is the cube diagonal.
    if (n == 2) {
        const double r = sqrt(3.0) / (2.0 * (1.0 + sqrt(3.0)));
        cout << r << ' ' << r << ' ' << r << '\n';
        cout << 1-r << ' ' << 1-r << ' ' << 1-r << '\n';
        return 0;
    }

    // Compare a conventional cubic grid with a finite FCC lattice.  In the
    // latter, integer triples of even coordinate sum have nearest-neighbour
    // distance sqrt(2)*a.  Giving them a margin sqrt(2)*a/2 makes both the
    // sphere contacts and the cube faces tight.
    int k = 1;
    while (1LL*k*k*k < n) ++k;
    double cubeR = 1.0 / (2.0*k);

    int m = 1;
    auto fccCount = [](long long s) { return (s*s*s + (s & 1)) / 2; };
    while (fccCount(m) < n) ++m;
    double a = 1.0 / (double(m - 1) + sqrt(2.0));
    double fccR = a / sqrt(2.0);

    if (cubeR >= fccR) {
        // A prefix of a grid remains a packing, so no assumptions on n are
        // needed beyond the number of requested points.
        int made = 0;
        for (int x = 0; x < k && made < n; ++x)
            for (int y = 0; y < k && made < n; ++y)
                for (int z = 0; z < k && made < n; ++z, ++made)
                    cout << (x + 0.5) / k << ' ' << (y + 0.5) / k << ' '
                         << (z + 0.5) / k << '\n';
    } else {
        int made = 0;
        for (int x = 0; x < m && made < n; ++x)
            for (int y = 0; y < m && made < n; ++y)
                for (int z = 0; z < m && made < n; ++z)
                    if (((x + y + z) & 1) == 0) {
                        double X = fccR + x*a, Y = fccR + y*a, Z = fccR + z*a;
                        cout << X << ' ' << Y << ' ' << Z << '\n';
                        ++made;
                    }
    }
    return 0;
}
