#include <bits/stdc++.h>
using namespace std;

struct Point { long double x, y, z; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    long long n;
    if (!(cin >> n)) return 0;

    vector<Point> out;
    long double chosen = -1;

    // The ordinary grid is also a certified fallback for the small capacity
    // jumps of a finite FCC block.
    long long m = 1;
    while (m * m * m < n) ++m;
    long double gridR = 1.0L / (2.0L * m);

    // An FCC lattice is the integer lattice points whose coordinate sum is
    // even.  Its shortest vector is (1,1,0), so lattice step sqrt(2)*r
    // realizes a sphere separation of exactly 2r.
    long long q = 1;
    auto fccCount = [](long long s) -> long long {
        return (s & 1) ? (s * s * s + 1) / 2 : s * s * s / 2;
    };
    while (fccCount(q) < n) ++q;
    long double fccR = 1.0L / (2.0L + (q - 1) * sqrtl(2.0L));

    if (fccR > gridR) {
        chosen = fccR;
        long double step = sqrtl(2.0L) * chosen;
        for (long long i = 0; i < q && (long long)out.size() < n; ++i)
            for (long long j = 0; j < q && (long long)out.size() < n; ++j)
                for (long long k = 0; k < q && (long long)out.size() < n; ++k)
                    if (((i + j + k) & 1) == 0)
                        out.push_back({chosen + i * step, chosen + j * step, chosen + k * step});
    } else {
        chosen = gridR;
        for (long long i = 0; i < m && (long long)out.size() < n; ++i)
            for (long long j = 0; j < m && (long long)out.size() < n; ++j)
                for (long long k = 0; k < m && (long long)out.size() < n; ++k)
                    out.push_back({(2 * i + 1) * chosen, (2 * j + 1) * chosen, (2 * k + 1) * chosen});
    }

    // Exact small packings improve the two first nontrivial cases without
    // changing the construction used at general n.
    if (n == 2) {
        long double r = sqrtl(3.0L) / (2.0L * (1.0L + sqrtl(3.0L)));
        out = {{r,r,r}, {1-r,1-r,1-r}};
    } else if (n == 3 || n == 4) {
        long double r = sqrtl(2.0L) / (2.0L * (1.0L + sqrtl(2.0L)));
        vector<Point> tet = {{r,r,r}, {r,1-r,1-r}, {1-r,r,1-r}, {1-r,1-r,r}};
        out.assign(tet.begin(), tet.begin() + n);
    }

    cout << setprecision(17);
    for (const auto &p : out)
        cout << (double)p.x << ' ' << (double)p.y << ' ' << (double)p.z << '\n';
    return 0;
}
