#include <bits/stdc++.h>
using namespace std;

// A finite lattice must leave a margin at the faces: the checker takes that
// margin into account, rather than using a radius supplied by the contestant.
enum Kind { SC, FCC, BCC };

static long long capacity(Kind k, int q) {
    long long a = (q + 1) / 2, b = q / 2;
    if (k == SC) return 1LL * q * q * q;
    if (k == FCC) return a*a*a + 3LL*a*b*b; // i+j+k even
    return a*a*a + b*b*b;                    // i,j,k have one common parity
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    const double rt2 = sqrt(2.0), rt3 = sqrt(3.0);
    double bestR = -1.0;
    int bestQ = 0;
    Kind best = SC;

    // q is the number of coordinate levels in each direction.  For each
    // lattice choose the smallest q whose finite cell contains n points.
    for (Kind k : {SC, FCC, BCC}) {
        int q = 1;
        while (capacity(k, q) < n) ++q;
        double r;
        if (k == SC) r = 1.0 / (2.0 * q);
        else if (k == FCC) r = 1.0 / (2.0 + rt2 * (q - 1));
        else r = rt3 / (2.0 * (q - 1 + rt3));
        if (r > bestR) bestR = r, bestQ = q, best = k;
    }

    // All listed coordinates lie in [r,1-r].  FCC uses integer triples of
    // even sum (nearest distance sqrt(2)*h); BCC uses triples of common parity
    // (nearest distance sqrt(3)*h).  Thus both nearest distance and face
    // margin certify the selected radius without an O(n^2) search.
    double h;
    if (best == SC) h = 2.0 * bestR;
    else if (best == FCC) h = rt2 * bestR;
    else h = 2.0 * bestR / rt3;

    int made = 0;
    cout << setprecision(17);
    for (int i = 0; i < bestQ && made < n; ++i)
        for (int j = 0; j < bestQ && made < n; ++j)
            for (int k = 0; k < bestQ && made < n; ++k) {
                bool take;
                if (best == SC) take = true;
                else if (best == FCC) take = ((i + j + k) & 1) == 0;
                else take = ((i & 1) == (j & 1) && (j & 1) == (k & 1));
                if (!take) continue;
                double x = bestR + h * i;
                double y = bestR + h * j;
                double z = bestR + h * k;
                cout << x << ' ' << y << ' ' << z << '\n';
                ++made;
            }
    return 0;
}
