#include <bits/stdc++.h>
using namespace std;

struct City { long long x, y; unsigned long long key; int id; };

// Hilbert index of a point in the 31-bit square.  Nearby points tend to have
// nearby indices, unlike the input's one-dimensional x ordering.
static unsigned long long hilbert(unsigned int x, unsigned int y) {
    unsigned long long d = 0;
    for (unsigned int s = 1u << 30; s; s >>= 1) {
        unsigned int rx = (x & s) != 0, ry = (y & s) != 0;
        d += (unsigned long long)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = (s - 1) - x; y = (s - 1) - y; }
            swap(x, y);
        }
    }
    return d;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<City> a(n);
    for (int i = 0; i < n; ++i) {
        cin >> a[i].x >> a[i].y;
        // Coordinates are in [-1e9,1e9], safely inside a 31-bit square.
        a[i].key = hilbert((unsigned int)(a[i].x + 1000000000LL),
                           (unsigned int)(a[i].y + 1000000000LL));
        a[i].id = i;
    }
    vector<char> prime(n, true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int p = 2; 1LL * p * p < n; ++p)
        if (prime[p]) for (int q = p * p; q < n; q += p) prime[q] = false;

    sort(a.begin(), a.end(), [](const City& u, const City& v) {
        return u.key != v.key ? u.key < v.key : u.id < v.id;
    });
    int z = 0;
    while (a[z].id != 0) ++z;
    vector<int> route(n + 1);
    route[0] = route[n] = 0;
    for (int k = 1; k < n; ++k) route[k] = a[(z + k) % n].id;

    // IDs are not sorted after the Hilbert sort; retain coordinate lookup.
    vector<long long> xs(n), ys(n);
    for (const auto &c : a) xs[c.id] = c.x, ys[c.id] = c.y;
    auto distance = [&](int u, int v) {
        return hypot((double)xs[u] - xs[v], (double)ys[u] - ys[v]);
    };
    auto weighted = [&](int step, int u, int v) {
        double w = (step % 10 == 0 && !prime[u]) ? 1.1 : 1.0;
        return w * distance(u, v);
    };

    // A bounded 2-opt pass repairs the occasional short-range crossing in the
    // curve order.  Every accepted move is checked with the actual step costs.
    const int WINDOW = 10;
    for (int pass = 0; pass < 1; ++pass) {
        for (int i = 1; i < n - 1; ++i) {
            int best = -1;
            double bestGain = 0;
            for (int j = i + 1; j < n && j <= i + WINDOW; ++j) {
                double oldEnds = distance(route[i-1], route[i]) + distance(route[j], route[j+1]);
                double newEnds = distance(route[i-1], route[j]) + distance(route[i], route[j+1]);
                if (oldEnds - newEnds > bestGain) bestGain = oldEnds - newEnds, best = j;
            }
            if (best < 0) continue;
            double oldCost = 0, newCost = 0;
            for (int t = i; t <= best + 1; ++t) {
                oldCost += weighted(t, route[t-1], route[t]);
                int uPos = t - 1, vPos = t;
                int u = (uPos >= i && uPos <= best) ? route[i + best - uPos] : route[uPos];
                int v = (vPos >= i && vPos <= best) ? route[i + best - vPos] : route[vPos];
                newCost += weighted(t, u, v);
            }
            if (newCost + 1e-9 < oldCost) reverse(route.begin() + i, route.begin() + best + 1);
        }
    }
    cout << n + 1 << '\n';
    for (int v : route) cout << v << '\n';
}
