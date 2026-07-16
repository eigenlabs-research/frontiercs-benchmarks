#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; unsigned long long key; int id; };

// Hilbert index in a 2^31 by 2^31 square.  Unlike an x-sort, consecutive
// Hilbert indices tend to be geographically close in both coordinates.
static unsigned long long hilbert(unsigned long long ux, unsigned long long uy) {
    long long x = (long long)ux, y = (long long)uy;
    unsigned long long d = 0;
    for (long long s = 1LL << 30; s; s >>= 1) {
        long long rx = (x & s) != 0, ry = (y & s) != 0;
        d += (unsigned long long)s * s * ((3 * rx) ^ ry);
        // This is the standard Hilbert quadrant rotation; its square has
        // side s, not the full coordinate-domain side.
        if (!ry) {
            // rot() reflects in the full 2^31 square, not the current
            // bit-sized quadrant.  Reflecting about s corrupts coordinates
            // after a lower-half rotation and destroys Hilbert locality.
            if (rx) { x = (1LL << 31) - 1 - x; y = (1LL << 31) - 1 - y; }
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
    vector<Point> a(n);
    long long minx = LLONG_MAX, miny = LLONG_MAX;
    for (int i = 0; i < n; ++i) {
        cin >> a[i].x >> a[i].y;
        a[i].id = i;
        minx = min(minx, a[i].x); miny = min(miny, a[i].y);
    }
    for (auto &v : a)
        v.key = hilbert((unsigned long long)(v.x - minx), (unsigned long long)(v.y - miny));
    sort(a.begin(), a.end(), [](const Point& u, const Point& v) {
        return u.key != v.key ? u.key < v.key : u.id < v.id;
    });

    vector<char> prime(n, true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    // Rotate the spatial ordering at city zero, so zero has its two Hilbert
    // neighbors rather than being artificially connected to an endpoint.
    int z = 0;
    while (a[z].id != 0) ++z;
    vector<int> p(n + 1);
    p[0] = p[n] = 0;
    for (int k = 1; k < n; ++k) p[k] = a[(z + k) % n].id;

    // Coordinates are indexed by city ID below; keep them separate from the
    // Hilbert-sorted array.
    vector<long long> X(n), Y(n);
    for (const auto &v : a) X[v.id] = v.x, Y[v.id] = v.y;
    auto edge = [&](const vector<int>& r, int t, int sa = -1, int sb = -1) {
        auto at = [&](int i) { return i == sa ? r[sb] : (i == sb ? r[sa] : r[i]); };
        int u = at(t - 1), v = at(t);
        double w = (t % 10 == 0 && !prime[u]) ? 1.1 : 1.0;
        return w * hypot((double)X[u] - X[v], (double)Y[u] - Y[v]);
    };
    auto swapGain = [&](int i, int j) {
        int ts[4] = {i, i + 1, j, j + 1};
        sort(ts, ts + 4);
        double before = 0, after = 0;
        for (int q = 0; q < 4; ++q) if ((q == 0 || ts[q] != ts[q-1]) && ts[q] >= 1 && ts[q] <= n) {
            before += edge(p, ts[q]);
            after += edge(p, ts[q], i, j);
        }
        return before - after;
    };
    auto improveSwap = [&](int i, int j) {
        if (i < 1 || j < 1 || i >= n || j >= n || i == j) return false;
        if (swapGain(i, j) > 1e-7) { swap(p[i], p[j]); return true; }
        return false;
    };

    // The only weighted positions are sources at 9,19,... .  Seek a nearby
    // prime for each one, accepting solely on the exact full objective.
    for (int pass = 0; pass < 2; ++pass) {
        for (int pos = 9; pos < n; pos += 10) if (!prime[p[pos]]) {
            int best = -1; double gain = 1e-7;
            for (int q = max(1, pos - 8); q <= min(n - 1, pos + 8); ++q) if (prime[p[q]]) {
                double g = swapGain(pos, q);
                if (g > gain) gain = g, best = q;
            }
            if (best != -1) swap(p[pos], p[best]);
        }
    }
    // A cheap exact local cleanup also fixes small Hilbert discontinuities.
    for (int pass = 0; pass < 2; ++pass)
        for (int i = 1; i + 1 < n; ++i) improveSwap(i, i + 1);

    // Test the missing local move from the incumbent: a 2-opt reversal.  Its
    // ordinary internal edge lengths cancel under reversal; only its two
    // boundary edges and the (few) tenth-step carrot surcharges must be
    // reconsidered.  Thus this is exact, rather than a geometric surrogate.
    auto dcity = [&](int u, int v) {
        return hypot((double)X[u] - X[v], (double)Y[u] - Y[v]);
    };
    auto reversalGain = [&](int l, int r) {
        double before = dcity(p[l - 1], p[l]) + dcity(p[r], p[r + 1]);
        double after = dcity(p[l - 1], p[r]) + dcity(p[l], p[r + 1]);
        // At a route index k in [l,r], the reversed city is old p[l+r-k].
        auto newAt = [&](int k) { return (k < l || k > r) ? p[k] : p[l + r - k]; };
        for (int t = ((l + 9) / 10) * 10; t <= r + 1; t += 10) {
            int ou = p[t - 1], ov = p[t];
            int nu = newAt(t - 1), nv = newAt(t);
            if (!prime[ou]) before += 0.1 * dcity(ou, ov);
            if (!prime[nu]) after += 0.1 * dcity(nu, nv);
        }
        return before - after;
    };
    const int extraLens[] = {40, 48, 64};
    for (int l = 1; l < n; ++l) {
        int bestR = -1;
        double bestGain = 1e-7;
        for (int len = 2; len <= 32 && l + len - 1 < n; ++len) {
            double g = reversalGain(l, l + len - 1);
            if (g > bestGain) bestGain = g, bestR = l + len - 1;
        }
        for (int len : extraLens) if (l + len - 1 < n) {
            double g = reversalGain(l, l + len - 1);
            if (g > bestGain) bestGain = g, bestR = l + len - 1;
        }
        if (bestR != -1) reverse(p.begin() + l, p.begin() + bestR + 1);
    }

    cout << n + 1 << '\n';
    for (int v : p) cout << v << '\n';
}
