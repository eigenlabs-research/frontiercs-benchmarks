#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

static uint64_t hilbert(uint32_t x, uint32_t y) {
    // The usual xy-to-Hilbert conversion on a 2^21 by 2^21 square.
    uint64_t d = 0;
    for (uint32_t s = 1u << 20; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = (1u << 21) - 1 - x; y = (1u << 21) - 1 - y; }
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
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &p : a) {
        cin >> p.x >> p.y;
        minx = min(minx, p.x); maxx = max(maxx, p.x);
        miny = min(miny, p.y); maxy = max(maxy, p.y);
    }
    vector<char> prime(n, true);
    if (n) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    auto edge = [&](int u, int v, int step) -> double {
        double dx = (double)a[u].x - (double)a[v].x;
        double dy = (double)a[u].y - (double)a[v].y;
        double z = hypot(dx, dy);
        return (step % 10 == 0 && !prime[u]) ? z * 1.1 : z;
    };
    auto score = [&](const vector<int>& p) {
        double ans = 0;
        for (int t = 1; t <= n; ++t) ans += edge(p[t-1], p[t], t);
        return ans;
    };
    auto fromOrder = [&](const vector<int>& q) {
        vector<int> p; p.reserve(n + 1); p.push_back(0);
        for (int v : q) if (v != 0) p.push_back(v);
        p.push_back(0);
        return p;
    };
    // Interpret q as a cyclic spatial order and cut that cycle at the depot.
    auto cutAtZero = [&](const vector<int>& q) {
        vector<int> p; p.reserve(n + 1); p.push_back(0);
        int z = (int)(find(q.begin(), q.end(), 0) - q.begin());
        for (int d = 1; d < n; ++d) p.push_back(q[(z + d) % n]);
        p.push_back(0);
        return p;
    };

    vector<int> best, order;
    order.reserve(max(0, n - 1));
    for (int i = 1; i < n; ++i) order.push_back(i);
    best = fromOrder(order);
    double bestScore = score(best);
    auto consider = [&](vector<int>&& p) {
        double s = score(p);
        if (s < bestScore) { bestScore = s; best = move(p); }
    };

    // Horizontal sweeps are retained as a robust alternative for layered inputs.
    const int strips[] = {2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128};
    for (int b : strips) {
        vector<vector<int>> bucket(b);
        long long span = maxy - miny;
        for (int i = 0; i < n; ++i) {
            int k = span ? (int)((__int128)(a[i].y - miny) * b / (span + 1)) : 0;
            if (k >= b) k = b - 1;
            bucket[k].push_back(i); // input IDs are already increasing in x
        }
        vector<int> q; q.reserve(n);
        for (int k = 0; k < b; ++k) {
            if (k & 1) for (auto it = bucket[k].rbegin(); it != bucket[k].rend(); ++it) q.push_back(*it);
            else q.insert(q.end(), bucket[k].begin(), bucket[k].end());
        }
        // This is the incumbent strip cycle; cutting it at city 0 keeps its local links.
        consider(cutAtZero(q));
        // Also allow the depot to act as an open-path endpoint rather than a cycle cut.
        vector<int> withoutZero; withoutZero.reserve(n - 1);
        for (int v : q) if (v) withoutZero.push_back(v);
        consider(fromOrder(withoutZero));
    }

    // Unlike a strip sweep, Hilbert order recursively keeps both x and y nearby.
    const uint32_t LIM = (1u << 21) - 1;
    vector<pair<uint64_t,int>> keyed; keyed.reserve(n);
    long long sx = maxx - minx, sy = maxy - miny;
    for (int i = 0; i < n; ++i) {
        uint32_t x = sx ? (uint32_t)((__int128)(a[i].x - minx) * LIM / sx) : 0;
        uint32_t y = sy ? (uint32_t)((__int128)(a[i].y - miny) * LIM / sy) : 0;
        keyed.push_back({hilbert(x, y), i});
    }
    sort(keyed.begin(), keyed.end());
    int zeroAt = 0;
    for (int i = 0; i < n; ++i) if (keyed[i].second == 0) { zeroAt = i; break; }
    vector<int> h; h.reserve(n - 1);
    for (int d = 1; d < n; ++d) h.push_back(keyed[(zeroAt + d) % n].second);
    consider(fromOrder(h));
    h.clear();
    for (int d = 1; d < n; ++d) h.push_back(keyed[(zeroAt - d + n) % n].second);
    consider(fromOrder(h));

    // Exact local cleanup: only the three position-weighted edges touched by a swap change.
    for (int pass = 0; pass < 4; ++pass) {
        bool changed = false;
        int lo = (pass & 1) ? n - 2 : 1;
        int hi = (pass & 1) ? 0 : n - 1;
        int inc = (pass & 1) ? -1 : 1;
        for (int i = lo; i != hi; i += inc) {
            double before = edge(best[i-1], best[i], i) + edge(best[i], best[i+1], i+1) + edge(best[i+1], best[i+2], i+2);
            double after  = edge(best[i-1], best[i+1], i) + edge(best[i+1], best[i], i+1) + edge(best[i], best[i+2], i+2);
            if (after < before) { swap(best[i], best[i+1]); changed = true; }
        }
        if (!changed) break;
    }
    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
}
