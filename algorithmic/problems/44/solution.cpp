#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> a(n);
    for (auto &p : a) cin >> p.x >> p.y;

    vector<char> prime(n, true);
    if (n > 0) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    auto distance = [&](int u, int v) {
        long double dx = (long double)a[u].x - a[v].x;
        long double dy = (long double)a[u].y - a[v].y;
        return sqrtl(dx * dx + dy * dy);
    };
    auto routeScore = [&](const vector<int>& p) {
        long double s = 0;
        for (int t = 1; t <= n; ++t) {
            long double d = distance(p[t - 1], p[t]);
            if (t % 10 == 0 && !prime[p[t - 1]]) d *= 1.1L;
            s += d;
        }
        return s;
    };

    // This is an exact local improvement: changes are kept only when their
    // complete effect on the weighted objective is beneficial.
    auto improveCarrots = [&](vector<int>& p) {
        auto edge = [&](int t) {
            long double d = distance(p[t - 1], p[t]);
            if (t % 10 == 0 && !prime[p[t - 1]]) d *= 1.1L;
            return d;
        };
        for (int pos = 9; pos < n; pos += 10) {
            if (prime[p[pos]]) continue;
            int best = -1;
            long double bestDelta = 0;
            // This runs for every spatial candidate; a short forward window
            // keeps the exact local search within the global time budget.
            int hi = min(n - 1, pos + 12);
            for (int q = pos + 1; q <= hi; ++q) if (prime[p[q]]) {
                int ts[4] = {pos, pos + 1, q, q + 1};
                sort(ts, ts + 4);
                long double before = 0;
                for (int k = 0; k < 4; ++k)
                    if (k == 0 || ts[k] != ts[k - 1]) before += edge(ts[k]);
                swap(p[pos], p[q]);
                long double after = 0;
                for (int k = 0; k < 4; ++k)
                    if (k == 0 || ts[k] != ts[k - 1]) after += edge(ts[k]);
                swap(p[pos], p[q]);
                if (after - before < bestDelta) bestDelta = after - before, best = q;
            }
            if (best != -1) swap(p[pos], p[best]);
        }
    };

    long long minx = a[0].x, maxx = a[0].x, miny = a[0].y, maxy = a[0].y;
    for (auto p : a) {
        minx = min(minx, p.x); maxx = max(maxx, p.x);
        miny = min(miny, p.y); maxy = max(maxy, p.y);
    }
    auto scaled = [](long long v, long long lo, long long hi) -> uint32_t {
        if (hi == lo) return 0;
        return (uint32_t)(((__int128)(v - lo) * ((1ULL << 31) - 1)) / (hi - lo));
    };
    auto morton = [&](int id, bool flip) {
        uint32_t x = scaled(a[id].x, minx, maxx), y = scaled(a[id].y, miny, maxy);
        if (flip) swap(x, y);
        uint64_t z = 0;
        for (int b = 30; b >= 0; --b) {
            z = (z << 1) | ((x >> b) & 1U);
            z = (z << 1) | ((y >> b) & 1U);
        }
        return z;
    };

    // Unlike Z-order, a Hilbert curve keeps successive recursive cells adjacent.
    // These variants distinguish whether the incumbent's Morton discontinuities,
    // rather than its carrot adjustment, are the limiting factor.
    auto hilbert = [](uint32_t x, uint32_t y) -> uint64_t {
        const uint32_t mask = (1U << 31) - 1;
        uint64_t d = 0;
        for (uint32_t s = 1U << 30; s; s >>= 1) {
            uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
            d += (uint64_t)s * s * ((3 * rx) ^ ry);
            if (!ry) {
                if (rx) { x = mask - x; y = mask - y; }
                swap(x, y);
            }
        }
        return d;
    };
    auto sortByKey = [&](auto key) {
        vector<pair<uint64_t, int>> keyed;
        keyed.reserve(n);
        for (int id = 0; id < n; ++id) keyed.push_back({key(id), id});
        sort(keyed.begin(), keyed.end());
        vector<int> order;
        order.reserve(n);
        for (auto [_, id] : keyed) order.push_back(id);
        return order;
    };

    vector<int> best;
    long double bestScore = numeric_limits<long double>::infinity();
    auto consider = [&](vector<int> order) {
        order.erase(remove(order.begin(), order.end(), 0), order.end());
        // A spatial ordering is an open curve; choose its less expensive end
        // to attach to the mandatory depot before making carrot adjustments.
        if (order.size() > 1 && distance(0, order.back()) < distance(0, order.front()))
            reverse(order.begin(), order.end());
        vector<int> p;
        p.reserve(n + 1);
        p.push_back(0);
        p.insert(p.end(), order.begin(), order.end());
        p.push_back(0);
        improveCarrots(p);
        long double s = routeScore(p);
        if (s < bestScore) bestScore = s, best.swap(p);
    };

    vector<int> ids(n);
    iota(ids.begin(), ids.end(), 0);
    consider(ids); // strengthened input-order baseline is always represented.

    ids = sortByKey([&](int id) { return morton(id, false); });
    consider(ids);
    ids = sortByKey([&](int id) { return morton(id, true); });
    consider(ids);
    ids = sortByKey([&](int id) {
        // x is already strictly ordered in the intended input, but retain an
        // explicit two-coordinate key for a deterministic fallback.
        uint64_t yy = (uint64_t)(a[id].y - miny);
        return (yy << 31) ^ (uint64_t)scaled(a[id].x, minx, maxx);
    });
    consider(ids);

    const uint32_t hmask = (1U << 31) - 1;
    ids = sortByKey([&](int id) {
        return hilbert(scaled(a[id].x, minx, maxx), scaled(a[id].y, miny, maxy));
    });
    consider(ids);
    ids = sortByKey([&](int id) {
        return hilbert(scaled(a[id].y, miny, maxy), scaled(a[id].x, minx, maxx));
    });
    consider(ids);
    ids = sortByKey([&](int id) {
        return hilbert(hmask - scaled(a[id].x, minx, maxx), scaled(a[id].y, miny, maxy));
    });
    consider(ids);
    ids = sortByKey([&](int id) {
        return hilbert(scaled(a[id].x, minx, maxx), hmask - scaled(a[id].y, miny, maxy));
    });
    consider(ids);

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
}
