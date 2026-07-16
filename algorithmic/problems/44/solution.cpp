#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

static uint64_t hilbert(uint32_t x, uint32_t y) {
    uint64_t d = 0;
    for (uint32_t s = 1u << 30; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * s * ((3 * rx) ^ ry);
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
    vector<Point> p(n);
    for (auto &q : p) cin >> q.x >> q.y;

    vector<char> prime(n, true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    // step is deliberately an argument: the checker charges the multiplier
    // from the source city of this global (not local) route step.
    auto edge = [&](int a, int b, int step) -> double {
        double dx = double(p[a].x) - double(p[b].x);
        double dy = double(p[a].y) - double(p[b].y);
        double z = sqrt(dx * dx + dy * dy);
        if (step % 10 == 0 && !prime[a]) z *= 1.1;
        return z;
    };
    auto cost = [&](const vector<int>& r) {
        double ans = 0;
        for (int t = 1; t <= n; ++t) ans += edge(r[t - 1], r[t], t);
        return ans;
    };

    vector<vector<int>> candidates;
    auto addCyclicOrder = [&](vector<pair<uint64_t, int>> v) {
        sort(v.begin(), v.end());
        int at = 0;
        while (v[at].second != 0) ++at;
        for (int dir : {1, -1}) {
            vector<int> r;
            r.reserve(n + 1);
            r.push_back(0);
            for (int k = 1; k < n; ++k) {
                int j = (at + dir * k) % n;
                if (j < 0) j += n;
                r.push_back(v[j].second);
            }
            r.push_back(0);
            candidates.push_back(move(r));
        }
    };

    vector<pair<uint64_t, int>> vx, vy, hxy, hyx, mxy, myx;
    vx.reserve(n); vy.reserve(n); hxy.reserve(n); hyx.reserve(n);
    mxy.reserve(n); myx.reserve(n);
    constexpr long long SHIFT = 1000000000LL;
    for (int i = 0; i < n; ++i) {
        uint32_t X = uint32_t(p[i].x + SHIFT), Y = uint32_t(p[i].y + SHIFT);
        vx.push_back({X, i}); vy.push_back({Y, i});
        hxy.push_back({hilbert(X, Y), i}); hyx.push_back({hilbert(Y, X), i});
        mxy.push_back({(uint64_t(X) << 31) | Y, i});
        myx.push_back({(uint64_t(Y) << 31) | X, i});
    }
    addCyclicOrder(move(vx)); addCyclicOrder(move(vy));
    addCyclicOrder(move(hxy)); addCyclicOrder(move(hyx));
    addCyclicOrder(move(mxy)); addCyclicOrder(move(myx));
    const int geometricCount = (int)candidates.size();

    // A raster is deliberately a separate family from the space-filling
    // curves.  Equal-height strips retain coherent rows even when many cities
    // have exactly the same y coordinate (where rank-quantile strips split a
    // row arbitrarily).  Alternating x direction avoids a return jump inside
    // every strip.
    long long ymin = p[0].y, ymax = p[0].y;
    for (const auto &q : p) ymin = min(ymin, q.y), ymax = max(ymax, q.y);
    auto addStripOrder = [&](int bands) {
        vector<pair<uint64_t, int>> v;
        v.reserve(n);
        unsigned long long span = (unsigned long long)(ymax - ymin) + 1;
        for (int i = 0; i < n; ++i) {
            int b = int((unsigned long long)(p[i].y - ymin) * bands / span);
            if (b >= bands) b = bands - 1;
            uint32_t X = uint32_t(p[i].x + SHIFT);
            uint32_t within = (b & 1) ? ((1u << 31) - 1 - X) : X;
            v.push_back({(uint64_t(b) << 31) | within, i});
        }
        addCyclicOrder(move(v));
    };
    addStripOrder(8);
    addStripOrder(32);

    auto bestIn = [&](int first, int last) {
        int best = first;
        double bestCost = cost(candidates[best]);
        for (int i = first + 1; i < last; ++i) {
            double v = cost(candidates[i]);
            if (v < bestCost) bestCost = v, best = i;
        }
        return best;
    };
    auto improve = [&](vector<int> route) {
        // Adjacent swaps are just length-two reversals. Testing all short
        // reversals is bounded exact 2-opt: only its <=9 changed steps need
        // recomputation, including their global carrot multipliers.
        for (int pass = 0; pass < 2; ++pass) {
            bool changed = false;
            for (int i = 1; i < n; ++i) {
                for (int len = 2; len <= 8 && i + len - 1 < n; ++len) {
                    int j = i + len - 1;
                    double before = 0, after = 0;
                    for (int t = i; t <= j + 1; ++t)
                        before += edge(route[t - 1], route[t], t);
                    after += edge(route[i - 1], route[j], i);
                    for (int k = j; k > i; --k)
                        after += edge(route[k], route[k - 1], i + 1 + j - k);
                    after += edge(route[i], route[j + 1], j + 1);
                    if (after + 1e-7 < before) {
                        reverse(route.begin() + i, route.begin() + j + 1);
                        changed = true;
                    }
                }
            }
            if (!changed) break;
        }
        return route;
    };

    // Raw family ranking is not enough: a raster may expose substantially
    // more local reversals. Refine its best representative independently.
    vector<int> route = improve(candidates[bestIn(0, geometricCount)]);
    vector<int> stripRoute = improve(candidates[bestIn(geometricCount, (int)candidates.size())]);
    if (cost(stripRoute) < cost(route)) route = move(stripRoute);

    cout << n + 1 << '\n';
    for (int v : route) cout << v << '\n';
}
