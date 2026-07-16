#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert distance of a point in a 2^bits by 2^bits square.
static uint64_t hilbert(uint32_t x, uint32_t y, int bits) {
    uint64_t d = 0;
    for (uint32_t s = 1u << (bits - 1); s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * uint64_t(s) * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) {
                x = (uint32_t)((1u << bits) - 1) - x;
                y = (uint32_t)((1u << bits) - 1) - y;
            }
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
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &q : p) {
        cin >> q.x >> q.y;
        minx = min(minx, q.x); maxx = max(maxx, q.x);
        miny = min(miny, q.y); maxy = max(maxy, q.y);
    }
    vector<char> prime(n, true);
    if (n > 0) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    constexpr int B = 21;
    constexpr uint32_t LIM = (1u << B) - 1;
    vector<uint32_t> qx(n), qy(n);
    long double dx = maxx - minx, dy = maxy - miny;
    for (int i = 0; i < n; ++i) {
        qx[i] = dx == 0 ? 0 : (uint32_t)((long double)(p[i].x - minx) * LIM / dx);
        qy[i] = dy == 0 ? 0 : (uint32_t)((long double)(p[i].y - miny) * LIM / dy);
    }
    auto distance = [&](int a, int b) {
        return hypot((double)p[a].x - p[b].x, (double)p[a].y - p[b].y);
    };
    auto cost = [&](const vector<int>& r) {
        long double ans = 0;
        for (int t = 1; t <= n; ++t) {
            double e = distance(r[t - 1], r[t]);
            if (t % 10 == 0 && !prime[r[t - 1]]) e *= 1.1;
            ans += e;
        }
        return ans;
    };
    vector<int> best;
    long double bestCost = numeric_limits<long double>::infinity();
    auto consider = [&](vector<int> ord) {
        int z = find(ord.begin(), ord.end(), 0) - ord.begin();
        vector<int> route; route.reserve(n + 1); route.push_back(0);
        for (int k = 1; k < n; ++k) route.push_back(ord[(z + k) % n]);
        route.push_back(0);
        long double v = cost(route);
        if (v < bestCost) bestCost = v, best.swap(route);
    };

    // All square symmetries choose favorable Hilbert endpoints and traversal direction.
    for (int mode = 0; mode < 8; ++mode) {
        vector<pair<uint64_t,int>> keyed; keyed.reserve(n);
        for (int i = 0; i < n; ++i) {
            uint32_t a = qx[i], b = qy[i];
            if (mode & 1) a = LIM - a;
            if (mode & 2) b = LIM - b;
            if (mode & 4) swap(a, b);
            keyed.push_back({hilbert(a, b, B), i});
        }
        sort(keyed.begin(), keyed.end());
        vector<int> ord; ord.reserve(n);
        for (auto [key, id] : keyed) ord.push_back(id);
        consider(ord);
        reverse(ord.begin(), ord.end());
        consider(ord);
    }
    // The supplied monotone order is also a safe candidate for nearly one-dimensional data.
    vector<int> monotone(n); iota(monotone.begin(), monotone.end(), 0);
    consider(monotone);
    reverse(monotone.begin(), monotone.end());
    consider(monotone);

    // Hilbert gives a good global tour, but its joins can still cross locally.  Improve
    // only short reversals and score them with the actual position-dependent objective.
    // This is deliberately bounded: it is a safe O(N) refinement after the O(N log N)
    // construction even at the largest input size.
    auto edgeCost = [&](const vector<int>& r, int t) -> long double {
        int a = r[t - 1], b = r[t];
        // Match the checker exactly: it uses double-precision hypot.
        double d = hypot((double)p[a].x - p[b].x, (double)p[a].y - p[b].y);
        if (t % 10 == 0 && !prime[a]) d *= 1.1;
        return d;
    };
    // Reverse positions [l,r].  Only steps l through r+1 change, including
    // the carrot multiplier because the source city at those positions changes.
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int l = 1; l < n - 1; ++l) {
            int last = min(n - 1, l + 6);
            for (int r = l + 1; r <= last; ++r) {
                long double before = 0;
                for (int t = l; t <= r + 1; ++t) before += edgeCost(best, t);
                reverse(best.begin() + l, best.begin() + r + 1);
                long double after = 0;
                for (int t = l; t <= r + 1; ++t) after += edgeCost(best, t);
                if (after + 1e-9L < before) {
                    changed = true;
                } else {
                    reverse(best.begin() + l, best.begin() + r + 1);
                }
            }
        }
        if (!changed) break;
    }

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
    return 0;
}
