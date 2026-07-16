#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert transpose algorithm for two 31-bit coordinates.
static uint64_t hilbert(uint32_t x, uint32_t y) {
    uint32_t a[2] = {x, y};
    const uint32_t M = 1u << 30;
    for (uint32_t q = M; q > 1; q >>= 1) {
        uint32_t p = q - 1;
        for (int i = 1; i < 2; ++i) {
            if (a[i] & q) a[0] ^= p;
            else {
                uint32_t t = (a[0] ^ a[i]) & p;
                a[0] ^= t; a[i] ^= t;
            }
        }
    }
    a[1] ^= a[0];
    uint32_t t = 0;
    for (uint32_t q = M; q > 1; q >>= 1) if (a[1] & q) t ^= q - 1;
    a[0] ^= t; a[1] ^= t;
    uint64_t ans = 0;
    for (int b = 30; b >= 0; --b) {
        ans = (ans << 1) | ((a[0] >> b) & 1u);
        ans = (ans << 1) | ((a[1] >> b) & 1u);
    }
    return ans;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> pt(n);
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &p : pt) {
        cin >> p.x >> p.y;
        minx = min(minx, p.x); maxx = max(maxx, p.x);
        miny = min(miny, p.y); maxy = max(maxy, p.y);
    }
    vector<char> prime(max(2, n), true);
    prime[0] = false; prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    auto edge = [&](const vector<int>& r, int i) -> long double {
        int j = (i + 1 == n ? 0 : i + 1);
        const Point &a = pt[r[i]], &b = pt[r[j]];
        long double dx = (long double)a.x - b.x, dy = (long double)a.y - b.y;
        long double v = sqrtl(dx * dx + dy * dy);
        if ((i + 1) % 10 == 0 && !prime[r[i]]) v *= 1.1L;
        return v;
    };
    auto cost = [&](const vector<int>& r) {
        long double z = 0;
        for (int i = 0; i < n; ++i) z += edge(r, i);
        return z;
    };
    auto improveCarrots = [&](vector<int>& r) {
        // A swap affects only four directed edges, so each proposed carrot repair is exact.
        for (int pos = 9; pos < n; pos += 10) {
            if (prime[r[pos]]) continue;
            int best = -1;
            long double bestDelta = 0;
            int lo = max(1, pos - 40), hi = min(n - 1, pos + 40);
            for (int q = lo; q <= hi; ++q) {
                if (!prime[r[q]] || ((q + 1) % 10 == 0)) continue;
                int ids[4] = {(pos - 1 + n) % n, pos, (q - 1 + n) % n, q};
                sort(ids, ids + 4);
                int m = int(unique(ids, ids + 4) - ids);
                long double before = 0, after = 0;
                for (int k = 0; k < m; ++k) before += edge(r, ids[k]);
                swap(r[pos], r[q]);
                for (int k = 0; k < m; ++k) after += edge(r, ids[k]);
                swap(r[pos], r[q]);
                if (after - before < bestDelta) bestDelta = after - before, best = q;
            }
            if (best != -1) swap(r[pos], r[best]);
        }
    };

    vector<int> best(n);
    iota(best.begin(), best.end(), 0); // always retain a valid baseline candidate
    long double bestCost = cost(best);
    const uint64_t LIM = (1ULL << 31) - 1;
    auto scaled = [&](long long v, long long mn, long long mx) -> uint32_t {
        if (mn == mx) return uint32_t(LIM / 2);
        __int128 num = (__int128)(v - mn) * LIM;
        return uint32_t(num / (mx - mn));
    };
    for (int fx = 0; fx < 2; ++fx) for (int fy = 0; fy < 2; ++fy) {
        vector<pair<uint64_t,int>> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i) {
            uint64_t xx = scaled(pt[i].x, minx, maxx), yy = scaled(pt[i].y, miny, maxy);
            if (fx) xx = LIM - xx;
            if (fy) yy = LIM - yy;
            order.push_back({hilbert(uint32_t(xx), uint32_t(yy)), i});
        }
        sort(order.begin(), order.end());
        for (int rev = 0; rev < 2; ++rev) {
            vector<int> r; r.reserve(n);
            int at = 0;
            while (order[at].second != 0) ++at;
            r.push_back(0);
            for (int k = 1; k < n; ++k) {
                int z = rev ? (at - k + n) % n : (at + k) % n;
                r.push_back(order[z].second);
            }
            long double c = cost(r);
            if (c < bestCost) bestCost = c, best.swap(r);
        }
    }
    improveCarrots(best);
    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
    cout << 0 << '\n';
}
