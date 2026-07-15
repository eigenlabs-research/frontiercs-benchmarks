#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

static uint64_t hilbert(uint32_t x, uint32_t y, int pw = 31, int rot = 0) {
    if (!pw) return 0;
    uint32_t h = 1u << (pw - 1);
    int seg = (x < h) ? ((y < h) ? 0 : 3) : ((y < h) ? 1 : 2);
    seg = (seg + rot) & 3;
    static const int delta[4] = {3, 0, 0, 1};
    uint32_t nx = x & (h - 1), ny = y & (h - 1);
    uint64_t sub = 1ULL << (2 * pw - 2);
    uint64_t add = hilbert(nx, ny, pw - 1, (rot + delta[seg]) & 3);
    return uint64_t(seg) * sub + ((seg == 1 || seg == 2) ? add : sub - add - 1);
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
    if (n) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    auto edge = [&](const vector<int>& r, int pos) { // edge r[pos] -> r[pos+1]
        const Point &a = p[r[pos]], &b = p[r[pos + 1]];
        double v = hypot((double)a.x - b.x, (double)a.y - b.y);
        if ((pos + 1) % 10 == 0 && !prime[r[pos]]) v *= 1.1;
        return v;
    };
    auto score = [&](const vector<int>& r) {
        double ans = 0;
        for (int i = 0; i < n; ++i) ans += edge(r, i);
        return ans;
    };
    // A cyclic spatial ordering can be cut at city zero without losing its locality.
    auto cutAtZero = [&](const vector<int>& cyc, bool backwards) {
        int z = find(cyc.begin(), cyc.end(), 0) - cyc.begin();
        vector<int> r; r.reserve(n + 1); r.push_back(0);
        for (int k = 1; k < n; ++k) {
            int at = backwards ? (z - k + n) % n : (z + k) % n;
            r.push_back(cyc[at]);
        }
        r.push_back(0);
        return r;
    };

    vector<vector<int>> candidates;
    vector<int> byx(n); iota(byx.begin(), byx.end(), 0);
    candidates.push_back(cutAtZero(byx, false));
    candidates.push_back(cutAtZero(byx, true));

    vector<int> horder(n); iota(horder.begin(), horder.end(), 0);
    auto norm = [](long long v, long long lo, long long hi) -> uint32_t {
        if (hi == lo) return 0;
        long double q = (long double)(v - lo) / (long double)(hi - lo);
        return (uint32_t)(q * 2147483647.0L + .5L);
    };
    vector<uint64_t> hkey(n);
    for (int i = 0; i < n; ++i) hkey[i] = hilbert(norm(p[i].x,minx,maxx), norm(p[i].y,miny,maxy));
    sort(horder.begin(), horder.end(), [&](int a, int b) { return hkey[a] < hkey[b]; });
    candidates.push_back(cutAtZero(horder, false));
    candidates.push_back(cutAtZero(horder, true));

    // Horizontal strip tours catch the common 'alternating rows' geometry that x order misses.
    vector<int> yorder(n); iota(yorder.begin(), yorder.end(), 0);
    sort(yorder.begin(), yorder.end(), [&](int a, int b) { return p[a].y == p[b].y ? p[a].x < p[b].x : p[a].y < p[b].y; });
    for (int bands : {2, 4, 8, 16, 32}) {
        if (bands >= n) break;
        vector<int> cyc; cyc.reserve(n);
        for (int b = 0; b < bands; ++b) {
            int l = (long long)b * n / bands, rr = (long long)(b + 1) * n / bands;
            vector<int> strip(yorder.begin() + l, yorder.begin() + rr);
            sort(strip.begin(), strip.end(), [&](int a, int c) { return p[a].x < p[c].x; });
            if (b & 1) reverse(strip.begin(), strip.end());
            cyc.insert(cyc.end(), strip.begin(), strip.end());
        }
        candidates.push_back(cutAtZero(cyc, false));
        candidates.push_back(cutAtZero(cyc, true));
    }

    vector<pair<double,int>> ranks;
    for (int i = 0; i < (int)candidates.size(); ++i) ranks.push_back({score(candidates[i]), i});
    sort(ranks.begin(), ranks.end());
    vector<int> best = candidates[ranks[0].second];
    double bestScore = ranks[0].first;
    auto started = chrono::steady_clock::now();
    // Exact short 2-opt moves also account for the position-dependent carrot multiplier.
    for (int which = 0; which < min(2, (int)ranks.size()); ++which) {
        vector<int> r = candidates[ranks[which].second];
        for (int pass = 0; pass < 3; ++pass) {
            bool changed = false;
            for (int l = 1; l < n - 1; ++l) {
                if (chrono::duration<double>(chrono::steady_clock::now() - started).count() > 2.05) break;
                for (int d = 1; d <= 7 && l + d < n; ++d) {
                    int rr = l + d;
                    double old = 0; for (int e = l - 1; e <= rr; ++e) old += edge(r, e);
                    reverse(r.begin() + l, r.begin() + rr + 1);
                    double nw = 0; for (int e = l - 1; e <= rr; ++e) nw += edge(r, e);
                    if (nw + 1e-7 < old) changed = true;
                    else reverse(r.begin() + l, r.begin() + rr + 1);
                }
            }
            if (!changed) break;
        }
        double s = score(r);
        if (s < bestScore) bestScore = s, best.swap(r);
        if (chrono::duration<double>(chrono::steady_clock::now() - started).count() > 2.1) break;
    }
    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
}
