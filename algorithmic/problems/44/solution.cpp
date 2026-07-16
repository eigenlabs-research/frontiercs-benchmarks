#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

static unsigned long long hilbert(unsigned x, unsigned y) {
    unsigned long long d = 0;
    for (unsigned s = 1u << 20; s; s >>= 1) {
        unsigned rx = (x & s) != 0, ry = (y & s) != 0;
        d += (unsigned long long)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = (s - 1) - x; y = (s - 1) - y; }
            swap(x, y);
        }
    }
    return d;
}

static unsigned long long morton(unsigned x, unsigned y) {
    unsigned long long z = 0;
    for (int b = 0; b < 21; ++b) {
        z |= (unsigned long long)((x >> b) & 1) << (2 * b);
        z |= (unsigned long long)((y >> b) & 1) << (2 * b + 1);
    }
    return z;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> p(n);
    long long xmin = LLONG_MAX, xmax = LLONG_MIN, ymin = LLONG_MAX, ymax = LLONG_MIN;
    for (auto &q : p) {
        cin >> q.x >> q.y;
        xmin = min(xmin, q.x); xmax = max(xmax, q.x);
        ymin = min(ymin, q.y); ymax = max(ymax, q.y);
    }
    vector<char> prime(max(2, n), true);
    prime[0] = false; prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    const unsigned LIM = (1u << 21) - 1;
    auto scale = [&](long long v, long long lo, long long hi) -> unsigned {
        if (hi == lo) return 0;
        long double q = (long double)(v - lo) * LIM / (long double)(hi - lo);
        if (q < 0) return 0;
        if (q > LIM) return LIM;
        return (unsigned)q;
    };
    vector<unsigned> nx(n), ny(n);
    for (int i = 0; i < n; ++i) {
        nx[i] = scale(p[i].x, xmin, xmax);
        ny[i] = scale(p[i].y, ymin, ymax);
    }
    auto edge = [&](int a, int b) {
        return hypot((double)p[a].x - p[b].x, (double)p[a].y - p[b].y);
    };
    auto cost = [&](const vector<int>& r) {
        long double ans = 0;
        for (int t = 1; t <= n; ++t) {
            double w = (t % 10 == 0 && !prime[r[t - 1]]) ? 1.1 : 1.0;
            ans += w * edge(r[t - 1], r[t]);
        }
        return ans;
    };

    vector<int> best;
    long double bestCost = numeric_limits<long double>::infinity();
    // Both space-filling curve families and all dihedral coordinate orientations.
    for (int family = 0; family < 2; ++family) for (int mode = 0; mode < 8; ++mode) {
        vector<pair<unsigned long long,int>> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i) {
            unsigned a = nx[i], b = ny[i];
            if (mode & 1) a = LIM - a;
            if (mode & 2) b = LIM - b;
            if (mode & 4) swap(a, b);
            order.push_back({family == 0 ? hilbert(a, b) : morton(a, b), i});
        }
        sort(order.begin(), order.end());
        int z = 0;
        while (order[z].second != 0) ++z;
        for (int backwards = 0; backwards < 2; ++backwards) {
            vector<int> r; r.reserve(n + 1); r.push_back(0);
            for (int k = 1; k < n; ++k) {
                int at = backwards ? (z - k + n) % n : (z + k) % n;
                r.push_back(order[at].second);
            }
            r.push_back(0);
            long double c = cost(r);
            if (c < bestCost) { bestCost = c; best.swap(r); }
        }
    }

    // Exact bounded 2-opt.  A reversal leaves every non-tenth internal edge
    // unchanged (only its direction changes), so its delta needs its two
    // boundary edges and the internal carrot steps, rather than all len edges.
    // This makes a substantially wider neighborhood affordable.
    for (int pass = 0; pass < 2; ++pass) {
        for (int l = 1; l < n - 1; ++l) {
            for (int len = 2; len <= 32 && l + len - 1 <= n - 1; ++len) {
                int r = l + len - 1;
                long double oldPart = 0, newPart = 0;
                auto weight = [&](int t, int source) {
                    return (t % 10 == 0 && !prime[source]) ? 1.1 : 1.0;
                };
                // The two boundary edges are always changed.
                oldPart += weight(l, best[l - 1]) * edge(best[l - 1], best[l]);
                newPart += weight(l, best[l - 1]) * edge(best[l - 1], best[r]);
                oldPart += weight(r + 1, best[r]) * edge(best[r], best[r + 1]);
                newPart += weight(r + 1, best[l]) * edge(best[l], best[r + 1]);
                // Interior distances are the same edges in reverse order; only
                // carrot multipliers at global tenth steps may differ.
                int first = ((l + 10) / 10) * 10;
                for (int t = first; t <= r; t += 10) {
                    int q = l + r - t + 1;
                    oldPart += weight(t, best[t - 1]) * edge(best[t - 1], best[t]);
                    newPart += weight(t, best[q]) * edge(best[q], best[q - 1]);
                }
                if (newPart + 1e-6L < oldPart) {
                    reverse(best.begin() + l, best.begin() + r + 1);
                    bestCost += newPart - oldPart;
                }
            }
        }
    }
    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
}
