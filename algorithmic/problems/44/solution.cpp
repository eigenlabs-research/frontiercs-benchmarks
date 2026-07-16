#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert rank in a 2^30 square.  Nearby ranks tend to be geographically nearby.
static uint64_t hilbert(int64_t x, int64_t y) {
    uint64_t d = 0;
    for (int64_t s = 1LL << 29; s; s >>= 1) {
        int rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * uint64_t(s) * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = s - 1 - x; y = s - 1 - y; }
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
    long long lo = LLONG_MAX, hi = LLONG_MIN;
    for (auto &q : p) {
        cin >> q.x >> q.y;
        lo = min(lo, min(q.x, q.y));
        hi = max(hi, max(q.x, q.y));
    }

    const uint32_t M = (1u << 30) - 1;
    vector<uint32_t> px(n), py(n);
    long double range = (long double)hi - lo;
    for (int i = 0; i < n; ++i) {
        if (range != 0) {
            px[i] = uint32_t((long double)(p[i].x - lo) * M / range);
            py[i] = uint32_t((long double)(p[i].y - lo) * M / range);
        }
    }
    vector<char> prime(n, true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    auto dist = [&](int a, int b) {
        return hypot((double)p[a].x - p[b].x, (double)p[a].y - p[b].y);
    };
    auto multiplier = [&](int step, int source) {
        return (step % 10 == 0 && !prime[source]) ? 1.1 : 1.0;
    };
    // r includes both fixed endpoint copies of city zero.
    auto cost = [&](const vector<int>& r) {
        double ans = 0;
        for (int e = 1; e <= n; ++e)
            ans += multiplier(e, r[e - 1]) * dist(r[e - 1], r[e]);
        return ans;
    };

    vector<int> best;
    double bestCost = numeric_limits<double>::infinity();
    for (int mode = 0; mode < 8; ++mode) {
        vector<pair<uint64_t, int>> key;
        key.reserve(n - 1);
        for (int id = 1; id < n; ++id) {
            uint32_t a = px[id], b = py[id];
            if (mode & 1) a = M - a;
            if (mode & 2) b = M - b;
            if (mode & 4) swap(a, b);
            key.push_back({hilbert(a, b), id});
        }
        sort(key.begin(), key.end());
        vector<int> r;
        r.reserve(n + 1);
        r.push_back(0);
        for (auto [_, id] : key) r.push_back(id);
        r.push_back(0);
        double v = cost(r);
        if (v < bestCost) { bestCost = v; best.swap(r); }
    }

    // Discriminating ablation: bounded 2-opt, evaluated against the exact carrot
    // objective.  Only penalty positions inside a reversal need special handling;
    // all other internal edge lengths are unchanged by reversal.
    const int WINDOW = 24, PASSES = 3;
    for (int pass = 0; pass < PASSES; ++pass) {
        bool changed = false;
        for (int i = 0; i + 2 < n; ++i) {
            int bestJ = -1;
            double bestDelta = 0;
            int lim = min(n - 1, i + WINDOW);
            for (int j = i + 2; j <= lim; ++j) {
                // Reverses vertices best[i+1..j], changing boundary edges.
                double delta = 0;
                delta += multiplier(i + 1, best[i]) * (dist(best[i], best[j]) - dist(best[i], best[i + 1]));
                delta += multiplier(j + 1, best[i + 1]) * dist(best[i + 1], best[j + 1]);
                delta -= multiplier(j + 1, best[j]) * dist(best[j], best[j + 1]);
                // An internal edge moves from old step q to new step e. Account
                // separately for penalties gained at e and lost at q.
                int first = ((i + 2 + 9) / 10) * 10;
                for (int e = first; e <= j; e += 10) {
                    int q = j - e + i + 2; // new source is best[q]
                    delta += 0.1 * dist(best[q - 1], best[q]) * (!prime[best[q]]);
                }
                for (int q = first; q <= j; q += 10) {
                    delta -= 0.1 * dist(best[q - 1], best[q]) * (!prime[best[q - 1]]);
                }
                if (delta < bestDelta) bestDelta = delta, bestJ = j;
            }
            if (bestJ != -1) {
                reverse(best.begin() + i + 1, best.begin() + bestJ + 1);
                bestCost += bestDelta;
                changed = true;
            }
        }
        if (!changed) break;
    }

    cout << n + 1 << '\n';
    for (int id : best) cout << id << '\n';
    return 0;
}
