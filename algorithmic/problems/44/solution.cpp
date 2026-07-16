#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

static uint64_t hilbert(uint32_t x, uint32_t y) {
    const uint32_t side = 1u << 21;
    uint64_t d = 0;
    for (uint32_t s = side >> 1; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = side - 1 - x; y = side - 1 - y; }
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
    long long minx = LLONG_MAX, miny = LLONG_MAX, maxx = LLONG_MIN, maxy = LLONG_MIN;
    for (auto &q : a) {
        cin >> q.x >> q.y;
        minx = min(minx, q.x); maxx = max(maxx, q.x);
        miny = min(miny, q.y); maxy = max(maxy, q.y);
    }
    vector<char> prime(n, true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    auto edge = [&](int t, int u, int v) {
        double z = hypot(double(a[u].x) - a[v].x, double(a[u].y) - a[v].y);
        return (t % 10 == 0 && !prime[u]) ? 1.1 * z : z;
    };
    auto cost = [&](const vector<int>& p) {
        double ans = 0;
        for (int t = 1; t <= n; ++t) ans += edge(t, p[t - 1], p[t]);
        return ans;
    };

    // Incumbent guard: a square-normalized Hilbert ordering.
    long long span = max(maxx - minx, maxy - miny);
    const long long SCALE = (1LL << 21) - 1;
    vector<pair<uint64_t,int>> key;
    key.reserve(n - 1);
    for (int i = 1; i < n; ++i) {
        uint32_t xx = uint32_t((a[i].x - minx) * SCALE / span);
        uint32_t yy = uint32_t((a[i].y - miny) * SCALE / span);
        key.push_back({hilbert(xx, yy), i});
    }
    sort(key.begin(), key.end());
    vector<int> best(n + 1);
    best[0] = best[n] = 0;
    for (int i = 1; i < n; ++i) best[i] = key[i - 1].second;
    double bestCost = cost(best);
    reverse(best.begin() + 1, best.begin() + n);
    double reverseCost = cost(best);
    if (reverseCost < bestCost) bestCost = reverseCost;
    else reverse(best.begin() + 1, best.begin() + n);

    // Different representation: y-stripes are traversed as alternating long
    // x sweeps.  It is a multiresolution decomposition, not a point-by-point
    // nearest-neighbour search, and is especially useful for lanes and grids.
    auto consider = [&](const vector<int>& order) {
        vector<int> q(n + 1);
        q[0] = q[n] = 0;
        for (int i = 1; i < n; ++i) q[i] = order[i - 1];
        // Evaluate both orientations before a winning route can be swapped out
        // of q; otherwise a forward winner suppresses its reverse orientation.
        double c = cost(q);
        if (c < bestCost) { bestCost = c; best = q; }
        reverse(q.begin() + 1, q.begin() + n);
        c = cost(q);
        if (c < bestCost) { bestCost = c; best = q; }
    };
    // B=1 also cost-checks the strengthened monotone-x baseline.  Bucketing
    // preserves x order inside a stripe, so each resolution is linear time.
    const int levels[] = {1, 2, 4, 8, 16, 32, 64, 128, 256};
    long long yspan = maxy - miny;
    for (int wanted : levels) {
        int bands = min(wanted, max(1, n - 1));
        vector<vector<int>> bucket(bands);
        for (int id = 1; id < n; ++id) {
            int b = int((a[id].y - miny) * bands / (yspan + 1));
            if (b >= bands) b = bands - 1;
            bucket[b].push_back(id);
        }
        vector<int> order;
        order.reserve(n - 1);
        for (int b = 0; b < bands; ++b) {
            if ((b & 1) == 0) order.insert(order.end(), bucket[b].begin(), bucket[b].end());
            else order.insert(order.end(), bucket[b].rbegin(), bucket[b].rend());
        }
        consider(order);
    }
    // Complementary x-stripes: cities are supplied in x order, so stripes are
    // contiguous ID ranges and each is swept by y.  This was the missing half
    // of the advertised axis-aligned decomposition; it handles vertical lanes.
    long long xspan = maxx - minx;
    for (int wanted : levels) {
        int bands = min(wanted, max(1, n - 1));
        vector<vector<int>> bucket(bands);
        for (int id = 1; id < n; ++id) {
            int b = int((a[id].x - minx) * bands / (xspan + 1));
            if (b >= bands) b = bands - 1;
            bucket[b].push_back(id);
        }
        vector<int> order;
        order.reserve(n - 1);
        for (int b = 0; b < bands; ++b) {
            auto &v = bucket[b];
            sort(v.begin(), v.end(), [&](int u, int v) {
                if (a[u].y != a[v].y) return a[u].y < a[v].y;
                return u < v;
            });
            if ((b & 1) == 0) order.insert(order.end(), v.begin(), v.end());
            else order.insert(order.end(), v.rbegin(), v.rend());
        }
        consider(order);
    }

    // Exact checker-aligned short 2-opt refinement of whichever construction
    // won.  Reversal evaluates all affected position-dependent carrot edges.
    int passes = (n <= 50000 ? 2 : 1);
    for (int pass = 0; pass < passes; ++pass) {
        bool changed = false;
        for (int len = 2; len <= 8; ++len) for (int l = 1; l + len - 1 < n; ++l) {
            int r = l + len - 1;
            double oldCost = 0, newCost = 0;
            for (int t = l; t <= r + 1; ++t) {
                oldCost += edge(t, best[t - 1], best[t]);
                int sp = t - 1, dp = t;
                int ns = (l <= sp && sp <= r) ? l + r - sp : sp;
                int nd = (l <= dp && dp <= r) ? l + r - dp : dp;
                newCost += edge(t, best[ns], best[nd]);
            }
            if (newCost + 1e-9 < oldCost) {
                reverse(best.begin() + l, best.begin() + r + 1);
                changed = true;
            }
        }
        if (!changed) break;
    }
    // A short reversal cannot repair an interleaving between distant positions.
    // Use Hilbert-near city pairs as a small geometric candidate set, and find
    // the best full 2-opt reversal among it.  Its carrot delta is O(1): normal
    // edge lengths are unchanged inside a reversal, while the two kinds of
    // penalty contributions are range sums grouped by position residue.
    vector<pair<int,int>> candidates;
    const int nearCount = 12;
    // A Hilbert curve is local but has discontinuities at cell boundaries.
    // Retain its candidate set, then independently add short x/y-order
    // neighborhoods.  These are only proposals: every accepted reversal is
    // still selected by the exact carrot-aware delta below.
    const int axisNear = 2;
    candidates.reserve(size_t(n) * (nearCount + 2 * axisNear));
    for (int d = 1; d <= nearCount; ++d)
        for (int k = 0; k + d < (int)key.size(); ++k)
            candidates.push_back({key[k].second, key[k + d].second});
    for (int d = 1; d <= axisNear; ++d)
        for (int u = 1; u + d < n; ++u)
            candidates.push_back({u, u + d});
    vector<int> yOrder(n - 1);
    iota(yOrder.begin(), yOrder.end(), 1);
    sort(yOrder.begin(), yOrder.end(), [&](int u, int v) {
        if (a[u].y != a[v].y) return a[u].y < a[v].y;
        return u < v;
    });
    for (int d = 1; d <= axisNear; ++d)
        for (int k = 0; k + d < n - 1; ++k)
            candidates.push_back({yOrder[k], yOrder[k + d]});

    int globalMoves = (n <= 50000 ? 12 : 4);
    for (int round = 0; round < globalMoves; ++round) {
        vector<int> pos(n);
        for (int i = 0; i <= n; ++i) pos[best[i]] = i;
        // oldPref: edge k has source best[k].  revPref: after reversal of
        // edge k, its source is best[k+1].
        vector<double> oldPref[10], revPref[10];
        for (int z = 0; z < 10; ++z) {
            oldPref[z].assign(n + 1, 0.0);
            revPref[z].assign(n + 1, 0.0);
        }
        for (int k = 0; k < n; ++k) {
            double dd = hypot(double(a[best[k]].x) - a[best[k + 1]].x,
                              double(a[best[k]].y) - a[best[k + 1]].y);
            for (int z = 0; z < 10; ++z) {
                oldPref[z][k + 1] = oldPref[z][k];
                revPref[z][k + 1] = revPref[z][k];
            }
            int z = k % 10;
            if (!prime[best[k]]) oldPref[z][k + 1] += dd;
            if (!prime[best[k + 1]]) revPref[z][k + 1] += dd;
        }
        auto range = [](const vector<double>& v, int l, int r) {
            return l > r ? 0.0 : v[r + 1] - v[l];
        };
        double bestDelta = -1e-7;
        int takeI = -1, takeJ = -1;
        for (auto [u, v] : candidates) {
            int i = pos[u], j = pos[v];
            if (i > j) swap(i, j);
            if (i < 1 || j >= n || j <= i + 1) continue;
            auto dd = [&](int x, int y) {
                return hypot(double(a[x].x) - a[y].x, double(a[x].y) - a[y].y);
            };
            // Replace (i,i+1),(j,j+1), reversing positions i+1..j.
            double delta = dd(best[i], best[j]) + dd(best[i + 1], best[j + 1])
                         - dd(best[i], best[i + 1]) - dd(best[j], best[j + 1]);
            double oldExtra = range(oldPref[9], i, j);
            int residue = (i + j + 1) % 10;
            double newExtra = range(revPref[residue], i + 1, j - 1);
            int t1 = i + 1, t2 = j + 1;
            if (t1 % 10 == 0 && !prime[best[i]]) newExtra += dd(best[i], best[j]);
            if (t2 % 10 == 0 && !prime[best[i + 1]]) newExtra += dd(best[i + 1], best[j + 1]);
            delta += 0.1 * (newExtra - oldExtra);
            if (delta < bestDelta) {
                bestDelta = delta;
                takeI = i; takeJ = j;
            }
        }
        if (takeI < 0) break;
        reverse(best.begin() + takeI + 1, best.begin() + takeJ + 1);
    }

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
}
