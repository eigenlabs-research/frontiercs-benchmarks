#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };
static vector<Point> p;
static vector<char> primeId;

static inline double edgeCost(int step, int a, int b) {
    double dx = (double)p[a].x - p[b].x;
    double dy = (double)p[a].y - p[b].y;
    double d = sqrt(dx * dx + dy * dy);
    return (step % 10 == 0 && !primeId[a]) ? 1.1 * d : d;
}

static double routeCost(const vector<int>& r) {
    double ans = 0;
    for (int t = 1; t < (int)r.size(); ++t) ans += edgeCost(t, r[t-1], r[t]);
    return ans;
}

// A reversal changes the city which starts every internal edge, so the usual
// position-independent 2-opt formula is invalid under the carrot rule.  Keep
// the move short and score every affected edge exactly instead.
static void shortReverseImprove(vector<int>& r) {
    const int n = (int)r.size() - 1;
    const int maxLen = 9;
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int i = 1; i < n; ++i) {
            int bestJ = -1;
            double bestDelta = -1e-7;
            for (int j = i + 1; j < n && j - i + 1 <= maxLen; ++j) {
                double oldv = 0, newv = 0;
                // Steps i..j+1 are precisely the affected edges.
                for (int t = i; t <= j + 1; ++t) {
                    oldv += edgeCost(t, r[t-1], r[t]);
                    int na = (t - 1 < i || t - 1 > j) ? r[t-1] : r[i + j - (t - 1)];
                    int nb = (t < i || t > j) ? r[t] : r[i + j - t];
                    newv += edgeCost(t, na, nb);
                }
                double delta = newv - oldv;
                if (delta < bestDelta) bestDelta = delta, bestJ = j;
            }
            if (bestJ != -1) {
                reverse(r.begin() + i, r.begin() + bestJ + 1);
                changed = true;
            }
        }
        if (!changed) break;
    }
}

// Build a cycle by sweeping contiguous x slabs.  Within a slab, a y sweep is
// reversed on alternating slabs; this is a strip decomposition rather than a
// point-wise space-filling curve and is particularly useful for layered data.
static void addSerpentineCandidates(vector<vector<int>>& candidates, int n, int slabs) {
    vector<vector<int>> bucket(slabs);
    for (int b = 0; b < slabs; ++b) {
        int l = (long long)b * n / slabs;
        int r = (long long)(b + 1) * n / slabs;
        bucket[b].reserve(r - l);
        for (int i = l; i < r; ++i) bucket[b].push_back(i);
        sort(bucket[b].begin(), bucket[b].end(), [](int a, int b) {
            if (p[a].y != p[b].y) return p[a].y < p[b].y;
            return p[a].x < p[b].x;
        });
    }
    for (int direction : {1, -1}) for (int phase : {0, 1}) {
        vector<int> cycle; cycle.reserve(n);
        for (int q = 0; q < slabs; ++q) {
            int b = direction == 1 ? q : slabs - 1 - q;
            bool up = ((q + phase) & 1) == 0;
            if (up) cycle.insert(cycle.end(), bucket[b].begin(), bucket[b].end());
            else cycle.insert(cycle.end(), bucket[b].rbegin(), bucket[b].rend());
        }
        int at = (int)(find(cycle.begin(), cycle.end(), 0) - cycle.begin());
        vector<int> r; r.reserve(n + 1); r.push_back(0);
        for (int k = 1; k < n; ++k) r.push_back(cycle[(at + k) % n]);
        r.push_back(0);
        candidates.push_back(move(r));
    }
}

static uint64_t hilbert(uint32_t x, uint32_t y) {
    const uint32_t mask = (1u << 31) - 1;
    uint64_t d = 0;
    for (uint32_t s = 1u << 30; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += (uint64_t)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = mask - x; y = mask - y; }
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
    p.resize(n);
    vector<uint32_t> ux(n), uy(n);
    // Coordinates are in [-1e9,1e9], so this is an exact embedding in 31 bits.
    for (int i = 0; i < n; ++i) {
        cin >> p[i].x >> p[i].y;
        ux[i] = (uint32_t)(p[i].x + 1000000000LL);
        uy[i] = (uint32_t)(p[i].y + 1000000000LL);
    }
    primeId.assign(n, true);
    if (n > 0) primeId[0] = false;
    if (n > 1) primeId[1] = false;
    for (int q = 2; 1LL * q * q < n; ++q)
        if (primeId[q]) for (int j = q*q; j < n; j += q) primeId[j] = false;

    vector<vector<int>> candidates;
    // Retain the monotone-x route as a guard for instances where it is already good.
    vector<int> base; base.reserve(n+1);
    for (int i = 0; i < n; ++i) base.push_back(i);
    base.push_back(0);
    candidates.push_back(move(base));

    const uint32_t M = (1u << 31) - 1;
    for (int mode = 0; mode < 4; ++mode) {
        vector<pair<uint64_t,int>> order; order.reserve(n);
        for (int i = 0; i < n; ++i) {
            uint32_t a = ux[i], b = uy[i];
            if (mode == 1) a = M - a;
            else if (mode == 2) b = M - b;
            else if (mode == 3) swap(a, b);
            order.push_back({hilbert(a,b), i});
        }
        sort(order.begin(), order.end());
        int z = 0;
        while (order[z].second != 0) ++z;
        for (int dir : {1, -1}) {
            vector<int> r; r.reserve(n+1); r.push_back(0);
            for (int k = 1; k < n; ++k) r.push_back(order[(z + dir*k + n) % n].second);
            r.push_back(0);
            candidates.push_back(move(r));
        }
    }

    // Multi-resolution x-strip sweeps provide a separate representation from
    // Hilbert: each slab is traversed vertically, with neighboring slabs in
    // opposite directions.  Equal-population slabs remain robust to uneven
    // coordinate scales and use the input's guaranteed x ordering.
    for (int slabs = 1; slabs <= min(n, 256); slabs <<= 1)
        addSerpentineCandidates(candidates, n, slabs);

    // Refine the best spatial cycle.  This is deliberately bounded: it tests
    // whether missed short topology, rather than the global ordering, is the
    // incumbent's limiting factor while remaining safe at 200000 cities.
    vector<pair<double,int>> rank;
    for (int i = 0; i < (int)candidates.size(); ++i) rank.push_back({routeCost(candidates[i]), i});
    sort(rank.begin(), rank.end());
    int take = min(1, (int)rank.size());
    for (int k = 0; k < take; ++k) shortReverseImprove(candidates[rank[k].second]);

    int best = 0;
    double bestCost = routeCost(candidates[0]);
    for (int i = 1; i < (int)candidates.size(); ++i) {
        double v = routeCost(candidates[i]);
        if (v < bestCost) bestCost = v, best = i;
    }
    cout << n + 1 << '\n';
    for (int v : candidates[best]) cout << v << '\n';
}
