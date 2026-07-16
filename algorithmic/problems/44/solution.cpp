#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };
static vector<Point> pt;
static vector<char> primeId;

static inline double edgeCost(const vector<int>& p, int k) {
    int a = p[k], b = p[k + 1];
    double dx = (double)pt[a].x - pt[b].x;
    double dy = (double)pt[a].y - pt[b].y;
    double w = ((k % 10 == 9) && !primeId[a]) ? 1.1 : 1.0;
    return w * hypot(dx, dy);
}
static double routeCost(const vector<int>& p) {
    double r = 0;
    for (int i = 0; i + 1 < (int)p.size(); ++i) r += edgeCost(p, i);
    return r;
}

// A 31-bit Morton (Z-order) key keeps points close in both coordinates,
// unlike the supplied monotone-x order on interleaved or clustered inputs.
static uint64_t morton(uint32_t x, uint32_t y) {
    uint64_t z = 0;
    for (int b = 0; b < 31; ++b) {
        z |= (uint64_t)((x >> b) & 1) << (2 * b + 1);
        z |= (uint64_t)((y >> b) & 1) << (2 * b);
    }
    return z;
}

// Unlike Z-order, a Hilbert traversal does not jump between diagonally touching
// quadrants.  It is a separate spatial representation, not a local mutation
// of the Morton candidate.
static uint64_t hilbert(uint32_t x0, uint32_t y0) {
    long long x = x0, y = y0;
    uint64_t d = 0;
    for (long long s = 1LL << 30; s; s >>= 1) {
        long long rx = (x & s) != 0, ry = (y & s) != 0;
        d += (uint64_t)(s * s) * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = s - 1 - x; y = s - 1 - y; }
            swap(x, y);
        }
    }
    return d;
}

// Greedy short 2-opt is deliberately bounded: it repairs crossings and small
// ordering errors without making a quadratic TSP search infeasible at 200k.
static void localImprove(vector<int>& p, int W = 8) {
    int n = (int)p.size() - 1;
    for (int i = 0; i < n - 1; ++i) {
        int bestj = -1;
        double bestDelta = 0.0;
        for (int j = i + 2; j <= min(n - 1, i + W); ++j) {
            double oldc = 0, newc = 0;
            // Reversing vertices [i+1,j] changes exactly these source edges,
            // including their step-dependent carrot multipliers.
            for (int k = i; k <= j; ++k) oldc += edgeCost(p, k);
            reverse(p.begin() + i + 1, p.begin() + j + 1);
            for (int k = i; k <= j; ++k) newc += edgeCost(p, k);
            reverse(p.begin() + i + 1, p.begin() + j + 1);
            if (newc - oldc < bestDelta) {
                bestDelta = newc - oldc;
                bestj = j;
            }
        }
        if (bestj != -1) reverse(p.begin() + i + 1, p.begin() + bestj + 1);
    }
}

// This is an explicit test of the carrot mechanism.  The old solver only
// incidentally moved prime IDs through 2-opt.  At each expensive source slot,
// try exchanging it with a nearby prime and accept only the exact improvement.
static void improvePenaltySlots(vector<int>& p) {
    const int n = (int)p.size() - 1;
    const int W = 18;
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int k = 9; k < n; k += 10) {
            if (primeId[p[k]]) continue;
            int best = -1;
            double bestDelta = 0.0;
            for (int j = max(1, k - W); j <= min(n - 1, k + W); ++j) {
                if (j == k || !primeId[p[j]]) continue;
                int e[4] = {k - 1, k, j - 1, j};
                sort(e, e + 4);
                double oldc = 0, newc = 0;
                for (int q = 0; q < 4; ++q) {
                    if (q && e[q] == e[q - 1]) continue;
                    oldc += edgeCost(p, e[q]);
                }
                swap(p[k], p[j]);
                for (int q = 0; q < 4; ++q) {
                    if (q && e[q] == e[q - 1]) continue;
                    newc += edgeCost(p, e[q]);
                }
                swap(p[k], p[j]);
                if (newc - oldc < bestDelta) bestDelta = newc - oldc, best = j;
            }
            if (best != -1) {
                swap(p[k], p[best]);
                changed = true;
            }
        }
        if (!changed) break;
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    pt.resize(n);
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &q : pt) {
        cin >> q.x >> q.y;
        minx = min(minx, q.x); maxx = max(maxx, q.x);
        miny = min(miny, q.y); maxy = max(maxy, q.y);
    }
    primeId.assign(n, true);
    primeId[0] = false; if (n > 1) primeId[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (primeId[i])
        for (int j = i * i; j < n; j += i) primeId[j] = false;

    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    const uint64_t M = (1ull << 31) - 1;
    long long rx = maxx - minx, ry = maxy - miny;
    vector<uint64_t> key(n), hkey(n);
    for (int i = 0; i < n; ++i) {
        uint32_t xx = (uint64_t)(pt[i].x - minx) * M / rx;
        uint32_t yy = ry ? (uint64_t)(pt[i].y - miny) * M / ry : 0;
        key[i] = morton(xx, yy);
        hkey[i] = hilbert(xx, yy);
    }
    sort(order.begin(), order.end(), [&](int a, int b) { return key[a] < key[b]; });
    vector<int> hilbertOrder(n);
    iota(hilbertOrder.begin(), hilbertOrder.end(), 0);
    sort(hilbertOrder.begin(), hilbertOrder.end(), [&](int a, int b) { return hkey[a] < hkey[b]; });

    auto makeRoute = [&](const vector<int>& v, bool backwards) {
        vector<int> p; p.reserve(n + 1); p.push_back(0);
        int z = find(v.begin(), v.end(), 0) - v.begin();
        for (int t = 1; t < n; ++t) {
            int at = backwards ? (z - t + n) % n : (z + t) % n;
            p.push_back(v[at]);
        }
        p.push_back(0);
        return p;
    };
    vector<int> byX(n); iota(byX.begin(), byX.end(), 0);
    vector<vector<int>> candidates;
    candidates.push_back(makeRoute(byX, false)); // exactly the supplied baseline
    candidates.push_back(makeRoute(byX, true));
    candidates.push_back(makeRoute(order, false));
    candidates.push_back(makeRoute(order, true));
    candidates.push_back(makeRoute(hilbertOrder, false));
    candidates.push_back(makeRoute(hilbertOrder, true));

    // Selecting on raw construction cost can discard a tour whose local geometry
    // repairs make it best. Give every construction the same inexpensive polish
    // before the selection gate; reserve the full pass for the selected tour.
    vector<int> best;
    double bestCost = numeric_limits<double>::infinity();
    for (auto& candidate : candidates) {
        localImprove(candidate, 4);
        double c = routeCost(candidate);
        if (c < bestCost) { bestCost = c; best = move(candidate); }
    }
    improvePenaltySlots(best);
    // A final local repair preserves the invariant that every accepted move is
    // beneficial under the full, position-dependent objective.
    localImprove(best);
    cout << n + 1 << '\n';
    for (int id : best) cout << id << '\n';
}
