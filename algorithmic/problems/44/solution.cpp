#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };
static vector<Point> pt;
static vector<char> primeId;
static int n;

static inline double edgeCost(const vector<int>& p, int pos) {
    // pos is the source position, hence its step is pos+1.
    int a = p[pos], b = p[pos + 1];
    double dx = (double)pt[a].x - pt[b].x;
    double dy = (double)pt[a].y - pt[b].y;
    double d = hypot(dx, dy);
    return (((pos + 1) % 10 == 0 && !primeId[a]) ? 1.1 : 1.0) * d;
}
static double score(const vector<int>& p) {
    double ans = 0;
    for (int i = 0; i < n; ++i) ans += edgeCost(p, i);
    return ans;
}

// Standard Hilbert-coordinate conversion.  The eight choices of axis/reflection
// below are important because city 0 is fixed as the beginning of the tour.
static uint64_t hilbert(uint32_t x, uint32_t y) {
    uint64_t d = 0;
    for (uint32_t s = 1u << 30; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += (uint64_t)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = (uint32_t)((1u << 31) - 1) - x; y = (uint32_t)((1u << 31) - 1) - y; }
            swap(x, y);
        }
    }
    return d;
}

// Exact bounded 2-opt descent.  On a reversal, ordinary internal edges keep
// the same length and multiplier; only the two boundary edges and internal
// tenth-step edges need to be evaluated.  This makes a substantially wider
// neighborhood practical even at the largest input size.
static void shortReverseDescent(vector<int>& p) {
    const int maxGap = 24; // reverse 2 through 25 consecutive cities
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int l = 1; l < n; ++l) {
            for (int gap = 1; gap <= maxGap && l + gap < n; ++gap) {
                int r = l + gap;
                // The changed boundary edges are l-1 and r.  The edges
                // strictly inside the reversal only matter at carrot slots.
                double before = edgeCost(p, l - 1) + edgeCost(p, r);
                for (int e = l; e < r; ++e)
                    if ((e + 1) % 10 == 0) before += edgeCost(p, e);

                auto cityAfter = [&](int pos) -> int {
                    return (l <= pos && pos <= r) ? p[l + r - pos] : p[pos];
                };
                auto newEdge = [&](int e) -> double {
                    int a = cityAfter(e), b = cityAfter(e + 1);
                    double dx = (double)pt[a].x - pt[b].x;
                    double dy = (double)pt[a].y - pt[b].y;
                    double d = hypot(dx, dy);
                    return (((e + 1) % 10 == 0 && !primeId[a]) ? 1.1 : 1.0) * d;
                };
                double after = newEdge(l - 1) + newEdge(r);
                for (int e = l; e < r; ++e)
                    if ((e + 1) % 10 == 0) after += newEdge(e);

                if (after + 1e-7 < before) {
                    reverse(p.begin() + l, p.begin() + r + 1);
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    if (!(cin >> n)) return 0;
    pt.resize(n);
    for (auto &q : pt) cin >> q.x >> q.y;

    primeId.assign(max(2, n), true);
    primeId[0] = false; primeId[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (primeId[i]) for (int j = i * i; j < n; j += i) primeId[j] = false;

    vector<int> best(n + 1);
    iota(best.begin(), best.end(), 0);
    best[n] = 0;
    double bestCost = score(best); // retaining this rules out a bad spatial order

    long long minx = pt[0].x, maxx = pt[0].x, miny = pt[0].y, maxy = pt[0].y;
    for (auto q : pt) {
        minx = min(minx, q.x); maxx = max(maxx, q.x);
        miny = min(miny, q.y); maxy = max(maxy, q.y);
    }
    const uint64_t LIM = (1ULL << 31) - 1;
    vector<pair<uint64_t,int>> keyed(n - 1);
    for (int orient = 0; orient < 8; ++orient) {
        for (int id = 1; id < n; ++id) {
            uint64_t ax = maxx == minx ? LIM / 2 : (uint64_t)((__int128)(pt[id].x - minx) * LIM / (maxx - minx));
            uint64_t ay = maxy == miny ? LIM / 2 : (uint64_t)((__int128)(pt[id].y - miny) * LIM / (maxy - miny));
            uint64_t u = ax, v = ay;
            if (orient & 1) u = LIM - u;
            if (orient & 2) v = LIM - v;
            if (orient & 4) swap(u, v);
            keyed[id - 1] = {hilbert((uint32_t)u, (uint32_t)v), id};
        }
        sort(keyed.begin(), keyed.end());
        vector<int> cand; cand.reserve(n + 1);
        cand.push_back(0);
        for (auto z : keyed) cand.push_back(z.second);
        cand.push_back(0);
        double c = score(cand);
        if (c < bestCost) { bestCost = c; best.swap(cand); }
    }

    // Discriminating repair: test whether the Hilbert route still has cheap
    // local crossings under the *actual* indexed carrot objective.
    vector<int> improved = best;
    shortReverseDescent(improved);
    double improvedCost = score(improved);
    if (improvedCost < bestCost) best.swap(improved);

    cout << n + 1 << '\n';
    for (int id : best) cout << id << '\n';
}
