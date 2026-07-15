#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Order a square by a Hilbert curve.  The recursion only has 21 levels here.
static uint64_t hilbert(uint32_t x, uint32_t y, int level, int rotation = 0) {
    if (level == 0) return 0;
    uint32_t half = 1u << (level - 1);
    int segment;
    if (x < half) segment = (y < half ? 0 : 3);
    else          segment = (y < half ? 1 : 2);
    segment = (segment + rotation) & 3;
    static const int delta[4] = {3, 0, 0, 1};
    uint64_t block = 1ULL << (2 * level - 2);
    uint64_t inside = hilbert(x & (half - 1), y & (half - 1), level - 1,
                              (rotation + delta[segment]) & 3);
    return uint64_t(segment) * block + ((segment == 1 || segment == 2) ? inside : block - 1 - inside);
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
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    auto normalized = [](long long v, long long lo, long long hi) -> uint32_t {
        static constexpr uint64_t LIM = (1u << 21) - 1;
        if (lo == hi) return 0;
        return (uint32_t)((__int128)(v - lo) * LIM / (hi - lo));
    };
    vector<uint32_t> xx(n), yy(n);
    for (int i = 0; i < n; ++i) {
        xx[i] = normalized(p[i].x, minx, maxx);
        yy[i] = normalized(p[i].y, miny, maxy);
    }

    auto cost = [&](const vector<int>& r) {
        double ans = 0;
        for (int i = 0; i < n; ++i) {
            int a = r[i], b = r[(i + 1) % n];
            double d = hypot((double)p[a].x - p[b].x, (double)p[a].y - p[b].y);
            if ((i + 1) % 10 == 0 && !prime[a]) d *= 1.1;
            ans += d;
        }
        return ans;
    };

    // Retain the supplied monotone order as a safe candidate.
    vector<int> best(n), candidate(n);
    iota(best.begin(), best.end(), 0);
    double bestCost = cost(best);
    candidate[0] = 0;
    for (int i = 1; i < n; ++i) candidate[i] = n - i;
    double reverseCost = cost(candidate);
    if (reverseCost < bestCost) { bestCost = reverseCost; best = candidate; }

    static constexpr uint32_t LIM = (1u << 21) - 1;
    // The eight reflections/axis swaps avoid committing to an arbitrary Hilbert orientation.
    for (int mode = 0; mode < 8; ++mode) {
        vector<pair<uint64_t,int>> keyed;
        keyed.reserve(n);
        for (int i = 0; i < n; ++i) {
            uint32_t a = xx[i], b = yy[i];
            if (mode & 1) swap(a, b);
            if (mode & 2) a = LIM - a;
            if (mode & 4) b = LIM - b;
            keyed.push_back({hilbert(a, b, 21), i});
        }
        sort(keyed.begin(), keyed.end());
        int at = 0;
        while (keyed[at].second != 0) ++at;
        for (int direction : {-1, 1}) {
            candidate[0] = 0;
            for (int k = 1; k < n; ++k) {
                int j = at + direction * k;
                if (j < 0) j += n;
                else if (j >= n) j -= n;
                candidate[k] = keyed[j].second;
            }
            double value = cost(candidate);
            if (value < bestCost) {
                bestCost = value;
                best = candidate;
            }
        }
    }

    // The space-filling order is globally good, but its cell boundaries can leave
    // adjacent inversions.  Exchange descent is deliberately evaluated with the
    // real step numbers: swapping positions i and i+1 changes three directed,
    // potentially carrot-weighted edges.
    auto edgeCost = [&](int i) {
        int a = best[i], b = best[(i + 1) % n];
        double d = hypot((double)p[a].x - p[b].x, (double)p[a].y - p[b].y);
        if ((i + 1) % 10 == 0 && !prime[a]) d *= 1.1;
        return d;
    };
    // A few passes capture local inversions without turning this into a
    // quadratic TSP search.  Strict acceptance means this cannot worsen the
    // chosen incumbent route under the checker objective.
    for (int pass = 0; pass < 4; ++pass) {
        bool changed = false;
        for (int i = 1; i + 1 < n; ++i) {
            double before = edgeCost(i - 1) + edgeCost(i) + edgeCost(i + 1);
            swap(best[i], best[i + 1]);
            double after = edgeCost(i - 1) + edgeCost(i) + edgeCost(i + 1);
            if (after < before) {
                changed = true;
            } else {
                swap(best[i], best[i + 1]);
            }
        }
        if (!changed) break;
    }

    // Adjacent exchanges cannot remove a crossing whose endpoints are separated
    // by a short run.  For a 2-opt reversal [i+1, j], ordinary interior edge
    // lengths are unchanged; only the two boundary edges and carrot surcharge
    // terms on interior tenth steps need re-evaluation.  This makes an exact
    // bounded 2-opt pass cheap even at the largest N.
    auto rawDist = [&](int a, int b) {
        return hypot((double)p[a].x - p[b].x, (double)p[a].y - p[b].y);
    };
    auto weighted = [&](int edge, int a, int b) {
        double d = rawDist(a, b);
        if ((edge + 1) % 10 == 0 && !prime[a]) d *= 1.1;
        return d;
    };
    constexpr int MAX_2OPT_GAP = 24;
    for (int i = 0; i + 2 < n; ++i) {
        int last = min(n - 1, i + MAX_2OPT_GAP);
        for (int j = i + 2; j <= last; ++j) {
            // Boundary edges are fully replaced by the reversal.
            double before = weighted(i, best[i], best[i + 1]) +
                            weighted(j, best[j], best[(j + 1) % n]);
            double after = weighted(i, best[i], best[j]) +
                           weighted(j, best[i + 1], best[(j + 1) % n]);
            // An interior edge retains its distance after reversal, so account
            // only for its possible 10% source-city surcharge.
            int firstPenalty = ((i + 1 + 9) / 10) * 10 - 1;
            for (int k = firstPenalty; k < j; k += 10) {
                if (k < i + 1) continue;
                if (!prime[best[k]]) before += 0.1 * rawDist(best[k], best[k + 1]);
                int sourceAfter = i + j + 1 - k;
                if (!prime[best[sourceAfter]])
                    after += 0.1 * rawDist(best[sourceAfter], best[sourceAfter - 1]);
            }
            if (after + 1e-9 < before)
                reverse(best.begin() + i + 1, best.begin() + j + 1);
        }
    }

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
    cout << 0 << '\n';
    return 0;
}
