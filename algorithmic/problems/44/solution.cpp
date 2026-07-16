#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert index on the fixed coordinate square [-1e9,1e9]^2.  Keeping the
// original scale (rather than independently normalizing axes) preserves the
// geometry which the curve is meant to approximate.
static unsigned long long hilbert(unsigned int x, unsigned int y) {
    unsigned long long d = 0;
    for (unsigned int s = 1u << 30; s; s >>= 1) {
        unsigned int rx = (x & s) != 0, ry = (y & s) != 0;
        d += (unsigned long long)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) {
                x = (1u << 31) - 1 - x;
                y = (1u << 31) - 1 - y;
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
    vector<Point> a(n);
    for (auto &p : a) cin >> p.x >> p.y;

    vector<char> prime(n, true);
    if (n > 0) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    // Match the checker's arithmetic: it uses double hypot and a double sum.
    auto edgeCost = [&](const vector<int>& r, int t) -> double {
        int u = r[t - 1], v = r[t];
        double dx = (double)a[u].x - (double)a[v].x;
        double dy = (double)a[u].y - (double)a[v].y;
        double z = hypot(dx, dy);
        if (t % 10 == 0 && !prime[u]) z *= 1.1;
        return z;
    };
    auto totalCost = [&](const vector<int>& r) {
        double ans = 0;
        for (int t = 1; t <= n; ++t) ans += edgeCost(r, t);
        return ans;
    };

    // Adjacent exchanges are cheap to assess exactly, including the changing
    // step multipliers.  A few sweeps remove local inversions/noise left by a
    // space filling curve without risking expensive global search.
    auto polish = [&](vector<int> r, int passes) {
        for (int pass = 0; pass < passes; ++pass) {
            bool backwards = pass & 1;
            for (int q = 1; q <= n - 2; ++q) {
                int i = backwards ? n - 1 - q : q;
                double before = 0;
                for (int t = i; t <= i + 2; ++t) before += edgeCost(r, t);
                swap(r[i], r[i + 1]);
                double after = 0;
                for (int t = i; t <= i + 2; ++t) after += edgeCost(r, t);
                if (after + 1e-7L >= before) swap(r[i], r[i + 1]);
            }
        }
        return r;
    };

    // A Hilbert curve is directional.  Its eight square symmetries give
    // different, equally geometric skeletons; this matters because city 0 is
    // fixed and the carrot positions make a reversed tour non-equivalent.
    // Keep all orientations as alternatives, rather than assuming the native
    // Hilbert orientation happens to align with the input distribution.
    const unsigned int LIM = 2000000000u;
    vector<int> ans, originalOrientation;
    bool ansFromOriginal = true;
    double bestCost = numeric_limits<double>::infinity();
    double originalCost = numeric_limits<double>::infinity();
    for (int orientation = 0; orientation < 8; ++orientation) {
        vector<pair<unsigned long long,int>> keyed;
        keyed.reserve(n - 1);
        for (int i = 1; i < n; ++i) {
            unsigned int x = (unsigned int)(a[i].x + 1000000000LL);
            unsigned int y = (unsigned int)(a[i].y + 1000000000LL);
            unsigned int u, v;
            switch (orientation) {
                case 0: u = x;       v = y;       break;
                case 1: u = y;       v = LIM - x; break;
                case 2: u = LIM - x; v = LIM - y; break;
                case 3: u = LIM - y; v = x;       break;
                case 4: u = x;       v = LIM - y; break;
                case 5: u = LIM - x; v = y;       break;
                case 6: u = y;       v = x;       break;
                default:u = LIM - y; v = LIM - x; break;
            }
            keyed.push_back({hilbert(u, v), i});
        }
        sort(keyed.begin(), keyed.end());
        vector<int> base(n + 1);
        base[0] = base[n] = 0;
        for (int i = 1; i < n; ++i) base[i] = keyed[i - 1].second;
        for (int direction = 0; direction < 2; ++direction) {
            vector<int> candidate(n + 1);
            candidate[0] = candidate[n] = 0;
            for (int i = 1; i < n; ++i)
                candidate[i] = direction ? base[n - i] : base[i];
            candidate = polish(move(candidate), orientation == 0 ? 4 : 2);
            double cost = totalCost(candidate);
            if (orientation == 0 && cost < originalCost) {
                originalCost = cost;
                originalOrientation = candidate;
            }
            if (cost < bestCost) {
                bestCost = cost;
                ansFromOriginal = (orientation == 0);
                ans = move(candidate);
            }
        }
    }

    // Discriminating neighborhood: an adjacent descent cannot remove a small
    // crossing when every intermediate adjacent swap is uphill.  Test exact
    // 2-opt reversals of short spans.  Reversing an edge preserves its length,
    // but changes both its source and step number, so all affected carrot
    // multipliers are explicitly re-evaluated.  Keeping edge lengths makes
    // this O(n * span) arithmetic rather than an expensive geometric search.
    auto shortTwoOpt = [&](vector<int> r) {
        vector<double> len(n + 1);
        auto rawDist = [&](int u, int v) {
            double dx = (double)a[u].x - (double)a[v].x;
            double dy = (double)a[u].y - (double)a[v].y;
            return hypot(dx, dy);
        };
        auto weighted = [&](int t, int u, double d) {
            return (t % 10 == 0 && !prime[u]) ? d * 1.1 : d;
        };
        for (int t = 1; t <= n; ++t) len[t] = rawDist(r[t - 1], r[t]);
        // Exhaust the very short neighborhood, then probe a few larger
        // lengths.  A curve-order discontinuity frequently needs more than
        // twelve vertices to uncross, but trying every such length would not
        // fit the largest instances.  The sparse lengths distinguish that
        // failure mode while retaining linear-size search.
        const int spans[] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                             16, 24, 32, 48};
        for (int pass = 0; pass < 2; ++pass) {
            for (int l = 1; l <= n - 3; ++l) {
                for (int span : spans) {
                    if (l + span > n - 1) continue;
                    int rr = l + span;
                    // Ordinary internal edge lengths have exactly the same
                    // sum after a reversal.  Only its two boundary edges and
                    // internal tenth-step surcharges can differ.  Evaluating
                    // only those terms keeps the sparse long probe viable.
                    double before = weighted(l, r[l - 1], len[l])
                                  + weighted(rr + 1, r[rr], len[rr + 1]);
                    double after = weighted(l, r[l - 1], rawDist(r[l - 1], r[rr]))
                                 + weighted(rr + 1, r[l], rawDist(r[l], r[rr + 1]));
                    // Boundary edge l was already included above; only
                    // internal reversed edges l+1..rr need surcharge deltas.
                    int firstPenalty = ((l + 1 + 9) / 10) * 10;
                    for (int t = firstPenalty; t <= rr; t += 10) {
                        int j = l + rr - t;
                        if (!prime[r[t - 1]]) before += 0.1 * len[t];
                        if (!prime[r[j + 1]]) after += 0.1 * len[j + 1];
                    }
                    if (after + 1e-7 < before) {
                        reverse(r.begin() + l, r.begin() + rr + 1);
                        for (int t = l; t <= rr + 1; ++t)
                            len[t] = rawDist(r[t - 1], r[t]);
                    }
                }
            }
        }
        return r;
    };
    vector<int> improved = shortTwoOpt(ans);
    // This final comparison is also a guard against any numerical ambiguity
    // in the local delta calculation: the incumbent route is never discarded
    // unless the checker-equivalent full objective is lower.
    if (totalCost(improved) < totalCost(ans)) ans = move(improved);
    // The original orientation with its two directions is exactly the prior
    // incumbent before short 2-opt.  Also test its completed incumbent path,
    // so selecting a promising alternate skeleton cannot regress it.
    if (!ansFromOriginal) {
        vector<int> priorImproved = shortTwoOpt(originalOrientation);
        if (totalCost(priorImproved) < totalCost(ans)) ans = move(priorImproved);
    }

    cout << n + 1 << '\n';
    for (int v : ans) cout << v << '\n';
    return 0;
}
