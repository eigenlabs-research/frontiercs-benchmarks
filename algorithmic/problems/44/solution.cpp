#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert index on a 2^31 by 2^31 grid.  Unlike sorting by x, this keeps
// nearby points together in both coordinates.
static uint64_t hilbert(uint32_t ux, uint32_t uy) {
    // Signed working coordinates are intentional: the standard Hilbert
    // rotation can temporarily reflect a coordinate below zero.
    int64_t x = ux, y = uy;
    uint64_t d = 0;
    for (int64_t s = 1LL << 30; s; s >>= 1) {
        int rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * uint64_t(s) * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) {
                x = s - 1 - x;
                y = s - 1 - y;
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
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &p : a) {
        cin >> p.x >> p.y;
        minx = min(minx, p.x); maxx = max(maxx, p.x);
        miny = min(miny, p.y); maxy = max(maxy, p.y);
    }

    vector<char> prime(n, true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    const long double SCALE = 2147483647.0L;
    long double rangex = (long double)maxx - minx, rangey = (long double)maxy - miny;
    vector<pair<uint64_t,int>> keyed;
    keyed.reserve(n);
    for (int i = 0; i < n; ++i) {
        uint32_t xx = rangex == 0 ? 0 : (uint32_t)(((long double)a[i].x - minx) * SCALE / rangex);
        uint32_t yy = rangey == 0 ? 0 : (uint32_t)(((long double)a[i].y - miny) * SCALE / rangey);
        keyed.push_back({hilbert(xx, yy), i});
    }
    sort(keyed.begin(), keyed.end());

    auto edge = [&](int u, int v) -> double {
        return hypot((double)a[u].x - a[v].x, (double)a[u].y - a[v].y);
    };
    auto cost = [&](const vector<int>& r) {
        double ans = 0;
        for (int t = 1; t <= n; ++t) {
            double w = (t % 10 == 0 && !prime[r[t - 1]]) ? 1.1 : 1.0;
            ans += w * edge(r[t - 1], r[t]);
        }
        return ans;
    };

    // Moving a nearby prime to a tenth-step source is accepted only when its
    // exact change to the complete objective is beneficial.
    auto carrotImprove = [&](vector<int>& r) {
        auto localCost = [&](const vector<int>& q, int u, int v) {
            int ts[4] = {u, u + 1, v, v + 1};
            double z = 0;
            for (int ii = 0; ii < 4; ++ii) {
                int t = ts[ii];
                bool seen = false;
                for (int jj = 0; jj < ii; ++jj) if (ts[jj] == t) seen = true;
                if (seen || t < 1 || t > n) continue;
                double w = (t % 10 == 0 && !prime[q[t - 1]]) ? 1.1 : 1.0;
                z += w * edge(q[t - 1], q[t]);
            }
            return z;
        };
        for (int t = 10; t <= n; t += 10) {
            int u = t - 1;
            if (prime[r[u]]) continue;
            int best = -1;
            double beforeBest = 0, afterBest = 0;
            for (int v = max(1, u - 6); v <= min(n - 1, u + 6); ++v) {
                if (v == u || !prime[r[v]]) continue;
                double before = localCost(r, u, v);
                swap(r[u], r[v]);
                double after = localCost(r, u, v);
                swap(r[u], r[v]);
                if (after + 1e-7 < before && (best < 0 || after - before < afterBest - beforeBest)) {
                    best = v; beforeBest = before; afterBest = after;
                }
            }
            if (best >= 0) swap(r[u], r[best]);
        }
    };

    // This is an exact bounded 2-opt search, including the positional carrot
    // multipliers.  It specifically tests whether the space-filling order
    // still contains short local crossings.  For a reversal [l,r], all
    // ordinary interior edges retain their lengths; only its two boundary
    // edges and the few tenth edges inside the window need to be rechecked.
    auto boundedTwoOpt = [&](vector<int>& r) {
        // Two orientations × three full scans × 28 candidates performs hundreds
        // of millions of hypot calls at N=200000.  Keep the refinement bounded
        // enough to remain a reliable improvement rather than a timeout risk.
        const int WINDOW = 12;
        for (int pass = 0; pass < 1; ++pass) {
            bool changed = false;
            for (int l = 1; l < n - 1; ++l) {
                int bestR = -1;
                double bestDelta = -1e-8;
                for (int len = 2; len <= WINDOW && l + len - 1 < n; ++len) {
                    int rr = l + len - 1;
                    double delta = edge(r[l - 1], r[rr]) + edge(r[l], r[rr + 1])
                                 - edge(r[l - 1], r[l]) - edge(r[rr], r[rr + 1]);
                    // Re-evaluate exactly those penalty edges whose endpoints
                    // or source identity can change under this reversal.
                    int first = ((l + 9) / 10) * 10;
                    for (int t = first; t <= rr + 1; t += 10) {
                        auto at = [&](int p) {
                            return (p >= l && p <= rr) ? r[l + rr - p] : r[p];
                        };
                        int oldA = r[t - 1], oldB = r[t];
                        int newA = at(t - 1), newB = at(t);
                        if (!prime[oldA]) delta -= 0.1 * edge(oldA, oldB);
                        if (!prime[newA]) delta += 0.1 * edge(newA, newB);
                    }
                    if (delta < bestDelta) bestDelta = delta, bestR = rr;
                }
                if (bestR >= 0) {
                    reverse(r.begin() + l, r.begin() + bestR + 1);
                    changed = true;
                }
            }
            if (!changed) break;
        }
    };

    // Adjacent exchanges are a cheap complementary descent: unlike 2-opt they
    // can repair a single bad ordering decision and deliberately choose which
    // city is the source of a carrot-sensitive edge.  Only three step costs
    // can change, so every accepted move has an exact objective delta.
    auto adjacentDescent = [&](vector<int>& r) {
        // One scan keeps this complementary pass safely below the expensive
        // bounded 2-opt budget on the 200k-city instances.
        for (int i = 1; i < n - 1; ++i) {
            int A = r[i - 1], B = r[i], C = r[i + 1], D = r[i + 2];
            auto weight = [&](int t, int source) {
                return (t % 10 == 0 && !prime[source]) ? 1.1 : 1.0;
            };
            double bc = edge(B, C);
            double before = weight(i, A) * edge(A, B)
                          + weight(i + 1, B) * bc
                          + weight(i + 2, C) * edge(C, D);
            double after = weight(i, A) * edge(A, C)
                         + weight(i + 1, C) * bc
                         + weight(i + 2, B) * edge(B, D);
            if (after + 1e-7 < before) swap(r[i], r[i + 1]);
        }
    };

    vector<int> order;
    order.reserve(n);
    for (auto [key, id] : keyed) order.push_back(id);
    int zero = find(order.begin(), order.end(), 0) - order.begin();

    vector<vector<int>> candidates;
    // The two directions have identical unpenalized cycle edges but put
    // different city IDs on the carrot-sensitive step positions.
    for (int direction : {1, -1}) {
        vector<int> r; r.reserve(n + 1); r.push_back(0);
        for (int k = 1; k < n; ++k) {
            int pos = (zero + direction * k) % n;
            if (pos < 0) pos += n;
            r.push_back(order[pos]);
        }
        r.push_back(0);
        carrotImprove(r);
        boundedTwoOpt(r);
        adjacentDescent(r);
        candidates.push_back(move(r));
    }
    // A deliberately different tour family: split the points into equal-population
    // horizontal ranks, then sweep each strip in the opposite x direction.  Unlike
    // Hilbert, this retains a long, coherent traversal across each layer and is
    // useful for data made of horizontal bands or highly anisotropic clouds.
    vector<int> byY(n), rankY(n);
    iota(byY.begin(), byY.end(), 0);
    sort(byY.begin(), byY.end(), [&](int u, int v) {
        if (a[u].y != a[v].y) return a[u].y < a[v].y;
        return a[u].x < a[v].x;
    });
    for (int k = 0; k < n; ++k) rankY[byY[k]] = k;
    vector<int> stripCounts = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512};
    long double aspect = rangex > 0 ? rangey / rangex : (long double)n;
    int idealStrips = max(1, min(n - 1, (int)llround(sqrt((long double)n * aspect))));
    stripCounts.push_back(idealStrips);
    sort(stripCounts.begin(), stripCounts.end());
    stripCounts.erase(unique(stripCounts.begin(), stripCounts.end()), stripCounts.end());
    for (int rows : stripCounts) {
        if (rows >= n) continue;
        vector<vector<int>> bucket(rows);
        // IDs are x-ordered by the input promise, so each bucket is already an
        // x-sorted strip without another sort.
        for (int id = 0; id < n; ++id)
            bucket[(long long)rankY[id] * rows / n].push_back(id);
        vector<int> sweep;
        sweep.reserve(n);
        for (int b = 0; b < rows; ++b) {
            if ((b & 1) == 0) {
                for (int id : bucket[b]) sweep.push_back(id);
            } else {
                for (auto it = bucket[b].rbegin(); it != bucket[b].rend(); ++it) sweep.push_back(*it);
            }
        }
        int atZero = find(sweep.begin(), sweep.end(), 0) - sweep.begin();
        for (int direction : {1, -1}) {
            vector<int> r;
            r.reserve(n + 1);
            r.push_back(0);
            for (int k = 1; k < n; ++k) {
                int pos = (atZero + direction * k) % n;
                if (pos < 0) pos += n;
                r.push_back(sweep[pos]);
            }
            r.push_back(0);
            candidates.push_back(move(r));
        }
    }

    // Never lose to the supplied monotone baseline on adversarial geometry.
    vector<int> baseline(n + 1);
    iota(baseline.begin(), baseline.end(), 0);
    baseline[n] = 0;
    candidates.push_back(move(baseline));

    int take = 0;
    double best = cost(candidates[0]);
    for (int i = 1; i < (int)candidates.size(); ++i) {
        double cur = cost(candidates[i]);
        if (cur < best) best = cur, take = i;
    }
    cout << n + 1 << '\n';
    for (int id : candidates[take]) cout << id << '\n';
}
