#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Standard xy-to-Hilbert conversion. 31 bits cover the complete input range.
static uint64_t hilbert(uint32_t x, uint32_t y) {
    uint64_t d = 0;
    for (uint32_t s = 1u << 30; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * s * ((3 * rx) ^ ry);
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
    for (auto &p : a) cin >> p.x >> p.y;

    vector<char> prime(max(2, n), true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    auto edge = [&](const vector<int>& r, int source) -> long double {
        if (source < 0 || source >= n) return 0;
        int u = r[source], v = r[source + 1];
        long double dx = (long double)a[u].x - a[v].x;
        long double dy = (long double)a[u].y - a[v].y;
        long double z = sqrtl(dx * dx + dy * dy);
        if ((source + 1) % 10 == 0 && !prime[u]) z *= 1.1L;
        return z;
    };
    auto cost = [&](const vector<int>& r) {
        long double ans = 0;
        for (int i = 0; i < n; ++i) ans += edge(r, i);
        return ans;
    };

    // A small exact local repair: at penalized positions, exchange with a nearby
    // prime only when the full weighted objective decreases.
    auto carrotRepair = [&](vector<int>& r) {
        const int W = 80;
        for (int p = 9; p < n; p += 10) {
            if (prime[r[p]]) continue;
            int best = -1;
            long double bestDelta = 0;
            for (int q = max(1, p - W); q <= min(n - 1, p + W); ++q) {
                if (!prime[r[q]] || (q + 1) % 10 == 0) continue;
                // A repair evaluates O(N*W) swaps; keep this exact affected-
                // edge set on the stack rather than heap-allocating a vector
                // for every trial.
                int ids[4] = {p - 1, p, q - 1, q};
                sort(ids, ids + 4);
                long double before = 0;
                for (int z = 0; z < 4; ++z)
                    if (z == 0 || ids[z] != ids[z - 1]) before += edge(r, ids[z]);
                swap(r[p], r[q]);
                long double after = 0;
                for (int z = 0; z < 4; ++z)
                    if (z == 0 || ids[z] != ids[z - 1]) after += edge(r, ids[z]);
                swap(r[p], r[q]);
                if (after - before < bestDelta) {
                    bestDelta = after - before;
                    best = q;
                }
            }
            if (best != -1) swap(r[p], r[best]);
        }
    };

    vector<pair<uint64_t, int>> order;
    order.reserve(n);
    constexpr uint64_t SCALE = 2147483647ULL;
    constexpr uint64_t RANGE = 2000000000ULL;
    for (int i = 0; i < n; ++i) {
        uint32_t xx = uint32_t((uint64_t)(a[i].x + 1000000000LL) * SCALE / RANGE);
        uint32_t yy = uint32_t((uint64_t)(a[i].y + 1000000000LL) * SCALE / RANGE);
        order.push_back({hilbert(xx, yy), i});
    }
    sort(order.begin(), order.end(), [](auto A, auto B) {
        return A.first != B.first ? A.first < B.first : A.second < B.second;
    });
    int at0 = 0;
    while (order[at0].second != 0) ++at0;

    auto makeRoute = [&](bool backwards) {
        vector<int> r;
        r.reserve(n + 1);
        r.push_back(0);
        for (int k = 1; k < n; ++k) {
            int j = backwards ? (at0 - k + n) % n : (at0 + k) % n;
            r.push_back(order[j].second);
        }
        r.push_back(0);
        carrotRepair(r);
        return r;
    };

    vector<int> best = makeRoute(false);
    long double bestCost = cost(best);
    vector<int> other = makeRoute(true);
    if (cost(other) < bestCost) { best = move(other); bestCost = cost(best); }

    // Diversification branch: partition the x-ordered input into vertical slabs,
    // sweep each slab in alternating y directions, and regard the result as a
    // cycle which can be cut at city zero.  This is deliberately unlike the
    // space-filling curve and is useful on long, horizontally banded instances.
    auto makeSlabCycle = [&](int slabs, bool backwards) {
        vector<int> cyc;
        cyc.reserve(n);
        for (int s = 0; s < slabs; ++s) {
            int lo = (long long)s * n / slabs;
            int hi = (long long)(s + 1) * n / slabs;
            vector<int> part;
            part.reserve(hi - lo);
            for (int i = lo; i < hi; ++i) part.push_back(i);
            sort(part.begin(), part.end(), [&](int u, int v) {
                if (a[u].y != a[v].y) return a[u].y < a[v].y;
                return u < v;
            });
            if ((s & 1) != 0) reverse(part.begin(), part.end());
            cyc.insert(cyc.end(), part.begin(), part.end());
        }
        int z = int(find(cyc.begin(), cyc.end(), 0) - cyc.begin());
        vector<int> r;
        r.reserve(n + 1);
        r.push_back(0);
        for (int k = 1; k < n; ++k) {
            int j = backwards ? (z - k + n) % n : (z + k) % n;
            r.push_back(cyc[j]);
        }
        r.push_back(0);
        return r;
    };
    // A transposed sweep is a separate candidate family: rank cities by y,
    // traverse each horizontal layer in x order, and reverse alternate layers.
    // It is particularly effective when the input consists of horizontal bands.
    vector<int> yRank(n);
    iota(yRank.begin(), yRank.end(), 0);
    sort(yRank.begin(), yRank.end(), [&](int u, int v) {
        if (a[u].y != a[v].y) return a[u].y < a[v].y;
        return u < v;
    });
    auto makeYSlabCycle = [&](int slabs, bool backwards) {
        vector<int> byY = yRank;
        vector<int> cyc;
        cyc.reserve(n);
        for (int s = 0; s < slabs; ++s) {
            int lo = (long long)s * n / slabs;
            int hi = (long long)(s + 1) * n / slabs;
            // IDs are already x-ordered, but sorting makes this independent
            // of that presentation detail and handles tied x defensively.
            sort(byY.begin() + lo, byY.begin() + hi, [&](int u, int v) {
                if (a[u].x != a[v].x) return a[u].x < a[v].x;
                return u < v;
            });
            if (s & 1) reverse(byY.begin() + lo, byY.begin() + hi);
            cyc.insert(cyc.end(), byY.begin() + lo, byY.begin() + hi);
        }
        int z = int(find(cyc.begin(), cyc.end(), 0) - cyc.begin());
        vector<int> r;
        r.reserve(n + 1);
        r.push_back(0);
        for (int k = 1; k < n; ++k) {
            int j = backwards ? (z - k + n) % n : (z + k) % n;
            r.push_back(cyc[j]);
        }
        r.push_back(0);
        return r;
    };

    // Screen using raw geometry, then spend the prime-aware repair budget only
    // on the strongest representatives of this separate route family.
    vector<pair<long double, vector<int>>> slabPool;
    for (int slabs : {2, 4, 8, 16, 32, 64}) {
        if (slabs >= n) continue;
        for (bool backwards : {false, true}) {
            vector<int> r = makeSlabCycle(slabs, backwards);
            slabPool.push_back({cost(r), move(r)});
        }
    }
    sort(slabPool.begin(), slabPool.end(), [](const auto& u, const auto& v) {
        return u.first < v.first;
    });
    for (int i = 0; i < (int)slabPool.size() && i < 3; ++i) {
        vector<int> candidate = move(slabPool[i].second);
        carrotRepair(candidate);
        long double candidateCost = cost(candidate);
        if (candidateCost < bestCost) {
            best = move(candidate);
            bestCost = candidateCost;
        }
    }

    // Unlike adding another always-repaired family, screen all y sweeps first
    // and repair only their strongest representative.  This keeps the large-N
    // search budget bounded while retaining a route that can win on layers.
    pair<long double, vector<int>> bestYSweep = {numeric_limits<long double>::infinity(), {}};
    for (int slabs : {2, 4, 8, 16, 32}) {
        if (slabs >= n) continue;
        for (bool backwards : {false, true}) {
            vector<int> r = makeYSlabCycle(slabs, backwards);
            long double raw = cost(r);
            if (raw < bestYSweep.first) bestYSweep = {raw, move(r)};
        }
    }
    if (!bestYSweep.second.empty() && bestYSweep.first < bestCost * 1.02L) {
        vector<int> candidate = move(bestYSweep.second);
        carrotRepair(candidate);
        long double candidateCost = cost(candidate);
        if (candidateCost < bestCost) {
            best = move(candidate);
            bestCost = candidateCost;
        }
    }

    // Retain the stated monotone tour as a safety net for degenerate geometry.
    vector<int> monotone(n + 1);
    iota(monotone.begin(), monotone.end() - 1, 0);
    monotone[n] = 0;
    if (cost(monotone) < bestCost) { best = move(monotone); bestCost = cost(best); }

    // Discriminating local search: a Hilbert traversal can still contain small
    // crossings.  Unlike ordinary 2-opt, reversal changes the city that starts
    // each tenth edge, so score every affected penalized edge exactly.
    auto twoOpt = [&](vector<int>& r) {
        const int W = 12;
        for (int i = 0; i <= n - 3; ++i) {
            int bestJ = -1;
            long double bestDelta = 0;
            for (int j = i + 2; j <= min(n - 1, i + W); ++j) {
                long double before = 0, after = 0;
                for (int k = i; k <= j; ++k) {
                    // The unpenalized interior-edge total is unchanged by
                    // reversal; only boundaries and carrot edges matter.
                    if (k != i && k != j && (k + 1) % 10 != 0) continue;
                    int ou = r[k], ov = r[k + 1];
                    int nu, nv;
                    if (k == i) { nu = r[i]; nv = r[j]; }
                    else if (k == j) { nu = r[i + 1]; nv = r[j + 1]; }
                    else {
                        // Position k maps to the mirrored old position inside
                        // the reversed interval.
                        nu = r[i + 1 + j - k];
                        nv = r[i + j - k];
                    }
                    auto weighted = [&](int u, int v) {
                        long double dx = (long double)a[u].x - a[v].x;
                        long double dy = (long double)a[u].y - a[v].y;
                        long double d = sqrtl(dx * dx + dy * dy);
                        if ((k + 1) % 10 == 0 && !prime[u]) d *= 1.1L;
                        return d;
                    };
                    before += weighted(ou, ov);
                    after += weighted(nu, nv);
                }
                if (after - before < bestDelta) {
                    bestDelta = after - before;
                    bestJ = j;
                }
            }
            if (bestJ != -1) reverse(r.begin() + i + 1, r.begin() + bestJ + 1);
        }
    };
    twoOpt(best);

    // A reversal cannot change the relative order of a short detour.  Probe the
    // complementary one-city relocation neighborhood, but score *all* edges
    // whose positions or endpoints change: moving a city also shifts carrot
    // positions in the intervening segment.
    auto relocate = [&](vector<int>& r) {
        const int W = 5;
        auto weighted = [&](int pos, int u, int v) {
            long double dx = (long double)a[u].x - a[v].x;
            long double dy = (long double)a[u].y - a[v].y;
            long double d = sqrtl(dx * dx + dy * dy);
            if ((pos + 1) % 10 == 0 && !prime[u]) d *= 1.1L;
            return d;
        };
        for (int p = 1; p < n; ++p) {
            int bestQ = -1;
            long double bestDelta = 0;
            for (int q = max(0, p - W); q <= min(n - 1, p + W); ++q) {
                if (q == p || q == p - 1) continue;
                int lo = (q < p ? q : p - 1);
                int hi = (q < p ? p : q);
                long double before = 0, after = 0;
                for (int k = lo; k <= hi; ++k) {
                    before += weighted(k, r[k], r[k + 1]);
                    int u, v;
                    if (q < p) {
                        auto at = [&](int z) {
                            if (z == q + 1) return r[p];
                            if (z >= q + 2 && z <= p) return r[z - 1];
                            return r[z];
                        };
                        u = at(k); v = at(k + 1);
                    } else {
                        auto at = [&](int z) {
                            if (z >= p && z < q) return r[z + 1];
                            if (z == q) return r[p];
                            return r[z];
                        };
                        u = at(k); v = at(k + 1);
                    }
                    after += weighted(k, u, v);
                }
                if (after - before < bestDelta) {
                    bestDelta = after - before;
                    bestQ = q;
                }
            }
            if (bestQ != -1) {
                if (bestQ < p)
                    rotate(r.begin() + bestQ + 1, r.begin() + p, r.begin() + p + 1);
                else
                    rotate(r.begin() + p, r.begin() + p + 1, r.begin() + bestQ + 1);
            }
        }
    };
    // Preserve the already selected route if numerical accumulation somehow
    // disagrees with the exact affected-edge deltas used by the sweep.
    vector<int> beforeRelocate = best;
    long double beforeRelocateCost = cost(best);
    relocate(best);
    if (cost(best) > beforeRelocateCost) best = move(beforeRelocate);

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
}
