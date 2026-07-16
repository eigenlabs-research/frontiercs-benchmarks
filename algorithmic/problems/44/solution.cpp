#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert order has the same O(n log n) construction cost as Morton order,
// but does not make a long jump between adjacent Morton quadrants.
static uint64_t hilbertKey(long long x, long long y) {
    uint32_t a = (uint32_t)(x + 1000000000LL);
    uint32_t b = (uint32_t)(y + 1000000000LL);
    const uint32_t TOP = 1U << 30; // shifted coordinates use at most 31 bits

    // Convert the two coordinates to the Hilbert transpose (Skilling's
    // in-place inverse Gray transform), then interleave its high bits first.
    for (uint32_t q = TOP; q > 1; q >>= 1) {
        uint32_t p = q - 1;
        if (b & q) {
            a ^= p;
        } else {
            uint32_t t = (a ^ b) & p;
            a ^= t;
            b ^= t;
        }
    }
    b ^= a;
    uint32_t t = 0;
    for (uint32_t q = TOP; q > 1; q >>= 1)
        if (b & q) t ^= q - 1;
    a ^= t;
    b ^= t;

    uint64_t z = 0;
    for (int bit = 30; bit >= 0; --bit) {
        z = (z << 1) | ((a >> bit) & 1U);
        z = (z << 1) | ((b >> bit) & 1U);
    }
    return z;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    vector<Point> pt(n);
    for (auto &p : pt) cin >> p.x >> p.y;

    vector<char> prime(n, true);
    if (n > 0) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i])
            for (int j = i * i; j < n; j += i) prime[j] = false;

    auto edgeCost = [&](const vector<int>& r, int t) {
        int a = r[t - 1], b = r[t];
        double dx = (double)pt[a].x - pt[b].x;
        double dy = (double)pt[a].y - pt[b].y;
        double d = hypot(dx, dy);
        return (t % 10 == 0 && !prime[a]) ? 1.1 * d : d;
    };
    auto routeCost = [&](const vector<int>& r) {
        double ans = 0;
        for (int t = 1; t <= n; ++t) ans += edgeCost(r, t);
        return ans;
    };

    // Keep the supplied monotone-x path as a safe fallback.
    vector<int> baseline(n + 1);
    iota(baseline.begin(), baseline.end() - 1, 0);
    baseline[n] = 0;
    double baseCost = routeCost(baseline);

    vector<uint64_t> key(n);
    for (int i = 0; i < n; ++i) key[i] = hilbertKey(pt[i].x, pt[i].y);
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int u, int v) {
        return key[u] != key[v] ? key[u] < key[v] : u < v;
    });
    int zero = (int)(find(order.begin(), order.end(), 0) - order.begin());
    auto makeRoute = [&](int dir) {
        vector<int> r;
        r.reserve(n + 1);
        r.push_back(0);
        for (int k = 1; k < n; ++k)
            r.push_back(order[(zero + dir * k + n * 2) % n]);
        r.push_back(0);
        return r;
    };
    vector<int> route = makeRoute(1), reverseRoute = makeRoute(-1);
    if (routeCost(reverseRoute) < routeCost(route)) route.swap(reverseRoute);

    // A directional sweep is deliberately a different construction family from
    // a space-filling curve.  It is particularly effective for long thin clouds
    // and parallel bands: points are swept in projection order rather than
    // repeatedly changing quadrant.  Keep only its exact-cost winner, so this
    // portfolio cannot displace the established initializer gratuitously.
    double chosenCost = routeCost(route);
    const int directions[][2] = {{0,1}, {1,1}, {1,-1}, {2,1}, {1,2}};
    for (auto &dir : directions) {
        vector<int> sweep;
        sweep.reserve(n - 1);
        for (int i = 1; i < n; ++i) sweep.push_back(i);
        long long dx = dir[0], dy = dir[1];
        sort(sweep.begin(), sweep.end(), [&](int a, int b) {
            __int128 pa = (__int128)pt[a].x * dx + (__int128)pt[a].y * dy;
            __int128 pb = (__int128)pt[b].x * dx + (__int128)pt[b].y * dy;
            return pa != pb ? pa < pb : a < b;
        });
        for (int backwards = 0; backwards < 2; ++backwards) {
            vector<int> candidate;
            candidate.reserve(n + 1);
            candidate.push_back(0);
            if (!backwards) candidate.insert(candidate.end(), sweep.begin(), sweep.end());
            else candidate.insert(candidate.end(), sweep.rbegin(), sweep.rend());
            candidate.push_back(0);
            double c = routeCost(candidate);
            if (c + 1e-6 < chosenCost) {
                chosenCost = c;
                route.swap(candidate);
            }
        }
    }

    // A strip snake is an ablation of a plain projection sweep.  Within every
    // primary-coordinate strip it reverses the secondary ordering, avoiding the
    // repeated band crossings which a one-dimensional sort makes on lattices.
    // Test both axes; each is scored under the complete carrot contract before
    // it is allowed to replace the established portfolio winner.
    for (int axis = 0; axis < 2; ++axis) {
        vector<long long> primary(n), secondary(n);
        long long lo = 0, hi = 0;
        for (int i = 0; i < n; ++i) {
            primary[i] = axis == 0 ? pt[i].x : pt[i].y;
            secondary[i] = axis == 0 ? pt[i].y : pt[i].x;
            if (i == 0) lo = hi = primary[i];
            else { lo = min(lo, primary[i]); hi = max(hi, primary[i]); }
        }
        int strips = max(2, min(512, (int)sqrt((double)n) + 1));
        auto stripOf = [&](int v) {
            if (hi == lo) return 0;
            return (int)((primary[v] - lo) * strips / (hi - lo + 1));
        };
        vector<int> snake;
        snake.reserve(n - 1);
        for (int i = 1; i < n; ++i) snake.push_back(i);
        sort(snake.begin(), snake.end(), [&](int a, int b) {
            int sa = stripOf(a), sb = stripOf(b);
            if (sa != sb) return sa < sb;
            if (secondary[a] != secondary[b])
                return (sa & 1) ? secondary[a] > secondary[b] : secondary[a] < secondary[b];
            return a < b;
        });
        for (int backwards = 0; backwards < 2; ++backwards) {
            vector<int> candidate;
            candidate.reserve(n + 1);
            candidate.push_back(0);
            if (!backwards) candidate.insert(candidate.end(), snake.begin(), snake.end());
            else candidate.insert(candidate.end(), snake.rbegin(), snake.rend());
            candidate.push_back(0);
            double c = routeCost(candidate);
            if (c + 1e-6 < chosenCost) {
                chosenCost = c;
                route.swap(candidate);
            }
        }
    }

    // Bounded 2-opt is evaluated with the actual position-dependent carrot cost.
    // Run it again after prime placement: swaps can expose a short crossing.
    const int WINDOW = 8;
    vector<double> ec(n + 1);
    auto runTwoOpt = [&]() {
        for (int t = 1; t <= n; ++t) ec[t] = edgeCost(route, t);
        for (int l = 1; l <= n - 2; ++l) {
            int bestR = -1;
            double bestDelta = 0.0;
            for (int r = l + 1; r <= min(n - 1, l + WINDOW); ++r) {
                double oldPart = 0.0, newPart = 0.0;
                for (int t = l; t <= r + 1; ++t) {
                    oldPart += ec[t];
                    int ia = t - 1, ib = t;
                    int a = (ia >= l && ia <= r) ? route[l + r - ia] : route[ia];
                    int b = (ib >= l && ib <= r) ? route[l + r - ib] : route[ib];
                    double dx = (double)pt[a].x - pt[b].x;
                    double dy = (double)pt[a].y - pt[b].y;
                    double d = hypot(dx, dy);
                    newPart += (t % 10 == 0 && !prime[a]) ? 1.1 * d : d;
                }
                if (newPart - oldPart < bestDelta - 1e-7) {
                    bestDelta = newPart - oldPart;
                    bestR = r;
                }
            }
            if (bestR != -1) {
                reverse(route.begin() + l, route.begin() + bestR + 1);
                for (int t = l; t <= bestR + 1; ++t) ec[t] = edgeCost(route, t);
            }
        }
    };
    runTwoOpt();

    // At each costly step, try nearby prime cities.  A swap is retained only if
    // all affected geometric and carrot terms improve exactly.
    for (int p = 9; p <= n - 1; p += 10) {
        if (prime[route[p]]) continue;
        int bestQ = -1;
        double bestDelta = 0.0;
        for (int q = max(1, p - WINDOW); q <= min(n - 1, p + WINDOW); ++q) {
            if (!prime[route[q]]) continue;
            vector<int> affected = {p, p + 1, q, q + 1};
            sort(affected.begin(), affected.end());
            affected.erase(unique(affected.begin(), affected.end()), affected.end());
            double oldPart = 0.0;
            for (int t : affected) oldPart += ec[t];
            swap(route[p], route[q]);
            double newPart = 0.0;
            for (int t : affected) newPart += edgeCost(route, t);
            swap(route[p], route[q]);
            if (newPart - oldPart < bestDelta - 1e-7) {
                bestDelta = newPart - oldPart;
                bestQ = q;
            }
        }
        if (bestQ != -1) {
            swap(route[p], route[bestQ]);
            for (int t : {p, p + 1, bestQ, bestQ + 1}) ec[t] = edgeCost(route, t);
        }
    }

    // Prime moves are exact improvements themselves, but can disrupt the local
    // geometry optimized above.  Restore the same bounded geometric optimum.
    runTwoOpt();

    if (routeCost(route) > baseCost + 1e-6) route.swap(baseline);
    cout << n + 1 << '\n';
    for (int v : route) cout << v << '\n';
    return 0;
}
