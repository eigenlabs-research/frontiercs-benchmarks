#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

static uint64_t hilbert(uint32_t x, uint32_t y) {
    uint64_t d = 0;
    const uint32_t mask = (1u << 31) - 1;
    for (uint32_t s = 1u << 30; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * s * ((3 * rx) ^ ry);
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
    vector<Point> a(n);
    for (auto &p : a) cin >> p.x >> p.y;

    vector<char> prime(max(2, n), true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    auto edge = [&](const vector<int>& r, int t) -> double { // step t, 1..n
        long double dx = (long double)a[r[t-1]].x - a[r[t]].x;
        long double dy = (long double)a[r[t-1]].y - a[r[t]].y;
        double z = hypot((double)dx, (double)dy);
        if (t % 10 == 0 && !prime[r[t-1]]) z *= 1.1;
        return z;
    };
    auto cost = [&](const vector<int>& r) {
        double ans = 0;
        for (int t = 1; t <= n; ++t) ans += edge(r, t);
        return ans;
    };

    // The cyclic x-order is retained as a safe candidate; reversing it changes carrot slots.
    vector<int> best(n + 1);
    iota(best.begin(), best.end() - 1, 0);
    best[n] = 0;
    double bestCost = cost(best);
    vector<int> rev(n + 1); rev[0] = rev[n] = 0;
    for (int i = 1; i < n; ++i) rev[i] = n - i;
    double c = cost(rev);
    if (c < bestCost) { bestCost = c; best.swap(rev); }

    const uint32_t LIM = (1u << 31) - 1;
    vector<pair<uint64_t,int>> order(n);
    // Dihedral variants avoid committing to an arbitrary Hilbert orientation.
    for (int type = 0; type < 8; ++type) {
        for (int i = 0; i < n; ++i) {
            uint32_t u = uint32_t(a[i].x + 1000000000LL);
            uint32_t v = uint32_t(a[i].y + 1000000000LL), X, Y;
            switch (type) {
                case 0: X=u; Y=v; break; case 1: X=LIM-u; Y=v; break;
                case 2: X=u; Y=LIM-v; break; case 3: X=LIM-u; Y=LIM-v; break;
                case 4: X=v; Y=u; break; case 5: X=LIM-v; Y=u; break;
                case 6: X=v; Y=LIM-u; break; default: X=LIM-v; Y=LIM-u;
            }
            order[i] = {hilbert(X,Y), i};
        }
        sort(order.begin(), order.end());
        int at = 0;
        while (order[at].second != 0) ++at;
        vector<int> r(n + 1); r[0] = r[n] = 0;
        for (int k = 1; k < n; ++k) r[k] = order[(at + k) % n].second;
        c = cost(r);
        if (c < bestCost) { bestCost = c; best.swap(r); }
    }

    // Exact local objective updates: adjacent moves repair small curve artifacts and account for
    // both geometry and the position-dependent multiplier.
    auto swapDelta = [&](vector<int>& r, int p, int q) {
        int ts[4] = {p, p+1, q, q+1};
        sort(ts, ts+4);
        double old = 0, neu = 0;
        int last = -1;
        for (int z : ts) if (z != last && z >= 1 && z <= n) { old += edge(r,z); last=z; }
        swap(r[p], r[q]);
        last = -1;
        for (int z : ts) if (z != last && z >= 1 && z <= n) { neu += edge(r,z); last=z; }
        swap(r[p], r[q]);
        return neu - old;
    };
    auto improveAdjacent = [&]() {
        bool changed = false;
        for (int p = 1; p + 1 < n; ++p) {
            double d = swapDelta(best, p, p+1);
            if (d < -1e-7) { swap(best[p], best[p+1]); bestCost += d; changed = true; }
        }
        return changed;
    };
    improveAdjacent(); improveAdjacent();

    // Adjacent exchanges cannot remove a small crossing without first taking a worse step.
    // Test that limitation with short 2-opt reversals.  Unlike ordinary TSP 2-opt, reversing
    // changes which city supplies every carrot slot inside the segment, so score all affected
    // edges explicitly.  The small bound keeps this linear-time after the Hilbert sorts.
    auto reverseDelta = [&](vector<int>& r, int p, int q) {
        double old = 0, neu = 0;
        for (int t = p; t <= q + 1; ++t) old += edge(r, t);
        reverse(r.begin() + p, r.begin() + q + 1);
        for (int t = p; t <= q + 1; ++t) neu += edge(r, t);
        reverse(r.begin() + p, r.begin() + q + 1);
        return neu - old;
    };
    for (int p = 1; p + 2 < n; ++p) {
        int take = -1;
        double gain = 0;
        for (int len = 3; len <= 10 && p + len - 1 < n; ++len) {
            int q = p + len - 1;
            double d = reverseDelta(best, p, q);
            if (d < gain) gain = d, take = q;
        }
        if (take != -1) {
            reverse(best.begin() + p, best.begin() + take + 1);
            bestCost += gain;
        }
    }

    // At a carrot slot, look for a nearby prime ID.  Restricting exchanges to a local window
    // keeps the spatial tour intact while allowing the exact objective to exploit prime sources.
    for (int p = 9; p < n; p += 10) {
        if (prime[best[p]]) continue;
        int lo = max(1, p - 30), hi = min(n - 1, p + 30), take = -1;
        double gain = 0;
        for (int q = lo; q <= hi; ++q) if (prime[best[q]]) {
            double d = swapDelta(best, p, q);
            if (d < gain) gain = d, take = q;
        }
        if (take != -1) { swap(best[p], best[take]); bestCost += gain; }
    }
    improveAdjacent(); improveAdjacent();

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
    return 0;
}
