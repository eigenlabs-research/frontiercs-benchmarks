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

    auto plainEdge = [&](int u, int v) -> double {
        long double dx = (long double)a[u].x - a[v].x;
        long double dy = (long double)a[u].y - a[v].y;
        return hypot((double)dx, (double)dy);
    };
    auto edge = [&](const vector<int>& r, int t) -> double { // step t, 1..n
        double z = plainEdge(r[t-1], r[t]);
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

    // A balanced kd traversal is deliberately a different representation from a space-filling
    // curve: it recursively separates alternating coordinate axes.  It is useful on
    // inputs made of a few broad clouds or long strips, where a single Hilbert orientation can
    // make an unlucky inter-cloud turn.  Coordinate signs give the possible child visit orders.
    auto tryKD = [&](int firstAxis, int sx, int sy) {
        vector<int> ids(n), seq;
        iota(ids.begin(), ids.end(), 0);
        seq.reserve(n);
        auto before = [&](int u, int v, int axis) {
            long long pu = axis ? a[u].y : a[u].x;
            long long pv = axis ? a[v].y : a[v].x;
            int sign = axis ? sy : sx;
            if (pu != pv) return sign > 0 ? pu < pv : pu > pv;
            // A deterministic secondary key keeps partitions reproducible on aligned data.
            long long qu = axis ? a[u].x : a[u].y;
            long long qv = axis ? a[v].x : a[v].y;
            return sign > 0 ? qu < qv : qu > qv;
        };
        auto build = [&](auto&& self, int l, int r, int depth) -> void {
            if (r - l <= 24) {
                int axis = (firstAxis + depth) & 1;
                sort(ids.begin() + l, ids.begin() + r,
                     [&](int u, int v) { return before(u, v, axis); });
                seq.insert(seq.end(), ids.begin() + l, ids.begin() + r);
                return;
            }
            int axis = (firstAxis + depth) & 1;
            int m = (l + r) >> 1;
            nth_element(ids.begin() + l, ids.begin() + m, ids.begin() + r,
                        [&](int u, int v) { return before(u, v, axis); });
            self(self, l, m, depth + 1);
            self(self, m, r, depth + 1);
        };
        build(build, 0, n, 0);
        for (int direction = 0; direction < 2; ++direction) {
            vector<int> r(n + 1); r[0] = r[n] = 0;
            int at = int(find(seq.begin(), seq.end(), 0) - seq.begin());
            for (int k = 1; k < n; ++k) {
                int z = direction ? (at - k + n) % n : (at + k) % n;
                r[k] = seq[z];
            }
            double kdCost = cost(r);
            if (kdCost < bestCost) { bestCost = kdCost; best.swap(r); }
        }
    };
    for (int axis = 0; axis < 2; ++axis)
        for (int sx : {-1, 1}) for (int sy : {-1, 1})
            tryKD(axis, sx, sy);

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

    // A reversal preserves a block's membership and cannot remove a one-city spur.  Relocate
    // one city within a short window instead.  The ordinary length delta changes only three
    // links; only carrot links crossed while the intervening block shifts need enumeration.
    // This is substantially cheaper than rescoring every edge in the moved interval.
    auto relocateDelta = [&](int p, int q) {
        int A = best[p - 1], X = best[p], B = best[p + 1];
        int C, D;
        if (p < q) { C = best[q]; D = best[q + 1]; }
        else       { C = best[q - 1]; D = best[q]; }
        double d = plainEdge(A, B) + plainEdge(C, X) + plainEdge(X, D)
                 - plainEdge(A, X) - plainEdge(X, B) - plainEdge(C, D);
        auto movedAt = [&](int i) {
            if (p < q) {
                if (i >= p && i < q) return best[i + 1];
                if (i == q) return best[p];
            } else {
                if (i == q) return best[p];
                if (i > q && i <= p) return best[i - 1];
            }
            return best[i];
        };
        int lo = min(p, q), hi = max(p, q) + 1;
        for (int t = ((lo + 9) / 10) * 10; t <= hi; t += 10) {
            int ou = best[t - 1], ov = best[t];
            int nu = movedAt(t - 1), nv = movedAt(t);
            if (!prime[ou]) d -= .1 * plainEdge(ou, ov);
            if (!prime[nu]) d += .1 * plainEdge(nu, nv);
        }
        return d;
    };
    const int shifts[] = {2, 3, 5, 8, 13, 21, 32};
    for (int p = 1; p < n; ++p) {
        int take = -1;
        double gain = 0;
        for (int s : shifts) for (int dir : {-1, 1}) {
            int q = p + dir * s;
            if (q < 1 || q >= n) continue;
            double d = relocateDelta(p, q);
            if (d < gain) gain = d, take = q;
        }
        if (take != -1) {
            int v = best[p];
            if (p < take) {
                for (int i = p; i < take; ++i) best[i] = best[i + 1];
            } else {
                for (int i = p; i > take; --i) best[i] = best[i - 1];
            }
            best[take] = v;
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
