#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert order is a locality-preserving ordering of the plane.  Unlike the
// input order, it does not assume that nearby x coordinates are nearby cities.
static unsigned long long hilbert(long long x, long long y) {
    const unsigned long long side = 1ULL << 31;
    unsigned long long d = 0;
    for (unsigned long long s = side >> 1; s; s >>= 1) {
        unsigned long long rx = (x & (long long)s) != 0, ry = (y & (long long)s) != 0;
        d += s * s * ((3 * rx) ^ ry);
        if (!ry) {
            // This is the standard Hilbert rotation at the current scale.
            // Coordinates can temporarily become negative; their low bits are
            // exactly what the following scales need.
            if (rx) { x = (long long)s - 1 - x; y = (long long)s - 1 - y; }
            swap(x, y);
        }
    }
    return d;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    const auto started = chrono::steady_clock::now();
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
    if (n > 0) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    const unsigned long long LIM = (1ULL << 31) - 1;
    vector<unsigned long long> xx(n), yy(n);
    long long dx = maxx - minx, dy = maxy - miny;
    for (int i = 0; i < n; ++i) {
        xx[i] = dx ? (unsigned long long)((__int128)(a[i].x - minx) * LIM / dx) : 0;
        yy[i] = dy ? (unsigned long long)((__int128)(a[i].y - miny) * LIM / dy) : 0;
    }
    auto distance = [&](int u, int v) {
        return hypot((double)a[u].x - a[v].x, (double)a[u].y - a[v].y);
    };
    auto cost = [&](const vector<int>& r) {
        double z = 0;
        for (int j = 0; j < n; ++j) {
            double w = ((j + 1) % 10 == 0 && !prime[r[j]]) ? 1.1 : 1.0;
            z += w * distance(r[j], r[j + 1]);
        }
        return z;
    };

    // The eight rotations/reflections are all the same construction, but a
    // fixed North Pole can cut different Hilbert cycles at different places.
    vector<int> best;
    double bestCost = numeric_limits<double>::infinity();
    for (int mode = 0; mode < 8; ++mode) {
        vector<pair<unsigned long long,int>> v;
        v.reserve(n);
        for (int i = 0; i < n; ++i) {
            unsigned long long x = xx[i], y = yy[i];
            if (mode & 1) x = LIM - x;
            if (mode & 2) y = LIM - y;
            if (mode & 4) swap(x, y);
            v.push_back({hilbert(x, y), i});
        }
        sort(v.begin(), v.end());
        int at = 0;
        while (v[at].second != 0) ++at;
        vector<int> r; r.reserve(n + 1);
        for (int k = 0; k < n; ++k) r.push_back(v[(at + k) % n].second);
        r.push_back(0);
        // Reversing a cycle is another valid traversal, and penalty positions
        // make its objective value different.
        double c = cost(r);
        if (c < bestCost) { bestCost = c; best = r; }
        reverse(r.begin() + 1, r.end() - 1);
        c = cost(r);
        if (c < bestCost) { bestCost = c; best = r; }
    }

    auto localCost = [&](const vector<int>& r, const vector<int>& pos) {
        double z = 0;
        for (int j : pos) if (0 <= j && j < n) {
            double w = ((j + 1) % 10 == 0 && !prime[r[j]]) ? 1.1 : 1.0;
            z += w * distance(r[j], r[j + 1]);
        }
        return z;
    };
    auto trySwap = [&](int u, int v) {
        if (u == v || u <= 0 || v <= 0 || u >= n || v >= n) return;
        vector<int> q = {u-1, u, v-1, v};
        sort(q.begin(), q.end()); q.erase(unique(q.begin(), q.end()), q.end());
        double before = localCost(best, q);
        swap(best[u], best[v]);
        double after = localCost(best, q);
        if (after >= before) swap(best[u], best[v]);
    };

    // At a tenth departure, cheaply move a nearby prime into that position
    // only when the complete affected local objective improves.
    for (int p = 9; p < n; p += 10) if (!prime[best[p]]) {
        int chosen = -1;
        double chosenDelta = 0;
        for (int d = -12; d <= 12; ++d) {
            int q = p + d;
            if (q <= 0 || q >= n || !prime[best[q]]) continue;
            vector<int> e = {p-1, p, q-1, q};
            sort(e.begin(), e.end()); e.erase(unique(e.begin(), e.end()), e.end());
            double before = localCost(best, e);
            swap(best[p], best[q]);
            double after = localCost(best, e);
            swap(best[p], best[q]);
            if (after - before < chosenDelta) { chosenDelta = after - before; chosen = q; }
        }
        if (chosen != -1) trySwap(p, chosen);
    }
    // Discriminating experiment: a geometrically shorter 2-opt move is not
    // necessarily better here, since reversing a segment changes which city
    // departs at its tenth positions.  Evaluate the *exact* delta.  The
    // unpenalized internal edges cancel under reversal, so only the two ends
    // and the (few) tenth positions inside the short segment need inspection.
    auto tryTwoOpt = [&](int i, int q) {
        // Replace (i-1,i),(q,q+1) by (i-1,q),(i,q+1), reversing [i,q].
        int a0 = best[i - 1], b0 = best[i], c0 = best[q], d0 = best[q + 1];
        auto one = [&](int j, int u, int v) {
            double w = ((j + 1) % 10 == 0 && !prime[u]) ? 1.1 : 1.0;
            return w * distance(u, v);
        };
        double before = one(i - 1, a0, b0) + one(q, c0, d0);
        double after  = one(i - 1, a0, c0) + one(q, b0, d0);
        int first = i + (9 - (i % 10) + 10) % 10;
        for (int j = first; j < q; j += 10) {
            before += one(j, best[j], best[j + 1]);
            // At position j after reversal the old edge indexed i+q-j-1
            // is traversed in the opposite direction.
            after += one(j, best[i + q - j], best[i + q - j - 1]);
        }
        if (after < before) reverse(best.begin() + i, best.begin() + q + 1);
    };

    // Discriminating runtime/quality experiment: the prior width of six may
    // leave ordinary short crossings behind, whereas every accepted wider
    // move is an exact improvement.  Use width twelve, but make the search
    // time-bounded: if spatial sorting was expensive on a largest instance,
    // immediately retain the already-valid Hilbert tour rather than risking
    // a timeout (the competing explanation for a lost score).
    const int WINDOW = 12;
    const auto refineDeadline = started + chrono::milliseconds(1650);
    bool outOfRefinementTime = false;
    for (int i = 1; i + 1 < n && !outOfRefinementTime; ++i) {
        if (chrono::steady_clock::now() >= refineDeadline) break;
        for (int q = i + 1; q <= min(n - 1, i + WINDOW); ++q) {
            tryTwoOpt(i, q);
            if ((q & 31) == 0 && chrono::steady_clock::now() >= refineDeadline) {
                outOfRefinementTime = true;
                break;
            }
        }
    }

    cout << n + 1 << '\n';
    for (int id : best) cout << id << '\n';
}
