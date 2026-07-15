#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    vector<Point> p(n);
    long long xmin = LLONG_MAX, xmax = LLONG_MIN, ymin = LLONG_MAX, ymax = LLONG_MIN;
    for (auto &q : p) {
        cin >> q.x >> q.y;
        xmin = min(xmin, q.x); xmax = max(xmax, q.x);
        ymin = min(ymin, q.y); ymax = max(ymax, q.y);
    }

    vector<char> prime(n, true);
    if (n > 0) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    const unsigned long long LIM = (1ULL << 31) - 1;
    auto scale = [&](long long v, long long lo, long long hi) -> unsigned int {
        if (hi == lo) return 0;
        unsigned long long d = (unsigned long long)(v - lo);
        unsigned long long range = (unsigned long long)(hi - lo);
        return (unsigned int)(d * LIM / range);
    };
    vector<unsigned int> bx(n), by(n);
    for (int i = 0; i < n; ++i) {
        bx[i] = scale(p[i].x, xmin, xmax);
        by[i] = scale(p[i].y, ymin, ymax);
    }

    // Index on a Hilbert space-filling curve.  Nearby points tend to remain
    // nearby in this order, unlike the given (x-only) ordering.
    auto hilbert = [](unsigned int x, unsigned int y) -> unsigned long long {
        unsigned long long ans = 0;
        for (int s = 30; s >= 0; --s) {
            unsigned int rx = (x >> s) & 1U, ry = (y >> s) & 1U;
            ans |= (unsigned long long)((3U * rx) ^ ry) << (2 * s);
            if (!ry) {
                if (rx) {
                    unsigned int mask = (1U << (s + 1)) - 1U;
                    x = mask - x; y = mask - y;
                }
                swap(x, y);
            }
        }
        return ans;
    };

    auto edge = [&](const vector<int>& r, int t) -> double {
        int a = r[t - 1], b = r[t];
        // Coordinate differences are bounded by 2e9, so this is well within
        // double's precise, safe range and is appreciably cheaper than hypot
        // in the many evaluations made by the local search.
        double dx = (double)p[a].x - p[b].x, dy = (double)p[a].y - p[b].y;
        double d = sqrt(dx * dx + dy * dy);
        if (t % 10 == 0 && !prime[a]) d *= 1.1;
        return d;
    };
    auto cost = [&](const vector<int>& r) -> double {
        double z = 0;
        for (int t = 1; t <= n; ++t) z += edge(r, t);
        return z;
    };

    vector<vector<int>> candidates;
    for (int mode = 0; mode < 8; ++mode) {
        vector<pair<unsigned long long,int>> a;
        a.reserve(n);
        bool sw = mode & 4, fx = mode & 1, fy = mode & 2;
        for (int i = 0; i < n; ++i) {
            unsigned int x = bx[i], y = by[i];
            if (sw) swap(x, y);
            if (fx) x = (unsigned int)LIM - x;
            if (fy) y = (unsigned int)LIM - y;
            a.push_back({hilbert(x, y), i});
        }
        sort(a.begin(), a.end());
        vector<int> ord(n);
        int at0 = 0;
        for (int i = 0; i < n; ++i) { ord[i] = a[i].second; if (ord[i] == 0) at0 = i; }

        auto add = [&](const vector<int>& middle) {
            vector<int> r; r.reserve(n + 1); r.push_back(0);
            for (int v : middle) if (v != 0) r.push_back(v);
            r.push_back(0); candidates.push_back(move(r));
        };
        add(ord);
        reverse(ord.begin(), ord.end());
        add(ord);
        // Also regard the curve as a cycle and cut it at the fixed depot.
        // This makes both edges incident to city 0 locally spatial whenever
        // city 0 is in the interior of the curve order.
        ord.resize(n);
        for (int i = 0; i < n; ++i) ord[i] = a[i].second;
        vector<int> cyc;
        for (int i = at0 + 1; i < n; ++i) cyc.push_back(ord[i]);
        for (int i = 0; i < at0; ++i) cyc.push_back(ord[i]);
        add(cyc);
        reverse(cyc.begin(), cyc.end());
        add(cyc);
    }

    vector<pair<double,int>> ranked;
    ranked.reserve(candidates.size());
    for (int i = 0; i < (int)candidates.size(); ++i) ranked.push_back({cost(candidates[i]), i});
    sort(ranked.begin(), ranked.end());

    // The carrot rule only affects sparse positions.  On the best geometric
    // tours, make inexpensive local exchanges to place a prime at such a
    // source position, accepting an exchange only when its exact cost drops.
    int take = min(4, (int)ranked.size());
    for (int q = 0; q < take; ++q) {
        vector<int>& r = candidates[ranked[q].second];
        for (int pos = 9; pos < n; pos += 10) {
            if (prime[r[pos]]) continue;
            int best = -1;
            double bestDelta = 0;
            for (int k = max(1, pos - 6); k <= min(n - 1, pos + 6); ++k) {
                if (!prime[r[k]]) continue;
                vector<int> ts = {pos, pos + 1, k, k + 1};
                sort(ts.begin(), ts.end()); ts.erase(unique(ts.begin(), ts.end()), ts.end());
                double before = 0, after = 0;
                for (int t : ts) if (t >= 1 && t <= n) before += edge(r, t);
                swap(r[pos], r[k]);
                for (int t : ts) if (t >= 1 && t <= n) after += edge(r, t);
                swap(r[pos], r[k]);
                if (after - before < bestDelta) { bestDelta = after - before; best = k; }
            }
            if (best != -1) swap(r[pos], r[best]);
        }
    }

    // Broader local 2-opt.  A reversal preserves every *unweighted* internal
    // edge (distances are symmetric), so its delta consists of two boundary
    // edges plus only the occasional carrot edges inside the block.  Computing
    // that delta directly permits a substantially wider search than repeatedly
    // scoring the whole reversed block.
    int polish = min(2, (int)ranked.size());
    const int WINDOW = 24;
    auto d = [&](int a, int b) -> double {
        double dx = (double)p[a].x - p[b].x, dy = (double)p[a].y - p[b].y;
        return sqrt(dx * dx + dy * dy);
    };
    for (int q = 0; q < polish; ++q) {
        vector<int>& r = candidates[ranked[q].second];
        // One pass can create a newly profitable reversal just behind the
        // changed block.  A second bounded sweep captures those dependencies
        // while retaining linear time for a fixed window.
        for (int pass = 0; pass < 2; ++pass) {
            bool changed = false;
            for (int i = 1; i < n - 1; ++i) {
                int bestj = -1;
                double bestDelta = 0.0;
                int hi = min(n - 1, i + WINDOW - 1);
                for (int j = i + 1; j <= hi; ++j) {
                    double oldL = d(r[i - 1], r[i]);
                    double oldR = d(r[j], r[j + 1]);
                    double newL = d(r[i - 1], r[j]);
                    double newR = d(r[i], r[j + 1]);
                    double delta = newL + newR - oldL - oldR;
                    // The boundary sources at i and j+1 may change their penalty.
                    if (i % 10 == 0 && !prime[r[i - 1]]) delta += .1 * (newL - oldL);
                    if ((j + 1) % 10 == 0) {
                        delta += .1 * ((!prime[r[i]]) * newR - (!prime[r[j]]) * oldR);
                    }
                    // At an internal step t, the reversed edge is old edge q,
                    // traversed from r[q] to r[q-1], q = i+j-t+1.
                    int first = ((i + 10) / 10) * 10;
                    for (int t = first; t <= j; t += 10) {
                        int k = i + j - t + 1;
                        double before = d(r[t - 1], r[t]);
                        double after = d(r[k - 1], r[k]);
                        delta += .1 * ((!prime[r[k]]) * after - (!prime[r[t - 1]]) * before);
                    }
                    if (delta < bestDelta) bestDelta = delta, bestj = j;
                }
                if (bestj != -1) {
                    reverse(r.begin() + i, r.begin() + bestj + 1);
                    changed = true;
                }
            }
            if (!changed) break;
        }
    }

    int answer = 0;
    double best = cost(candidates[0]);
    for (int i = 1; i < (int)candidates.size(); ++i) {
        double v = cost(candidates[i]);
        if (v < best) best = v, answer = i;
    }
    cout << n + 1 << '\n';
    for (int v : candidates[answer]) cout << v << '\n';
}
