#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; int id; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> a(n);
    vector<long long> X(n), Y(n);
    for (int i = 0; i < n; ++i) {
        cin >> X[i] >> Y[i];
        a[i] = {X[i], Y[i], i};
    }
    vector<char> prime(n, true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    // An adaptive kd traversal is deliberately used instead of a fixed global
    // space-filling curve.  Splits follow the longer local dimension, so a
    // thin cluster is not repeatedly cut across its short dimension.  At each
    // join either child can be reversed; choose the closest pair of exposed
    // endpoints, producing a path rather than merely a set of spatial bins.
    vector<int> order;
    order.reserve(n);
    auto dist = [&](int u, int v) {
        return hypot((double)X[u] - X[v], (double)Y[u] - Y[v]);
    };
    function<pair<int,int>(int,int)> build = [&](int lo, int hi) -> pair<int,int> {
        int m = hi - lo;
        if (m <= 24) {
            long long minx = a[lo].x, maxx = a[lo].x, miny = a[lo].y, maxy = a[lo].y;
            for (int i = lo + 1; i < hi; ++i) {
                minx = min(minx, a[i].x); maxx = max(maxx, a[i].x);
                miny = min(miny, a[i].y); maxy = max(maxy, a[i].y);
            }
            bool byX = maxx - minx >= maxy - miny;
            sort(a.begin() + lo, a.begin() + hi, [&](const Point& u, const Point& v) {
                long long pu = byX ? u.x : u.y, pv = byX ? v.x : v.y;
                return pu != pv ? pu < pv : u.id < v.id;
            });
            int begin = order.size();
            for (int i = lo; i < hi; ++i) order.push_back(a[i].id);
            return {order[begin], order.back()};
        }
        long long minx = a[lo].x, maxx = a[lo].x, miny = a[lo].y, maxy = a[lo].y;
        for (int i = lo + 1; i < hi; ++i) {
            minx = min(minx, a[i].x); maxx = max(maxx, a[i].x);
            miny = min(miny, a[i].y); maxy = max(maxy, a[i].y);
        }
        bool byX = maxx - minx >= maxy - miny;
        int mid = lo + m / 2;
        nth_element(a.begin() + lo, a.begin() + mid, a.begin() + hi, [&](const Point& u, const Point& v) {
            long long pu = byX ? u.x : u.y, pv = byX ? v.x : v.y;
            return pu != pv ? pu < pv : u.id < v.id;
        });
        int start = order.size();
        auto L = build(lo, mid);
        int cut = order.size();
        auto R = build(mid, hi);
        // Test orientations; reversing a completed subtree is exact and cheap
        // in total because each level reverses O(n) IDs.
        double best = 1e300; bool revL = false, revR = false;
        for (int x = 0; x < 2; ++x) for (int y = 0; y < 2; ++y) {
            int le = x ? L.first : L.second;
            int rs = y ? R.second : R.first;
            double d = dist(le, rs);
            if (d < best) best = d, revL = x, revR = y;
        }
        if (revL) reverse(order.begin() + start, order.begin() + cut);
        if (revR) reverse(order.begin() + cut, order.end());
        int first = revL ? L.second : L.first;
        int last = revR ? R.first : R.second;
        return {first, last};
    };
    build(0, n);

    int z = find(order.begin(), order.end(), 0) - order.begin();
    vector<int> p(n + 1);
    p[0] = p[n] = 0;
    for (int k = 1; k < n; ++k) p[k] = order[(z + k) % n];

    auto edge = [&](int t, int sa = -1, int sb = -1) {
        auto at = [&](int i) { return i == sa ? p[sb] : (i == sb ? p[sa] : p[i]); };
        int u = at(t - 1), v = at(t);
        return (t % 10 == 0 && !prime[u] ? 1.1 : 1.0) * dist(u, v);
    };
    auto swapGain = [&](int i, int j) {
        int ts[4] = {i, i + 1, j, j + 1};
        sort(ts, ts + 4);
        double before = 0, after = 0;
        for (int q = 0; q < 4; ++q) if ((!q || ts[q] != ts[q-1]) && ts[q] >= 1 && ts[q] <= n) {
            before += edge(ts[q]); after += edge(ts[q], i, j);
        }
        return before - after;
    };
    for (int pass = 0; pass < 2; ++pass)
        for (int pos = 9; pos < n; pos += 10) if (!prime[p[pos]]) {
            int best = -1; double gain = 1e-7;
            for (int q = max(1, pos - 8); q <= min(n - 1, pos + 8); ++q) if (prime[p[q]]) {
                double g = swapGain(pos, q);
                if (g > gain) gain = g, best = q;
            }
            if (best >= 0) swap(p[pos], p[best]);
        }
    for (int pass = 0; pass < 2; ++pass)
        for (int i = 1; i + 1 < n; ++i) if (swapGain(i, i + 1) > 1e-7) swap(p[i], p[i + 1]);

    // Long reversal checks can otherwise call hypot O(length/10) times each.
    // Cache current edge lengths; an internal reversed edge has the same length.
    vector<double> edgeLen(n + 1);
    for (int t = 1; t <= n; ++t) edgeLen[t] = dist(p[t-1], p[t]);
    auto reversalGain = [&](int l, int r) {
        double before = edgeLen[l] + edgeLen[r+1];
        double after = dist(p[l-1], p[r]) + dist(p[l], p[r+1]);
        for (int t = ((l + 9) / 10) * 10; t <= r + 1; t += 10) {
            int ou = p[t-1];
            int nu = (t == l ? p[l-1] : t == r+1 ? p[l] : p[l+r-t+1]);
            double oldLen = edgeLen[t];
            double newLen = (t == l ? dist(p[l-1], p[r]) :
                             t == r+1 ? dist(p[l], p[r+1]) : edgeLen[l+r-t+1]);
            if (!prime[ou]) before += .1 * oldLen;
            if (!prime[nu]) after += .1 * newLen;
        }
        return before - after;
    };
    // Keep the cheap short cleanup from the incumbent before considering
    // non-local repairs.
    const int lengths[] = {2,3,4,5,6,7,8,10,12,16,20,24,28,32,40,48,64};
    for (int l = 1; l < n; ++l) {
        int best = -1; double gain = 1e-7;
        for (int len : lengths) if (l + len - 1 < n) {
            int r = l + len - 1; double g = reversalGain(l, r);
            if (g > gain) gain = g, best = r;
        }
        if (best >= 0) {
            reverse(p.begin() + l, p.begin() + best + 1);
            for (int t = l; t <= best + 1; ++t) edgeLen[t] = dist(p[t-1], p[t]);
        }
    }

    // The construction's only non-local decision is the closing edge of the
    // kd path.  Look for long 2-opt moves whose new endpoint is close in the
    // input's x order.  Unlike the former fixed-length pass, these moves can
    // repair a crossing between distant portions of the tour.  A small exact
    // evaluation budget keeps the carrot-aware reversal calculation bounded.
    vector<int> where(n);
    for (int i = 1; i < n; ++i) where[p[i]] = i;
    int exactChecks = 0, accepted = 0;
    const int CHECK_LIMIT = 2400, ACCEPT_LIMIT = 120, WINDOW = 24;
    for (int l = 1; l < n && exactChecks < CHECK_LIMIT && accepted < ACCEPT_LIMIT; ++l) {
        int u = p[l - 1];
        for (int delta = 1; delta <= WINDOW && exactChecks < CHECK_LIMIT && accepted < ACCEPT_LIMIT; ++delta) {
            int cand[2] = {u - delta, u + delta};
            for (int q = 0; q < 2 && exactChecks < CHECK_LIMIT && accepted < ACCEPT_LIMIT; ++q) {
                int v = cand[q];
                if (v <= 0 || v >= n) continue;
                int r = where[v];
                if (r <= l + 63 || r >= n) continue;
                // Most candidate pairs cannot improve even the unpenalized
                // tour, so reject them before the O(length/10) exact check.
                double geometric = dist(p[l-1], p[l]) + dist(p[r], p[r+1])
                                 - dist(p[l-1], p[r]) - dist(p[l], p[r+1]);
                if (geometric <= 1e-7) continue;
                ++exactChecks;
                if (reversalGain(l, r) > 1e-7) {
                    reverse(p.begin() + l, p.begin() + r + 1);
                    for (int k = l; k <= r; ++k) where[p[k]] = k;
                    for (int t = l; t <= r + 1; ++t) edgeLen[t] = dist(p[t-1], p[t]);
                    ++accepted;
                }
            }
        }
    }
    cout << n + 1 << '\n';
    for (int v : p) cout << v << '\n';
}
