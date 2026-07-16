#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };
static vector<Point> pt;
static vector<char> primeCity;
static int n;

// A deletion-aware kd tree lets a greedy tour use actual geometric proximity,
// rather than assuming that a space-filling-curve neighbor is always useful.
struct KDNode {
    int city, left = -1, right = -1, parent = -1, cnt = 1;
    long long lx, rx, ly, ry;
    bool alive = true;
};
static vector<KDNode> kd;

static int buildKD(vector<int>& a, int l, int r, int par) {
    if (l >= r) return -1;
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (int i = l; i < r; ++i) {
        minx = min(minx, pt[a[i]].x); maxx = max(maxx, pt[a[i]].x);
        miny = min(miny, pt[a[i]].y); maxy = max(maxy, pt[a[i]].y);
    }
    int axis = (maxx - minx >= maxy - miny) ? 0 : 1;
    int m = (l + r) / 2;
    nth_element(a.begin() + l, a.begin() + m, a.begin() + r, [&](int u, int v) {
        return axis ? pt[u].y < pt[v].y : pt[u].x < pt[v].x;
    });
    int v = (int)kd.size();
    kd.push_back({a[m], -1, -1, par, r - l, minx, maxx, miny, maxy, true});
    int left = buildKD(a, l, m, v), right = buildKD(a, m + 1, r, v);
    kd[v].left = left; kd[v].right = right;
    return v;
}
static inline double boxDist2(const KDNode& q, const Point& p) {
    double dx = p.x < q.lx ? (double)q.lx - p.x : (p.x > q.rx ? (double)p.x - q.rx : 0.0);
    double dy = p.y < q.ly ? (double)q.ly - p.y : (p.y > q.ry ? (double)p.y - q.ry : 0.0);
    return dx * dx + dy * dy;
}
static void nearestKD(int v, const Point& from, int& best, double& bestD) {
    if (v < 0 || kd[v].cnt == 0 || boxDist2(kd[v], from) >= bestD) return;
    KDNode& q = kd[v];
    if (q.alive) {
        double dx = (double)pt[q.city].x - from.x, dy = (double)pt[q.city].y - from.y;
        double d = dx * dx + dy * dy;
        if (d < bestD) bestD = d, best = q.city;
    }
    int a = q.left, b = q.right;
    if (a >= 0 && b >= 0 && boxDist2(kd[b], from) < boxDist2(kd[a], from)) swap(a, b);
    nearestKD(a, from, best, bestD);
    nearestKD(b, from, best, bestD);
}
static void eraseKD(int v) {
    for (; v >= 0; v = kd[v].parent) --kd[v].cnt;
}
static vector<int> greedyNearestTour() {
    vector<int> ids(n), where(n);
    iota(ids.begin(), ids.end(), 0);
    kd.clear(); kd.reserve(n);
    int root = buildKD(ids, 0, n, -1);
    for (int i = 0; i < n; ++i) where[kd[i].city] = i;
    vector<int> r(n + 1); r[0] = r[n] = 0;
    int cur = 0;
    kd[where[0]].alive = false; eraseKD(where[0]);
    for (int pos = 1; pos < n; ++pos) {
        int nxt = -1; double d = numeric_limits<double>::infinity();
        nearestKD(root, pt[cur], nxt, d);
        r[pos] = cur = nxt;
        kd[where[cur]].alive = false; eraseKD(where[cur]);
    }
    return r;
}

static inline double edgeCost(const vector<int>& r, int t) {
    int a = r[t - 1], b = r[t];
    double d = hypot((double)pt[a].x - pt[b].x, (double)pt[a].y - pt[b].y);
    return (t % 10 == 0 && !primeCity[a]) ? 1.1 * d : d;
}
static double routeCost(const vector<int>& r) {
    double ans = 0;
    for (int t = 1; t <= n; ++t) ans += edgeCost(r, t);
    return ans;
}

// Sum just the edges whose endpoints or penalty source can change after swapping positions a,b.
static double nearbyCost(const vector<int>& r, int a, int b) {
    int ts[4] = {a, a + 1, b, b + 1};
    double ans = 0;
    for (int z = 0; z < 4; ++z) {
        int t = ts[z];
        if (t < 1 || t > n) continue;
        bool seen = false;
        for (int q = 0; q < z; ++q) if (ts[q] == t) seen = true;
        if (!seen) ans += edgeCost(r, t);
    }
    return ans;
}

static void improve(vector<int>& r) {
    // Put a nearby prime at a penalized source only when its complete local effect is beneficial.
    for (int t = 10; t <= n; t += 10) {
        int k = t - 1;
        if (primeCity[r[k]]) continue;
        int best = -1;
        double bestDelta = 0.0;
        int lo = max(1, k - 30), hi = min(n - 1, k + 30);
        for (int j = lo; j <= hi; ++j) if (primeCity[r[j]]) {
            double before = nearbyCost(r, k, j);
            swap(r[k], r[j]);
            double after = nearbyCost(r, k, j);
            swap(r[k], r[j]);
            if (after - before < bestDelta) bestDelta = after - before, best = j;
        }
        if (best != -1) swap(r[k], r[best]);
    }
    // Remove short Hilbert discontinuities and retain only strictly improving moves.
    for (int pass = 0; pass < 2; ++pass) {
        for (int i = 1; i + 1 < n; ++i) {
            double before = nearbyCost(r, i, i + 1);
            swap(r[i], r[i + 1]);
            double after = nearbyCost(r, i, i + 1);
            if (after >= before) swap(r[i], r[i + 1]);
        }
    }
}

// An adjacent exchange cannot remove a crossing that needs several consecutive
// cities reversed.  Unlike ordinary TSP 2-opt, carrot multipliers make every
// internal edge relevant, so evaluate each short reversal exactly.
static void boundedTwoOpt(vector<int>& r) {
    for (int len = 3; len <= 8; ++len) {
        for (int i = 1; i + len - 1 < n; ++i) {
            int j = i + len - 1;
            double before = 0.0;
            for (int t = i; t <= j + 1; ++t) before += edgeCost(r, t);
            reverse(r.begin() + i, r.begin() + j + 1);
            double after = 0.0;
            for (int t = i; t <= j + 1; ++t) after += edgeCost(r, t);
            if (after >= before) reverse(r.begin() + i, r.begin() + j + 1);
        }
    }
}

static uint64_t hilbert(uint32_t x, uint32_t y) {
    const uint32_t M = (1u << 21) - 1;
    uint64_t d = 0;
    for (uint32_t s = 1u << 20; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) x = M - x, y = M - y;
            swap(x, y);
        }
    }
    return d;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    if (!(cin >> n)) return 0;
    pt.resize(n);
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto& p : pt) {
        cin >> p.x >> p.y;
        minx = min(minx, p.x); maxx = max(maxx, p.x);
        miny = min(miny, p.y); maxy = max(maxy, p.y);
    }
    primeCity.assign(n, true);
    primeCity[0] = false;
    if (n > 1) primeCity[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (primeCity[i])
        for (int j = i * i; j < n; j += i) primeCity[j] = false;

    // Compare every construction only after the same exact repair.  Otherwise
    // a raw Hilbert/monotone loser can be discarded despite becoming best once
    // its penalty positions and short discontinuities are repaired.
    vector<int> monotone(n + 1);
    monotone[0] = monotone[n] = 0;
    for (int i = 1; i < n; ++i) monotone[i] = i;
    vector<int> best = monotone;
    improve(best);
    double bestCost = routeCost(best);
    vector<int> rev = monotone;
    reverse(rev.begin() + 1, rev.begin() + n);
    improve(rev);
    double c = routeCost(rev);
    if (c < bestCost) bestCost = c, best = rev;

    const uint32_t M = (1u << 21) - 1;
    auto norm = [&](long long v, long long lo, long long hi) -> uint32_t {
        if (hi == lo) return M / 2;
        return (uint32_t)((__int128)(v - lo) * M / (hi - lo));
    };
    vector<uint32_t> xx(n), yy(n);
    for (int i = 0; i < n; ++i) xx[i] = norm(pt[i].x, minx, maxx), yy[i] = norm(pt[i].y, miny, maxy);

    for (int mode = 0; mode < 8; ++mode) {
        vector<pair<uint64_t,int>> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i) {
            uint32_t x = xx[i], y = yy[i], a, b;
            switch (mode) {
                case 0: a=x; b=y; break; case 1: a=y; b=x; break;
                case 2: a=M-x; b=y; break; case 3: a=x; b=M-y; break;
                case 4: a=M-x; b=M-y; break; case 5: a=M-y; b=x; break;
                case 6: a=y; b=M-x; break; default: a=M-y; b=M-x; break;
            }
            order.push_back({hilbert(a,b), i});
        }
        sort(order.begin(), order.end());
        int at = 0;
        while (order[at].second != 0) ++at;
        for (int direction : {-1, 1}) {
            vector<int> r(n + 1); r[0] = r[n] = 0;
            for (int k = 1; k < n; ++k) r[k] = order[(at + direction * k + n) % n].second;
            improve(r);
            double value = routeCost(r);
            if (value < bestCost) bestCost = value, best = move(r);
        }
    }
    // Discriminating candidate: this is selected only by the same exact
    // objective, so a poor greedy path cannot displace the incumbent.
    vector<int> greedy = greedyNearestTour();
    // A raw greedy route can lose to Hilbert even when the same exact local
    // repair makes it the best candidate, so repair before the comparison.
    improve(greedy);
    c = routeCost(greedy);
    if (c < bestCost) bestCost = c, best = move(greedy);

    improve(best);
    // Final-only exact short 2-opt is monotone: it cannot displace a portfolio
    // winner unless it improves that same checker objective.
    boundedTwoOpt(best);
    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
}
