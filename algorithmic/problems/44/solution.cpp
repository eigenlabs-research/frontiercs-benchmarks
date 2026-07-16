#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert index in a 2^31 by 2^31 square.
// The Hilbert rotation temporarily creates negative coordinates.  Keep its
// working state signed: unsigned reflection would wrap and corrupt later bits.
static unsigned long long hilbert(long long x, long long y) {
    unsigned long long d = 0;
    for (long long s = 1LL << 30; s; s >>= 1) {
        unsigned int rx = (x & s) != 0, ry = (y & s) != 0;
        d += (unsigned long long)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = s - 1 - x; y = s - 1 - y; }
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
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    const unsigned int LIM = 0x7fffffffu;
    vector<unsigned int> xx(n), yy(n);
    long double sx = (long double)maxx - minx, sy = (long double)maxy - miny;
    for (int i = 0; i < n; ++i) {
        xx[i] = sx == 0 ? 0 : (unsigned int)((long double)(a[i].x - minx) * LIM / sx);
        yy[i] = sy == 0 ? 0 : (unsigned int)((long double)(a[i].y - miny) * LIM / sy);
    }
    auto edge = [&](const vector<int>& r, int t) {
        int u = r[t - 1], v = r[t];
        double z = hypot((double)a[u].x - a[v].x, (double)a[u].y - a[v].y);
        if (t % 10 == 0 && !prime[u]) z *= 1.1;
        return z;
    };
    auto cost = [&](const vector<int>& r) {
        double s = 0;
        for (int t = 1; t <= n; ++t) s += edge(r, t);
        return s;
    };

    // The incumbent curve is retained as a safe competing representation.
    vector<int> best;
    double bestCost = numeric_limits<double>::infinity();
    // Reflections and axis exchanges give the curve different entry and exit regions.
    for (int mode = 0; mode < 8; ++mode) {
        vector<pair<unsigned long long, int>> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i) {
            unsigned int x = xx[i], y = yy[i];
            if (mode & 1) x = LIM - x;
            if (mode & 2) y = LIM - y;
            if (mode & 4) swap(x, y);
            order.push_back({hilbert(x, y), i});
        }
        sort(order.begin(), order.end());
        int at = 0;
        while (order[at].second != 0) ++at;
        vector<int> r;
        r.reserve(n + 1);
        r.push_back(0);
        for (int k = 1; k < n; ++k) r.push_back(order[(at + k) % n].second);
        r.push_back(0);
        double c = cost(r);
        if (c < bestCost) { bestCost = c; best.swap(r); }
    }

    // Exact small local moves are shared by both global representations.
    auto improve = [&](vector<int>& route) {
        for (int l = 1; l < n - 1; ++l) {
            int take = -1;
            double gain = 0;
            for (int r = l + 1; r <= min(n - 1, l + 8); ++r) {
                double before = 0, after = 0;
                for (int t = l; t <= r + 1; ++t) before += edge(route, t);
                reverse(route.begin() + l, route.begin() + r + 1);
                for (int t = l; t <= r + 1; ++t) after += edge(route, t);
                reverse(route.begin() + l, route.begin() + r + 1);
                if (after - before < gain) gain = after - before, take = r;
            }
            if (take != -1) reverse(route.begin() + l, route.begin() + take + 1);
        }
        // A reversal cannot express moving one out-of-place point between two
        // nearby neighbours.  Test such Or-opt(1) moves with the full global
        // carrot objective, since every shifted position can change a multiplier.
        for (int p = 1; p < n; ++p) {
            int take = -1;
            double gain = 0;
            for (int q = max(1, p - 4); q <= min(n - 1, p + 4); ++q) {
                if (q == p) continue;
                int lo = min(p, q), hi = max(p, q);
                double before = 0, after = 0;
                for (int t = lo; t <= hi + 1; ++t) before += edge(route, t);
                if (q < p)
                    rotate(route.begin() + q, route.begin() + p, route.begin() + p + 1);
                else
                    rotate(route.begin() + p, route.begin() + p + 1, route.begin() + q + 1);
                for (int t = lo; t <= hi + 1; ++t) after += edge(route, t);
                // Undo by moving the same city from its new position back to p.
                if (q < p)
                    rotate(route.begin() + q, route.begin() + q + 1, route.begin() + p + 1);
                else
                    rotate(route.begin() + p, route.begin() + q, route.begin() + q + 1);
                if (after - before < gain) gain = after - before, take = q;
            }
            if (take != -1) {
                if (take < p)
                    rotate(route.begin() + take, route.begin() + p, route.begin() + p + 1);
                else
                    rotate(route.begin() + p, route.begin() + p + 1, route.begin() + take + 1);
            }
        }
        auto swapDelta = [&](int p, int q) {
            vector<int> ts = {p, p + 1, q, q + 1};
            sort(ts.begin(), ts.end()); ts.erase(unique(ts.begin(), ts.end()), ts.end());
            double before = 0;
            for (int t : ts) if (t >= 1 && t <= n) before += edge(route, t);
            swap(route[p], route[q]);
            double after = 0;
            for (int t : ts) if (t >= 1 && t <= n) after += edge(route, t);
            swap(route[p], route[q]);
            return after - before;
        };
        for (int p = 9; p < n; p += 10) if (!prime[route[p]]) {
            int choose = -1; double gain = 0;
            for (int q = max(1, p - 40); q <= min(n - 1, p + 40); ++q) {
                if (q % 10 == 9 || !prime[route[q]]) continue;
                double d = swapDelta(p, q);
                if (d < gain) gain = d, choose = q;
            }
            if (choose != -1) swap(route[p], route[choose]);
        }
    };
    improve(best);
    bestCost = cost(best);

    // Diversified construction: a balanced kd tree supports deletion-aware exact
    // nearest-neighbour queries. Unlike a curve order, every next city is chosen
    // from the remaining spatial set, which is especially useful for separated clusters.
    struct Node { int city, left = -1, right = -1, parent = -1, cnt = 1; long long lx, rx, ly, ry; };
    vector<Node> tree; tree.reserve(n);
    vector<int> ids(n), nodeOf(n);
    iota(ids.begin(), ids.end(), 0);
    function<int(int,int,int,int)> build = [&](int lo, int hi, int depth, int parent) -> int {
        if (lo >= hi) return -1;
        int axis = depth & 1, mid = (lo + hi) / 2;
        nth_element(ids.begin()+lo, ids.begin()+mid, ids.begin()+hi, [&](int u, int v) {
            return axis ? a[u].y < a[v].y : a[u].x < a[v].x;
        });
        int z = (int)tree.size(), c = ids[mid];
        tree.push_back({c, -1, -1, parent, 1, a[c].x, a[c].x, a[c].y, a[c].y});
        nodeOf[c] = z;
        int L = build(lo, mid, depth+1, z), R = build(mid+1, hi, depth+1, z);
        tree[z].left = L; tree[z].right = R;
        for (int q : {L, R}) if (q != -1) {
            tree[z].cnt += tree[q].cnt;
            tree[z].lx = min(tree[z].lx, tree[q].lx); tree[z].rx = max(tree[z].rx, tree[q].rx);
            tree[z].ly = min(tree[z].ly, tree[q].ly); tree[z].ry = max(tree[z].ry, tree[q].ry);
        }
        return z;
    };
    int root = build(0, n, 0, -1);
    vector<char> alive(n, true);
    auto eraseCity = [&](int c) {
        alive[c] = false;
        for (int z = nodeOf[c]; z != -1; z = tree[z].parent) {
            Node &v = tree[z]; --v.cnt;
            // Keep bounds of only live descendants; this avoids late-stage scans
            // through boxes whose geometrically close points were already removed.
            v.lx = LLONG_MAX; v.rx = LLONG_MIN; v.ly = LLONG_MAX; v.ry = LLONG_MIN;
            auto absorb = [&](long long lx, long long rx, long long ly, long long ry) {
                v.lx=min(v.lx,lx); v.rx=max(v.rx,rx); v.ly=min(v.ly,ly); v.ry=max(v.ry,ry);
            };
            if (alive[v.city]) absorb(a[v.city].x,a[v.city].x,a[v.city].y,a[v.city].y);
            for (int q : {v.left, v.right}) if (q != -1 && tree[q].cnt)
                absorb(tree[q].lx,tree[q].rx,tree[q].ly,tree[q].ry);
        }
    };
    auto boxDist = [&](int z, int c) {
        long double dx = 0, dy = 0;
        if (a[c].x < tree[z].lx) dx = (long double)tree[z].lx - a[c].x;
        if (a[c].x > tree[z].rx) dx = (long double)a[c].x - tree[z].rx;
        if (a[c].y < tree[z].ly) dy = (long double)tree[z].ly - a[c].y;
        if (a[c].y > tree[z].ry) dy = (long double)a[c].y - tree[z].ry;
        return dx*dx + dy*dy;
    };
    vector<int> nn; nn.reserve(n+1); nn.push_back(0); eraseCity(0);
    for (int step = 1; step < n; ++step) {
        int from = nn.back(), chosen = -1;
        long double bd = numeric_limits<long double>::infinity();
        function<void(int)> query = [&](int z) {
            if (z == -1 || tree[z].cnt == 0 || boxDist(z, from) >= bd) return;
            int c = tree[z].city;
            if (alive[c]) { long double dx = (long double)a[c].x-a[from].x, dy = (long double)a[c].y-a[from].y, d=dx*dx+dy*dy; if (d < bd) bd=d, chosen=c; }
            int L=tree[z].left, R=tree[z].right;
            if (L != -1 && R != -1 && boxDist(R,from) < boxDist(L,from)) swap(L,R);
            query(L); query(R);
        };
        query(root); eraseCity(chosen); nn.push_back(chosen);
    }
    nn.push_back(0);
    improve(nn);
    double nnCost = cost(nn);
    if (nnCost < bestCost) best.swap(nn);

    cout << n + 1 << '\n';
    for (int id : best) cout << id << '\n';
    return 0;
}
