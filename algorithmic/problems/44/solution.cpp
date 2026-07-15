#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// A deletion-capable kd tree.  Splitting prime and non-prime cities lets the
// constructor ask for a nearby city of the required type without scanning.
struct KDTree {
    struct Node {
        long long lx, rx, ly, ry;
        int left = -1, right = -1, parent = -1, id = -1, live = 0;
    };
    const vector<Point>& p;
    vector<Node> tr;
    vector<int> leaf;

    KDTree(const vector<Point>& points, vector<int> ids) : p(points), leaf(points.size(), -1) {
        tr.reserve(ids.size() * 2 + 1);
        if (!ids.empty()) build(ids, 0, (int)ids.size(), -1);
    }

    int build(vector<int>& a, int l, int r, int par) {
        long long lx = p[a[l]].x, rx = lx, ly = p[a[l]].y, ry = ly;
        for (int i = l + 1; i < r; ++i) {
            lx = min(lx, p[a[i]].x); rx = max(rx, p[a[i]].x);
            ly = min(ly, p[a[i]].y); ry = max(ry, p[a[i]].y);
        }
        int v = (int)tr.size();
        tr.push_back({lx, rx, ly, ry, -1, -1, par, -1, r-l});
        if (r - l == 1) {
            tr[v].id = a[l];
            leaf[a[l]] = v;
            return v;
        }
        int mid = (l + r) / 2;
        bool byX = (rx - lx >= ry - ly);
        nth_element(a.begin()+l, a.begin()+mid, a.begin()+r, [&](int u, int w) {
            return byX ? p[u].x < p[w].x : p[u].y < p[w].y;
        });
        tr[v].left = build(a, l, mid, v);
        tr[v].right = build(a, mid, r, v);
        return v;
    }

    static long double sq(long long x) { return (long double)x * x; }
    long double pointDist(int a, int b) const {
        return sq(p[a].x-p[b].x) + sq(p[a].y-p[b].y);
    }
    long double boxDist(int v, int q) const {
        long double dx = 0, dy = 0;
        if (p[q].x < tr[v].lx) dx = (long double)tr[v].lx - p[q].x;
        else if (p[q].x > tr[v].rx) dx = (long double)p[q].x - tr[v].rx;
        if (p[q].y < tr[v].ly) dy = (long double)tr[v].ly - p[q].y;
        else if (p[q].y > tr[v].ry) dy = (long double)p[q].y - tr[v].ry;
        return dx*dx + dy*dy;
    }
    void search(int v, int q, int& ans, long double& best) const {
        if (v < 0 || tr[v].live == 0 || boxDist(v, q) >= best) return;
        if (tr[v].id >= 0) {
            long double d = pointDist(q, tr[v].id);
            if (d < best) best = d, ans = tr[v].id;
            return;
        }
        int a = tr[v].left, b = tr[v].right;
        if (a >= 0 && b >= 0 && boxDist(b, q) < boxDist(a, q)) swap(a, b);
        search(a, q, ans, best);
        search(b, q, ans, best);
    }
    int nearest(int q) const {
        int ans = -1;
        long double best = numeric_limits<long double>::infinity();
        if (!tr.empty()) search(0, q, ans, best);
        return ans;
    }
    void erase(int id) {
        int v = leaf[id];
        if (v < 0 || tr[v].live == 0) return;
        for (; v >= 0; v = tr[v].parent) --tr[v].live;
    }
};

static vector<char> sieve(int n) {
    vector<char> prime(n, true);
    if (n) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL*i*i < n; ++i)
        if (prime[i]) for (int j = i*i; j < n; j += i) prime[j] = false;
    return prime;
}

static long double routeCost(const vector<int>& r, const vector<Point>& p, const vector<char>& prime) {
    long double total = 0;
    int n = (int)p.size();
    for (int t = 1; t <= n; ++t) {
        int a = r[t-1], b = r[t];
        long double dx = (long double)p[a].x-p[b].x, dy = (long double)p[a].y-p[b].y;
        long double d = sqrtl(dx*dx + dy*dy);
        if (t % 10 == 0 && !prime[a]) d *= 1.1L;
        total += d;
    }
    return total;
}

// A bounded 2-opt pass is a cheap, general repair for the crossings left by a
// greedy construction.  Unlike ordinary 2-opt, every affected edge is priced
// at its actual global step, since reversal changes penalty sources as well.
static void shortTwoOpt(vector<int>& r, const vector<Point>& p, const vector<char>& prime) {
    const int n = (int)p.size();
    auto edge = [&](int a, int b, int t) -> long double {
        long double dx = (long double)p[a].x - p[b].x;
        long double dy = (long double)p[a].y - p[b].y;
        long double d = sqrtl(dx * dx + dy * dy);
        if (t % 10 == 0 && !prime[a]) d *= 1.1L;
        return d;
    };
    // At most ten vertices are reversed.  This keeps the pass linear while
    // still covering the short crossings characteristic of NN tours.
    for (int i = 1; i < n - 1; ++i) {
        int bestJ = -1;
        long double bestDelta = 0;
        for (int j = i + 1; j < n && j <= i + 9; ++j) {
            long double before = 0, after = 0;
            for (int t = i; t <= j + 1; ++t) {
                before += edge(r[t - 1], r[t], t);
                int u = (t - 1 >= i && t - 1 <= j) ? r[i + j - (t - 1)] : r[t - 1];
                int v = (t >= i && t <= j) ? r[i + j - t] : r[t];
                after += edge(u, v, t);
            }
            long double delta = after - before;
            if (delta < bestDelta) bestDelta = delta, bestJ = j;
        }
        if (bestJ >= 0) reverse(r.begin() + i, r.begin() + bestJ + 1);
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> p(n);
    for (auto& q : p) cin >> q.x >> q.y;
    vector<char> prime = sieve(n);
    vector<int> ps, np;
    int primeLeft = 0;
    for (int i = 1; i < n; ++i) {
        if (prime[i]) ps.push_back(i), ++primeLeft;
        else np.push_back(i);
    }
    // Compare the old hard constraint against its key ablation.  A prime source
    // avoids a tenth-edge surcharge, but forcing it can cost much more than 10%.
    auto buildNearest = [&](bool reservePrimes) {
        KDTree primeTree(p, ps), nonprimeTree(p, np);
        int primesRemaining = primeLeft, specialLeft = 0;
        for (int pos = 1; pos < n; ++pos) if ((pos + 1) % 10 == 0) ++specialLeft;
        vector<int> route;
        route.reserve(n + 1);
        route.push_back(0);
        int cur = 0;
        for (int pos = 1; pos < n; ++pos) {
            bool special = ((pos + 1) % 10 == 0);
            int a = primeTree.nearest(cur), b = nonprimeTree.nearest(cur), take = -1;
            if (reservePrimes) {
                if (special) take = a;
                else if (primesRemaining <= specialLeft) take = b;
                else if (a < 0) take = b;
                else if (b < 0) take = a;
                else take = primeTree.pointDist(cur, a) < nonprimeTree.pointDist(cur, b) ? a : b;
            } else {
                // This is the contract ablation: price the actual 1.1 factor
                // instead of treating a prime at a tenth position as mandatory.
                if (a < 0) take = b;
                else if (b < 0) take = a;
                else if (special)
                    take = primeTree.pointDist(cur, a) <= 1.21L * nonprimeTree.pointDist(cur, b) ? a : b;
                else take = primeTree.pointDist(cur, a) < nonprimeTree.pointDist(cur, b) ? a : b;
            }
            if (take < 0) take = (a >= 0 ? a : b);
            route.push_back(take);
            if (prime[take]) {
                primeTree.erase(take);
                --primesRemaining;
            } else nonprimeTree.erase(take);
            if (special) --specialLeft;
            cur = take;
        }
        route.push_back(0);
        return route;
    };
    vector<int> strictNN = buildNearest(true);
    vector<int> relaxedNN = buildNearest(false);

    // Runtime evaluation makes this a safe experiment: retain the incumbent
    // and monotone candidates whenever the relaxed policy loses.
    vector<int> monotone;
    monotone.reserve(n + 1);
    for (int i = 0; i < n; ++i) monotone.push_back(i);
    monotone.push_back(0);
    const vector<int>* answer = &strictNN;
    if (routeCost(relaxedNN, p, prime) < routeCost(*answer, p, prime)) answer = &relaxedNN;
    if (routeCost(monotone, p, prime) < routeCost(*answer, p, prime)) answer = &monotone;

    // Keep an untouched incumbent candidate: the local pass is accepted only
    // after an independent full-route evaluation.
    vector<int> repaired = *answer;
    shortTwoOpt(repaired, p, prime);
    if (routeCost(repaired, p, prime) < routeCost(*answer, p, prime)) answer = &repaired;
    cout << n + 1 << '\n';
    for (int v : *answer) cout << v << '\n';
}
