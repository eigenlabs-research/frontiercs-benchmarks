#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

static unsigned long long hilbert(unsigned x, unsigned y, unsigned n) {
    unsigned long long d = 0;
    for (unsigned s = n >> 1; s; s >>= 1) {
        unsigned rx = (x & s) != 0, ry = (y & s) != 0;
        d += (unsigned long long)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = n - 1 - x; y = n - 1 - y; }
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
    vector<Point> p(n);
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &a : p) {
        cin >> a.x >> a.y;
        minx = min(minx, a.x); maxx = max(maxx, a.x);
        miny = min(miny, a.y); maxy = max(maxy, a.y);
    }
    vector<char> prime(n, true);
    if (n) prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    auto dis = [&](int a, int b) -> double {
        return hypot((double)p[a].x - p[b].x, (double)p[a].y - p[b].y);
    };
    auto cost = [&](const vector<int>& q) {
        // q contains precisely the non-depot cities, in visiting order.
        double ans = 0;
        int prev = 0;
        for (int i = 0; i < (int)q.size(); ++i) {
            int t = i + 1;
            ans += ((t % 10 == 0 && !prime[prev]) ? 1.1 : 1.0) * dis(prev, q[i]);
            prev = q[i];
        }
        ans += ((n % 10 == 0 && !prime[prev]) ? 1.1 : 1.0) * dis(prev, 0);
        return ans;
    };

    vector<vector<int>> candidates;
    vector<int> byx(n - 1);
    iota(byx.begin(), byx.end(), 1);
    candidates.push_back(byx);
    reverse(byx.begin(), byx.end());
    candidates.push_back(byx);

    vector<int> byy(n - 1);
    iota(byy.begin(), byy.end(), 1);
    sort(byy.begin(), byy.end(), [&](int a, int b) {
        if (p[a].y != p[b].y) return p[a].y < p[b].y;
        return p[a].x < p[b].x;
    });
    candidates.push_back(byy);
    reverse(byy.begin(), byy.end());
    candidates.push_back(byy);

    const unsigned SIDE = 1u << 21;
    unsigned long long rx = maxx - minx, ry = maxy - miny;
    vector<pair<unsigned long long,int>> keys;
    keys.reserve(n - 1);
    for (int i = 1; i < n; ++i) {
        unsigned xx = rx ? (unsigned)(((unsigned long long)(p[i].x - minx) * (SIDE - 1)) / rx) : SIDE / 2;
        unsigned yy = ry ? (unsigned)(((unsigned long long)(p[i].y - miny) * (SIDE - 1)) / ry) : SIDE / 2;
        keys.push_back({hilbert(xx, yy, SIDE), i});
    }
    sort(keys.begin(), keys.end());
    vector<int> h;
    h.reserve(n - 1);
    for (auto z : keys) h.push_back(z.second);
    candidates.push_back(h);
    reverse(h.begin(), h.end());
    candidates.push_back(h);

    vector<int> best = candidates[0];
    double bestCost = cost(best);
    for (const auto &q : candidates) {
        double v = cost(q);
        if (v < bestCost) bestCost = v, best = q;
    }

    // Control branch: the incumbent's geometric short 2-opt cleanup.
    // It is retained so the experiment cannot discard a known useful route.
    int width = (n <= 5000 ? 32 : 10);
    int passes = (n <= 5000 ? 3 : 2);
    vector<int> geometric = best;
    for (int pass = 0; pass < passes; ++pass) {
        bool changed = false;
        for (int i = 0; i + 2 < n; ++i) {
            int a = (i == 0 ? 0 : geometric[i - 1]);
            int b = geometric[i];
            for (int j = i + 1; j < n - 1 && j <= i + width; ++j) {
                int c = geometric[j], d = (j + 1 == n - 1 ? 0 : geometric[j + 1]);
                if (dis(a,c) + dis(b,d) + 1e-9 < dis(a,b) + dis(c,d)) {
                    reverse(geometric.begin() + i, geometric.begin() + j + 1);
                    changed = true;
                    b = geometric[i];
                }
            }
        }
        if (!changed) break;
    }
    double v = cost(geometric);
    if (v < bestCost) bestCost = v, best = geometric;

    // Ablation branch: the same neighbourhood, but accept a reversal only after
    // evaluating every affected edge with its actual step multiplier.  Reversing
    // a segment preserves most distances but can change which city starts steps
    // 10,20,..., precisely what an unweighted 2-opt test cannot see.
    vector<int> carrot = best;
    auto routeCityAfterReverse = [&](const vector<int>& q, int i, int j, int pos) {
        if (pos == 0 || pos == n) return 0;
        int k = pos - 1;
        return (i <= k && k <= j) ? q[i + j - k] : q[k];
    };
    for (int pass = 0; pass < passes; ++pass) {
        bool changed = false;
        for (int i = 0; i + 2 < n; ++i) {
            for (int j = i + 1; j < n - 1 && j <= i + width; ++j) {
                double before = 0, after = 0;
                auto addEdge = [&](int t) {
                    int oa = (t == 1 ? 0 : carrot[t - 2]);
                    int ob = (t == n ? 0 : carrot[t - 1]);
                    double mult = (t % 10 == 0 && !prime[oa]) ? 1.1 : 1.0;
                    before += mult * dis(oa, ob);
                    int na = routeCityAfterReverse(carrot, i, j, t - 1);
                    int nb = routeCityAfterReverse(carrot, i, j, t);
                    mult = (t % 10 == 0 && !prime[na]) ? 1.1 : 1.0;
                    after += mult * dis(na, nb);
                };
                // Ordinary internal edges have equal length in either direction.
                // Recompute just the two boundary edges and internal penalty steps.
                addEdge(i + 1);
                addEdge(j + 2);
                int firstPenalty = ((i + 2 + 9) / 10) * 10;
                for (int t = firstPenalty; t <= j + 1; t += 10) addEdge(t);
                if (after + 1e-9 < before) {
                    reverse(carrot.begin() + i, carrot.begin() + j + 1);
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
    v = cost(carrot);
    if (v < bestCost) bestCost = v, best = move(carrot);

    cout << n + 1 << '\n' << 0 << '\n';
    for (int x : best) cout << x << '\n';
    cout << 0 << '\n';
}
