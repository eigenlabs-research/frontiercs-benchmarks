#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

static unsigned long long hilbert(unsigned x, unsigned y) {
    // A 21-bit Hilbert index.  Unlike sorting just by x, this keeps both
    // coordinates local on inputs made of several spatial clusters.
    unsigned long long d = 0;
    for (unsigned s = 1u << 20; s; s >>= 1) {
        unsigned rx = (x & s) != 0, ry = (y & s) != 0;
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
    vector<Point> p(n);
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &q : p) {
        cin >> q.x >> q.y;
        minx = min(minx, q.x); maxx = max(maxx, q.x);
        miny = min(miny, q.y); maxy = max(maxy, q.y);
    }
    vector<char> prime(n, true);
    prime[0] = false; if (n > 1) prime[1] = false;
    for (int i = 2; 1LL*i*i < n; ++i) if (prime[i])
        for (int j = i*i; j < n; j += i) prime[j] = false;

    auto edge = [&](int a, int b, int t) {
        double dx = (double)p[a].x - p[b].x, dy = (double)p[a].y - p[b].y;
        double v = sqrt(dx*dx + dy*dy);
        return (t % 10 == 0 && !prime[a]) ? 1.1*v : v;
    };
    auto cost = [&](const vector<int>& r) {
        double z = 0;
        for (int t = 1; t <= n; ++t) z += edge(r[t-1], r[t], t);
        return z;
    };
    // Turn any cyclic ordering into the required ordering beginning at 0.
    auto fromCycle = [&](const vector<int>& a, bool backwards) {
        int at = find(a.begin(), a.end(), 0) - a.begin();
        vector<int> r; r.reserve(n+1); r.push_back(0);
        for (int k = 1; k < n; ++k) {
            int z = backwards ? (at - k + n) % n : (at + k) % n;
            r.push_back(a[z]);
        }
        r.push_back(0); return r;
    };

    vector<vector<int>> candidates;
    vector<int> byx(n), byy(n), byh(n);
    iota(byx.begin(), byx.end(), 0);
    iota(byy.begin(), byy.end(), 0);
    iota(byh.begin(), byh.end(), 0);
    sort(byy.begin(), byy.end(), [&](int a, int b) {
        if (p[a].y != p[b].y) return p[a].y < p[b].y;
        return p[a].x < p[b].x;
    });
    long double rx = maxx - minx, ry = maxy - miny;
    vector<unsigned long long> key(n);
    for (int i = 0; i < n; ++i) {
        unsigned X = rx == 0 ? 0 : (unsigned)((p[i].x-minx) * ((1LL<<21)-1) / rx);
        unsigned Y = ry == 0 ? 0 : (unsigned)((p[i].y-miny) * ((1LL<<21)-1) / ry);
        key[i] = hilbert(X, Y);
    }
    sort(byh.begin(), byh.end(), [&](int a, int b) { return key[a] < key[b]; });
    for (auto *q : {&byx, &byy, &byh}) {
        candidates.push_back(fromCycle(*q, false));
        candidates.push_back(fromCycle(*q, true));
    }
    int best = 0;
    double bestCost = cost(candidates[0]);
    for (int i = 1; i < (int)candidates.size(); ++i) {
        double z = cost(candidates[i]);
        if (z < bestCost) bestCost = z, best = i;
    }
    vector<int> r = move(candidates[best]);

    // Short 2-opt moves repair crossings and, because the exact weighted
    // objective is evaluated, also place prime IDs favorably at carrot steps.
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int l = 1; l < n; ++l) {
            for (int len = 2; len <= 7 && l + len - 1 < n; ++len) {
                int rr = l + len - 1;
                double oldv = 0, newv = 0;
                for (int t = l; t <= rr + 1; ++t) {
                    oldv += edge(r[t-1], r[t], t);
                    int aidx = t-1, bidx = t;
                    int a = (aidx >= l && aidx <= rr) ? r[l+rr-aidx] : r[aidx];
                    int b = (bidx >= l && bidx <= rr) ? r[l+rr-bidx] : r[bidx];
                    newv += edge(a, b, t);
                }
                if (newv + 1e-7 < oldv) {
                    reverse(r.begin()+l, r.begin()+rr+1);
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
    cout << n+1 << '\n';
    for (int v : r) cout << v << '\n';
}
