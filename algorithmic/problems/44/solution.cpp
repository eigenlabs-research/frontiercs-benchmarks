#include <bits/stdc++.h>
using namespace std;

struct City {
    long long x, y;
    uint64_t key;
    int id;
};

static uint64_t hilbert(uint32_t x, uint32_t y) {
    // Coordinates are in a 2^31 by 2^31 square.  A Hilbert order keeps
    // consecutive points spatially close while still being cheap to sort.
    uint64_t d = 0;
    for (uint32_t s = 1u << 30; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * uint64_t(s) * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) {
                x = s - 1 - x;
                y = s - 1 - y;
            }
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
    vector<City> a(n);
    for (int i = 0; i < n; ++i) {
        cin >> a[i].x >> a[i].y;
        a[i].id = i;
        // Input bounds make this conversion fit in 31 bits.
        a[i].key = hilbert(uint32_t(a[i].x + 1000000000LL),
                           uint32_t(a[i].y + 1000000000LL));
    }

    vector<char> prime(n, true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int p = 2; 1LL * p * p < n; ++p)
        if (prime[p]) for (int q = p * p; q < n; q += p) prime[q] = false;

    auto edgeCost = [&](const vector<int>& r, int i) {
        // Edge i is route step i+1.
        int u = r[i], v = r[i + 1];
        double dx = double(a[u].x) - double(a[v].x);
        double dy = double(a[u].y) - double(a[v].y);
        double z = hypot(dx, dy);
        if ((i + 1) % 10 == 0 && !prime[u]) z *= 1.1;
        return z;
    };
    auto routeCost = [&](const vector<int>& r) {
        double ans = 0;
        for (int i = 0; i < n; ++i) ans += edgeCost(r, i);
        return ans;
    };

    // In each ten-step block, try placing a prime at the city which starts
    // the penalized edge.  Only a constant-size neighborhood changes, so the
    // best such swap can be selected using the real weighted objective.
    auto repairCarrots = [&](vector<int>& r) {
        for (int pos = 9; pos < n; pos += 10) {
            int lo = max(1, pos - 9);
            int best = -1;
            double bestDelta = 0.0; // retaining the current order is allowed
            for (int q = lo; q <= pos; ++q) {
                if (!prime[r[q]]) continue;
                vector<int> es;
                for (int e : {q - 1, q, pos - 1, pos})
                    if (e >= 0 && e < n && find(es.begin(), es.end(), e) == es.end()) es.push_back(e);
                double before = 0, after = 0;
                for (int e : es) before += edgeCost(r, e);
                swap(r[q], r[pos]);
                for (int e : es) after += edgeCost(r, e);
                swap(r[q], r[pos]);
                if (after - before < bestDelta) {
                    bestDelta = after - before;
                    best = q;
                }
            }
            if (best != -1) swap(r[best], r[pos]);
        }
    };

    // Exact short-range 2-opt is a discriminating refinement: Hilbert gives
    // locality, but it does not guarantee that its short consecutive arcs do
    // not cross.  Unlike ordinary TSP 2-opt, reversing an arc also relocates
    // carrot sources, so every proposed move is evaluated with those weights.
    auto shortTwoOpt = [&](vector<int>& r) {
        vector<double> d(n);
        auto rebuild = [&](int lo, int hi) {
            lo = max(lo, 0); hi = min(hi, n - 1);
            for (int e = lo; e <= hi; ++e) {
                double dx = double(a[r[e]].x) - double(a[r[e + 1]].x);
                double dy = double(a[r[e]].y) - double(a[r[e + 1]].y);
                d[e] = hypot(dx, dy);
            }
        };
        rebuild(0, n - 1);
        auto weight = [&](int e, int u) {
            return ((e + 1) % 10 == 0 && !prime[u]) ? 1.1 : 1.0;
        };
        const int W = 12;
        // Two sweeps permit improvements exposed by an earlier reversal,
        // while retaining a linear-size bounded search for N=200000.
        for (int pass = 0; pass < 2; ++pass) {
            for (int i = 0; i + 2 < n; ++i) {
                int top = min(n - 1, i + W);
                for (int j = i + 2; j <= top; ++j) {
                    double before = d[i] * weight(i, r[i]) +
                                    d[j] * weight(j, r[j]);
                    double after = 0;
                    double dx = double(a[r[i]].x) - double(a[r[j]].x);
                    double dy = double(a[r[i]].y) - double(a[r[j]].y);
                    after += hypot(dx, dy) * weight(i, r[i]);
                    dx = double(a[r[i + 1]].x) - double(a[r[j + 1]].x);
                    dy = double(a[r[i + 1]].y) - double(a[r[j + 1]].y);
                    after += hypot(dx, dy) * weight(j, r[i + 1]);
                    for (int t = i + 1; t < j; ++t) {
                        before += d[t] * weight(t, r[t]);
                        // Old edge t is traversed backwards at this new index.
                        int newEdge = i + j - t;
                        after += d[t] * weight(newEdge, r[t + 1]);
                    }
                    if (after < before - 1e-7) {
                        reverse(r.begin() + i + 1, r.begin() + j + 1);
                        rebuild(i, j);
                    }
                }
            }
        }
    };

    vector<City> h = a;
    sort(h.begin(), h.end(), [](const City& u, const City& v) {
        if (u.key != v.key) return u.key < v.key;
        return u.id < v.id;
    });
    int zero = 0;
    while (h[zero].id != 0) ++zero;

    vector<int> answer;
    double answerCost = numeric_limits<double>::infinity();
    auto consider = [&](const vector<int>& cycle, int start, int dir) {
        vector<int> r(n + 1);
        r[0] = r[n] = 0;
        for (int p = 1; p < n; ++p) {
            int j = (start + dir * p) % n;
            if (j < 0) j += n;
            r[p] = cycle[j];
        }
        repairCarrots(r);
        shortTwoOpt(r);
        // A reversal can put a non-prime back at a carrot source; restore
        // the same exact-cost local invariant before comparing candidates.
        repairCarrots(r);
        double c = routeCost(r);
        if (c < answerCost) {
            answerCost = c;
            answer.swap(r);
        }
    };

    vector<int> hc(n), xc(n);
    for (int i = 0; i < n; ++i) hc[i] = h[i].id, xc[i] = i;
    // Both orientations are worthwhile because the fixed North Pole and the
    // step-number-dependent carrot cost make an otherwise identical cycle
    // asymmetric.  The monotone order is retained as a safe candidate.
    consider(hc, zero, 1);
    consider(hc, zero, -1);
    consider(xc, 0, 1);
    consider(xc, 0, -1);

    cout << n + 1 << '\n';
    for (int v : answer) cout << v << '\n';
}
