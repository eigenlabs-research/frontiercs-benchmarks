#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

// Hilbert index on the fixed square containing all legal coordinates.
static unsigned long long hilbert(unsigned x, unsigned y) {
    unsigned long long ans = 0;
    for (unsigned s = 1u << 30; s; s >>= 1) {
        unsigned rx = (x & s) != 0, ry = (y & s) != 0;
        ans += (unsigned long long)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) {
                x = (1u << 31) - 1 - x;
                y = (1u << 31) - 1 - y;
            }
            swap(x, y);
        }
    }
    return ans;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> a(n);
    for (auto &p : a) cin >> p.x >> p.y;

    vector<char> prime(n, true);
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i)
        if (prime[i]) for (int j = i * i; j < n; j += i) prime[j] = false;

    auto d = [&](int u, int v) -> double {
        return hypot((double)a[u].x - a[v].x, (double)a[u].y - a[v].y);
    };
    auto edge = [&](int step, int u, int v) -> double {
        return (step % 10 == 0 && !prime[u] ? 1.1 : 1.0) * d(u, v);
    };
    auto cost = [&](const vector<int>& p) {
        double z = 0;
        for (int t = 1; t <= n; ++t) z += edge(t, p[t - 1], p[t]);
        return z;
    };

    vector<pair<unsigned long long,int>> keyed;
    keyed.reserve(n);
    for (int i = 0; i < n; ++i) {
        unsigned x = (unsigned)(a[i].x + 1000000000LL);
        unsigned y = (unsigned)(a[i].y + 1000000000LL);
        keyed.push_back({hilbert(x, y), i});
    }
    sort(keyed.begin(), keyed.end());

    vector<int> cyclic;
    cyclic.reserve(n);
    for (auto [key, id] : keyed) cyclic.push_back(id);
    int at = find(cyclic.begin(), cyclic.end(), 0) - cyclic.begin();
    vector<int> forward(n + 1), backward(n + 1), base(n + 1);
    forward[0] = backward[0] = base[0] = 0;
    for (int k = 1; k < n; ++k) {
        forward[k] = cyclic[(at + k) % n];
        backward[k] = cyclic[(at - k + n) % n];
        base[k] = k;
    }
    forward[n] = backward[n] = base[n] = 0;

    vector<int> route = forward;
    double best = cost(forward), cb = cost(backward), cx = cost(base);
    if (cb < best) route = backward, best = cb;
    if (cx < best) route = base, best = cx;

    // A second, deliberately different representation is a strip tour.  The
    // input is x-ordered, so equal-population x strips need no global sort;
    // within each strip we sweep y, reversing direction in successive strips.
    // This is especially useful when consecutive x values alternate between
    // distant horizontal bands, a geometry where a space-filling curve may
    // still make many unnecessary band crossings.
    auto considerStrips = [&](int strips) {
        // The two orientations differ only in each already-sorted stripe's
        // direction.  Building both from one sort avoids duplicate O(N log N)
        // work for every resolution.
        vector<int> up(n + 1), down(n + 1);
        up[0] = up[n] = down[0] = down[n] = 0;
        int m = n - 1;
        for (int s = 0; s < strips; ++s) {
            int l = 1 + (long long)m * s / strips;
            int r = 1 + (long long)m * (s + 1) / strips;
            vector<int> part;
            part.reserve(r - l);
            for (int id = l; id < r; ++id) part.push_back(id);
            sort(part.begin(), part.end(), [&](int u, int v) {
                if (a[u].y != a[v].y) return a[u].y < a[v].y;
                return u < v;
            });
            for (int k = 0; k < (int)part.size(); ++k) {
                up[l + k] = part[(s & 1) ? (int)part.size() - 1 - k : k];
                down[l + k] = part[(s & 1) ? k : (int)part.size() - 1 - k];
            }
        }
        double zu = cost(up), zd = cost(down);
        if (zu < best) route = move(up), best = zu;
        if (zd < best) route = move(down), best = zd;
    };
    // A logarithmic family covers different aspect ratios without turning the
    // construction into a sample-specific fixed grid.
    for (int b = 2; b < n; b <<= 1) {
        considerStrips(b);
        if (b >= 512) break;
    }

    // The transposed construction is deliberately separate from the x-strip
    // family: it keeps horizontal bands intact and sweeps along x inside each
    // band.  Equal-population bands avoid relying on any coordinate scale.
    vector<int> byY;
    byY.reserve(n - 1);
    for (int id = 1; id < n; ++id) byY.push_back(id);
    sort(byY.begin(), byY.end(), [&](int u, int v) {
        if (a[u].y != a[v].y) return a[u].y < a[v].y;
        return u < v;
    });
    auto considerHorizontalStrips = [&](int strips) {
        vector<int> right(n + 1), left(n + 1);
        right[0] = right[n] = left[0] = left[n] = 0;
        int m = n - 1;
        for (int s = 0; s < strips; ++s) {
            int l = (long long)m * s / strips;
            int r = (long long)m * (s + 1) / strips;
            vector<int> part(byY.begin() + l, byY.begin() + r);
            // IDs are in strictly increasing x order by the input contract.
            sort(part.begin(), part.end());
            for (int k = 0; k < (int)part.size(); ++k) {
                int pos = 1 + l + k;
                right[pos] = part[(s & 1) ? (int)part.size() - 1 - k : k];
                left[pos] = part[(s & 1) ? k : (int)part.size() - 1 - k];
            }
        }
        double zr = cost(right), zl = cost(left);
        if (zr < best) route = move(right), best = zr;
        if (zl < best) route = move(left), best = zl;
    };
    for (int b = 2; b < n; b <<= 1) {
        considerHorizontalStrips(b);
        // A few resolutions cover banded data without duplicating the full
        // logarithmic x-strip sorting budget.
        if (b >= 32) break;
    }

    // A short 2-opt pass repairs local discontinuities of the curve.  Its
    // delta is exact: the two boundary edges are evaluated in full below;
    // reversed internal edges retain their ordinary distances, so only their
    // every-tenth-step surcharges need reconsidering.
    const int span = 12;
    // A reversal can make an earlier start position profitable again.  Repeat
    // the same bounded, exact neighborhood once to test that missed-backtrack
    // effect, while retaining a strict fixed upper bound for large instances.
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int i = 1; i < n; ++i) {
            for (int len = 2; len <= span && i + len - 1 < n; ++len) {
                int j = i + len - 1;
                double oldc = edge(i, route[i-1], route[i]) + edge(j+1, route[j], route[j+1]);
                double newc = edge(i, route[i-1], route[j]) + edge(j+1, route[i], route[j+1]);
                // Internal base distances are merely reversed.  At a penalty step,
                // the source city and hence its 10% surcharge may change.
                // Step i is the left boundary, not an internal reversed edge.
                // Starting at i+1 is essential when i itself is divisible by 10.
                int first = ((i + 10) / 10) * 10;
                for (int t = first; t <= j; t += 10) {
                    oldc += edge(t, route[t-1], route[t]) - d(route[t-1], route[t]);
                    int ns = i + j - t + 1, nt = i + j - t;
                    newc += edge(t, route[ns], route[nt]) - d(route[ns], route[nt]);
                }
                if (newc + 1e-9 < oldc) {
                    reverse(route.begin() + i, route.begin() + j + 1);
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }

    cout << n + 1 << '\n';
    for (int v : route) cout << v << '\n';
    return 0;
}
