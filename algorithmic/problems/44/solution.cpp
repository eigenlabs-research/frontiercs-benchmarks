#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };

static uint64_t hilbert(uint32_t x, uint32_t y) {
    const uint32_t mask = 0x7fffffffu;
    uint64_t d = 0;
    for (uint32_t s = 1u << 30; s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += uint64_t(s) * uint64_t(s) * ((3 * rx) ^ ry);
        if (!ry) {
            if (rx) { x = mask - x; y = mask - y; }
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
    prime[0] = false;
    if (n > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < n; ++i) if (prime[i])
        for (int j = i * i; j < n; j += i) prime[j] = false;

    const uint32_t M = 0x7fffffffu;
    vector<uint32_t> X(n), Y(n);
    auto scale = [&](long long v, long long lo, long long hi) -> uint32_t {
        if (hi == lo) return M / 2;
        return (uint64_t)(v - lo) * M / (uint64_t)(hi - lo);
    };
    for (int i = 0; i < n; ++i) {
        X[i] = scale(p[i].x, minx, maxx);
        Y[i] = scale(p[i].y, miny, maxy);
    }
    auto edge = [&](int u, int v, int pos) {
        double z = hypot((double)p[u].x - p[v].x, (double)p[u].y - p[v].y);
        if ((pos + 1) % 10 == 0 && !prime[u]) z *= 1.1;
        return z;
    };
    auto score = [&](const vector<int>& a) {
        double ans = 0;
        for (int i = 0; i < n; ++i) ans += edge(a[i], a[(i + 1) % n], i);
        return ans;
    };
    auto makeRoute = [&](vector<pair<uint64_t,int>>& ord) {
        sort(ord.begin(), ord.end());
        vector<int> a; a.reserve(n); a.push_back(0);
        int cut = 0;
        while (ord[cut].second != 0) ++cut;
        for (int k = 1; k < n; ++k) a.push_back(ord[(cut + k) % n].second);
        return a;
    };
    auto tryCandidate = [&](vector<int> candidate, vector<int>& best, double &bestCost) {
        // A cycle can be cut at 0 in either direction.  Geometry is unchanged,
        // but the carrot sources are not, so score both exact directed tours.
        vector<int> reversed = candidate;
        reverse(reversed.begin() + 1, reversed.end());
        double c = score(candidate), rc = score(reversed);
        if (rc < c) { c = rc; candidate.swap(reversed); }
        if (c < bestCost) { bestCost = c; best.swap(candidate); }
    };

    vector<int> best;
    double bestCost = numeric_limits<double>::infinity();
    // The input-order cycle also has two distinct carrot assignments once cut at 0.
    vector<int> baseline(n);
    iota(baseline.begin(), baseline.end(), 0);
    tryCandidate(baseline, best, bestCost);
    // Space-filling-curve incumbent candidates.
    for (int mode = 0; mode < 5; ++mode) {
        vector<pair<uint64_t,int>> ord; ord.reserve(n);
        for (int i = 0; i < n; ++i) {
            uint32_t x = X[i], y = Y[i];
            if (mode == 1) y = M - y;
            else if (mode == 2) x = M - x;
            else if (mode == 3) swap(x, y);
            uint64_t key;
            if (mode == 4) {
                key = 0;
                for (int b = 0; b < 31; ++b)
                    key |= (uint64_t((x >> b) & 1) << (2*b)) | (uint64_t((y >> b) & 1) << (2*b+1));
            } else key = hilbert(x, y);
            ord.push_back({key, i});
        }
        tryCandidate(makeRoute(ord), best, bestCost);
    }

    // A deliberately different representation: partition the plane into wide
    // strips, then sweep each strip monotonically in alternating directions.
    // Unlike a Hilbert curve this preserves long runs within horizontal/vertical
    // layers, which is useful for anisotropic or banded point distributions.
    for (int axis = 0; axis < 2; ++axis) for (int bits : {2, 4, 6, 8}) {
        uint32_t bands = 1u << bits;
        vector<pair<uint64_t,int>> ord; ord.reserve(n);
        for (int i = 0; i < n; ++i) {
            uint32_t across = axis ? X[i] : Y[i];
            uint32_t along  = axis ? Y[i] : X[i];
            uint32_t b = across >> (31 - bits);
            uint32_t walk = (b & 1) ? (M - along) : along;
            ord.push_back({(uint64_t(b) << 31) | walk, i});
        }
        tryCandidate(makeRoute(ord), best, bestCost);
    }

    // Rotated strips are a separate sweep family: they keep diagonal layers
    // contiguous when x and y are strongly correlated or anti-correlated.
    // Only two resolutions are used to retain the incumbent's large-N budget.
    vector<long long> Dplus(n), Dminus(n);
    long long minPlus = LLONG_MAX, maxPlus = LLONG_MIN;
    long long minMinus = LLONG_MAX, maxMinus = LLONG_MIN;
    for (int i = 0; i < n; ++i) {
        Dplus[i] = p[i].x + p[i].y;
        Dminus[i] = p[i].x - p[i].y;
        minPlus = min(minPlus, Dplus[i]); maxPlus = max(maxPlus, Dplus[i]);
        minMinus = min(minMinus, Dminus[i]); maxMinus = max(maxMinus, Dminus[i]);
    }
    auto scaleDiagonal = [&](long long v, long long lo, long long hi) -> uint32_t {
        if (lo == hi) return M / 2;
        return (uint64_t)(v - lo) * M / (uint64_t)(hi - lo);
    };
    for (int orient = 0; orient < 2; ++orient) for (int bits : {3, 6}) {
        vector<pair<uint64_t,int>> ord; ord.reserve(n);
        for (int i = 0; i < n; ++i) {
            uint32_t across, along;
            if (!orient) {
                across = scaleDiagonal(Dplus[i], minPlus, maxPlus);
                along = scaleDiagonal(Dminus[i], minMinus, maxMinus);
            } else {
                across = scaleDiagonal(Dminus[i], minMinus, maxMinus);
                along = scaleDiagonal(Dplus[i], minPlus, maxPlus);
            }
            uint32_t b = across >> (31 - bits);
            uint32_t walk = (b & 1) ? (M - along) : along;
            ord.push_back({(uint64_t(b) << 31) | walk, i});
        }
        tryCandidate(makeRoute(ord), best, bestCost);
    }

    // Exact adjacent-swap descent, including the position-dependent carrot term.
    auto adjacentDescent = [&] {
        for (int pass = 0; pass < 4; ++pass) {
            bool changed = false;
            for (int l = 1; l + 1 < n; ++l) {
                double oldc = edge(best[l - 1], best[l], l - 1)
                            + edge(best[l], best[l + 1], l)
                            + edge(best[l + 1], best[(l + 2) % n], l + 1);
                double newc = edge(best[l - 1], best[l + 1], l - 1)
                            + edge(best[l + 1], best[l], l)
                            + edge(best[l], best[(l + 2) % n], l + 1);
                if (newc + 1e-8 < oldc) {
                    swap(best[l], best[l + 1]);
                    changed = true;
                }
            }
            if (!changed) break;
        }
    };
    adjacentDescent();
    // Small exact reversals repair local crossings without quadratic work.
    int maxDepth = (n <= 100000 ? 12 : 5);
    for (int l = 1; l < n; ++l) {
        for (int d = 1; d <= maxDepth && l + d < n; ++d) {
            int r = l + d;
            double oldc = 0, newc = 0;
            for (int pos = l - 1; pos <= r; ++pos) {
                oldc += edge(best[pos], best[(pos + 1) % n], pos);
                auto at = [&](int k) { return (k >= l && k <= r) ? best[l + r - k] : best[k]; };
                newc += edge(at(pos), at((pos + 1) % n), pos);
            }
            if (newc + 1e-8 < oldc) reverse(best.begin() + l, best.begin() + r + 1);
        }
    }
    // Reversals change neighboring pairs, so restore adjacent optimality afterward.
    adjacentDescent();

    cout << n + 1 << '\n';
    for (int v : best) cout << v << '\n';
    cout << 0 << '\n';
    return 0;
}
