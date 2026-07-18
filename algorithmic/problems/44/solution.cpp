#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x, y; };

static inline double euclid(const Pt& a, const Pt& b) {
    return hypot((double)a.x - (double)b.x, (double)a.y - (double)b.y);
}

static vector<char> sievePrime(int n) {
    vector<char> p(max(2, n), true);
    p[0] = false; p[1] = false;
    for (long long i = 2; i * i < n; ++i) if (p[(int)i])
        for (long long j = i * i; j < n; j += i) p[(int)j] = false;
    return p;
}

// Hilbert order for a square [0,2^bits).  bits<=31, result fits in uint64_t.
static uint64_t hilbertOrder(uint32_t x, uint32_t y, int bits) {
    uint64_t d = 0;
    for (int s = bits - 1; s >= 0; --s) {
        uint32_t rx = (x >> s) & 1U;
        uint32_t ry = (y >> s) & 1U;
        d += (uint64_t(1) << (2 * s)) * ((3U * rx) ^ ry);
        if (ry == 0) {
            if (rx == 1) {
                uint32_t mask = (bits == 32 ? 0xffffffffU : ((1U << bits) - 1U));
                x = mask - x;
                y = mask - y;
            }
            swap(x, y);
        }
    }
    return d;
}

static double routeCost(const vector<int>& r, const vector<Pt>& pts, const vector<char>& prime) {
    int N = (int)pts.size();
    double total = 0.0;
    for (int t = 1; t <= N; ++t) {
        double m = (t % 10 == 0 && !prime[r[t-1]]) ? 1.1 : 1.0;
        total += m * euclid(pts[r[t-1]], pts[r[t]]);
    }
    return total;
}

static inline double stepCost(int t, int a, int b, const vector<Pt>& pts, const vector<char>& prime) {
    double m = (t % 10 == 0 && !prime[a]) ? 1.1 : 1.0;
    return m * euclid(pts[a], pts[b]);
}

static bool trySwapDelta(vector<int>& r, int i, const vector<Pt>& pts, const vector<char>& prime) {
    // swap internal adjacent positions i and i+1, where 1 <= i < N-1
    int N = (int)pts.size();
    int a0 = r[i-1], a = r[i], b = r[i+1], c = r[i+2];
    double oldc = stepCost(i, a0, a, pts, prime)
                + stepCost(i+1, a, b, pts, prime)
                + stepCost(i+2, b, c, pts, prime);
    double newc = stepCost(i, a0, b, pts, prime)
                + stepCost(i+1, b, a, pts, prime)
                + stepCost(i+2, a, c, pts, prime);
    if (newc + 1e-9 < oldc) {
        swap(r[i], r[i+1]);
        return true;
    }
    return false;
}

static void adjacentImprove(vector<int>& r, const vector<Pt>& pts, const vector<char>& prime, int passes) {
    int N = (int)pts.size();
    if (N <= 3) return;
    for (int pass = 0; pass < passes; ++pass) {
        bool any = false;
        if ((pass & 1) == 0) {
            for (int i = 1; i <= N - 2; ++i) any |= trySwapDelta(r, i, pts, prime);
        } else {
            for (int i = N - 2; i >= 1; --i) any |= trySwapDelta(r, i, pts, prime);
        }
        if (!any) break;
    }
}

static inline int nodeAfterReversal(const vector<int>& r, int pos, int l, int rr) {
    if (pos >= l && pos <= rr) return r[l + rr - pos];
    return r[pos];
}

static bool tryReverseDelta(vector<int>& r, int l, int rr,
                            const vector<Pt>& pts, const vector<char>& prime) {
    // Reverse internal positions [l, rr].  This is a bounded exact 2-opt move:
    // besides the two boundary edges, carrot penalties can change on every
    // 10th step inside the reversed segment, so score the whole small window.
    double oldc = 0.0, newc = 0.0;
    for (int t = l; t <= rr + 1; ++t) {
        oldc += stepCost(t, r[t-1], r[t], pts, prime);
        int a = nodeAfterReversal(r, t - 1, l, rr);
        int b = nodeAfterReversal(r, t,     l, rr);
        newc += stepCost(t, a, b, pts, prime);
    }
    if (newc + 1e-9 < oldc) {
        reverse(r.begin() + l, r.begin() + rr + 1);
        return true;
    }
    return false;
}

static void boundedTwoOptImprove(vector<int>& r, const vector<Pt>& pts,
                                  const vector<char>& prime) {
    int N = (int)pts.size();
    if (N <= 4) return;
    int maxLen;
    int passes;
    if (N <= 2000) { maxLen = 55; passes = 3; }
    else if (N <= 10000) { maxLen = 32; passes = 2; }
    else if (N <= 50000) { maxLen = 18; passes = 1; }
    else { maxLen = 10; passes = 1; }

    for (int pass = 0; pass < passes; ++pass) {
        bool any = false;
        for (int l = 1; l <= N - 2; ++l) {
            int hi = min(N - 1, l + maxLen - 1);
            for (int rr = l + 2; rr <= hi; ++rr) {
                if (tryReverseDelta(r, l, rr, pts, prime)) {
                    any = true;
                    break; // restart choices from the updated local order
                }
            }
        }
        if (!any) break;
    }
}

static vector<int> makeRouteFromOrder(const vector<int>& cyc, int pos0, bool forward) {
    int N = (int)cyc.size();
    vector<int> r; r.reserve(N + 1);
    r.push_back(0);
    if (forward) {
        for (int k = 1; k < N; ++k) {
            int id = cyc[(pos0 + k) % N];
            if (id) r.push_back(id);
        }
    } else {
        for (int k = 1; k < N; ++k) {
            int id = cyc[(pos0 - k + N) % N];
            if (id) r.push_back(id);
        }
    }
    r.push_back(0);
    return r;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    if (!(cin >> N)) return 0;
    vector<Pt> pts(N);
    for (int i = 0; i < N; ++i) cin >> pts[i].x >> pts[i].y;

    vector<char> prime = sievePrime(N);

    vector<vector<int>> candidates;
    candidates.reserve(20);

    // Baseline direction and its reverse are cheap safety candidates.
    vector<int> base; base.reserve(N+1); base.push_back(0);
    for (int i = 1; i < N; ++i) base.push_back(i);
    base.push_back(0);
    candidates.push_back(base);
    vector<int> rev; rev.reserve(N+1); rev.push_back(0);
    for (int i = N-1; i >= 1; --i) rev.push_back(i);
    rev.push_back(0);
    candidates.push_back(rev);

    const int BITS = 21; // enough resolution for 200k cities; keeps Hilbert arithmetic simple.
    const uint32_t MAXC = (1U << BITS) - 1U;

    for (int mode = 0; mode < 8; ++mode) {
        vector<long long> ax(N), ay(N);
        long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
        for (int i = 0; i < N; ++i) {
            long long u = pts[i].x, v = pts[i].y;
            if (mode & 4) swap(u, v);
            if (mode & 1) u = -u;
            if (mode & 2) v = -v;
            ax[i] = u; ay[i] = v;
            minx = min(minx, u); maxx = max(maxx, u);
            miny = min(miny, v); maxy = max(maxy, v);
        }
        long long rangex = max(1LL, maxx - minx), rangey = max(1LL, maxy - miny);
        vector<pair<uint64_t,int>> ord; ord.reserve(N);
        for (int i = 0; i < N; ++i) {
            uint32_t ux = (uint32_t)((__int128)(ax[i] - minx) * MAXC / rangex);
            uint32_t uy = (uint32_t)((__int128)(ay[i] - miny) * MAXC / rangey);
            uint64_t h = hilbertOrder(ux, uy, BITS);
            ord.push_back({h, i});
        }
        sort(ord.begin(), ord.end(), [](const auto& a, const auto& b){
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });
        vector<int> cyc(N);
        int pos0 = 0;
        for (int i = 0; i < N; ++i) { cyc[i] = ord[i].second; if (cyc[i] == 0) pos0 = i; }
        candidates.push_back(makeRouteFromOrder(cyc, pos0, true));
        candidates.push_back(makeRouteFromOrder(cyc, pos0, false));
    }

    int bestIdx = 0;
    double best = routeCost(candidates[0], pts, prime);
    for (int i = 1; i < (int)candidates.size(); ++i) {
        double c = routeCost(candidates[i], pts, prime);
        if (c < best) { best = c; bestIdx = i; }
    }

    vector<int> ans = std::move(candidates[bestIdx]);
    adjacentImprove(ans, pts, prime, N <= 5000 ? 8 : 3);
    boundedTwoOptImprove(ans, pts, prime);
    adjacentImprove(ans, pts, prime, N <= 5000 ? 4 : 1);

    cout << N + 1 << '\n';
    for (int v : ans) cout << v << '\n';
    return 0;
}
