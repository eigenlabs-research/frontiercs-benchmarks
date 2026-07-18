#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x, y; };

static vector<char> sievePrime(int n) {
    vector<char> p(max(2, n), true);
    p[0] = false;
    if (n > 1) p[1] = false;
    for (long long i = 2; i * i < n; ++i) if (p[(int)i])
        for (long long j = i * i; j < n; j += i) p[(int)j] = false;
    return p;
}

static inline double euclid(const Pt& a, const Pt& b) {
    return hypot((double)a.x - (double)b.x, (double)a.y - (double)b.y);
}

static inline uint64_t hilbertIndex(uint32_t x, uint32_t y, int bits = 21) {
    uint64_t d = 0;
    uint32_t xx = x, yy = y;
    for (uint32_t s = 1u << (bits - 1); s; s >>= 1) {
        uint32_t rx = (xx & s) ? 1u : 0u;
        uint32_t ry = (yy & s) ? 1u : 0u;
        d += (uint64_t)s * (uint64_t)s * ((3u * rx) ^ ry);
        if (ry == 0) {
            if (rx == 1) {
                xx = (s - 1) - xx;
                yy = (s - 1) - yy;
            }
            swap(xx, yy);
        }
    }
    return d;
}

static inline pair<long long,long long> transformPoint(const Pt& p, int v) {
    long long x = p.x, y = p.y;
    switch (v) {
        case 0: return { x,  y};
        case 1: return { x, -y};
        case 2: return {-x,  y};
        case 3: return {-x, -y};
        case 4: return { y,  x};
        case 5: return { y, -x};
        case 6: return {-y,  x};
        default:return {-y, -x};
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    if (!(cin >> N)) return 0;
    vector<Pt> a(N);
    for (int i = 0; i < N; ++i) cin >> a[i].x >> a[i].y;

    vector<char> isPrime = sievePrime(N);

    auto edgeCost = [&](int step, int u, int v) -> double {
        double m = (step % 10 == 0 && !isPrime[u]) ? 1.1 : 1.0;
        return m * euclid(a[u], a[v]);
    };
    auto routeCost = [&](const vector<int>& r) -> double {
        double s = 0.0;
        for (int t = 1; t <= N; ++t) s += edgeCost(t, r[t-1], r[t]);
        return s;
    };

    vector<int> best(N + 1);
    best[0] = 0;
    for (int i = 1; i < N; ++i) best[i] = i;
    best[N] = 0;
    double bestCost = routeCost(best);
    auto consider = [&](const vector<int>& cand) {
        double c = routeCost(cand);
        if (c < bestCost) { bestCost = c; best = cand; }
    };

    const int BITS = 21;
    const long long GRID = (1LL << BITS) - 1;
    vector<pair<uint64_t,int>> keys;
    keys.reserve(max(0, N - 1));

    for (int v = 0; v < 8; ++v) {
        long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
        for (int i = 1; i < N; ++i) {
            auto q = transformPoint(a[i], v);
            minx = min(minx, q.first); maxx = max(maxx, q.first);
            miny = min(miny, q.second); maxy = max(maxy, q.second);
        }
        long long rangex = max(1LL, maxx - minx);
        long long rangey = max(1LL, maxy - miny);
        keys.clear();
        for (int i = 1; i < N; ++i) {
            auto q = transformPoint(a[i], v);
            uint32_t gx = (uint32_t)((__int128)(q.first - minx) * GRID / rangex);
            uint32_t gy = (uint32_t)((__int128)(q.second - miny) * GRID / rangey);
            keys.push_back({hilbertIndex(gx, gy, BITS), i});
        }
        sort(keys.begin(), keys.end(), [](const auto& A, const auto& B){
            if (A.first != B.first) return A.first < B.first;
            return A.second < B.second;
        });
        vector<int> cand(N + 1);
        cand[0] = 0; cand[N] = 0;
        for (int i = 1; i < N; ++i) cand[i] = keys[i-1].second;
        consider(cand);
        reverse(cand.begin() + 1, cand.begin() + N);
        consider(cand);
    }

    // Additional initial family: angular sweeps around natural centers.  This is cheap
    // and complements Hilbert order on ring/star or few-radial-cluster instances where
    // preserving angle avoids repeated long chords.
    auto addAngular = [&](long double cx, long double cy) {
        int M = N - 1;
        if (M <= 0) return;
        vector<pair<pair<double,double>,int>> ord;
        ord.reserve(M);
        for (int i = 1; i < N; ++i) {
            double ang = atan2((long double)a[i].y - cy, (long double)a[i].x - cx);
            double rr = ((double)a[i].x - (double)cx) * ((double)a[i].x - (double)cx)
                      + ((double)a[i].y - (double)cy) * ((double)a[i].y - (double)cy);
            ord.push_back({{ang, rr}, i});
        }
        sort(ord.begin(), ord.end(), [](const auto& A, const auto& B){
            if (A.first.first != B.first.first) return A.first.first < B.first.first;
            if (A.first.second != B.first.second) return A.first.second < B.first.second;
            return A.second < B.second;
        });
        vector<int> cyc(M);
        for (int i = 0; i < M; ++i) cyc[i] = ord[i].second;
        auto emitBestBreak = [&](const vector<int>& c) {
            int start = 0;
            double bestBreak = 1e300;
            for (int s = 0; s < M; ++s) {
                int first = c[s], last = c[(s + M - 1) % M];
                double v = euclid(a[0], a[first])
                         + ((N % 10 == 0 && !isPrime[last]) ? 1.1 : 1.0) * euclid(a[last], a[0])
                         - euclid(a[last], a[first]);
                if (v < bestBreak) { bestBreak = v; start = s; }
            }
            vector<int> cand(N + 1);
            cand[0] = 0; cand[N] = 0;
            for (int j = 0; j < M; ++j) cand[j + 1] = c[(start + j) % M];
            consider(cand);
        };
        emitBestBreak(cyc);
        reverse(cyc.begin(), cyc.end());
        emitBestBreak(cyc);
    };
    long double sx = 0, sy = 0;
    long long minx0 = LLONG_MAX, maxx0 = LLONG_MIN, miny0 = LLONG_MAX, maxy0 = LLONG_MIN;
    for (int i = 0; i < N; ++i) {
        sx += a[i].x; sy += a[i].y;
        minx0 = min(minx0, a[i].x); maxx0 = max(maxx0, a[i].x);
        miny0 = min(miny0, a[i].y); maxy0 = max(maxy0, a[i].y);
    }
    addAngular((long double)a[0].x, (long double)a[0].y);
    addAngular(sx / N, sy / N);
    addAngular(((long double)minx0 + maxx0) / 2.0L, ((long double)miny0 + maxy0) / 2.0L);

    // Cheap penalty-aware local polishing: adjacent exchanges accepted only by exact cost delta.
    for (int pass = 0; pass < 4; ++pass) {
        bool any = false;
        for (int i = 1; i <= N - 2; ++i) {
            int A = best[i-1], B = best[i], C = best[i+1], D = best[i+2];
            double before = edgeCost(i, A, B) + edgeCost(i+1, B, C) + edgeCost(i+2, C, D);
            double after  = edgeCost(i, A, C) + edgeCost(i+1, C, B) + edgeCost(i+2, B, D);
            if (after + 1e-7 < before) {
                swap(best[i], best[i+1]);
                bestCost += after - before;
                any = true;
                ++i;
            }
        }
        if (!any) break;
    }

    // Slightly wider bounded-window reordering.  This keeps the Hilbert skeleton but
    // can remove small kinks and place prime/non-prime sources on better modulo-10
    // steps.  The delta is evaluated exactly for the affected penalized edges.
    auto blockCost = [&](const vector<int>& r, int l, int w) -> double {
        double s = 0.0;
        for (int step = l; step <= l + w; ++step) s += edgeCost(step, r[step-1], r[step]);
        return s;
    };
    auto windowPolish = [&](int w, int passes) {
        if (N - 1 < w) return;
        vector<int> idx(w), cur(w), bestPerm(w);
        iota(idx.begin(), idx.end(), 0);
        for (int pass = 0; pass < passes; ++pass) {
            bool any = false;
            for (int l = 1; l + w - 1 <= N - 1; ++l) {
                for (int j = 0; j < w; ++j) cur[j] = best[l + j];
                double old = blockCost(best, l, w);
                double valBest = old;
                bestPerm = cur;
                sort(idx.begin(), idx.end());
                do {
                    bool identity = true;
                    for (int j = 0; j < w; ++j) if (idx[j] != j) { identity = false; break; }
                    if (identity) continue; // the current order was already costed
                    for (int j = 0; j < w; ++j) best[l + j] = cur[idx[j]];
                    double val = blockCost(best, l, w);
                    if (val + 1e-7 < valBest) {
                        valBest = val;
                        for (int j = 0; j < w; ++j) bestPerm[j] = best[l + j];
                    }
                } while (next_permutation(idx.begin(), idx.end()));
                for (int j = 0; j < w; ++j) best[l + j] = bestPerm[j];
                if (valBest + 1e-7 < old) {
                    bestCost += valBest - old;
                    any = true;
                    l += w - 1;
                }
            }
            if (!any) break;
        }
    };
    windowPolish(3, 2);
    if (N <= 50000) windowPolish(4, 1);

    cout << N + 1 << '\n';
    for (int id : best) cout << id << '\n';
    return 0;
}
