#include <bits/stdc++.h>
using namespace std;

struct City {
    long long x, y;
    int id;
};

static inline uint64_t hilbertOrder(uint32_t x, uint32_t y, int bits = 31) {
    uint64_t d = 0;
    for (int s = bits - 1; s >= 0; --s) {
        uint32_t rx = (x >> s) & 1U;
        uint32_t ry = (y >> s) & 1U;
        d += uint64_t((3U * rx) ^ ry) << (2 * s);
        if (ry == 0) {
            if (rx == 1) {
                uint32_t mask = (uint32_t(1) << bits) - 1U;
                x = mask - x;
                y = mask - y;
            }
            swap(x, y);
        }
    }
    return d;
}

static inline double dist2(const City& a, const City& b) {
    return hypot((double)a.x - (double)b.x, (double)a.y - (double)b.y);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    if (!(cin >> N)) return 0;
    vector<City> c(N);
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (int i = 0; i < N; ++i) {
        cin >> c[i].x >> c[i].y;
        c[i].id = i;
        minx = min(minx, c[i].x); maxx = max(maxx, c[i].x);
        miny = min(miny, c[i].y); maxy = max(maxy, c[i].y);
    }

    vector<char> prime(max(2, N), true);
    prime[0] = false; if (N > 1) prime[1] = false;
    for (int i = 2; 1LL * i * i < N; ++i) if (prime[i])
        for (long long j = 1LL * i * i; j < N; j += i) prime[(int)j] = false;

    auto cost = [&](const vector<int>& route) {
        double total = 0.0;
        for (int t = 1; t <= N; ++t) {
            int a = route[t - 1], b = route[t];
            double m = (t % 10 == 0 && !prime[a]) ? 1.1 : 1.0;
            total += m * dist2(c[a], c[b]);
        }
        return total;
    };

    vector<int> best;
    double bestCost = numeric_limits<double>::infinity();
    auto tryRoute = [&](vector<int> mid) {
        vector<int> r;
        r.reserve(N + 1);
        r.push_back(0);
        for (int id : mid) if (id != 0) r.push_back(id);
        r.push_back(0);
        if ((int)r.size() != N + 1) return;
        double v = cost(r);
        if (v < bestCost) { bestCost = v; best = move(r); }
    };

    // Always include the strengthened baseline as a safety net.
    vector<int> baseline;
    baseline.reserve(N - 1);
    for (int i = 1; i < N; ++i) baseline.push_back(i);
    tryRoute(baseline);

    const uint64_t M = (uint64_t(1) << 31) - 1;
    auto scaleCoord = [&](long long v, long long lo, long long hi) -> uint32_t {
        if (hi == lo) return 0;
        __int128 num = (__int128)(v - lo) * M;
        return (uint32_t)(num / (hi - lo));
    };

    vector<pair<uint64_t,int>> keys;
    keys.reserve(N - 1);
    for (int variant = 0; variant < 8; ++variant) {
        keys.clear();
        bool sw = variant & 1, invx = variant & 2, invy = variant & 4;
        for (int i = 1; i < N; ++i) {
            uint32_t sx = scaleCoord(c[i].x, minx, maxx);
            uint32_t sy = scaleCoord(c[i].y, miny, maxy);
            if (sw) swap(sx, sy);
            if (invx) sx = (uint32_t)M - sx;
            if (invy) sy = (uint32_t)M - sy;
            keys.emplace_back(hilbertOrder(sx, sy), i);
        }
        sort(keys.begin(), keys.end(), [](const auto& a, const auto& b){
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });
        vector<int> mid;
        mid.reserve(N - 1);
        for (auto &kv : keys) mid.push_back(kv.second);
        tryRoute(mid);
        reverse(mid.begin(), mid.end());
        tryRoute(mid);
    }

    // Cheap carrot-aware local refinement: accept short relocations that improve the
    // exact penalized cost on the affected span.  This generalizes adjacent swaps and
    // can repair small Hilbert-order glitches without risking O(N^2) search.
    if (N <= 200000 && !best.empty()) {
        auto edge = [&](int t, int u, int v) {
            double m = (t % 10 == 0 && !prime[u]) ? 1.1 : 1.0;
            return m * dist2(c[u], c[v]);
        };
        const int W = 5;
        bool changed = true;
        for (int pass = 0; pass < 2 && changed; ++pass) {
            changed = false;
            for (int i = 1; i < N; ++i) {
                bool did = false;
                for (int d = -W; d <= W && !did; ++d) {
                    if (d == 0) continue;
                    int k = i + d;
                    if (k < 1 || k >= N) continue;
                    int L = min(i, k), R = min(N, max(i, k) + 1);
                    int moved = best[i];
                    auto atNew = [&](int pos) -> int {
                        if (pos <= 0 || pos >= N) return best[pos];
                        if (d > 0) {
                            if (pos >= i && pos < k) return best[pos + 1];
                            if (pos == k) return moved;
                        } else {
                            if (pos == k) return moved;
                            if (pos > k && pos <= i) return best[pos - 1];
                        }
                        return best[pos];
                    };
                    double oldv = 0.0, newv = 0.0;
                    for (int t = L; t <= R; ++t) {
                        oldv += edge(t, best[t - 1], best[t]);
                        newv += edge(t, atNew(t - 1), atNew(t));
                    }
                    if (newv + 1e-9 < oldv) {
                        if (d > 0) {
                            for (int p = i; p < k; ++p) best[p] = best[p + 1];
                            best[k] = moved;
                        } else {
                            for (int p = i; p > k; --p) best[p] = best[p - 1];
                            best[k] = moved;
                        }
                        changed = did = true;
                    }
                }
            }
        }
    }

    cout << N + 1 << '\n';
    for (int id : best) cout << id << '\n';
    return 0;
}
