#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x, y; };

static uint64_t hilbertIndex(uint32_t x, uint32_t y, int bits) {
    uint64_t d = 0;
    for (int s = bits - 1; s >= 0; --s) {
        uint32_t rx = (x >> s) & 1u;
        uint32_t ry = (y >> s) & 1u;
        d += uint64_t((3u * rx) ^ ry) << (2 * s);
        if (ry == 0) {
            if (rx == 1) {
                uint32_t n = (1u << bits) - 1u;
                x = n - x;
                y = n - y;
            }
            swap(x, y);
        }
    }
    return d;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    if (!(cin >> N)) return 0;
    vector<Pt> p(N);
    for (int i = 0; i < N; ++i) cin >> p[i].x >> p[i].y;

    vector<char> prime(max(2, N), true);
    prime[0] = false; prime[1] = false;
    for (long long i = 2; i * i < N; ++i) if (prime[i])
        for (long long j = i * i; j < N; j += i) prime[(int)j] = false;

    auto dist = [&](int a, int b) -> double {
        long double dx = (long double)p[a].x - p[b].x;
        long double dy = (long double)p[a].y - p[b].y;
        return (double)sqrt(dx * dx + dy * dy);
    };
    auto edgeCost = [&](const vector<int>& r, int k) -> double { // edge k: r[k] -> r[k+1]
        if (k < 0 || k >= N) return 0.0;
        double m = ((k + 1) % 10 == 0 && !prime[r[k]]) ? 1.1 : 1.0;
        return m * dist(r[k], r[k + 1]);
    };
    auto routeCost = [&](const vector<int>& r) -> double {
        double s = 0.0;
        for (int k = 0; k < N; ++k) s += edgeCost(r, k);
        return s;
    };

    auto makeCycleRoute = [&](const vector<int>& ord, bool rev) {
        vector<int> r;
        r.reserve(N + 1);
        r.push_back(0);
        int z = int(find(ord.begin(), ord.end(), 0) - ord.begin());
        if (!rev) {
            for (int i = z + 1; i < N; ++i) r.push_back(ord[i]);
            for (int i = 0; i < z; ++i) r.push_back(ord[i]);
        } else {
            for (int i = z - 1; i >= 0; --i) r.push_back(ord[i]);
            for (int i = N - 1; i > z; --i) r.push_back(ord[i]);
        }
        r.push_back(0);
        return r;
    };

    vector<int> best(N + 1);
    for (int i = 0; i < N; ++i) best[i] = i;
    best[N] = 0;
    double bestCost = routeCost(best);

    long long minx = p[0].x, maxx = p[0].x, miny = p[0].y, maxy = p[0].y;
    for (auto &q : p) {
        minx = min(minx, q.x); maxx = max(maxx, q.x);
        miny = min(miny, q.y); maxy = max(maxy, q.y);
    }
    const int BITS = 21;
    const long long SCALE = (1LL << BITS) - 1;
    long long rangex = max(1LL, maxx - minx), rangey = max(1LL, maxy - miny);
    vector<uint32_t> nx(N), ny(N);
    for (int i = 0; i < N; ++i) {
        nx[i] = (uint32_t)((__int128)(p[i].x - minx) * SCALE / rangex);
        ny[i] = (uint32_t)((__int128)(p[i].y - miny) * SCALE / rangey);
    }

    // Try all eight symmetries of the bounding box before taking the Hilbert
    // order.  The curve is directional; completing this small dihedral set is a
    // cheap way to avoid unlucky long jumps caused only by orientation.
    for (int variant = 0; variant < 8; ++variant) {
        vector<pair<uint64_t,int>> keyed;
        keyed.reserve(N);
        for (int i = 0; i < N; ++i) {
            uint32_t x0 = nx[i], y0 = ny[i];
            uint32_t a, b;
            if (variant < 4) {
                a = (variant & 1) ? (uint32_t)SCALE - x0 : x0;
                b = (variant & 2) ? (uint32_t)SCALE - y0 : y0;
            } else {
                a = (variant & 1) ? (uint32_t)SCALE - y0 : y0;
                b = (variant & 2) ? (uint32_t)SCALE - x0 : x0;
            }
            keyed.push_back({hilbertIndex(a, b, BITS), i});
        }
        sort(keyed.begin(), keyed.end());
        vector<int> ord(N);
        for (int i = 0; i < N; ++i) ord[i] = keyed[i].second;
        for (int rev = 0; rev < 2; ++rev) {
            vector<int> r = makeCycleRoute(ord, rev);
            double c = routeCost(r);
            if (c < bestCost) { bestCost = c; best.swap(r); }
        }
    }

    auto deltaSwap = [&](const vector<int>& r, int a, int b) -> double {
        if (a == b || a <= 0 || b <= 0 || a >= N || b >= N) return 0.0;
        vector<int> ks = {a - 1, a, b - 1, b};
        sort(ks.begin(), ks.end());
        ks.erase(unique(ks.begin(), ks.end()), ks.end());
        double before = 0.0, after = 0.0;
        auto at = [&](int idx) -> int {
            if (idx == a) return r[b];
            if (idx == b) return r[a];
            return r[idx];
        };
        for (int k : ks) {
            before += edgeCost(r, k);
            if (k >= 0 && k < N) {
                double m = ((k + 1) % 10 == 0 && !prime[at(k)]) ? 1.1 : 1.0;
                after += m * dist(at(k), at(k + 1));
            }
        }
        return after - before;
    };

    auto deltaReverse = [&](const vector<int>& r, int a, int b) -> double {
        if (a <= 0 || b <= a || b >= N) return 0.0;
        double before = 0.0, after = 0.0;
        auto at = [&](int idx) -> int {
            if (idx >= a && idx <= b) return r[a + b - idx];
            return r[idx];
        };
        for (int k = a - 1; k <= b; ++k) {
            before += edgeCost(r, k);
            double m = ((k + 1) % 10 == 0 && !prime[at(k)]) ? 1.1 : 1.0;
            after += m * dist(at(k), at(k + 1));
        }
        return after - before;
    };

    // Bounded exact 2-opt polishing.  Hilbert order is usually close, so most
    // useful uncrossings are short in route order; evaluating short reversals
    // exactly also accounts for the 10th-step prime penalty inside the segment.
    const int REV_W = 14;
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int i = 1; i + 1 < N; ++i) {
            int bj = -1; double bd = -1e-7;
            int hi = min(N - 1, i + REV_W);
            for (int j = i + 1; j <= hi; ++j) {
                double d = deltaReverse(best, i, j);
                if (d < bd) { bd = d; bj = j; }
            }
            if (bj != -1) {
                reverse(best.begin() + i, best.begin() + bj + 1);
                bestCost += bd;
                changed = true;
                i = bj;
            }
        }
        if (!changed) break;
    }

    // Local carrot-aware polishing: put nearby prime city IDs at every tenth source
    // only when the exact penalized tour length decreases.
    const int W = 35;
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int pos = 9; pos <= N - 1; pos += 10) {
            if (prime[best[pos]]) continue;
            int bj = -1; double bd = -1e-9;
            int lo = max(1, pos - W), hi = min(N - 1, pos + W);
            for (int j = lo; j <= hi; ++j) {
                if (j == pos || !prime[best[j]]) continue;
                double d = deltaSwap(best, pos, j);
                if (d < bd) { bd = d; bj = j; }
            }
            if (bj != -1) {
                swap(best[pos], best[bj]);
                bestCost += bd;
                changed = true;
            }
        }
        if (!changed) break;
    }

    // One cheap adjacent-swap descent catches small ordering mistakes from ties and penalties.
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int i = 1; i + 1 < N; ++i) {
            double d = deltaSwap(best, i, i + 1);
            if (d < -1e-7) {
                swap(best[i], best[i + 1]);
                bestCost += d;
                changed = true;
                ++i;
            }
        }
        if (!changed) break;
    }

    cout << N + 1 << '\n';
    for (int v : best) cout << v << '\n';
    return 0;
}
