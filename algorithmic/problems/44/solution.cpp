#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };
static vector<Point> pt;
static vector<char> primeId;

static inline double edgeCost(int a, int b, int pos) { // pos is the zero-based source position
    double dx = (double)pt[a].x - pt[b].x;
    double dy = (double)pt[a].y - pt[b].y;
    double w = ((pos + 1) % 10 == 0 && !primeId[a]) ? 1.1 : 1.0;
    return w * hypot(dx, dy);
}

static double value(const vector<int>& q) {
    const int n = (int)q.size();
    double ans = 0;
    for (int i = 0; i < n; ++i) ans += edgeCost(q[i], q[(i + 1) % n], i);
    return ans;
}

// Exact delta for reversing the inclusive non-root interval [l,r].
static double reversalDelta(const vector<int>& q, int l, int r) {
    const int n = (int)q.size();
    double oldv = edgeCost(q[l - 1], q[l], l - 1);
    double newv = edgeCost(q[l - 1], q[r], l - 1);
    for (int p = l; p < r; ++p) {
        oldv += edgeCost(q[p], q[p + 1], p);
        int src = q[l + r - p];
        int dst = q[l + r - p - 1];
        newv += edgeCost(src, dst, p);
    }
    int after = (r + 1) % n;
    oldv += edgeCost(q[r], q[after], r);
    newv += edgeCost(q[l], q[after], r);
    return newv - oldv;
}

// A deliberately small neighbourhood keeps this robust at 200k cities while
// still removing local discontinuities introduced by a space filling order.
static void polish(vector<int>& q, int passes) {
    const int n = (int)q.size();
    const int span = 8;
    for (int pass = 0; pass < passes; ++pass) {
        bool changed = false;
        for (int l = 1; l < n - 1; ++l) {
            int hi = min(n - 1, l + span - 1);
            for (int r = l + 1; r <= hi; ++r) {
                if (reversalDelta(q, l, r) < -1e-7) {
                    reverse(q.begin() + l, q.begin() + r + 1);
                    changed = true;
                    break;
                }
            }
        }
        if (!changed) break;
    }
}

// Delta for moving q[i] immediately after q[j].  Unlike a conventional
// Or-opt delta, every affected source position is scored because its carrot
// multiplier may change when the intervening block shifts.
static double relocationDelta(const vector<int>& q, int i, int j) {
    const int n = (int)q.size();
    if (i == j || i == j + 1) return 0.0;
    const int lo = (i < j ? i - 1 : j);
    const int hi = (i < j ? j : i);
    auto movedAt = [&](int p) {
        if (i < j) {
            if (p >= i && p < j) return q[p + 1];
            if (p == j) return q[i];
        } else {
            if (p == j + 1) return q[i];
            if (p >= j + 2 && p <= i) return q[p - 1];
        }
        return q[p];
    };
    double oldv = 0, newv = 0;
    for (int p = lo; p <= hi; ++p) {
        oldv += edgeCost(q[p], q[(p + 1) % n], p);
        newv += edgeCost(movedAt(p), movedAt((p + 1) % n), p);
    }
    return newv - oldv;
}

static void relocatePolish(vector<int>& q) {
    const int n = (int)q.size();
    // This is intentionally a distinct, small Or-opt neighborhood.  Moving
    // farther is both costly and already well represented by the space order.
    for (int i = 1; i < n; ++i) {
        bool changed = false;
        for (int d = 2; d <= 6 && !changed; ++d) {
            for (int j : {i - d, i + d}) {
                if (j < 0 || j >= n) continue;
                if (relocationDelta(q, i, j) < -1e-7) {
                    // The displacement is bounded, so rotate applies this
                    // move in O(6), rather than an O(n) vector erase/insert.
                    if (i < j)
                        rotate(q.begin() + i, q.begin() + i + 1, q.begin() + j + 1);
                    else
                        rotate(q.begin() + j + 1, q.begin() + i, q.begin() + i + 1);
                    changed = true;
                    break;
                }
            }
        }
    }
}

static uint64_t hilbertKey(uint32_t x, uint32_t y) {
    constexpr uint32_t B = 20;
    uint64_t d = 0;
    for (uint32_t s = 1u << (B - 1); s; s >>= 1) {
        uint32_t rx = (x & s) != 0, ry = (y & s) != 0;
        d += (uint64_t)s * s * ((3 * rx) ^ ry);
        if (!ry) {
            // Hilbert rotation reflects inside this level's active square,
            // not across the full 20-bit coordinate range.
            if (rx) {
                uint32_t side = s << 1;
                x = side - 1 - x;
                y = side - 1 - y;
            }
            swap(x, y);
        }
    }
    return d;
}

static uint64_t mortonKey(uint32_t x, uint32_t y) {
    uint64_t z = 0;
    for (int b = 0; b < 20; ++b) {
        z |= (uint64_t)((x >> b) & 1) << (2 * b + 1);
        z |= (uint64_t)((y >> b) & 1) << (2 * b);
    }
    return z;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    pt.resize(n);
    long long minx = LLONG_MAX, maxx = LLONG_MIN, miny = LLONG_MAX, maxy = LLONG_MIN;
    for (auto &p : pt) {
        cin >> p.x >> p.y;
        minx = min(minx, p.x); maxx = max(maxx, p.x);
        miny = min(miny, p.y); maxy = max(maxy, p.y);
    }
    primeId.assign(n, true);
    if (n > 0) primeId[0] = false;
    if (n > 1) primeId[1] = false;
    for (int i = 2; (long long)i * i < n; ++i) if (primeId[i])
        for (int j = i * i; j < n; j += i) primeId[j] = false;

    const uint64_t SCALE = (1u << 20) - 1;
    unsigned long long rx = (unsigned long long)(maxx - minx);
    unsigned long long ry = (unsigned long long)(maxy - miny);
    vector<pair<uint64_t,int>> h, z;
    h.reserve(max(0, n - 1)); z.reserve(max(0, n - 1));
    for (int i = 1; i < n; ++i) {
        uint32_t x = rx ? (uint64_t)(pt[i].x - minx) * SCALE / rx : 0;
        uint32_t y = ry ? (uint64_t)(pt[i].y - miny) * SCALE / ry : 0;
        h.push_back({hilbertKey(x, y), i});
        z.push_back({mortonKey(x, y), i});
    }
    sort(h.begin(), h.end());
    sort(z.begin(), z.end());

    vector<int> best;
    double bestVal = numeric_limits<double>::infinity();
    auto consider = [&](const vector<pair<uint64_t,int>>& keys, bool backwards, int passes) {
        vector<int> q; q.reserve(n); q.push_back(0);
        if (!backwards) for (auto e : keys) q.push_back(e.second);
        else for (auto it = keys.rbegin(); it != keys.rend(); ++it) q.push_back(it->second);
        polish(q, passes);
        relocatePolish(q);
        double v = value(q);
        if (v < bestVal) { bestVal = v; best.swap(q); }
    };
    // Fully retain the incumbent's two Hilbert passes; each candidate also
    // gets a single exact relocation sweep, while Morton remains independent.
    consider(h, false, 2); consider(h, true, 2);
    consider(z, false, 2); consider(z, true, 2);

    cout << n + 1 << '\n';
    for (int id : best) cout << id << '\n';
    cout << 0 << '\n';
    return 0;
}
