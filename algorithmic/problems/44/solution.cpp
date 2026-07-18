#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x, y; };

static inline double d2d(const Pt& a, const Pt& b) {
    return hypot((double)a.x - (double)b.x, (double)a.y - (double)b.y);
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
    for (long long i = 2; i * i < N; ++i) if (prime[(int)i])
        for (long long j = i * i; j < N; j += i) prime[(int)j] = false;

    auto cost = [&](const vector<int>& r) -> double {
        double s = 0.0;
        for (int t = 1; t <= N; ++t) {
            int a = r[t-1], b = r[t];
            double m = (t % 10 == 0 && !prime[a]) ? 1.1 : 1.0;
            s += m * d2d(p[a], p[b]);
        }
        return s;
    };

    vector<int> best;
    double bestCost = numeric_limits<double>::infinity();

    auto consider = [&](vector<int>& mid) {
        vector<int> r;
        r.reserve(N + 1);
        r.push_back(0);
        for (int id : mid) if (id != 0) r.push_back(id);
        r.push_back(0);
        if ((int)r.size() != N + 1) return;
        double c = cost(r);
        if (c < bestCost) { bestCost = c; best.swap(r); }
    };

    // Baseline/order-by-input candidate: always valid and strong on nearly monotone instances.
    vector<int> mid;
    mid.reserve(max(0, N-1));
    for (int i = 1; i < N; ++i) mid.push_back(i);
    consider(mid);
    reverse(mid.begin(), mid.end());
    consider(mid);

    vector<int> ids;
    ids.reserve(max(0, N-1));
    for (int i = 1; i < N; ++i) ids.push_back(i);

    vector<int> bucketSizes;
    for (int b = 1; b < N; b <<= 1) bucketSizes.push_back(b);
    int extras[] = {3, 6, 10, 20, 40, 80, 160, 320, 640, 1280, 5000, 20000, 60000};
    for (int b : extras) if (b > 0 && b < N) bucketSizes.push_back(b);
    bucketSizes.push_back(max(1, N-1));
    sort(bucketSizes.begin(), bucketSizes.end());
    bucketSizes.erase(unique(bucketSizes.begin(), bucketSizes.end()), bucketSizes.end());

    // One general mechanism: spatial strip serpentine tours.  Sort by a primary
    // coordinate, cut into strips, and traverse each strip by the other coordinate,
    // alternating direction.  Trying both axes and parities handles rows, columns,
    // clusters and monotone cases without using sample-specific facts.
    for (int axis = 0; axis < 2; ++axis) {
        vector<int> base = ids;
        sort(base.begin(), base.end(), [&](int a, int b) {
            long long pa = axis == 0 ? p[a].x : p[a].y;
            long long pb = axis == 0 ? p[b].x : p[b].y;
            long long sa = axis == 0 ? p[a].y : p[a].x;
            long long sb = axis == 0 ? p[b].y : p[b].x;
            if (pa != pb) return pa < pb;
            if (sa != sb) return sa < sb;
            return a < b;
        });
        for (int revPrimary = 0; revPrimary < 2; ++revPrimary) {
            vector<int> prim = base;
            if (revPrimary) reverse(prim.begin(), prim.end());
            for (int B : bucketSizes) {
                for (int flip = 0; flip < 2; ++flip) {
                    vector<int> route = prim;
                    for (int l = 0, k = 0; l < (int)route.size(); l += B, ++k) {
                        int rr = min<int>(route.size(), l + B);
                        bool asc = ((k & 1) ^ flip) == 0;
                        sort(route.begin() + l, route.begin() + rr, [&](int a, int b) {
                            long long va = axis == 0 ? p[a].y : p[a].x;
                            long long vb = axis == 0 ? p[b].y : p[b].x;
                            long long ta = axis == 0 ? p[a].x : p[a].y;
                            long long tb = axis == 0 ? p[b].x : p[b].y;
                            if (va != vb) return asc ? va < vb : va > vb;
                            if (ta != tb) return asc ? ta < tb : ta > tb;
                            return asc ? a < b : a > b;
                        });
                    }
                    consider(route);
                }
            }
        }
    }

    // Cheap deterministic swap polishing of the selected tour.  Besides adjacent
    // swaps, try a short forward window; this often fixes local crossings and can
    // move prime city IDs onto carrot-sensitive source positions, but every move
    // is accepted only by exact penalized edge delta.
    auto edgeCost = [&](int step, int a, int b) -> double {
        double m = (step % 10 == 0 && !prime[a]) ? 1.1 : 1.0;
        return m * d2d(p[a], p[b]);
    };
    if (!best.empty() && N >= 4) {
        auto swapDelta = [&](int i, int j) -> double {
            int steps[4] = {i, i + 1, j, j + 1};
            double oldc = 0.0, newc = 0.0;
            auto val = [&](int pos) -> int {
                if (pos == i) return best[j];
                if (pos == j) return best[i];
                return best[pos];
            };
            for (int a = 0; a < 4; ++a) {
                int s = steps[a];
                if (s < 1 || s > N) continue;
                bool seen = false;
                for (int b = 0; b < a; ++b) if (steps[b] == s) seen = true;
                if (seen) continue;
                oldc += edgeCost(s, best[s-1], best[s]);
                newc += edgeCost(s, val(s-1), val(s));
            }
            return newc - oldc;
        };

        const int W = (N <= 5000 ? 32 : (N <= 50000 ? 16 : 8));
        bool changed = true;
        for (int pass = 0; pass < 2 && changed; ++pass) {
            changed = false;
            for (int i = 1; i < N; ++i) {
                int bestJ = -1;
                double bestD = -1e-9;
                int lim = min(N - 1, i + W);
                for (int j = i + 1; j <= lim; ++j) {
                    double d = swapDelta(i, j);
                    if (d < bestD) { bestD = d; bestJ = j; }
                }
                if (bestJ != -1) {
                    swap(best[i], best[bestJ]);
                    changed = true;
                }
            }
        }

        // Wider carrot-aware polishing: only positions that source a 10th step
        // have multiplier discontinuities, so spend extra comparisons there.
        // This can pull nearby prime IDs onto expensive steps (or reject it) by
        // the exact same local delta while keeping runtime linear.
        const int CW = (N <= 5000 ? 220 : (N <= 50000 ? 128 : 48));
        const int CP = (N <= 50000 ? 2 : 1);
        for (int pass = 0; pass < CP; ++pass) {
            bool any = false;
            for (int t = 10; t <= N; t += 10) {
                int i = t - 1;
                if (i <= 0 || i >= N) continue;
                int bestJ = -1;
                double bestD = -1e-9;
                int lo = max(1, i - CW), hi = min(N - 1, i + CW);
                for (int j = lo; j <= hi; ++j) if (j != i) {
                    // Prefer testing swaps that alter primality at the penalized
                    // source, but allow any exact improvement in the same window.
                    double d = swapDelta(min(i, j), max(i, j));
                    if (d < bestD) { bestD = d; bestJ = j; }
                }
                if (bestJ != -1) {
                    swap(best[i], best[bestJ]);
                    any = true;
                }
            }
            if (!any) break;
        }
    }

    cout << N + 1 << '\n';
    for (int v : best) cout << v << '\n';
    return 0;
}
