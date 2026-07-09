// Traveling Santa with carrot (prime every-10th-step) constraint.
// Strategy: x-order seed + multi-start windowed 2-opt + prime-aware swaps under 2.4s budget.
#include <bits/stdc++.h>
using namespace std;

static inline double dist2D(long long xa, long long ya, long long xb, long long yb) {
    double dx = (double)xa - (double)xb;
    double dy = (double)ya - (double)yb;
    return hypot(dx, dy);
}

struct Solver {
    int N;
    vector<long long> x, y;
    vector<char> isPrime;
    chrono::steady_clock::time_point t0;
    double timeLimit = 2.35;

    bool timeUp() const {
        return chrono::duration<double>(chrono::steady_clock::now() - t0).count() > timeLimit;
    }

    double edgeCost(int a, int b, int step /*1..N*/) const {
        double m = 1.0;
        if (step % 10 == 0 && !isPrime[a]) m = 1.1;
        return m * dist2D(x[a], y[a], x[b], y[b]);
    }

    // route has N+1 entries, route[0]=route[N]=0
    double tourCost(const vector<int>& route) const {
        double total = 0.0;
        for (int t = 1; t <= N; t++) {
            total += edgeCost(route[t - 1], route[t], t);
        }
        return total;
    }

    // 2-opt reverse segment [i..j] inclusive (1-based indices into middle cities)
    // route indices 0..N; cities at 1..N-1 are internal.
    // Reverse path between positions L and R (0-based in route), 0 < L <= R < N
    void reverseSeg(vector<int>& route, int L, int R) {
        while (L < R) {
            swap(route[L], route[R]);
            ++L;
            --R;
        }
    }

    // Delta for 2-opt that reconnects (a-b) and (c-d) into (a-c) and (b-d) by reversing b..c
    // positions: i-1 -- i ... j -- j+1, reverse i..j
    double delta2opt(const vector<int>& route, int i, int j) const {
        // steps involving edges (i-1,i) at step i, (j,j+1) at step j+1
        int a = route[i - 1], b = route[i];
        int c = route[j], d = route[j + 1];
        // After reverse: a-c, ..., b-d; internal step multipliers for reversed segment shift
        // Exact recompute of affected range is safer for prime penalties.
        // Affected steps: i .. j+1 (edges). Also all steps between change their source cities.
        // So recompute cost of steps i..j+1 before/after.
        double before = 0.0, after = 0.0;
        // before
        for (int t = i; t <= j + 1; t++) {
            before += edgeCost(route[t - 1], route[t], t);
        }
        // after: reversed i..j
        // new sequence: ... a, route[j], route[j-1], ..., route[i], d ...
        vector<int> mid;
        mid.reserve(j - i + 1);
        for (int k = j; k >= i; --k) mid.push_back(route[k]);
        // steps i: a->mid[0]=route[j]
        after += edgeCost(a, mid[0], i);
        for (int t = 1; t < (int)mid.size(); t++) {
            after += edgeCost(mid[t - 1], mid[t], i + t);
        }
        // step j+1: mid.back()=route[i] -> d
        after += edgeCost(mid.back(), d, j + 1);
        return after - before;
    }

    void twoOptWindow(vector<int>& route, int W, int passes) {
        int n = N; // last index of route is N (city 0)
        for (int pass = 0; pass < passes && !timeUp(); ++pass) {
            bool improved = false;
            for (int i = 1; i < n && !timeUp(); ++i) {
                int jmax = min(n - 1, i + W);
                for (int j = i + 1; j <= jmax; ++j) {
                    // skip tiny
                    if (j == i) continue;
                    double dlt = delta2opt(route, i, j);
                    if (dlt < -1e-9) {
                        reverseSeg(route, i, j);
                        improved = true;
                    }
                }
            }
            if (!improved) break;
        }
    }

    // Random neighborhood 2-opt: sample pairs
    void twoOptRandom(vector<int>& route, int trials) {
        int n = N;
        mt19937 rng((uint32_t)chrono::steady_clock::now().time_since_epoch().count() ^ 0xC0FFEE);
        uniform_int_distribution<int> distI(1, max(1, n - 2));
        for (int t = 0; t < trials && !timeUp(); ++t) {
            int i = distI(rng);
            int j = i + 1 + (rng() % min(80, n - 1 - i));
            if (j >= n) continue;
            double dlt = delta2opt(route, i, j);
            if (dlt < -1e-9) {
                reverseSeg(route, i, j);
            }
        }
    }

    // Try place primes at penalty source positions P[9], P[19], ...
    void primeAlign(vector<int>& route) {
        // positions pos = 9,19,... if pos < N
        vector<int> need;
        for (int pos = 9; pos < N; pos += 10) need.push_back(pos);
        if (need.empty()) return;

        // map city -> position in route
        vector<int> posOf(N, -1);
        for (int i = 0; i <= N; i++) posOf[route[i]] = i;

        // list primes among 2..N-1
        vector<int> primes;
        for (int c = 2; c < N; c++) if (isPrime[c]) primes.push_back(c);

        // For each needed pos, if current not prime, try swap with a nearby prime
        for (int pos : need) {
            if (timeUp()) break;
            int cur = route[pos];
            if (isPrime[cur]) continue;
            // search primes whose current position is close in index or space
            double bestDelta = 0.0;
            int bestP = -1;
            for (int pr : primes) {
                int ppos = posOf[pr];
                if (ppos <= 0 || ppos >= N) continue; // don't move fixed 0 ends
                if (isPrime[route[pos]]) break;
                // swap pos and ppos
                // evaluate cost delta by recomputing affected edges
                auto costAround = [&](int a, int b) {
                    // positions a,b: recompute steps involving neighbors of both
                    double s = 0;
                    auto addPos = [&](int p) {
                        if (p > 0) s += edgeCost(route[p - 1], route[p], p);
                        if (p + 1 <= N) s += edgeCost(route[p], route[p + 1], p + 1);
                    };
                    addPos(a);
                    if (b != a && b != a - 1 && b != a + 1) addPos(b);
                    else if (b != a) {
                        // overlapping; just full local band
                        int lo = min(a, b) - 1, hi = max(a, b) + 1;
                        lo = max(lo, 1);
                        hi = min(hi, N);
                        s = 0;
                        for (int t = lo; t <= hi; t++) s += edgeCost(route[t - 1], route[t], t);
                    }
                    return s;
                };
                // simpler: swap and recompute full? too slow. local band:
                int lo = max(1, min(pos, ppos) - 1);
                int hi = min(N, max(pos, ppos) + 1);
                double before = 0;
                for (int t = lo; t <= hi; t++) before += edgeCost(route[t - 1], route[t], t);
                swap(route[pos], route[ppos]);
                double after = 0;
                for (int t = lo; t <= hi; t++) after += edgeCost(route[t - 1], route[t], t);
                double dlt = after - before;
                if (dlt < bestDelta - 1e-12) {
                    bestDelta = dlt;
                    bestP = ppos;
                }
                // revert
                swap(route[pos], route[ppos]);
            }
            if (bestP >= 0) {
                swap(route[pos], route[bestP]);
                // fix posOf
                posOf[route[pos]] = pos;
                posOf[route[bestP]] = bestP;
            }
        }
    }

    // Or-opt / relocate single city to better slot (limited)
    void relocate(vector<int>& route, int trials) {
        mt19937 rng(1234567u ^ (uint32_t)N);
        if (N < 4) return;
        for (int t = 0; t < trials && !timeUp(); ++t) {
            int i = 1 + (rng() % (N - 1)); // city position to move
            int j = 1 + (rng() % (N - 1)); // insert before j
            if (i == j || i + 1 == j) continue;
            // remove i, insert before j
            int city = route[i];
            vector<int> nr;
            nr.reserve(N + 1);
            for (int k = 0; k <= N; k++) if (k != i) nr.push_back(route[k]);
            // after removal, adjust j
            int jj = j;
            if (j > i) jj--;
            // insert city at position jj
            nr.insert(nr.begin() + jj, city);
            if ((int)nr.size() != N + 1) continue;
            if (nr.front() != 0 || nr.back() != 0) continue;
            double c0 = tourCost(route);
            double c1 = tourCost(nr);
            if (c1 < c0 - 1e-9) route.swap(nr);
        }
    }

    vector<int> solve() {
        t0 = chrono::steady_clock::now();
        // seed: natural order 0,1,2,...,N-1,0 (x-sorted baseline)
        vector<int> best(N + 1);
        for (int i = 0; i < N; i++) best[i] = i;
        best[N] = 0;
        double bestCost = tourCost(best);

        // y-sorted alternate seed (among cities 1..N-1)
        {
            vector<int> mid(N - 1);
            iota(mid.begin(), mid.end(), 1);
            sort(mid.begin(), mid.end(), [&](int a, int b) {
                if (y[a] != y[b]) return y[a] < y[b];
                return x[a] < x[b];
            });
            vector<int> r(N + 1);
            r[0] = 0;
            for (int i = 0; i < N - 1; i++) r[i + 1] = mid[i];
            r[N] = 0;
            twoOptWindow(r, 40, 2);
            primeAlign(r);
            double c = tourCost(r);
            if (c < bestCost) {
                bestCost = c;
                best.swap(r);
            }
        }

        // zigzag by x-bands (even/odd strips)
        if (N >= 20 && !timeUp()) {
            vector<int> mid(N - 1);
            iota(mid.begin(), mid.end(), 1);
            // already x-sorted ids ≈ x order
            // reverse alternate blocks of size B
            for (int B : {8, 16, 32, 64}) {
                if (timeUp()) break;
                vector<int> r(N + 1);
                r[0] = 0;
                int idx = 1;
                bool rev = false;
                for (int start = 1; start < N; start += B) {
                    int end = min(N - 1, start + B - 1);
                    if (!rev) {
                        for (int c = start; c <= end; c++) r[idx++] = c;
                    } else {
                        for (int c = end; c >= start; c--) r[idx++] = c;
                    }
                    rev = !rev;
                }
                r[N] = 0;
                twoOptWindow(r, 60, 3);
                twoOptRandom(r, 5000);
                primeAlign(r);
                twoOptRandom(r, 3000);
                double c = tourCost(r);
                if (c < bestCost) {
                    bestCost = c;
                    best.swap(r);
                }
            }
        }

        // Improve best
        twoOptWindow(best, 100, 5);
        twoOptRandom(best, 20000);
        primeAlign(best);
        twoOptRandom(best, 10000);
        relocate(best, 2000);
        twoOptWindow(best, 50, 3);

        // Multi-start: random shuffles of small blocks then polish
        mt19937 rng(42u);
        int restarts = 0;
        while (!timeUp() && restarts < 40) {
            ++restarts;
            vector<int> r = best;
            // scramble a random window
            int L = 1 + (int)(rng() % max(1, N - 10));
            int R = min(N - 1, L + 5 + (int)(rng() % 40));
            shuffle(r.begin() + L, r.begin() + R + 1, rng);
            twoOptWindow(r, 80, 2);
            twoOptRandom(r, 4000);
            primeAlign(r);
            double c = tourCost(r);
            if (c < bestCost - 1e-9) {
                bestCost = c;
                best.swap(r);
            }
        }

        // final polish
        twoOptRandom(best, 15000);
        primeAlign(best);
        twoOptWindow(best, 120, 4);
        return best;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int N;
    if (!(cin >> N)) return 0;
    Solver S;
    S.N = N;
    S.x.resize(N);
    S.y.resize(N);
    for (int i = 0; i < N; i++) cin >> S.x[i] >> S.y[i];
    S.isPrime.assign(max(2, N), true);
    S.isPrime[0] = false;
    if (N > 1) S.isPrime[1] = false;
    for (int i = 2; 1LL * i * i < N; i++) {
        if (S.isPrime[i]) {
            for (int j = i * i; j < N; j += i) S.isPrime[j] = false;
        }
    }
    auto route = S.solve();
    cout << (N + 1) << "\n";
    for (int v : route) cout << v << "\n";
    return 0;
}
