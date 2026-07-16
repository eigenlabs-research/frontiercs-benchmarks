#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x, y; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    if (!(cin >> N)) return 0;
    vector<Pt> p(N);
    for (int i = 0; i < N; ++i) cin >> p[i].x >> p[i].y;

    vector<char> prime(max(2, N), true);
    prime[0] = false;
    if (N > 1) prime[1] = false;
    for (long long i = 2; i * i < N; ++i) if (prime[(int)i])
        for (long long j = i * i; j < N; j += i) prime[(int)j] = false;

    auto dist = [&](int a, int b) -> double {
        double dx = (double)p[a].x - (double)p[b].x;
        double dy = (double)p[a].y - (double)p[b].y;
        return hypot(dx, dy);
    };
    auto stepCost = [&](const vector<int>& r, int t) -> double {
        double m = (t % 10 == 0 && !prime[r[t-1]]) ? 1.1 : 1.0;
        return m * dist(r[t-1], r[t]);
    };
    auto totalCost = [&](const vector<int>& r) -> double {
        double s = 0.0;
        for (int t = 1; t <= N; ++t) s += stepCost(r, t);
        return s;
    };

    vector<int> best(N + 1);
    best[0] = 0; best[N] = 0;
    for (int i = 1; i < N; ++i) best[i] = i;
    double bestCost = totalCost(best);

    auto consider = [&](vector<int>& ord) {
        vector<int> r(N + 1);
        r[0] = 0; r[N] = 0;
        for (int i = 1; i < N; ++i) r[i] = ord[i-1];
        double c = totalCost(r);
        if (c < bestCost) { bestCost = c; best.swap(r); }
    };

    vector<int> ids;
    ids.reserve(max(0, N-1));
    for (int i = 1; i < N; ++i) ids.push_back(i);

    // Pure y-order, useful when input alternates between distant horizontal rows.
    vector<int> yord = ids;
    stable_sort(yord.begin(), yord.end(), [&](int a, int b){
        if (p[a].y != p[b].y) return p[a].y < p[b].y;
        return p[a].x < p[b].x;
    });
    consider(yord);
    reverse(yord.begin(), yord.end());
    consider(yord);
    reverse(yord.begin(), yord.end());

    vector<int> blockSizes;
    for (int b : {2,4,8,16,32,64,128,256,512,1024,2048,4096,8192})
        if (b < N) blockSizes.push_back(b);
    int sq = max(2, (int)(sqrt((double)max(1, N))));
    blockSizes.push_back(sq);
    blockSizes.push_back(max(2, sq * 2));
    sort(blockSizes.begin(), blockSizes.end());
    blockSizes.erase(unique(blockSizes.begin(), blockSizes.end()), blockSizes.end());

    // Horizontal bands: sort by y into bands, traverse each band by x, alternating direction.
    for (int B : blockSizes) {
        int M = (int)yord.size();
        int nb = (M + B - 1) / B;
        vector<vector<int>> bands(nb);
        for (int bi = 0; bi < nb; ++bi) {
            int l = bi * B, r = min(M, l + B);
            bands[bi].assign(yord.begin() + l, yord.begin() + r);
            sort(bands[bi].begin(), bands[bi].end()); // city id order is x order
        }
        for (int revBands = 0; revBands < 2; ++revBands) for (int startRev = 0; startRev < 2; ++startRev) {
            vector<int> ord; ord.reserve(M);
            for (int q = 0; q < nb; ++q) {
                int bi = revBands ? (nb - 1 - q) : q;
                bool rev = ((q & 1) ^ startRev);
                if (!rev) ord.insert(ord.end(), bands[bi].begin(), bands[bi].end());
                else ord.insert(ord.end(), bands[bi].rbegin(), bands[bi].rend());
            }
            consider(ord);
        }
    }

    // Vertical strips: keep x-bands, traverse each strip by y, alternating direction.
    for (int B : blockSizes) {
        int M = (int)ids.size();
        int nb = (M + B - 1) / B;
        vector<vector<int>> strips(nb);
        for (int bi = 0; bi < nb; ++bi) {
            int l = bi * B, r = min(M, l + B);
            strips[bi].assign(ids.begin() + l, ids.begin() + r);
            sort(strips[bi].begin(), strips[bi].end(), [&](int a, int b){
                if (p[a].y != p[b].y) return p[a].y < p[b].y;
                return a < b;
            });
        }
        for (int revStrips = 0; revStrips < 2; ++revStrips) for (int startRev = 0; startRev < 2; ++startRev) {
            vector<int> ord; ord.reserve(M);
            for (int q = 0; q < nb; ++q) {
                int bi = revStrips ? (nb - 1 - q) : q;
                bool rev = ((q & 1) ^ startRev);
                if (!rev) ord.insert(ord.end(), strips[bi].begin(), strips[bi].end());
                else ord.insert(ord.end(), strips[bi].rbegin(), strips[bi].rend());
            }
            consider(ord);
        }
    }

    auto swapDelta = [&](const vector<int>& r, int i, int j) -> double {
        if (i == j) return 0.0;
        vector<int> ts = {i, i+1, j, j+1};
        sort(ts.begin(), ts.end());
        ts.erase(unique(ts.begin(), ts.end()), ts.end());
        double before = 0.0, after = 0.0;
        for (int t : ts) if (1 <= t && t <= N) before += stepCost(r, t);
        auto val = [&](int idx) -> int {
            if (idx == i) return r[j];
            if (idx == j) return r[i];
            return r[idx];
        };
        for (int t : ts) if (1 <= t && t <= N) {
            int a = val(t-1), b = val(t);
            double m = (t % 10 == 0 && !prime[a]) ? 1.1 : 1.0;
            after += m * dist(a, b);
        }
        return after - before;
    };

    // A small deterministic local polish: adjacent swaps, then swaps that try to put primes
    // at every 10th step source city when it is locally worthwhile.
    for (int pass = 0; pass < 2; ++pass) {
        bool changed = false;
        for (int i = 1; i + 1 < N; ++i) {
            double d = swapDelta(best, i, i+1);
            if (d < -1e-7) { swap(best[i], best[i+1]); bestCost += d; changed = true; }
        }
        if (!changed) break;
    }

    auto isTargetPos = [&](int pos) { return pos >= 9 && pos < N && pos % 10 == 9; };
    const int R = 80;
    for (int pos = 9; pos < N; pos += 10) {
        if (prime[best[pos]]) continue;
        int bj = -1; double bd = 0.0;
        int lo = max(1, pos - R), hi = min(N - 1, pos + R);
        for (int j = lo; j <= hi; ++j) {
            if (j == pos || isTargetPos(j) || !prime[best[j]]) continue;
            double d = swapDelta(best, pos, j);
            if (d < bd) { bd = d; bj = j; }
        }
        if (bj != -1) { swap(best[pos], best[bj]); bestCost += bd; }
    }

    cout << N + 1 << '\n';
    for (int v : best) cout << v << '\n';
    return 0;
}
