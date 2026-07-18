#include <bits/stdc++.h>
using namespace std;

struct Solver {
    int J, M, N;
    vector<vector<int>> mach;
    vector<vector<long long>> dur;
    vector<vector<int>> posOfMach;
    vector<vector<long long>> rem; // remaining work including operation k
    vector<long long> machineLoad, jobLoad;
    mt19937_64 rng;
    chrono::steady_clock::time_point startTime;
    double timeLimit = 0.92;

    bool timeUp() const {
        return chrono::duration<double>(chrono::steady_clock::now() - startTime).count() > timeLimit;
    }

    long long eval(const vector<vector<int>>& ord) const {
        vector<int> indeg(N, 0);
        vector<vector<int>> adj(N);
        adj.reserve(N);
        auto id = [&](int j, int k){ return j * M + k; };
        for (int j = 0; j < J; ++j) {
            for (int k = 0; k + 1 < M; ++k) {
                adj[id(j,k)].push_back(id(j,k+1));
                indeg[id(j,k+1)]++;
            }
        }
        for (int m = 0; m < M; ++m) {
            for (int t = 0; t + 1 < J; ++t) {
                int a = ord[m][t], b = ord[m][t+1];
                int ka = posOfMach[a][m], kb = posOfMach[b][m];
                adj[id(a,ka)].push_back(id(b,kb));
                indeg[id(b,kb)]++;
            }
        }
        deque<int> q;
        vector<long long> st(N, 0);
        for (int i = 0; i < N; ++i) if (!indeg[i]) q.push_back(i);
        int seen = 0;
        long long best = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop_front(); seen++;
            int j = u / M, k = u % M;
            long long fin = st[u] + dur[j][k];
            if (fin > best) best = fin;
            for (int v : adj[u]) {
                if (st[v] < fin) st[v] = fin;
                if (--indeg[v] == 0) q.push_back(v);
            }
        }
        if (seen != N) return LLONG_MAX / 4;
        return best;
    }

    struct Result { long long ms; vector<vector<int>> ord; };

    double priorityScore(int rule, int j, int k, long long s, long long c,
                         const array<double,6>& w) {
        long long p = dur[j][k];
        long long r = rem[j][k];
        long long tail = r - p;
        switch (rule) {
            case 0: return (double)p;                         // SPT
            case 1: return -(double)p;                        // LPT
            case 2: return -(double)r;                        // most work remaining
            case 3: return (double)r;                         // least work remaining
            case 4: return (double)c;                         // earliest completion
            case 5: return (double)s * 4.0 - (double)tail;    // early, long tail
            case 6: return -(double)tail;                     // largest downstream tail
            case 7: return (double)tail;                      // smallest downstream tail
            case 8: return -(double)jobLoad[j];
            case 9: return (double)p - 0.35 * (double)tail;
            case 11: return (double)s - 0.55 * (double)tail - 0.10 * (double)p; // global long-tail urgency
            case 12: return (double)c - 0.75 * (double)tail + 0.15 * (double)p; // earliest finish with tail bias
            case 13: return -(double)machineLoad[mach[j][k]] - 0.25 * (double)tail; // bottleneck first
            default:
                return w[0]*(double)s + w[1]*(double)c + w[2]*(double)p
                     + w[3]*(double)r + w[4]*(double)tail + w[5]*(double)jobLoad[j];
        }
    }

    Result serialDispatch(int rule, const array<double,6>& w, bool randomTie) {
        // Serial schedule generation: choose one globally best next operation at a time.
        // This is a complementary representation to GT's conflict-set choice and is
        // especially useful when one bottleneck machine should dominate dispatching.
        vector<vector<int>> ord(M);
        for (int m = 0; m < M; ++m) ord[m].reserve(J);
        vector<int> nxt(J, 0);
        vector<long long> jr(J, 0), mr(M, 0);
        uniform_real_distribution<double> tiny(0.0, 1e-7);
        for (int done = 0; done < N; ++done) {
            int choose = -1;
            double bestScore = 1e300;
            for (int j = 0; j < J; ++j) if (nxt[j] < M) {
                int k = nxt[j], m = mach[j][k];
                long long s = max(jr[j], mr[m]);
                long long c = s + dur[j][k];
                double sc = priorityScore(rule, j, k, s, c, w);
                // Keep serial schedules from idling too aggressively unless the
                // rule strongly asks for it.
                sc += 1e-6 * (double)s;
                if (randomTie) sc += tiny(rng);
                else sc += j * 1e-9;
                if (sc < bestScore) { bestScore = sc; choose = j; }
            }
            int j = choose, k = nxt[j], m = mach[j][k];
            long long f = max(jr[j], mr[m]) + dur[j][k];
            ord[m].push_back(j);
            jr[j] = f; mr[m] = f; nxt[j]++;
        }
        return {eval(ord), move(ord)};
    }

    Result gifflerThompson(int rule, const array<double,6>& w, bool randomTie) {
        vector<vector<int>> ord(M);
        for (int m = 0; m < M; ++m) ord[m].reserve(J);
        vector<int> nxt(J, 0);
        vector<long long> jr(J, 0), mr(M, 0);
        uniform_real_distribution<double> tiny(0.0, 1e-7);
        for (int done = 0; done < N; ++done) {
            long long bestC = LLONG_MAX;
            int bestM = -1;
            for (int j = 0; j < J; ++j) if (nxt[j] < M) {
                int k = nxt[j], m = mach[j][k];
                long long s = max(jr[j], mr[m]);
                long long c = s + dur[j][k];
                if (c < bestC) { bestC = c; bestM = m; }
            }
            int choose = -1;
            double bestScore = 1e300;
            for (int j = 0; j < J; ++j) if (nxt[j] < M && mach[j][nxt[j]] == bestM) {
                int k = nxt[j];
                long long s = max(jr[j], mr[bestM]);
                if (s < bestC) {
                    long long c = s + dur[j][k];
                    double sc = priorityScore(rule, j, k, s, c, w);
                    if (randomTie) sc += tiny(rng);
                    else sc += j * 1e-9;
                    if (sc < bestScore) { bestScore = sc; choose = j; }
                }
            }
            if (choose < 0) { // Should not happen; fallback to any available operation on bestM.
                for (int j = 0; j < J; ++j) if (nxt[j] < M && mach[j][nxt[j]] == bestM) { choose = j; break; }
            }
            int j = choose, k = nxt[j], m = mach[j][k];
            long long s = max(jr[j], mr[m]);
            long long f = s + dur[j][k];
            ord[m].push_back(j);
            jr[j] = f; mr[m] = f; nxt[j]++;
        }
        return {eval(ord), move(ord)};
    }

    void tryAdjacentImprove(Result &best) {
        // A small, safe first-improvement pass over adjacent machine-order swaps.
        bool improved = true;
        while (improved && !timeUp()) {
            improved = false;
            for (int m = 0; m < M && !timeUp(); ++m) {
                for (int i = 0; i + 1 < J && !timeUp(); ++i) {
                    swap(best.ord[m][i], best.ord[m][i+1]);
                    long long v = eval(best.ord);
                    if (v < best.ms) {
                        best.ms = v;
                        improved = true;
                    } else {
                        swap(best.ord[m][i], best.ord[m][i+1]);
                    }
                }
            }
        }
    }

    void solve() {
        cin >> J >> M;
        N = J * M;
        mach.assign(J, vector<int>(M));
        dur.assign(J, vector<long long>(M));
        posOfMach.assign(J, vector<int>(M));
        machineLoad.assign(M, 0);
        jobLoad.assign(J, 0);
        for (int j = 0; j < J; ++j) {
            for (int k = 0; k < M; ++k) {
                cin >> mach[j][k] >> dur[j][k];
                posOfMach[j][mach[j][k]] = k;
                machineLoad[mach[j][k]] += dur[j][k];
                jobLoad[j] += dur[j][k];
            }
        }
        rem.assign(J, vector<long long>(M+1, 0));
        for (int j = 0; j < J; ++j)
            for (int k = M-1; k >= 0; --k) rem[j][k] = rem[j][k+1] + dur[j][k];

        startTime = chrono::steady_clock::now();
        rng.seed(chrono::high_resolution_clock::now().time_since_epoch().count() ^ (uint64_t)J<<32 ^ M);

        Result best{LLONG_MAX/4, {}};
        array<double,6> zero{0,0,0,0,0,0};
        for (int r = 0; r < 10 && !timeUp(); ++r) {
            Result res = gifflerThompson(r, zero, false);
            if (res.ms < best.ms) best = move(res);
        }
        for (int r : {2,4,5,6,9,11,12,13}) if (!timeUp()) {
            Result res = serialDispatch(r, zero, false);
            if (res.ms < best.ms) best = move(res);
        }

        uniform_real_distribution<double> coef(-1.0, 1.0);
        int iter = 0;
        while (!timeUp()) {
            array<double,6> w;
            // Scaled random linear priorities; signs are allowed to choose either direction.
            w[0] = coef(rng) * 0.6;   // start time
            w[1] = coef(rng) * 0.4;   // completion time
            w[2] = coef(rng) * 1.2;   // processing time
            w[3] = coef(rng) * 0.9;   // remaining work
            w[4] = coef(rng) * 0.9;   // downstream tail
            w[5] = coef(rng) * 0.3;   // total job size
            Result res = (iter % 3 == 2) ? serialDispatch(10, w, true) : gifflerThompson(10, w, true);
            if (res.ms < best.ms) best = move(res);
            if (++iter % 32 == 0 && chrono::duration<double>(chrono::steady_clock::now() - startTime).count() > 0.72) break;
        }
        tryAdjacentImprove(best);

        if (best.ord.empty()) {
            best.ord.assign(M, vector<int>());
            for (int m = 0; m < M; ++m) for (int j = 0; j < J; ++j) best.ord[m].push_back(j);
        }
        for (int m = 0; m < M; ++m) {
            for (int i = 0; i < J; ++i) {
                if (i) cout << ' ';
                cout << best.ord[m][i];
            }
            cout << '\n';
        }
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    Solver s;
    s.solve();
    return 0;
}
