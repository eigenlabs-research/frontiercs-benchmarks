#include <bits/stdc++.h>
using namespace std;

struct Solver {
    int J, M, N;
    vector<vector<int>> mach, opMach;
    vector<vector<long long>> dur, rem;
    vector<long long> mload;
    mt19937_64 rng;
    chrono::steady_clock::time_point deadline;

    struct Result {
        long long makespan = LLONG_MAX;
        vector<vector<int>> seq;
    };

    double rnd01() { return (rng() >> 11) * (1.0 / 9007199254740992.0); }

    long long evalExact(const vector<vector<int>>& seq) {
        const long long INF = LLONG_MAX/4;
        vector<int> indeg(N, 0), head(N, -1), to;
        to.reserve(2*N);
        auto addEdge = [&](int a, int b){ to.push_back(b); indeg[b]++; };
        for (int j=0;j<J;j++) for (int k=0;k+1<M;k++) addEdge(j*M+k, j*M+k+1);
        for (int m=0;m<M;m++) {
            if ((int)seq[m].size() != J) return INF;
            vector<int> seen(J,0);
            int prev = -1;
            for (int job: seq[m]) {
                if (job < 0 || job >= J || seen[job]++) return INF;
                int node = job*M + opMach[job][m];
                if (prev != -1) addEdge(prev, node);
                prev = node;
            }
        }
        vector<vector<int>> adj(N);
        adj.reserve(N);
        int e = 0;
        for (int j=0;j<J;j++) for (int k=0;k+1<M;k++,e++) adj[j*M+k].push_back(to[e]);
        for (int m=0;m<M;m++) for (int p=0;p+1<J;p++,e++) {
            int aJob = seq[m][p], bJob = seq[m][p+1];
            adj[aJob*M + opMach[aJob][m]].push_back(bJob*M + opMach[bJob][m]);
        }
        queue<int> q;
        vector<long long> dist(N, 0);
        for (int i=0;i<N;i++) if (!indeg[i]) q.push(i);
        int cnt = 0;
        long long best = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop(); cnt++;
            int j = u / M, k = u % M;
            long long fin = dist[u] + dur[j][k];
            best = max(best, fin);
            for (int v: adj[u]) {
                if (dist[v] < fin) dist[v] = fin;
                if (--indeg[v] == 0) q.push(v);
            }
        }
        return cnt == N ? best : INF;
    }

    long long evalCritical(const vector<vector<int>>& seq, vector<int>& pred, int& endNode) {
        const long long INF = LLONG_MAX/4;
        vector<int> indeg(N, 0);
        vector<vector<int>> adj(N);
        auto addEdge = [&](int a, int b){ adj[a].push_back(b); indeg[b]++; };
        for (int j=0;j<J;j++) for (int k=0;k+1<M;k++) addEdge(j*M+k, j*M+k+1);
        for (int m=0;m<M;m++) {
            if ((int)seq[m].size() != J) return INF;
            vector<int> seen(J,0);
            int prev = -1;
            for (int job: seq[m]) {
                if (job < 0 || job >= J || seen[job]++) return INF;
                int node = job*M + opMach[job][m];
                if (prev != -1) addEdge(prev, node);
                prev = node;
            }
        }
        queue<int> q;
        vector<long long> dist(N, 0);
        pred.assign(N, -1);
        for (int i=0;i<N;i++) if (!indeg[i]) q.push(i);
        int cnt = 0;
        long long bestSpan = 0;
        endNode = -1;
        while (!q.empty()) {
            int u = q.front(); q.pop(); cnt++;
            int j = u / M, k = u % M;
            long long fin = dist[u] + dur[j][k];
            if (fin > bestSpan) { bestSpan = fin; endNode = u; }
            for (int v: adj[u]) {
                if (dist[v] < fin) { dist[v] = fin; pred[v] = u; }
                if (--indeg[v] == 0) q.push(v);
            }
        }
        return cnt == N ? bestSpan : INF;
    }

    static void movePos(vector<int>& v, int from, int to) {
        if (from == to) return;
        int x = v[from];
        if (from < to) {
            for (int i=from;i<to;i++) v[i] = v[i+1];
            v[to] = x;
        } else {
            for (int i=from;i>to;i--) v[i] = v[i-1];
            v[to] = x;
        }
    }

    void localImprove(Result &best) {
        best.makespan = evalExact(best.seq);
        auto adjDeadline = deadline - chrono::milliseconds(120);
        bool improved = true;
        while (improved && chrono::steady_clock::now() < adjDeadline) {
            improved = false;
            long long passBest = best.makespan;
            int bm = -1, bp = -1;
            for (int m=0;m<M && chrono::steady_clock::now() < adjDeadline;m++) {
                for (int p=0;p+1<J;p++) {
                    swap(best.seq[m][p], best.seq[m][p+1]);
                    long long val = evalExact(best.seq);
                    swap(best.seq[m][p], best.seq[m][p+1]);
                    if (val < passBest) { passBest = val; bm = m; bp = p; }
                }
            }
            if (bm != -1) {
                swap(best.seq[bm][bp], best.seq[bm][bp+1]);
                best.makespan = passBest;
                improved = true;
            }
        }

        // Adjacent swaps are a small N5 neighborhood; after it stalls, try a
        // critical-path insertion neighborhood.  Only jobs bordering machine
        // arcs on one current critical path are moved, keeping the exact
        // acyclicity/makespan evaluator affordable under the 1s limit.
        while (chrono::steady_clock::now() < deadline) {
            vector<int> pred;
            int endNode = -1;
            long long cur = evalCritical(best.seq, pred, endNode);
            if (cur < best.makespan) best.makespan = cur;
            vector<vector<int>> pos(M, vector<int>(J));
            for (int m=0;m<M;m++) for (int p=0;p<J;p++) pos[m][best.seq[m][p]] = p;
            vector<tuple<int,int,int>> moves;
            for (int v=endNode; v!=-1 && pred[v]!=-1; v=pred[v]) {
                int u = pred[v];
                int ju = u / M, ku = u % M;
                int jv = v / M, kv = v % M;
                int mu = mach[ju][ku], mv = mach[jv][kv];
                if (mu == mv) {
                    int p = pos[mu][ju], q = pos[mu][jv];
                    if (abs(p-q) == 1) {
                        int lo = max(0, min(p,q)-4), hi = min(J-1, max(p,q)+4);
                        for (int t=lo;t<=hi;t++) {
                            if (t != p) moves.emplace_back(mu, p, t);
                            if (t != q) moves.emplace_back(mu, q, t);
                        }
                    }
                }
            }
            sort(moves.begin(), moves.end());
            moves.erase(unique(moves.begin(), moves.end()), moves.end());
            long long passBest = best.makespan;
            int bm=-1, bf=-1, bt=-1;
            for (auto [m,from,to]: moves) {
                if (chrono::steady_clock::now() >= deadline) break;
                movePos(best.seq[m], from, to);
                long long val = evalExact(best.seq);
                movePos(best.seq[m], to, from);
                if (val < passBest) { passBest = val; bm=m; bf=from; bt=to; }
            }
            if (bm == -1) break;
            movePos(best.seq[bm], bf, bt);
            best.makespan = passBest;
        }
    }

    Result build(const array<double,8>& w, bool gt) {
        vector<int> idx(J, 0);
        vector<long long> jr(J, 0), mr(M, 0);
        vector<vector<int>> seq(M);
        for (int m=0;m<M;m++) seq[m].reserve(J);
        long long cmax = 0;

        for (int done=0; done<N; ++done) {
            vector<int> cand;
            long long bestC = LLONG_MAX;
            int focusM = -1;
            if (gt) {
                for (int j=0;j<J;j++) if (idx[j] < M) {
                    int k = idx[j], m = mach[j][k];
                    long long st = max(jr[j], mr[m]);
                    long long ct = st + dur[j][k];
                    if (ct < bestC) { bestC = ct; focusM = m; }
                }
                for (int j=0;j<J;j++) if (idx[j] < M) {
                    int k = idx[j], m = mach[j][k];
                    long long st = max(jr[j], mr[m]);
                    if (m == focusM && st < bestC) cand.push_back(j);
                }
            } else {
                for (int j=0;j<J;j++) if (idx[j] < M) cand.push_back(j);
            }

            int chosen = cand[0];
            double bestScore = -1e300;
            for (int j: cand) {
                int k = idx[j], m = mach[j][k];
                long long p = dur[j][k];
                long long st = max(jr[j], mr[m]);
                long long waitJ = max(0LL, mr[m] - jr[j]);
                long long waitM = max(0LL, jr[j] - mr[m]);
                double slackToLB = (double)st;
                // Higher score is selected.  Features deliberately mix classic dispatching
                // priorities (MWKR/LWKR, SPT/LPT, bottleneck focus, earliest start) through
                // one active-schedule generator.
                double sc = 0.0;
                sc += w[0] * (double)rem[j][k];             // remaining job work
                sc += w[1] * (double)p;                     // current processing time
                sc += w[2] * (double)mload[m];              // global bottleneck weight
                sc += w[3] * (double)(M - k);               // remaining operation count
                sc -= w[4] * (double)st;                    // earliest start preference
                sc -= w[5] * (double)waitJ;                 // avoid jobs waiting for machine
                sc -= w[6] * (double)waitM;                 // avoid machines waiting for job
                sc += w[7] * (double)(rem[j][k] - p);       // tail after this operation
                sc += 1e-7 * rnd01();
                if (sc > bestScore) { bestScore = sc; chosen = j; }
            }

            int k = idx[chosen], m = mach[chosen][k];
            long long st = max(jr[chosen], mr[m]);
            long long ct = st + dur[chosen][k];
            jr[chosen] = mr[m] = ct;
            cmax = max(cmax, ct);
            seq[m].push_back(chosen);
            idx[chosen]++;
        }
        return {cmax, seq};
    }

    void run() {
        cin >> J >> M;
        N = J * M;
        mach.assign(J, vector<int>(M));
        opMach.assign(J, vector<int>(M));
        dur.assign(J, vector<long long>(M));
        mload.assign(M, 0);
        uint64_t seed = 1469598103934665603ULL;
        for (int j=0;j<J;j++) for (int k=0;k<M;k++) {
            cin >> mach[j][k] >> dur[j][k];
            opMach[j][mach[j][k]] = k;
            mload[mach[j][k]] += dur[j][k];
            seed ^= (uint64_t)(mach[j][k] + 1009) + 0x9e3779b97f4a7c15ULL + (seed<<6) + (seed>>2);
            seed ^= (uint64_t)(dur[j][k] + 9176) + 0x9e3779b97f4a7c15ULL + (seed<<6) + (seed>>2);
        }
        rng.seed(seed);
        rem.assign(J, vector<long long>(M+1, 0));
        for (int j=0;j<J;j++) for (int k=M-1;k>=0;k--) rem[j][k] = rem[j][k+1] + dur[j][k];

        auto startTime = chrono::steady_clock::now();
        deadline = startTime + chrono::milliseconds(740);
        vector<array<double,8>> rules = {
            array<double,8>{ 1, 0, 0, 0, 0, 0, 0, 0},      // MWKR
            array<double,8>{-1, 0, 0, 0, 0, 0, 0, 0},      // LWKR variant
            array<double,8>{ 0,-1, 0, 0, 0, 0, 0, 0},      // SPT
            array<double,8>{ 0, 1, 0, 0, 0, 0, 0, 0},      // LPT
            array<double,8>{ 1,-2, 0, 0, 0, 0, 0, 1},
            array<double,8>{ 1, 0, 1, 0, 1, 0, 0, 0},
            array<double,8>{ 1, 0, 2, 0, 0, 1, 0, 1},
            array<double,8>{ 0, 0, 1, 1, 1, 0, 0, 0},
            array<double,8>{ 2,-1, 1, 0, 1, 1, 0, 1},
            array<double,8>{ 2, 1, 1, 0, 1, 0, 1, 2},
            array<double,8>{ 0,-1, 2, 1, 0, 1, 0, 0},
            array<double,8>{ 1, 1, 0, 1, 1, 1, 1, 0}
        };

        Result best;
        // Deterministic, robust slate first.
        for (auto &r: rules) {
            Result a = build(r, true);
            if (a.makespan < best.makespan) best = std::move(a);
            Result b = build(r, false);
            if (b.makespan < best.makespan) best = std::move(b);
        }

        // Time-bounded randomized weight search around the same dispatching mechanism.
        int iter = 0;
        while (chrono::steady_clock::now() < deadline) {
            array<double,8> w{};
            // Bias toward known useful signs but allow diversity across instance families.
            w[0] = (rnd01()*4.0 - 1.0);       // rem usually positive
            w[1] = (rnd01()*4.0 - 2.0);       // p may be SPT or LPT
            w[2] = (rnd01()*3.0);             // bottleneck nonnegative
            w[3] = (rnd01()*2.0 - 0.5);
            w[4] = (rnd01()*2.0);             // prefer earlier starts
            w[5] = (rnd01()*2.0 - 0.5);
            w[6] = (rnd01()*2.0 - 0.5);
            w[7] = (rnd01()*4.0 - 1.0);
            Result r = build(w, (iter & 3) != 0);
            if (r.makespan < best.makespan) best = std::move(r);
            ++iter;
        }

        // Targeted improvement over the dispatch portfolio: polish the selected
        // machine permutations by exact adjacent-swap hill climbing.  Invalid
        // swaps are rejected by the topological evaluator.
        deadline = startTime + chrono::milliseconds(965);
        localImprove(best);

        for (int m=0;m<M;m++) {
            // Defensive completion should never be needed for schedules generated above.
            if ((int)best.seq[m].size() != J) {
                vector<int> seen(J,0);
                for (int x: best.seq[m]) if (0 <= x && x < J) seen[x]=1;
                for (int j=0;j<J;j++) if (!seen[j]) best.seq[m].push_back(j);
                best.seq[m].resize(J);
            }
            for (int i=0;i<J;i++) {
                if (i) cout << ' ';
                cout << best.seq[m][i];
            }
            cout << '\n';
        }
    }
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    Solver s;
    if (!(cin.peek(), cin.good())) return 0;
    s.run();
    return 0;
}
