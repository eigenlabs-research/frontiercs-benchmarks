#include <bits/stdc++.h>
using namespace std;

struct Solver {
    int J, M, N;
    vector<vector<int>> mach;
    vector<vector<long long>> p;
    vector<vector<long long>> rem;
    vector<long long> mload, jload;
    mt19937_64 rng;

    struct Result {
        long long makespan;
        vector<vector<int>> order;
    };

    uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    double urand() {
        return (rng() >> 11) * (1.0 / 9007199254740992.0);
    }

    double score_candidate(int rule, int j, int k, long long est, const vector<double>& w, double noise) {
        long long pr = p[j][k];
        long long r = rem[j][k];
        int opsLeft = M - k;
        switch (rule) {
            case 0: return (double)pr;                              // SPT among conflicts
            case 1: return -(double)pr;                             // LPT
            case 2: return -(double)r;                              // most work remaining
            case 3: return (double)r;                               // least work remaining
            case 4: return (double)est + 0.15 * pr;                 // earliest start, short tie
            case 5: return -(double)(r - pr) + 0.05 * pr;           // largest tail after this op
            case 6: return -(double)opsLeft * 1000000.0 - (double)r * 1e-3;
            case 7: return (double)opsLeft * 1000000.0 - (double)r * 1e-3;
            default:
                return w[0] * (double)pr + w[1] * (double)r +
                       w[2] * (double)est + w[3] * (double)(r - pr) +
                       w[4] * (double)opsLeft + noise * (urand() - 0.5);
        }
    }

    Result build(int rule, vector<double> w = {}, double noise = 0.0) {
        vector<int> nextOp(J, 0);
        vector<long long> jobReady(J, 0), machineReady(M, 0);
        vector<vector<int>> order(M);
        for (int m = 0; m < M; ++m) order[m].reserve(J);
        long long makespan = 0;

        for (int done = 0; done < N; ++done) {
            long long bestC = LLONG_MAX;
            int chosenMachine = -1;
            // Giffler-Thompson: find the operation with minimum earliest completion.
            for (int j = 0; j < J; ++j) if (nextOp[j] < M) {
                int k = nextOp[j], m = mach[j][k];
                long long est = max(jobReady[j], machineReady[m]);
                long long c = est + p[j][k];
                if (c < bestC || (c == bestC && m < chosenMachine)) {
                    bestC = c;
                    chosenMachine = m;
                }
            }

            int bestJ = -1;
            double bestS = numeric_limits<double>::infinity();
            // Resolve only operations competing for the critical machine before bestC.
            for (int j = 0; j < J; ++j) if (nextOp[j] < M) {
                int k = nextOp[j], m = mach[j][k];
                if (m != chosenMachine) continue;
                long long est = max(jobReady[j], machineReady[m]);
                if (est >= bestC) continue;
                double s = score_candidate(rule, j, k, est, w, noise);
                if (s < bestS || (s == bestS && j < bestJ)) {
                    bestS = s;
                    bestJ = j;
                }
            }
            if (bestJ < 0) { // Should not happen, but keep output valid.
                for (int j = 0; j < J; ++j) if (nextOp[j] < M && mach[j][nextOp[j]] == chosenMachine) { bestJ = j; break; }
            }

            int k = nextOp[bestJ], m = mach[bestJ][k];
            long long st = max(jobReady[bestJ], machineReady[m]);
            long long ft = st + p[bestJ][k];
            jobReady[bestJ] = ft;
            machineReady[m] = ft;
            makespan = max(makespan, ft);
            order[m].push_back(bestJ);
            nextOp[bestJ]++;
        }
        return {makespan, order};
    }

    void solve() {
        cin >> J >> M;
        N = J * M;
        mach.assign(J, vector<int>(M));
        p.assign(J, vector<long long>(M));
        rem.assign(J, vector<long long>(M + 1, 0));
        mload.assign(M, 0);
        jload.assign(J, 0);
        uint64_t h = splitmix64((uint64_t)J * 1000003ULL + M);
        for (int j = 0; j < J; ++j) {
            for (int k = 0; k < M; ++k) {
                cin >> mach[j][k] >> p[j][k];
                mload[mach[j][k]] += p[j][k];
                jload[j] += p[j][k];
                h = splitmix64(h ^ ((uint64_t)(mach[j][k] + 1) * 1000000007ULL) ^ (uint64_t)p[j][k]);
            }
        }
        for (int j = 0; j < J; ++j) {
            for (int k = M - 1; k >= 0; --k) rem[j][k] = rem[j][k + 1] + p[j][k];
        }
        rng.seed(h);

        auto start = chrono::steady_clock::now();
        const double LIMIT = 0.88; // leave compile/runtime margin for a 1s TL

        Result best{LLONG_MAX, {}};
        for (int r = 0; r < 8; ++r) {
            Result cur = build(r);
            if (cur.makespan < best.makespan) best = move(cur);
        }

        int iter = 0;
        while (true) {
            double elapsed = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            if (elapsed > LIMIT) break;
            vector<double> w(5);
            // Randomized members of one family: active schedules with linear dispatch priorities.
            // Scales are chosen so no single unit dominates solely due to magnitude.
            w[0] = (urand() * 4.0 - 2.0);             // processing time
            w[1] = (urand() * 4.0 - 2.0);             // remaining work
            w[2] = (urand() * 1.5 - 0.25);            // earliest start
            w[3] = (urand() * 4.0 - 2.0);             // tail work
            w[4] = (urand() * 200000.0 - 100000.0);   // operations remaining
            double noise = (iter % 3 == 0) ? (urand() * 1e6) : 0.0;
            Result cur = build(8, w, noise);
            if (cur.makespan < best.makespan) best = move(cur);
            ++iter;
        }

        // All generated schedules append every job exactly once on every machine.
        // In case of malformed input fallback to identity-sized lines.
        if ((int)best.order.size() != M) {
            best.order.assign(M, vector<int>());
            for (int m = 0; m < M; ++m) for (int j = 0; j < J; ++j) best.order[m].push_back(j);
        }
        for (int m = 0; m < M; ++m) {
            if ((int)best.order[m].size() != J) {
                vector<int> seen(J, 0), fixed;
                for (int x : best.order[m]) if (0 <= x && x < J && !seen[x]) seen[x] = 1, fixed.push_back(x);
                for (int j = 0; j < J; ++j) if (!seen[j]) fixed.push_back(j);
                best.order[m] = fixed;
            }
            for (int i = 0; i < J; ++i) {
                if (i) cout << ' ';
                cout << best.order[m][i];
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
