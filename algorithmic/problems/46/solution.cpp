#include <bits/stdc++.h>
using namespace std;

struct Op { int m; long long p; };
struct Instance {
    int J, M;
    vector<vector<Op>> ops;
    vector<vector<int>> posOnMachine;
    vector<long long> jobSum, machSum;
};

struct Schedule {
    long long C = LLONG_MAX;
    vector<vector<int>> seq;
};

static uint64_t rng_state = 88172645463325252ull;
static inline uint64_t xrng() {
    rng_state ^= rng_state << 7;
    rng_state ^= rng_state >> 9;
    return rng_state;
}
static inline double rnd01() { return (xrng() >> 11) * (1.0 / 9007199254740992.0); }

long long evaluate_sequence(const Instance& I, const vector<vector<int>>& seq) {
    int J = I.J, M = I.M, N = J * M;
    vector<int> indeg(N, 0);
    vector<vector<int>> adj(N);
    for (int u=0; u<N; ++u) adj[u].reserve(2);
    auto id = [M](int j, int k){ return j * M + k; };
    for (int j=0; j<J; ++j) {
        for (int k=0; k+1<M; ++k) {
            adj[id(j,k)].push_back(id(j,k+1));
            indeg[id(j,k+1)]++;
        }
    }
    for (int m=0; m<M; ++m) {
        if ((int)seq[m].size() != J) return LLONG_MAX/4;
        for (int a=0; a+1<J; ++a) {
            int j1 = seq[m][a], j2 = seq[m][a+1];
            if (j1<0 || j1>=J || j2<0 || j2>=J) return LLONG_MAX/4;
            int u = id(j1, I.posOnMachine[j1][m]);
            int v = id(j2, I.posOnMachine[j2][m]);
            adj[u].push_back(v);
            indeg[v]++;
        }
    }
    deque<int> q;
    vector<long long> dist(N, 0);
    for (int u=0; u<N; ++u) if (!indeg[u]) q.push_back(u);
    int seen = 0;
    long long C = 0;
    while (!q.empty()) {
        int u = q.front(); q.pop_front(); seen++;
        int j = u / M, k = u % M;
        long long fin = dist[u] + I.ops[j][k].p;
        C = max(C, fin);
        for (int v: adj[u]) {
            if (dist[v] < fin) dist[v] = fin;
            if (--indeg[v] == 0) q.push_back(v);
        }
    }
    return seen == N ? C : LLONG_MAX/4;
}

Schedule build_schedule(const Instance& I, int rule, const vector<double>& jw, const vector<double>& mw) {
    int J = I.J, M = I.M, N = J * M;
    vector<int> next(J, 0);
    vector<long long> jr(J, 0), mr(M, 0), rem = I.jobSum;
    vector<vector<int>> seq(M); for (int m=0;m<M;++m) seq[m].reserve(J);
    long long C = 0;

    for (int done = 0; done < N; ++done) {
        int best = -1;
        long double bestKey = 0;
        long long bestStart = 0;
        for (int j = 0; j < J; ++j) if (next[j] < M) {
            const Op &o = I.ops[j][next[j]];
            long long st = max(jr[j], mr[o.m]);
            long long fin = st + o.p;
            long double key;
            switch (rule) {
                case 0: // earliest finish, long remaining tie
                    key = (long double)fin * 1000000.0L - (long double)rem[j] - jw[j]; break;
                case 1: // earliest start, then longest processing on the ready frontier
                    key = (long double)st * 1000000.0L - (long double)o.p * 4.0L - (long double)rem[j] * 0.2L - jw[j]; break;
                case 2: // most work remaining (keeps long jobs moving)
                    key = (long double)st * 3.0L + (long double)fin * 0.4L - (long double)rem[j] * 2.0L - jw[j]; break;
                case 3: // shortest current operation among early starters
                    key = (long double)st * 1000000.0L + (long double)o.p * 10.0L - (long double)rem[j] * 0.05L - jw[j]; break;
                case 4: // bottleneck-oriented: don't idle heavy machines, favor their long jobs
                    key = (long double)st * 1000000.0L + (long double)fin - (long double)I.machSum[o.m] * 0.01L - (long double)o.p * mw[o.m] - jw[j]; break;
                case 5: // apparent tardiness/cost style using remaining tail
                    key = (long double)fin + (long double)max(0LL, st - jr[j]) * 2.0L - (long double)rem[j] * 0.6L - jw[j]; break;
                default: // random-key active schedule
                    key = (long double)st * (1.0L + mw[o.m]) + (long double)o.p * 0.35L - (long double)rem[j] * jw[j]; break;
            }
            if (best < 0 || key < bestKey || (key == bestKey && fin < bestStart)) {
                best = j; bestKey = key; bestStart = fin;
            }
        }
        Op o = I.ops[best][next[best]];
        long long st = max(jr[best], mr[o.m]);
        long long fin = st + o.p;
        jr[best] = mr[o.m] = fin;
        C = max(C, fin);
        seq[o.m].push_back(best);
        rem[best] -= o.p;
        next[best]++;
    }
    return {C, seq};
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Instance I;
    if (!(cin >> I.J >> I.M)) return 0;
    I.ops.assign(I.J, vector<Op>(I.M));
    I.posOnMachine.assign(I.J, vector<int>(I.M, -1));
    I.jobSum.assign(I.J, 0);
    I.machSum.assign(I.M, 0);
    for (int j=0;j<I.J;++j) {
        for (int k=0;k<I.M;++k) {
            int m; long long p; cin >> m >> p;
            I.ops[j][k] = {m,p};
            I.posOnMachine[j][m] = k;
            I.jobSum[j] += p;
            I.machSum[m] += p;
        }
    }

    // Deterministic seed from the instance, so repeated judging is reproducible.
    rng_state = 1469598103934665603ull;
    for (int j=0;j<I.J;++j) for (int k=0;k<I.M;++k) {
        rng_state ^= (uint64_t)(I.ops[j][k].m + 1000003) + (uint64_t)I.ops[j][k].p * 1000000007ull;
        rng_state *= 1099511628211ull;
    }

    auto startClock = chrono::steady_clock::now();
    auto elapsed = [&](){ return chrono::duration<double>(chrono::steady_clock::now() - startClock).count(); };

    Schedule best;
    vector<double> jw(I.J, 0.0), mw(I.M, 1.0);
    for (int r=0; r<6; ++r) {
        fill(jw.begin(), jw.end(), 0.0);
        fill(mw.begin(), mw.end(), 1.0);
        Schedule s = build_schedule(I, r, jw, mw);
        if (s.C < best.C) best = std::move(s);
    }

    long long maxMach = 1;
    for (long long x: I.machSum) maxMach = max(maxMach, x);
    int iter = 0;
    const double SEARCH_LIMIT = 0.74; // leave time for exact sequence improvement and runtime margin
    while (elapsed() < SEARCH_LIMIT) {
        int rule = 6 + (xrng() % 3);
        // Randomized but mechanism-consistent priority perturbations around active dispatch rules.
        for (int j=0;j<I.J;++j) {
            double a = 0.25 + 2.5 * rnd01();
            if (iter & 1) jw[j] = a * sqrt((double)max(1LL, I.jobSum[j]));
            else jw[j] = a;
        }
        for (int m=0;m<I.M;++m) {
            double load = (double)I.machSum[m] / (double)maxMach;
            mw[m] = 0.05 + 2.0 * load + 0.5 * rnd01();
        }
        Schedule s;
        if (rule == 6) s = build_schedule(I, xrng()%6, jw, mw);
        else s = build_schedule(I, rule, jw, mw);
        if (s.C < best.C) best = std::move(s);
        ++iter;
    }

    // Targeted local improvement: adjacent exchanges in the emitted machine permutations,
    // accepted only when an exact longest-path evaluation proves acyclicity and lower Cmax.
    if (!best.seq.empty()) best.C = evaluate_sequence(I, best.seq);
    const double LIMIT = 0.88;
    bool changed = true;
    while (changed && elapsed() < LIMIT) {
        changed = false;
        for (int m=0; m<I.M && elapsed() < LIMIT; ++m) {
            for (int a=0; a+1<I.J && elapsed() < LIMIT; ++a) {
                swap(best.seq[m][a], best.seq[m][a+1]);
                long long c = evaluate_sequence(I, best.seq);
                if (c < best.C) {
                    best.C = c;
                    changed = true;
                } else {
                    swap(best.seq[m][a], best.seq[m][a+1]);
                }
            }
        }
    }

    if (best.seq.empty()) { // should not happen, but produce a valid trivial order
        best.seq.assign(I.M, {});
        for (int m=0;m<I.M;++m) for (int j=0;j<I.J;++j) best.seq[m].push_back(j);
    }
    for (int m=0;m<I.M;++m) {
        for (int i=0;i<I.J;++i) {
            if (i) cout << ' ';
            cout << best.seq[m][i];
        }
        cout << '\n';
    }
    return 0;
}
