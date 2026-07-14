// Job Shop Scheduling (JSSP) — minimize makespan.
//
// Output: M lines, line m = a permutation of jobs {0..J-1} giving the order in
// which machine m processes the jobs. The judge builds the earliest-feasible
// schedule (longest path in the disjunctive graph) and scores clamp(1 - P/B, 0, 1).
//
// Approach:
//   1) Giffler-Thompson active-schedule construction under several priority rules;
//      keep the best makespan as the seed.
//   2) Critical-path neighborhood tabu search: swap two adjacent jobs on a machine
//      when the arc between them lies on a critical path. Reversing a critical arc
//      is always feasibility-preserving, but we still guard every candidate with a
//      full Kahn topo-count so any infeasible orientation is rejected.
//   3) Diversification (random adjacent swaps) on stagnation.
// A hard wall-clock budget keeps us under the 1s time limit.

#include <cstdio>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstdint>
using namespace std;

static int J, M, N;                      // jobs, machines, operations (= J*M)
static vector<int>       mach;           // mach[j*M+k]  = machine of k-th op of job j
static vector<long long> dur;            // dur[j*M+k]   = processing time
static vector<int>       posInJob;       // posInJob[j*M+m] = route index k where job j uses machine m
static vector<long long> rem;            // rem[j*M+k]   = remaining work of job j from op k onward

static inline int OP(int j, int k) { return j * M + k; }

// machNext/machPrev describe current machine orders; job route arcs are implicit
// (op k -> op k+1 within a job when k%M != M-1).
static vector<int>       machNext, machPrev;
static vector<long long> startT, tailT;
static vector<int>       indeg, order_;

// Compute makespan via longest path (Kahn). Returns -1 if the orientation is cyclic.
// If wantTail, also fills tailT[op] = longest path from op to a sink.
static long long computeSchedule(bool wantTail) {
    for (int op = 0; op < N; ++op) {
        int d = 0;
        if (op % M != 0) ++d;            // has job predecessor when k>0
        if (machPrev[op] >= 0) ++d;
        indeg[op] = d;
        startT[op] = 0;
    }
    int head = 0, tail = 0;
    for (int op = 0; op < N; ++op) if (indeg[op] == 0) order_[tail++] = op;
    long long makespan = 0;
    int processed = 0;
    while (head < tail) {
        int u = order_[head++];
        ++processed;
        long long fin = startT[u] + dur[u];
        if (fin > makespan) makespan = fin;
        int k = u % M;
        if (k + 1 < M) {                 // job successor
            int v = u + 1;
            if (startT[v] < fin) startT[v] = fin;
            if (--indeg[v] == 0) order_[tail++] = v;
        }
        int w = machNext[u];             // machine successor
        if (w >= 0) {
            if (startT[w] < fin) startT[w] = fin;
            if (--indeg[w] == 0) order_[tail++] = w;
        }
    }
    if (processed != N) return -1;       // cyclic -> infeasible
    if (wantTail) {
        for (int i = N - 1; i >= 0; --i) {
            int u = order_[i];
            long long t = dur[u];
            int k = u % M;
            if (k + 1 < M) t = max(t, dur[u] + tailT[u + 1]);
            int w = machNext[u];
            if (w >= 0) t = max(t, dur[u] + tailT[w]);
            tailT[u] = t;
        }
    }
    return makespan;
}

// Rebuild machNext/machPrev from machine job-orders `seq`.
static void buildLinks(const vector<vector<int>>& seq) {
    for (int op = 0; op < N; ++op) { machNext[op] = -1; machPrev[op] = -1; }
    for (int m = 0; m < M; ++m) {
        const vector<int>& s = seq[m];
        for (int i = 0; i + 1 < J; ++i) {
            int a = OP(s[i],   posInJob[s[i]   * M + m]);
            int b = OP(s[i+1], posInJob[s[i+1] * M + m]);
            machNext[a] = b;
            machPrev[b] = a;
        }
    }
}

// Giffler-Thompson active schedule generation.
// rule: 0=MWKR 1=SPT 2=LPT 3=LWKR 4=MOPNR 5=FCFS
static long long gifflerThompson(int rule, vector<vector<int>>& seq) {
    for (int m = 0; m < M; ++m) seq[m].clear();
    vector<int> routePos(J, 0);
    vector<long long> jobReady(J, 0), machReady(M, 0);
    int scheduled = 0;
    while (scheduled < N) {
        long long bestEC = -1; int bestMach = -1;
        for (int j = 0; j < J; ++j) {
            int k = routePos[j];
            if (k >= M) continue;
            int op = OP(j, k), m = mach[op];
            long long EC = max(jobReady[j], machReady[m]) + dur[op];
            if (bestEC < 0 || EC < bestEC) { bestEC = EC; bestMach = m; }
        }
        int m = bestMach;
        int chosenJob = -1; double chosenKey = 0;
        for (int j = 0; j < J; ++j) {
            int k = routePos[j];
            if (k >= M) continue;
            int op = OP(j, k);
            if (mach[op] != m) continue;
            long long ES = max(jobReady[j], machReady[m]);
            if (ES >= bestEC) continue;  // not in the conflict set at machine m
            double key;
            switch (rule) {
                case 0: key = -(double)rem[op];  break; // MWKR
                case 1: key =  (double)dur[op];  break; // SPT
                case 2: key = -(double)dur[op];  break; // LPT
                case 3: key =  (double)rem[op];  break; // LWKR
                case 4: key = -(double)(M - k);  break; // MOPNR
                default:key =  (double)j;        break; // FCFS
            }
            if (chosenJob < 0 || key < chosenKey) { chosenKey = key; chosenJob = j; }
        }
        int j = chosenJob, k = routePos[j], op = OP(j, k);
        long long fin = max(jobReady[j], machReady[m]) + dur[op];
        seq[m].push_back(j);
        jobReady[j] = fin; machReady[m] = fin;
        routePos[j] = k + 1;
        ++scheduled;
    }
    buildLinks(seq);
    return computeSchedule(false);
}

// deterministic xorshift RNG (reproducible scoring, no wall-clock/PID dependence)
static uint64_t rngState = 0x9E3779B97F4A7C15ULL;
static inline uint64_t rnd() {
    rngState ^= rngState << 13; rngState ^= rngState >> 7; rngState ^= rngState << 17;
    return rngState;
}

int main() {
    auto t0 = chrono::steady_clock::now();
    if (scanf("%d %d", &J, &M) != 2) return 0;
    N = J * M;
    mach.assign(N, 0); dur.assign(N, 0);
    posInJob.assign(N, 0); rem.assign(N, 0);
    for (int j = 0; j < J; ++j) {
        for (int k = 0; k < M; ++k) {
            int m; long long p;
            if (scanf("%d %lld", &m, &p) != 2) return 0;
            mach[OP(j, k)] = m; dur[OP(j, k)] = p;
            posInJob[j * M + m] = k;
        }
        long long acc = 0;
        for (int k = M - 1; k >= 0; --k) { acc += dur[OP(j, k)]; rem[OP(j, k)] = acc; }
    }

    machNext.assign(N, -1); machPrev.assign(N, -1);
    startT.assign(N, 0);   tailT.assign(N, 0);
    indeg.assign(N, 0);    order_.assign(N, 0);

    rngState = 0x9E3779B97F4A7C15ULL ^ ((uint64_t)J << 32) ^ (uint64_t)(M * 2654435761u);
    for (int op = 0; op < N && op < 64; ++op) rngState ^= (uint64_t)dur[op] * 0x100000001B3ULL + op;
    if (rngState == 0) rngState = 1;

    // ---- construction ----
    vector<vector<int>> seq(M), best(M), cur(M);
    long long bestMk = -1;
    for (int rule = 0; rule < 6; ++rule) {
        long long mk = gifflerThompson(rule, seq);
        if (mk > 0 && (bestMk < 0 || mk < bestMk)) { bestMk = mk; best = seq; }
    }
    if (bestMk < 0) {                    // fallback: route-order sequences (always feasible)
        for (int m = 0; m < M; ++m) { best[m].resize(J); for (int j = 0; j < J; ++j) best[m][j] = j; }
        buildLinks(best);
        bestMk = computeSchedule(false);
    }
    cur = best;
    buildLinks(cur);
    long long curMk = computeSchedule(true);   // fills startT/tailT for cur

    // ---- critical-path tabu search ----
    const double BUDGET = 0.78;          // seconds, comfortably under the 1s limit
    vector<int> tabuUntil(N, 0);
    const int TENURE = max(6, N / 12);
    int iter = 0;

    // Iterated Local Search: run critical-path tabu until it stagnates, then kick from
    // the global best and re-optimize. Kick strength grows when kicks keep failing.
    vector<vector<int>> localBest = best;      // best of the current ILS restart
    long long localBestMk = bestMk;
    int stagn = 0;                              // consecutive non-improving tabu iters
    const int MAX_STAGN = max(150, 2 * N);      // stagnation before a kick
    int kickStrength = max(2, N / 40);
    int failedKicks = 0;

    vector<int> candA, candB, candM;

    auto timeUp = [&]() {
        return chrono::duration<double>(chrono::steady_clock::now() - t0).count() > BUDGET;
    };
    // Apply `k` random adjacent machine swaps to `cur` (feasibility guarded by caller).
    auto kick = [&](int k) {
        for (int t = 0; t < k; ++t) {
            int m = rnd() % M, i = rnd() % (J - 1);
            swap(cur[m][i], cur[m][i + 1]);
        }
    };

    while (J >= 2) {
        if (timeUp()) break;
        ++iter;

        // collect critical machine arcs a->b (both endpoints on a critical path, arc tight)
        candA.clear(); candB.clear(); candM.clear();
        for (int op = 0; op < N; ++op) {
            int b = machNext[op];
            if (b < 0) continue;
            if (startT[op] + dur[op] == startT[b] &&
                startT[op] + tailT[op] == curMk &&
                startT[b]  + tailT[b]  == curMk) {
                candA.push_back(op); candB.push_back(b); candM.push_back(mach[op]);
            }
        }

        bool improvedLocal = false;
        if (!candA.empty()) {
            // evaluate each candidate: swap the two adjacent jobs in cur[m], eval, revert.
            long long bestMoveMk = -1; int bestIdx = -1;
            for (size_t c = 0; c < candA.size(); ++c) {
                if ((c & 63) == 0 && timeUp()) break;   // bail mid-iteration near the budget
                int a = candA[c], b = candB[c];          // machine arc a->b (machNext[a]==b)
                // reverse the arc via local link surgery: pa->a->b->nb  =>  pa->b->a->nb
                int pa = machPrev[a], nb = machNext[b];
                machPrev[b] = pa; machNext[b] = a;
                machPrev[a] = b;  machNext[a] = nb;
                if (pa >= 0) machNext[pa] = b;
                if (nb >= 0) machPrev[nb] = a;
                long long mk = computeSchedule(false);
                // revert
                machPrev[a] = pa; machNext[a] = b;
                machPrev[b] = a;  machNext[b] = nb;
                if (pa >= 0) machNext[pa] = a;
                if (nb >= 0) machPrev[nb] = b;
                if (mk < 0) continue;
                bool aspire = (mk < localBestMk);
                bool tabu = (tabuUntil[a] > iter || tabuUntil[b] > iter);
                if (tabu && !aspire) continue;
                if (bestIdx < 0 || mk < bestMoveMk) { bestMoveMk = mk; bestIdx = (int)c; }
            }
            if (bestIdx < 0) bestIdx = 0;    // all moves tabu -> pick first candidate anyway

            int m = candM[bestIdx];
            int ja = candA[bestIdx] / M, jb = candB[bestIdx] / M;
            {
                vector<int>& s = cur[m];
                int ia = -1, ib = -1;
                for (int i = 0; i < J; ++i) { if (s[i] == ja) ia = i; else if (s[i] == jb) ib = i; }
                if (ia >= 0 && ib >= 0) swap(s[ia], s[ib]);
            }
            buildLinks(cur);
            curMk = computeSchedule(true);
            tabuUntil[candA[bestIdx]] = iter + TENURE;
            tabuUntil[candB[bestIdx]] = iter + TENURE;

            if (curMk >= 0 && curMk < localBestMk) { localBestMk = curMk; localBest = cur; improvedLocal = true; }
            if (curMk >= 0 && curMk < bestMk)      { bestMk = curMk; best = cur; }
        }

        if (improvedLocal) stagn = 0; else ++stagn;

        // ILS kick when the local search stalls (or there were no critical arcs to move)
        if (stagn >= MAX_STAGN || candA.empty()) {
            if (localBestMk <= bestMk) failedKicks++; else failedKicks = 0;
            // grow perturbation the longer we fail to find something new
            int k = kickStrength + min(3 * failedKicks, 4 * kickStrength);
            cur = best;                       // always kick from the global best (better-walk ILS)
            kick(k);
            buildLinks(cur);
            long long mk = computeSchedule(true);
            if (mk < 0) { cur = best; buildLinks(cur); mk = computeSchedule(true); }  // undo cyclic kick
            curMk = mk;
            localBest = cur; localBestMk = mk;
            for (int i = 0; i < N; ++i) tabuUntil[i] = 0;   // fresh tabu for the new restart
            stagn = 0;
            if (curMk >= 0 && curMk < bestMk) { bestMk = curMk; best = cur; }
        }
    }

    // ---- output ----
    for (int m = 0; m < M; ++m)
        for (int j = 0; j < J; ++j) printf(j + 1 < J ? "%d " : "%d\n", best[m][j]);
    return 0;
}
