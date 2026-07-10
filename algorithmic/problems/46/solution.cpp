#include <cstdio>
#include <vector>
#include <queue>
#include <cstring>
#include <chrono>
#include <random>
#include <algorithm>

// Job Shop Scheduling (JSSP) solver.
//
// Input : J jobs and M machines. Each job visits every machine exactly once in
//         a job-specific order, with a positive processing time per operation.
// Output: for each machine m (0..M-1) a permutation of 0..J-1 giving the order
//         in which that machine processes the jobs. The judge builds the
//         earliest-feasible schedule from these machine orders together with
//         the job chains; its makespan is the longest path in the disjunctive
//         graph and is exactly what we minimize.
//
// Algorithm:
//   1. Construct a feasible incumbent with a critical-path list scheduler.
//      Each job's operations are released (in order of decreasing remaining
//      job-chain length) only after the previous operation of that job has
//      finished, and each operation starts as early as its machine and its
//      job predecessor allow. This greedily derives machine orders that are
//      always acyclic, so the emitted output is always valid.
//   2. Score any candidate with the exact objective the judge uses: the longest
//      path in the disjunctive graph, computed by a Kahn topological pass
//      carrying earliest-finish times. Orders that would create a cycle are
//      detected and rejected, so an invalid schedule is never produced.
//   3. Refine the incumbent with first-improvement adjacent transpositions on
//      each machine (the classic 2-opt neighborhood for this representation),
//      then escape local optima via iterated local search with random
//      perturbations, retaining the best makespan found. All work is bounded
//      by a hard evaluation budget and a wall-clock deadline.

using namespace std;

static int J = 0, M = 0, N = 0;
static vector<int> machine_of;       // machine index of operation (j*M + k)
static vector<long long> proc;       // processing time of operation
static vector<int> job_pred, job_next;  // job-chain predecessor / successor
static vector<int> mac_pred, mac_next;  // machine-order predecessor / successor
static vector<int> inDeg, qbuf;      // scratch: in-degree and Kahn queue
static vector<long long> maxPred;    // scratch: best predecessor finish
static vector<long long> finish;     // scratch: finish time of each operation

static long long gEvals = 0;
// Termination is dominated by a hard evaluation budget (not just wall-clock:
// a count bound guarantees the time cap is never exceeded). Each makespan call
// is O(N) with no allocation, so this budget stays safely inside the per-case
// time limit across machines while still permitting substantial local search.
static const long long EVAL_CAP = 60000;          // hard bound on evaluations
static const long long TIME_BUDGET_MS = 500;       // wall-clock safety bound

static const long long INF = (long long)4e18;

// Earliest-feasible makespan implied by the machine orders in `order` (which
// store operation ids). Returns INF if the orders induce a cycle. Flat arrays
// are reused so there is no per-call allocation in this hot path.
static long long makespan(const vector<vector<int>>& order) {
    ++gEvals;

    // Machine linkage from the orders. Every operation appears on exactly one
    // machine, so mac_pred is overwritten for all ops; mac_next is only set for
    // non-last entries, hence clear it first.
    memset(mac_next.data(), -1, sizeof(int) * (size_t)N);
    for (int m = 0; m < M; ++m) {
        const vector<int>& om = order[m];
        int prev = -1;
        for (size_t p = 0; p < om.size(); ++p) {
            int o = om[p];
            mac_pred[o] = prev;
            if (prev != -1) mac_next[prev] = o;
            prev = o;
        }
    }

    // In-degree from job chain + machine order; seed the queue with sources.
    int processed = 0, qh = 0, qt = 0;
    for (int o = 0; o < N; ++o) {
        int d = (job_pred[o] != -1) + (mac_pred[o] != -1);
        inDeg[o] = d;
        maxPred[o] = 0;
        if (d == 0) { finish[o] = proc[o]; qbuf[qt++] = o; }
    }

    // Kahn topological pass carrying longest finish times.
    while (qh < qt) {
        int u = qbuf[qh++];
        ++processed;
        long long f = finish[u];
        int s = job_next[u];
        if (s != -1) {
            if (f > maxPred[s]) maxPred[s] = f;
            if (--inDeg[s] == 0) { finish[s] = proc[s] + maxPred[s]; qbuf[qt++] = s; }
        }
        s = mac_next[u];
        if (s != -1) {
            if (f > maxPred[s]) maxPred[s] = f;
            if (--inDeg[s] == 0) { finish[s] = proc[s] + maxPred[s]; qbuf[qt++] = s; }
        }
    }
    if (processed < N) return INF;  // cycle -> orders would be invalid

    long long best = 0;
    for (int o = 0; o < N; ++o) if (finish[o] > best) best = finish[o];
    return best;
}

// Critical-path level: remaining processing time along the job chain to the
// end. Higher level => more urgent to schedule first.
static vector<long long> comp_level() {
    vector<long long> lvl(N, 0);
    for (int j = J - 1; j >= 0; --j)
        for (int k = M - 1; k >= 0; --k) {
            int o = j * M + k;
            lvl[o] = proc[o] + (k + 1 < M ? lvl[j * M + k + 1] : 0);
        }
    return lvl;
}

// Build machine orders with a priority-dispatched list schedule. Operation id
// is j*M + k, so the job predecessor is op-1 and the job id is op/M.
static vector<vector<int>> construct(const vector<long long>& prio) {
    vector<vector<int>> mo(M);
    vector<long long> machineFree(M, 0), jobFinish(N, 0);
    using P = pair<long long, int>;
    priority_queue<P> pq;
    for (int j = 0; j < J; ++j) pq.push({prio[j * M], j * M});
    while (!pq.empty()) {
        int op = pq.top().second; pq.pop();
        int m = machine_of[op];
        long long pred = (op % M == 0) ? 0 : jobFinish[op - 1];       // job-chain finish
        long long start = machineFree[m] > pred ? machineFree[m] : pred;
        long long fin = start + proc[op];
        machineFree[m] = fin;
        jobFinish[op] = fin;
        mo[m].push_back(op);                                          // record operation id on machine
        if (op % M + 1 < M) pq.push({prio[op + 1], op + 1});          // release next op of job
    }
    return mo;
}

// First-improvement adjacent-transposition local search on every machine.
// Keeps an improving swap and continues scanning; repeats passes until no
// improvement in a full pass, the pass cap, or the evaluation/time budget.
static void local_search(vector<vector<int>>& mo, long long& cur,
                         chrono::steady_clock::time_point deadline, int maxPasses) {
    bool improved = true;
    int passes = 0;
    while (improved && passes < maxPasses) {
        improved = false;
        ++passes;
        for (int m = 0; m < M; ++m) {
            vector<int>& o = mo[m];
            int sz = (int)o.size();
            if (sz < 2) continue;
            for (int i = 0; i + 1 < sz; ++i) {
                swap(o[i], o[i + 1]);
                long long ns = makespan(mo);
                if (ns < cur) {
                    cur = ns;
                    improved = true;            // accept and advance
                } else {
                    swap(o[i], o[i + 1]);       // reject: restore
                }
                if ((gEvals & 31) == 0) {      // periodic budget guard
                    if (gEvals >= EVAL_CAP || chrono::steady_clock::now() >= deadline)
                        return;
                }
            }
        }
    }
}

int main() {
    if (scanf("%d%d", &J, &M) != 2) return 0;
    N = J * M;
    machine_of.assign(N, 0);            proc.assign(N, 0);
    job_pred.assign(N, -1);             job_next.assign(N, -1);
    mac_pred.assign(N, -1);             mac_next.assign(N, -1);
    inDeg.assign(N, 0);                 maxPred.assign(N, 0);
    finish.assign(N, 0);                qbuf.assign(N, 0);

    for (int j = 0; j < J; ++j)
        for (int k = 0; k < M; ++k) {
            int op = j * M + k, machine;
            long long p;
            if (scanf("%d%lld", &machine, &p) != 2) return 1;
            machine_of[op] = machine;
            proc[op] = p;
            job_pred[op] = (k > 0) ? op - 1 : -1;
            job_next[op] = (k + 1 < M) ? op + 1 : -1;
        }

    // Several dispatch rules give robust starts across instance types
    // (random, flow-shop-like, bottleneck-heavy); keep the best by makespan.
    vector<long long> p_level = comp_level();
    vector<long long> p_lpt(N), p_spt(N);
    for (int o = 0; o < N; ++o) { p_lpt[o] = proc[o]; p_spt[o] = -proc[o]; }

    vector<vector<int>> o_level = construct(p_level);
    vector<vector<int>> o_lpt = construct(p_lpt);
    vector<vector<int>> o_spt = construct(p_spt);
    long long m_level = makespan(o_level);
    long long m_lpt = makespan(o_lpt);
    long long m_spt = makespan(o_spt);

    vector<vector<int>> best_order;
    long long best_ms;
    if (m_level <= m_lpt && m_level <= m_spt) { best_order = o_level; best_ms = m_level; }
    else if (m_lpt <= m_spt)                  { best_order = o_lpt;  best_ms = m_lpt;  }
    else                                      { best_order = o_spt;  best_ms = m_spt;  }

    // Per-instance deterministic seed (input-hashed) for reproducible results.
    uint64_t seed = 1469598103934665603ULL;
    auto mix = [&](uint64_t x) { seed = (seed ^ x) * 1099511628211ULL; };
    mix((uint64_t)J); mix((uint64_t)M);
    for (int o = 0; o < N; ++o) { mix((uint64_t)machine_of[o]); mix((uint64_t)proc[o]); }
    mt19937_64 rng(seed);

    auto deadline = chrono::steady_clock::now() + chrono::milliseconds(TIME_BUDGET_MS);

    // Local optimum from the best construction (allow a thorough first pass).
    local_search(best_order, best_ms, deadline, 24);

    // Iterated local search: perturb and re-optimize, retaining the global best.
    while (chrono::steady_clock::now() < deadline && gEvals < EVAL_CAP) {
        vector<vector<int>> cand = best_order;
        int K = 2 + (int)(rng() % 3);  // 2..4 random adjacent transpositions
        for (int t = 0; t < K; ++t) {
            int m = (int)(rng() % M);
            vector<int>& o = cand[m];
            int sz = (int)o.size();
            if (sz < 2) continue;
            int i = (int)(rng() % (sz - 1));
            swap(o[i], o[i + 1]);
        }
        long long cur = makespan(cand);
        if (cur >= INF) continue;                 // reject infeasible orders
        local_search(cand, cur, deadline, 5);
        if (cur < best_ms) { best_ms = cur; best_order = cand; }
    }

    for (int m = 0; m < M; ++m) {
        const vector<int>& om = best_order[m];
        for (int j = 0; j < J; ++j)
            printf("%d%s", om[j] / M, j + 1 < J ? " " : "\n");  // op id -> job id
    }
    return 0;
}
