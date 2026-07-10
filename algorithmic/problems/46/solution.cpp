// Job Shop Scheduling (JSSP) solver.
//
// Input: J jobs, M machines.  Each job visits every machine exactly once in a
//   job-specific order, with a positive integer processing time per operation.
// Output: for each machine m (0..M-1) a permutation of 0..J-1 giving the order
//   in which that machine processes the jobs.  The judge derives the
//   earliest-feasible schedule (longest path in the disjunctive graph) and
//   scores clamp(1 - makespan / baseline, 0, 1).
//
// Algorithm:
//   1. Construct feasible schedules with a critical-path list scheduler using
//      several priority rules (remaining job-chain length, LPT, SPT).  List
//      scheduling always yields acyclic machine orders, so output is always
//      valid.
//   2. Evaluate any candidate by the exact judge objective: the longest path
//      in the disjunctive graph (job-chain arcs + machine-order arcs), via
//      Kahn's topological sort carrying earliest-finish times.  Cyclic orders
//      are detected and rejected.
//   3. Refine the best construction with first-improvement adjacent
//      transpositions on each machine, scanning machines in decreasing order
//      of total processing load (bottleneck machines first for faster
//      convergence).
//   4. Escape local optima via iterated local search with current-solution
//      regrowth: perturb the current incumbent with random adjacent swaps,
//      re-optimize, and always accept the result as the new incumbent while
//      tracking the global best.  Periodic restarts from fresh randomized
//      constructions provide strong diversification.  All compute is bounded
//      by a wall-clock deadline and an evaluation cap.

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <chrono>
using namespace std;

static int J = 0, M = 0, N = 0;

// Per-operation data.  Operation id = j*M + k (job j, position k in route).
static vector<int>        opMach;    // machine index of each operation
static vector<long long>  opTime;    // processing time of each operation
static vector<long long>  opLevel;   // critical-path level (remaining job-chain time)
static vector<int>        machOrder; // machines sorted by decreasing total load

// Reusable graph buffers for makespan evaluation (no per-call allocation).
static vector<int>        head;      // adjacency-list head (linked by index)
static vector<int>        nxt;       // next-arc link
static vector<int>        to;        // arc destination
static vector<int>        indeg;     // in-degree counter
static vector<long long>  dist;      // longest-path finish-to-start time
static vector<int>        qbuf;      // Kahn queue
static int                nArcs;
static long long          gEvals;

static const long long INF = (long long)4e18;
static const long long EVAL_CAP = 500000;       // safety cap on evaluations
static const int       TIME_BUDGET_MS = 400;   // wall-clock budget (well under 1.1s)
static const int       INIT_BUDGET_MS = 100;   // initial local search budget

static inline void addArc(int u, int v) {
    nxt[nArcs] = head[u];
    to[nArcs]  = v;
    head[u]    = nArcs++;
}

// Makespan implied by machine orders `mo`: longest path in the disjunctive
// graph.  Returns -1 if the orders induce a cycle (output would be invalid).
// Flat buffers are reused so there is no per-call allocation in this hot path.
static long long makespan(const vector<vector<int>>& mo) {
    ++gEvals;

    // Reset adjacency and in-degree.
    memset(head.data(),  -1, sizeof(int) * (size_t)N);
    memset(indeg.data(),  0, sizeof(int) * (size_t)N);
    nArcs = 0;

    // Job-chain arcs: op(j,k) -> op(j,k+1).
    for (int j = 0; j < J; j++) {
        int base = j * M;
        for (int k = 0; k + 1 < M; k++) {
            addArc(base + k, base + k + 1);
            ++indeg[base + k + 1];
        }
    }

    // Machine-order arcs: consecutive operations on each machine.
    for (int m = 0; m < M; m++) {
        const vector<int>& o = mo[m];
        for (size_t i = 1; i < o.size(); i++) {
            addArc(o[i - 1], o[i]);
            ++indeg[o[i]];
        }
    }

    // Kahn topological pass carrying longest finish times.
    memset(dist.data(), 0, sizeof(long long) * (size_t)N);
    int qh = 0, qt = 0;
    for (int i = 0; i < N; i++)
        if (!indeg[i]) qbuf[qt++] = i;

    int done = 0;
    long long best = 0;
    while (qh < qt) {
        int u = qbuf[qh++];
        ++done;
        long long f = dist[u] + opTime[u];        // earliest finish of u
        if (f > best) best = f;
        for (int e = head[u]; e != -1; e = nxt[e]) {
            int v = to[e];
            if (f > dist[v]) dist[v] = f;          // relax v's start time
            if (--indeg[v] == 0) qbuf[qt++] = v;
        }
    }
    if (done < N) return -1;                        // cycle: orders infeasible
    return best;
}

// Critical-path list scheduling.  Operations are dispatched in order of the
// given static priority; an operation is released only after its job
// predecessor finishes.  Always produces acyclic (valid) machine orders.
// When `rng` is non-null, a random tie-break tag is mixed into the priority
// to diversify the dispatch order and yield different feasible schedules.
static vector<vector<int>> construct(const vector<long long>& prio,
                                     mt19937_64* rng = nullptr) {
    vector<vector<int>> mo(M);
    vector<long long> mf(M, 0), jf(N, 0);           // machine-free, job-finish
    using P = pair<pair<long long, uint32_t>, int>; // ((priority, tiebreak), op_id)
    priority_queue<P> pq;
    auto tag = [&]() { return rng ? (uint32_t)((*rng)() >> 33) : 0u; };
    for (int j = 0; j < J; j++)
        pq.push({{prio[j * M], tag()}, j * M});
    while (!pq.empty()) {
        int op = pq.top().second; pq.pop();         // op_id stored correctly
        int m = opMach[op];
        long long pred = (op % M == 0) ? 0 : jf[op - 1];  // job-chain finish
        long long fin  = max(mf[m], pred) + opTime[op];
        mf[m] = fin;
        jf[op] = fin;
        mo[m].push_back(op);
        if (op % M + 1 < M)                         // release next operation
            pq.push({{prio[op + 1], tag()}, op + 1});
    }
    return mo;
}

// First-improvement adjacent-transposition local search on every machine.
// Machines are scanned in decreasing order of total load (bottleneck first)
// for faster convergence.  Accepts the first improving swap found, rescans,
// and repeats until no improvement in a full pass, the pass cap, or the
// budget is exhausted.
static void localSearch(vector<vector<int>>& mo, long long& cur,
                        chrono::steady_clock::time_point dl, int maxPasses) {
    bool improved = true;
    int  passes   = 0;
    while (improved && passes < maxPasses) {
        improved = false;
        ++passes;
        for (int mi = 0; mi < M; mi++) {
            int m = machOrder[mi];
            vector<int>& o = mo[m];
            int sz = (int)o.size();
            if (sz < 2) continue;
            for (int i = 0; i + 1 < sz; i++) {
                swap(o[i], o[i + 1]);
                long long ns = makespan(mo);
                if (ns >= 0 && ns < cur) {
                    cur = ns;
                    improved = true;                 // accept and continue
                } else {
                    swap(o[i], o[i + 1]);            // reject: restore
                }
                // Periodic budget guard (cheap time check every 64 evals).
                if ((gEvals & 63) == 0 &&
                    (gEvals >= EVAL_CAP || chrono::steady_clock::now() >= dl))
                    return;
            }
        }
    }
}

int main() {
    if (scanf("%d%d", &J, &M) != 2) return 0;
    N = J * M;

    opMach.assign(N, 0);
    opTime.assign(N, 0);
    opLevel.assign(N, 0);
    head.assign(N, -1);
    nxt.assign(2 * N + 8, 0);
    to.assign(2 * N + 8, 0);
    indeg.assign(N, 0);
    dist.assign(N, 0);
    qbuf.assign(N, 0);
    gEvals = 0;

    for (int j = 0; j < J; j++)
        for (int k = 0; k < M; k++) {
            int op = j * M + k;
            if (scanf("%d%lld", &opMach[op], &opTime[op]) != 2) return 1;
        }

    // Critical-path levels: remaining processing time along the job chain to
    // the end.  Higher level = more urgent to schedule.
    for (int j = 0; j < J; j++) {
        int base = j * M;
        opLevel[base + M - 1] = opTime[base + M - 1];
        for (int k = M - 2; k >= 0; k--)
            opLevel[base + k] = opTime[base + k] + opLevel[base + k + 1];
    }

    // Machine order by decreasing total load (bottleneck machines first in
    // local search for faster convergence).
    machOrder.resize(M);
    for (int i = 0; i < M; i++) machOrder[i] = i;
    vector<long long> load(M, 0);
    for (int i = 0; i < N; i++) load[opMach[i]] += opTime[i];
    sort(machOrder.begin(), machOrder.end(),
         [&](int a, int b) { return load[a] > load[b]; });

    // Deterministic per-input seed (FNV-1a hash of the instance).
    uint64_t seed = 1469598103934665603ULL;
    auto mix = [&](uint64_t x) { seed = (seed ^ x) * 1099511628211ULL; };
    mix((uint64_t)J); mix((uint64_t)M);
    for (int i = 0; i < N; i++) {
        mix((uint64_t)opMach[i]);
        mix((uint64_t)opTime[i]);
    }
    mt19937_64 rng(seed);

    // Priority vectors for several dispatch rules.
    vector<long long> pLevel(N), pLpt(N), pSpt(N);
    for (int i = 0; i < N; i++) {
        pLevel[i] = opLevel[i];
        pLpt[i]   = opTime[i];
        pSpt[i]   = -opTime[i];
    }

    // Try multiple constructions; keep the best by makespan.
    vector<vector<int>> bestMo;
    long long best = INF;

    auto tryConstruct = [&](const vector<long long>& p, bool randomize) {
        auto mo = construct(p, randomize ? &rng : nullptr);
        long long ms = makespan(mo);
        if (ms >= 0 && ms < best) {
            best = ms;
            bestMo = std::move(mo);
        }
    };

    tryConstruct(pLevel, false);   // critical-path level (remaining job-chain)
    tryConstruct(pLpt,   false);   // longest processing time first
    tryConstruct(pSpt,   false);   // shortest processing time first
    tryConstruct(pLevel, true);    // randomized level (tie-break)

    // Deadlines: initial search gets at most INIT_BUDGET_MS; the total
    // computation (ILS + polish) is bounded by TIME_BUDGET_MS from here.
    auto t0     = chrono::steady_clock::now();
    auto dlInit = t0 + chrono::milliseconds(INIT_BUDGET_MS);
    auto dlTotal = t0 + chrono::milliseconds(TIME_BUDGET_MS);

    // Initial local search from the best construction.
    localSearch(bestMo, best, dlInit, 1000);

    // Iterated local search with current-solution regrowth (standard ILS):
    // perturb the current incumbent, re-optimize, and always accept as the
    // new current.  The global best is tracked separately.  Every 8 iterations
    // a fresh randomized construction provides a strong diversification kick.
    vector<vector<int>> curMo = bestMo;
    int iter = 0;
    while (chrono::steady_clock::now() < dlTotal && gEvals < EVAL_CAP) {
        if ((iter & 7) == 7) {
            // Diversification: restart from a new randomized construction.
            vector<vector<int>> cand = construct(pLevel, &rng);
            long long ns = makespan(cand);
            if (ns >= 0) {
                localSearch(cand, ns, dlTotal, 30);
                curMo = cand;
                if (ns < best) { best = ns; bestMo = cand; }
            } else {
                curMo = bestMo;
            }
        } else {
            // Perturbation: 5..15 random adjacent transpositions on current.
            vector<vector<int>> cand = curMo;
            int K = 5 + (int)(rng() % 11);
            for (int t = 0; t < K; t++) {
                int m = (int)(rng() % M);
                vector<int>& o = cand[m];
                int sz = (int)o.size();
                if (sz < 2) continue;
                int i = (int)(rng() % (sz - 1));
                swap(o[i], o[i + 1]);
            }
            long long ns = makespan(cand);
            if (ns >= 0) {
                localSearch(cand, ns, dlTotal, 4);
                curMo = cand;
                if (ns < best) { best = ns; bestMo = cand; }
            }
        }
        ++iter;
    }

    // Final polish from the best incumbent.
    localSearch(bestMo, best, dlTotal, 1000);

    // Emit one permutation per machine.  Operation id j*M+k maps to job j.
    for (int m = 0; m < M; m++) {
        for (int j = 0; j < J; j++) {
            if (j) putchar(' ');
            printf("%d", bestMo[m][j] / M);
        }
        putchar('\n');
    }
    return 0;
}
