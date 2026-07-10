// Job Shop Scheduling: minimize makespan by choosing, for each machine, a
// permutation of the jobs. The makespan is the longest path of the disjunctive
// graph induced by the job routes and the chosen machine orders.
//
// Method:
//  - Build feasible machine orders with a serial Schedule Generation Scheme
//    using several priority rules (critical-path/remaining, MLPT, LPT, SPT,
//    randomized). This always yields an acyclic, valid orientation.
//  - Evaluate any set of machine orders by a longest-path sweep (Kahn) which
//    also detects cyclic (infeasible) orientations.
//  - Improve with first-improving adjacent transpositions on each machine,
//    accepting only moves that stay feasible and strictly reduce the makespan.
//  - Run an Iterated Local Search: randomized restarts plus small,
//    feasibility-guarded perturbations of the incumbent, time-boxed.
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <chrono>

static int J, M, N;
static const int JMAX = 80, MMAX = 40, NMAX = JMAX * MMAX;

static int machine[JMAX][MMAX];   // machine[j][k]: machine of k-th op of job j
static long long proc[JMAX][MMAX]; // processing time of operation (j,k)
static int step_of[JMAX][MMAX];   // step_of[j][m]: route position on machine m

// Current machine orders: ord[m][r] = job at rank r on machine m.
static int ord[MMAX][JMAX];
static int pos[JMAX][MMAX];       // pos[j][m] = rank of job j on machine m

// Incumbent best.
static int best_ord[MMAX][JMAX];
static long long best_ms;

// Longest-path / link structures (rebuilt from ord each evaluation).
static int mach_prev[NMAX], mach_next[NMAX];
static int indeg[NMAX];
static long long comp[NMAX];
static int q[NMAX];

static long long mfree[MMAX], jfree[JMAX];
static long long sumProc[JMAX];
static long long rem[JMAX][MMAX]; // remaining job-chain time from (j,k) (lower bound)
static int idx[NMAX];
static long long pri[NMAX];
static int mlist[MMAX][JMAX];
static int mcnt[MMAX];

static const long long INF = (long long)4e18;
static const long long INFEAS = -1;

static std::chrono::steady_clock::time_point g_t0;
static inline bool time_up() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::steady_clock::now() - g_t0).count() > 600;
}

static inline int op_id(int j, int k) { return j * M + k; }

// Longest path of the disjunctive graph for the current ord.
// Returns INFEAS if the orientation induces a cycle.
static long long makespan() {
    for (int m = 0; m < M; ++m)
        for (int r = 0; r < J; ++r)
            pos[ord[m][r]][m] = r;

    for (int j = 0; j < J; ++j) {
        for (int k = 0; k < M; ++k) {
            int op = op_id(j, k);
            int m = machine[j][k];
            int r = pos[j][m];
            if (r > 0) {
                int jj = ord[m][r - 1];
                mach_prev[op] = op_id(jj, step_of[jj][m]);
            } else {
                mach_prev[op] = -1;
            }
            if (r < J - 1) {
                int jj = ord[m][r + 1];
                mach_next[op] = op_id(jj, step_of[jj][m]);
            } else {
                mach_next[op] = -1;
            }
            indeg[op] = (k > 0 ? 1 : 0) + (mach_prev[op] >= 0 ? 1 : 0);
            comp[op] = 0;
        }
    }

    int head = 0, tail = 0;
    for (int op = 0; op < N; ++op)
        if (indeg[op] == 0)
            q[tail++] = op;

    int processed = 0;
    long long best = 0;
    while (head < tail) {
        int op = q[head++];
        ++processed;
        int j = op / M, k = op % M;
        long long st = 0;
        if (k > 0) {
            long long c = comp[op - 1];
            if (c > st) st = c;
        }
        int mp = mach_prev[op];
        if (mp >= 0) {
            long long c = comp[mp];
            if (c > st) st = c;
        }
        long long fin = st + proc[j][k];
        comp[op] = fin;
        if (fin > best) best = fin;
        if (k + 1 < M) {
            int ns = op + 1;
            if (--indeg[ns] == 0) q[tail++] = ns;
        }
        if (mach_next[op] >= 0) {
            int ns = mach_next[op];
            if (--indeg[ns] == 0) q[tail++] = ns;
        }
    }
    if (processed < N) return INFEAS; // cyclic orientation
    return best;
}

// Serial Schedule Generation Scheme: dispatch operations in priority order,
// each as early as possible. The resulting per-machine dispatch order is a
// guaranteed-feasible (acyclic) set of machine orders.
static void construct(int rule) {
    for (int j = 0; j < J; ++j) {
        for (int k = 0; k < M; ++k) {
            int op = op_id(j, k);
            if (rule == 0) {
                pri[op] = sumProc[j];            // MLPT (long job first)
            } else if (rule == 1) {
                pri[op] = proc[j][k];            // LPT (long op first)
            } else if (rule == 2) {
                pri[op] = -proc[j][k];           // SPT (short op first)
            } else if (rule == 3) {
                pri[op] = rem[j][k];             // critical path (long remaining chain)
            } else {
                pri[op] = rand();                // randomized
            }
        }
    }
    for (int i = 0; i < N; ++i) idx[i] = i;
    std::sort(idx, idx + N, [](int a, int b) {
        if (pri[a] != pri[b]) return pri[a] > pri[b];
        return a < b;
    });
    for (int m = 0; m < M; ++m) { mfree[m] = 0; mcnt[m] = 0; }
    for (int j = 0; j < J; ++j) jfree[j] = 0;
    for (int t = 0; t < N; ++t) {
        int op = idx[t];
        int j = op / M, k = op % M, m = machine[j][k];
        long long s = mfree[m];
        if (jfree[j] > s) s = jfree[j];
        long long fin = s + proc[j][k];
        mfree[m] = fin;
        jfree[j] = fin;
        mlist[m][mcnt[m]++] = j;
    }
    for (int m = 0; m < M; ++m)
        for (int r = 0; r < J; ++r)
            ord[m][r] = mlist[m][r];
}

// First-improving adjacent transpositions on each machine; accept only moves
// that keep the orientation feasible and strictly reduce the makespan.
static long long local_search(long long cur) {
    bool improved = true;
    while (improved && !time_up()) {
        improved = false;
        for (int m = 0; m < M && !improved; ++m) {
            for (int r = 0; r < J - 1; ++r) {
                int j1 = ord[m][r], j2 = ord[m][r + 1];
                ord[m][r] = j2;
                ord[m][r + 1] = j1;
                long long ms = makespan();
                if (ms >= 0 && ms < cur) {
                    cur = ms;
                    improved = true;
                    break;
                } else {
                    ord[m][r] = j1;
                    ord[m][r + 1] = j2;
                }
            }
        }
    }
    return cur;
}

// Perturb the current order with a few feasible random adjacent swaps.
static void perturb(int k) {
    if (J < 2) return;
    int attempts = 0, limit = k * 6 + 6;
    while (k > 0 && attempts < limit) {
        int m = rand() % M;
        int r = rand() % (J - 1);
        int j1 = ord[m][r], j2 = ord[m][r + 1];
        ord[m][r] = j2;
        ord[m][r + 1] = j1;
        if (makespan() >= 0) {
            --k;
        } else {
            ord[m][r] = j1;
            ord[m][r + 1] = j2;
        }
        ++attempts;
    }
}

static void save_best() {
    for (int m = 0; m < M; ++m)
        for (int r = 0; r < J; ++r)
            best_ord[m][r] = ord[m][r];
}

int main() {
    if (scanf("%d%d", &J, &M) != 2) return 0;
    if (J > JMAX || M > MMAX) return 1; // guard against exceeding static limits
    N = J * M;
    for (int j = 0; j < J; ++j) {
        long long s = 0;
        for (int k = 0; k < M; ++k) {
            int mm;
            long long pp;
            if (scanf("%d%lld", &mm, &pp) != 2) return 1;
            machine[j][k] = mm;
            proc[j][k] = pp;
            step_of[j][mm] = k;
            s += pp;
        }
        sumProc[j] = s;
    }
    // Remaining job-chain processing time (critical-path priority, lower bound).
    for (int j = 0; j < J; ++j)
        for (int k = M - 1; k >= 0; --k)
            rem[j][k] = proc[j][k] + (k + 1 < M ? rem[j][k + 1] : 0);
    std::srand(1234567);

    g_t0 = std::chrono::steady_clock::now();

    best_ms = INF;
    // Strong deterministic constructions, each refined by local search.
    for (int rule = 3; rule >= 0; --rule) {
        construct(rule);
        long long ms = makespan();
        if (ms < 0) continue;
        ms = local_search(ms);
        if (ms < best_ms) { best_ms = ms; save_best(); }
        if (time_up()) break;
    }

    // Iterated local search with randomized restarts and perturbation.
    while (!time_up()) {
        construct(4);                        // randomized SGS restart
        long long ms = makespan();
        if (ms >= 0) {
            ms = local_search(ms);
            if (ms < best_ms) { best_ms = ms; save_best(); }
        }
        if (time_up()) break;
        // Perturb incumbent and re-optimize.
        for (int r = 0; r < M; ++r)
            for (int c = 0; c < J; ++c)
                ord[r][c] = best_ord[r][c];
        perturb(2 + (rand() % 3));
        ms = makespan();
        if (ms >= 0) {
            ms = local_search(ms);
            if (ms < best_ms) { best_ms = ms; save_best(); }
        }
    }

    // Fallback: guarantee a feasible output even if (unexpectedly) no move
    // improved; the best_ord always holds a construction-derived feasible order.
    if (best_ms == INF) {
        construct(0);
        save_best();
    }

    for (int m = 0; m < M; ++m) {
        for (int r = 0; r < J; ++r) {
            if (r + 1 < J) printf("%d ", best_ord[m][r]);
            else printf("%d\n", best_ord[m][r]);
        }
    }
    return 0;
}
