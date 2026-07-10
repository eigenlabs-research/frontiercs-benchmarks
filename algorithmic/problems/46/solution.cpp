#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <chrono>

// Job Shop Scheduling heuristic solver.
// Emits one permutation of the jobs for each machine. The judge turns these
// machine orders into the earliest-feasible schedule (longest path in the
// disjunctive graph formed by the machine arcs and the job precedence chains)
// and measures its makespan.
//
// Two complementary strategies are used:
//  * Flow-shop instances (every job visits machines 0..M-1 in natural order)
//    are tackled with the NEH insertion heuristic followed by an iterated
//    adjacent-transposition search over a single permutation, evaluated with
//    the standard O(J*M) flow-shop makespan dynamic program.
//  * General job-shop instances use several serial schedule generation priority
//    rules (critical-remaining-work, shortest-processing-time,
//    longest-processing-time), each refined by first-improvement adjacent
//    transposition local search inside an iterated-local-search loop. Feasible
//    (cycle-free) orders are required and verified by a longest-path check.

#define MAXDIM 96
#define MAXOP  (MAXDIM * MAXDIM)

static int J, M, n;                       // jobs, machines, total operations
static long long procOf[MAXOP];           // processing time of operation (j*M+k)
static int8_t machOf[MAXOP];              // machine of operation
static long long remWork[MAXOP];          // remaining job processing from this op to job end
static int jobPosMachine[MAXDIM][MAXDIM]; // route position k of machine m in job j

static int machOrder[MAXDIM][MAXDIM];     // machOrder[m][i] = job processed i-th on machine m
static int machLen[MAXDIM];               // jobs assigned to each machine (== J)
static int bestOrder[MAXDIM][MAXDIM];     // best machine order found so far

// Reusable workspace for the longest-path evaluation (general instances).
static int indeg[MAXOP], outcnt[MAXOP], s0[MAXOP], s1[MAXOP], qbuf[MAXOP];
static long long startT[MAXOP];

// Longest path (makespan) in the disjunctive graph defined by the current
// machine orders and the job precedence chains. Returns -1 if the orders induce
// a cycle (infeasible); otherwise returns the makespan.
static inline long long makespan() {
    for (int i = 0; i < n; ++i) {
        indeg[i] = outcnt[i] = 0;
        s0[i] = s1[i] = -1;
        startT[i] = 0;
    }
    // Machine arcs: consecutive jobs on each machine define the orientation.
    for (int m = 0; m < M; ++m) {
        for (int i = 0; i + 1 < machLen[m]; ++i) {
            int ja = machOrder[m][i], jb = machOrder[m][i + 1];
            int u = ja * M + jobPosMachine[ja][m];
            int v = jb * M + jobPosMachine[jb][m];
            if (outcnt[u] == 0) s0[u] = v; else s1[u] = v;
            ++outcnt[u]; ++indeg[v];
        }
    }
    // Job arcs: consecutive operations along each job's route.
    for (int j = 0; j < J; ++j) {
        for (int k = 0; k + 1 < M; ++k) {
            int u = j * M + k, v = j * M + (k + 1);
            if (outcnt[u] == 0) s0[u] = v; else s1[u] = v;
            ++outcnt[u]; ++indeg[v];
        }
    }
    // Kahn topological traversal carrying earliest start times.
    int head = 0, tail = 0;
    for (int i = 0; i < n; ++i) if (indeg[i] == 0) qbuf[tail++] = i;
    int processed = 0;
    long long res = 0;
    while (head < tail) {
        int u = qbuf[head++];
        long long comp = startT[u] + procOf[u];
        if (comp > res) res = comp;
        for (int t = 0; t < outcnt[u]; ++t) {
            int v = (t == 0) ? s0[u] : s1[u];
            if (comp > startT[v]) startT[v] = comp;
            if (--indeg[v] == 0) qbuf[tail++] = v;
        }
        ++processed;
    }
    if (processed < n) return -1; // cyclic orientation -> infeasible
    return res;
}

static void copyOrders(int dst[MAXDIM][MAXDIM], const int src[MAXDIM][MAXDIM]) {
    for (int m = 0; m < M; ++m)
        for (int i = 0; i < J; ++i) dst[m][i] = src[m][i];
}

// Serial schedule generation with a selectable priority rule among eligible
// operations (whose job predecessor is already scheduled). The chosen operation
// is placed as early as possible on its machine, which always yields a feasible
// (cycle-free) machine ordering.
//   mode 0: critical path (largest remaining job work), ties by longer processing
//   mode 1: shortest processing time first
//   mode 2: longest processing time first
static void serialSGS(int mode) {
    static long long machAvail[MAXDIM], jobDone[MAXDIM];
    static int schedCount[MAXDIM];
    static bool done[MAXOP];
    for (int i = 0; i < n; ++i) done[i] = false;
    for (int m = 0; m < M; ++m) { machAvail[m] = 0; machLen[m] = 0; }
    for (int j = 0; j < J; ++j) { jobDone[j] = 0; schedCount[j] = 0; }
    int scheduled = 0;
    while (scheduled < n) {
        int best = -1;
        for (int op = 0; op < n; ++op) {
            if (done[op]) continue;
            int j = op / M, k = op % M;
            if (schedCount[j] != k) continue; // job predecessor not scheduled yet
            if (best < 0) { best = op; continue; }
            bool take = false;
            if (mode == 0) {
                if (remWork[op] != remWork[best]) take = remWork[op] > remWork[best];
                else if (procOf[op] != procOf[best]) take = procOf[op] > procOf[best];
                else take = rand() & 1;
            } else if (mode == 1) {
                if (procOf[op] != procOf[best]) take = procOf[op] < procOf[best];
                else take = rand() & 1;
            } else {
                if (procOf[op] != procOf[best]) take = procOf[op] > procOf[best];
                else take = rand() & 1;
            }
            if (take) best = op;
        }
        if (best < 0) break; // no eligible operation (should not occur for valid input)
        int j = best / M, m = machOf[best];
        long long st = machAvail[m] > jobDone[j] ? machAvail[m] : jobDone[j];
        long long en = st + procOf[best];
        machAvail[m] = en;
        jobDone[j] = en;
        done[best] = true;
        schedCount[j]++;
        machOrder[m][machLen[m]++] = j;
        ++scheduled;
    }
}

// Iterated local search for general job-shop instances. Descends from the best
// multi-start order by first-improving adjacent transpositions per machine,
// escapes local optima by random-walk perturbation, and restarts from a fresh
// schedule (cycling priority rules) when no global progress is made. Only
// feasible (cycle-free) orders are accepted; the best order found is restored.
static void optimize() {
    int modes[3] = {0, 1, 2};
    long long best = -1;
    int restartMode = 0;
    // Multi-start: keep the lowest-makespan order among several priority rules.
    for (int mi = 0; mi < 3; ++mi) {
        serialSGS(modes[mi]);
        long long ms = makespan();
        if (ms >= 0 && (best < 0 || ms < best)) {
            best = ms;
            copyOrders(bestOrder, machOrder);
        }
    }
    copyOrders(machOrder, bestOrder);
    long long cur = best;

    using namespace std::chrono;
    auto t0 = steady_clock::now();
    const long long LIMIT_US = 500000LL; // wall-clock budget for refinement
    auto over = [&]() {
        return duration_cast<microseconds>(steady_clock::now() - t0).count() > LIMIT_US;
    };
    auto nowUs = [&]() {
        return duration_cast<microseconds>(steady_clock::now() - t0).count();
    };
    long long lastImprove = nowUs();

    static int perm[MAXDIM];
    auto firstImp = [&]() -> bool {
        for (int a = 0; a < M; ++a) {
            int m = perm[a];
            if (machLen[m] < 2) continue;
            for (int i = 0; i + 1 < machLen[m]; ++i) {
                std::swap(machOrder[m][i], machOrder[m][i + 1]);
                long long nm = makespan();
                if (nm >= 0 && nm < cur) {
                    cur = nm;
                    if (cur < best) {
                        best = cur;
                        copyOrders(bestOrder, machOrder);
                        lastImprove = nowUs();
                    }
                    return true;
                }
                std::swap(machOrder[m][i], machOrder[m][i + 1]);
            }
        }
        return false;
    };

    for (;;) {
        if (over()) break;
        // Restart from a fresh schedule when no progress for a while.
        if (nowUs() - lastImprove > 120000) {
            serialSGS(modes[restartMode]);
            restartMode = (restartMode + 1) % 3;
            long long nm = makespan();
            if (nm >= 0 && (nm < best || best < 0)) { best = nm; copyOrders(bestOrder, machOrder); }
            cur = nm;
            lastImprove = nowUs();
            continue;
        }
        // Randomized machine scan order for first-improvement descent.
        for (int mm = 0; mm < M; ++mm) perm[mm] = mm;
        for (int mm = M - 1; mm > 0; --mm) {
            int r = rand() % (mm + 1);
            std::swap(perm[mm], perm[r]);
        }
        if (!firstImp()) {
            // Local optimum reached: restore the best known order, then perturb.
            if (best < cur) { copyOrders(machOrder, bestOrder); cur = best; }
            for (int it = 0; it < 32 && !over(); ++it) {
                int m = rand() % M;
                if (machLen[m] < 2) continue;
                int i = rand() % (machLen[m] - 1);
                std::swap(machOrder[m][i], machOrder[m][i + 1]);
                long long nm = makespan();
                if (nm < 0) { std::swap(machOrder[m][i], machOrder[m][i + 1]); continue; }
                cur = nm;
            }
        }
    }
    copyOrders(machOrder, bestOrder);
}

// True when every job visits machines in the natural order 0,1,...,M-1, i.e. a
// (permutation) flow shop. Permutation schedules are always feasible here.
static bool isFlowShop() {
    for (int j = 0; j < J; ++j)
        for (int k = 0; k < M; ++k)
            if (machOf[j * M + k] != k) return false;
    return true;
}

// Makespan of a permutation flow-shop schedule (same job order on every
// machine) covering the first cnt jobs. Standard O(cnt*M) longest-path DP.
static long long flowMs(const int *order, int cnt) {
    static long long prevRow[MAXDIM], curRow[MAXDIM];
    long long acc = 0;
    for (int m = 0; m < M; ++m) { acc += procOf[order[0] * M + m]; prevRow[m] = acc; }
    for (int i = 1; i < cnt; ++i) {
        curRow[0] = prevRow[0] + procOf[order[i] * M];
        for (int m = 1; m < M; ++m) {
            long long up = prevRow[m];
            long long left = curRow[m - 1];
            long long mx = up > left ? up : left;
            curRow[m] = mx + procOf[order[i] * M + m];
        }
        for (int m = 0; m < M; ++m) prevRow[m] = curRow[m];
    }
    return prevRow[M - 1];
}

// NEH heuristic: order jobs by descending total processing time, then insert
// each job at the position minimizing the flow-shop makespan.
static void neh(int *order) {
    static int idx[MAXDIM], tmp[MAXDIM];
    static long long tot[MAXDIM];
    for (int j = 0; j < J; ++j) {
        long long s = 0;
        for (int k = 0; k < M; ++k) s += procOf[j * M + k];
        tot[j] = s;
        idx[j] = j;
    }
    // Selection sort: descending total processing time.
    for (int i = 0; i < J; ++i) {
        int p = i;
        for (int l = i + 1; l < J; ++l)
            if (tot[idx[l]] > tot[idx[p]]) p = l;
        if (p != i) std::swap(idx[i], idx[p]);
    }
    int len = 1;
    order[0] = idx[0];
    for (int t = 1; t < J; ++t) {
        int j = idx[t];
        long long bestMs = -1;
        int bestPos = 0;
        for (int pos = 0; pos <= len; ++pos) {
            for (int a = 0; a < pos; ++a) tmp[a] = order[a];
            tmp[pos] = j;
            for (int a = pos; a < len; ++a) tmp[a + 1] = order[a];
            long long ms = flowMs(tmp, len + 1);
            if (bestMs < 0 || ms < bestMs) { bestMs = ms; bestPos = pos; }
        }
        for (int a = len; a > bestPos; --a) order[a] = order[a - 1];
        order[bestPos] = j;
        ++len;
    }
}

// Iterated local search over a single permutation for flow-shop instances,
// evaluated with the fast flow-shop makespan. First-improving adjacent
// transpositions plus random-walk perturbation, restoring the best found.
static void flowOptimize() {
    static int order[MAXDIM], bestPerm[MAXDIM];
    neh(order);
    if (J <= 1) {
        for (int m = 0; m < M; ++m) {
            for (int i = 0; i < J; ++i) machOrder[m][i] = order[i];
            machLen[m] = J;
        }
        return;
    }
    long long cur = flowMs(order, J);
    long long best = cur;
    for (int i = 0; i < J; ++i) bestPerm[i] = order[i];

    using namespace std::chrono;
    auto t0 = steady_clock::now();
    const long long LIMIT_US = 500000LL;
    auto over = [&]() {
        return duration_cast<microseconds>(steady_clock::now() - t0).count() > LIMIT_US;
    };

    auto firstImp = [&]() -> bool {
        for (int i = 0; i + 1 < J; ++i) {
            std::swap(order[i], order[i + 1]);
            long long nm = flowMs(order, J);
            if (nm < cur) {
                cur = nm;
                if (cur < best) {
                    best = cur;
                    for (int z = 0; z < J; ++z) bestPerm[z] = order[z];
                }
                return true;
            }
            std::swap(order[i], order[i + 1]);
        }
        return false;
    };

    while (!over()) {
        if (!firstImp()) {
            // Local optimum: restore best, then perturb with adjacent swaps.
            for (int z = 0; z < J; ++z) order[z] = bestPerm[z];
            cur = best;
            for (int it = 0; it < 60 && !over(); ++it) {
                int i = rand() % (J - 1);
                std::swap(order[i], order[i + 1]);
                cur = flowMs(order, J);
                if (cur < best) {
                    best = cur;
                    for (int z = 0; z < J; ++z) bestPerm[z] = order[z];
                }
            }
        }
    }
    for (int m = 0; m < M; ++m) {
        for (int i = 0; i < J; ++i) machOrder[m][i] = bestPerm[i];
        machLen[m] = J;
    }
}

int main() {
    if (scanf("%d %d", &J, &M) != 2) return 0;
    n = J * M;
    for (int j = 0; j < J; ++j) {
        for (int m = 0; m < M; ++m) jobPosMachine[j][m] = -1;
        for (int k = 0; k < M; ++k) {
            int mm;
            long long pp;
            if (scanf("%d %lld", &mm, &pp) != 2) return 1;
            int op = j * M + k;
            machOf[op] = (int8_t)mm;
            procOf[op] = pp;
            jobPosMachine[j][mm] = k;
        }
    }
    for (int j = 0; j < J; ++j) {
        long long acc = 0;
        for (int k = M - 1; k >= 0; --k) {
            acc += procOf[j * M + k];
            remWork[j * M + k] = acc;
        }
    }
    std::srand(20240710);
    if (isFlowShop()) flowOptimize(); else optimize();
    for (int m = 0; m < M; ++m) {
        for (int i = 0; i < J; ++i) {
            if (i) std::putchar(' ');
            std::printf("%d", machOrder[m][i]);
        }
        std::putchar('\n');
    }
    return 0;
}
