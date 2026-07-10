#include <cstdio>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <chrono>
using namespace std;

// Job Shop Scheduling solver.
//
// Constructs feasible schedules with a serial schedule generation scheme
// using priority dispatching rules (remaining-route-time and shortest-
// processing-time variants). Machine processing orders are extracted directly
// from the generated schedule, which guarantees the disjunctive graph is
// acyclic and the output is always valid.
//
// The machine orders are then refined by an adjacent-transposition local
// search (per machine, first-improvement passes) plus an iterated local
// search that perturbs the incumbent and re-optimizes. A wall-clock budget
// keeps every case well inside the time limit.

static int J, M;
static vector<vector<int>> mmat;   // mmat[j][k] : machine of the k-th operation of job j
static vector<vector<int>> pmat;   // pmat[j][k] : processing time of that operation
static vector<vector<int>> inv;    // inv[j][m]   : step index where job j visits machine m

// Scratch buffers reused across makespan evaluations to avoid allocations.
static vector<int> indeg;    // size J*M
static vector<int> pos;      // size M*J : pos[m*J + j] = position of job j on machine m
static vector<int> qbuf;     // topological queue of flattened node ids
static vector<long long> ef; // earliest finish time per operation

// Build a feasible schedule by dispatching operations in priority order
// (larger priority value is dispatched earlier). Returns the makespan and
// fills 'orders' with the resulting machine processing orders.
static long long serial_sgs(const vector<vector<long long>>& prio,
                            vector<vector<int>>& orders) {
    vector<long long> mfree(M, 0), jfin(J, 0);
    using S = pair<int, int>;
    auto cmp = [&](const S& a, const S& b) {
        return prio[a.first][a.second] < prio[b.first][b.second];
    };
    priority_queue<S, vector<S>, decltype(cmp)> pq(cmp);
    for (int j = 0; j < J; ++j) pq.push({j, 0});
    orders.assign(M, vector<int>());
    for (int m = 0; m < M; ++m) orders[m].reserve(J);
    long long makespan = 0;
    while (!pq.empty()) {
        auto [j, k] = pq.top(); pq.pop();
        int m = mmat[j][k];
        long long s = max(mfree[m], jfin[j]);
        long long f = s + pmat[j][k];
        if (f > makespan) makespan = f;
        mfree[m] = f;
        jfin[j] = f;
        orders[m].push_back(j);
        if (k + 1 < M) pq.push({j, k + 1});
    }
    return makespan;
}

// Earliest-feasible makespan for given machine orders, computed as the
// longest path in the disjunctive graph (job chain arcs + oriented machine
// arcs) via Kahn topological relaxation. Returns -1 if a cycle is detected.
static long long evaluate(const vector<vector<int>>& orders) {
    const int N = J * M;
    fill(indeg.begin(), indeg.end(), 0);
    fill(ef.begin(), ef.end(), 0);
    fill(pos.begin(), pos.end(), -1);

    for (int j = 0; j < J; ++j)
        for (int k = 1; k < M; ++k)
            indeg[j * M + k]++;                       // job-chain predecessor
    for (int m = 0; m < M; ++m)
        for (int i = 0; i < J; ++i)
            pos[m * J + orders[m][i]] = i;            // position of each job on the machine
    for (int m = 0; m < M; ++m)
        for (int i = 1; i < J; ++i) {                 // machine-arc predecessors
            int jb = orders[m][i];
            int kb = inv[jb][m];
            indeg[jb * M + kb]++;
        }

    qbuf.clear();
    int qhead = 0;
    for (int j = 0; j < J; ++j)
        for (int k = 0; k < M; ++k)
            if (!indeg[j * M + k]) {
                ef[j * M + k] = pmat[j][k];
                qbuf.push_back(j * M + k);
            }

    int done = 0;
    long long best = 0;
    while (qhead < (int)qbuf.size()) {
        int id = qbuf[qhead++];
        ++done;
        int j = id / M, k = id % M;
        long long f = ef[id];
        if (f > best) best = f;
        if (k + 1 < M) {                              // job successor
            int t = j * M + (k + 1);
            if (f > ef[t]) ef[t] = f;
            if (--indeg[t] == 0) {
                ef[t] += pmat[j][k + 1];
                qbuf.push_back(t);
            }
        }
        int m = mmat[j][k];
        int pp = pos[m * J + j];
        if (pp + 1 < J) {                             // machine successor
            int nb = orders[m][pp + 1];
            int kb2 = inv[nb][m];
            int t = nb * M + kb2;
            if (f > ef[t]) ef[t] = f;
            if (--indeg[t] == 0) {
                ef[t] += pmat[nb][kb2];
                qbuf.push_back(t);
            }
        }
    }
    return done < N ? -1 : best;
}

using Deadline = chrono::steady_clock::time_point;

static const long long INF = (long long)4e18;
static long long gBest = INF;
static vector<vector<int>> gOrders;

// First-improvement adjacent-transposition local search per machine.
// Updates the global incumbent whenever a better solution is found.
static void local_search(vector<vector<int>>& orders, Deadline deadline) {
    long long cur = evaluate(orders);
    if (cur < 0) return;
    if (cur < gBest) { gBest = cur; gOrders = orders; }
    bool improved = true;
    while (improved && chrono::steady_clock::now() < deadline) {
        improved = false;
        for (int m = 0; m < M; ++m) {
            if (chrono::steady_clock::now() >= deadline) break;
            for (int i = 0; i + 1 < J; ++i) {
                if (chrono::steady_clock::now() >= deadline) break;
                swap(orders[m][i], orders[m][i + 1]);
                long long ns = evaluate(orders);
                if (ns >= 0 && ns < cur) {
                    cur = ns;
                    improved = true;
                    if (ns < gBest) { gBest = ns; gOrders = orders; }
                } else {
                    swap(orders[m][i], orders[m][i + 1]); // revert
                }
            }
        }
    }
}

int main() {
    if (scanf("%d%d", &J, &M) != 2) return 0;
    if (J <= 0 || M <= 0) return 0;

    mmat.assign(J, vector<int>(M));
    pmat.assign(J, vector<int>(M));
    inv.assign(J, vector<int>(M, -1));
    for (int j = 0; j < J; ++j)
        for (int k = 0; k < M; ++k) {
            int m, p;
            if (scanf("%d%d", &m, &p) != 2) return 1;
            mmat[j][k] = m;
            pmat[j][k] = p;
            inv[j][m] = k;
        }

    indeg.assign((size_t)J * M, 0);
    ef.assign((size_t)J * M, 0);
    pos.assign((size_t)M * J, -1);

    Deadline start = chrono::steady_clock::now();
    Deadline deadline = start + chrono::milliseconds(700);
    mt19937 rng(20240710);

    auto run_rule = [&](const vector<vector<long long>>& prio) {
        vector<vector<int>> orders;
        serial_sgs(prio, orders);
        local_search(orders, deadline);
    };

    // Rule 1: longest remaining processing time along the route (LRPT).
    {
        vector<vector<long long>> prio(J, vector<long long>(M));
        for (int j = 0; j < J; ++j) {
            long long acc = 0;
            for (int k = M - 1; k >= 0; --k) {
                acc += pmat[j][k];
                prio[j][k] = acc;
            }
        }
        run_rule(prio);
    }
    // Rule 2: shortest processing time on the current operation (SPT).
    {
        vector<vector<long long>> prio(J, vector<long long>(M));
        for (int j = 0; j < J; ++j)
            for (int k = 0; k < M; ++k)
                prio[j][k] = -(long long)pmat[j][k];
        run_rule(prio);
    }

    // Iterated local search: perturb the incumbent with a few random
    // adjacent transpositions, then re-optimize. Only valid (acyclic)
    // perturbations are accepted.
    if (J >= 2) {
        for (int it = 0; it < 24 && chrono::steady_clock::now() < deadline; ++it) {
            vector<vector<int>> cand = gOrders;
            for (int t = 0; t < 5; ++t) {
                int m = (int)(rng() % (unsigned)M);
                int i = (int)(rng() % (unsigned)(J - 1));
                swap(cand[m][i], cand[m][i + 1]);
                if (evaluate(cand) < 0) swap(cand[m][i], cand[m][i + 1]); // revert if cyclic
            }
            local_search(cand, deadline);
        }
    }

    // Fallback: a feasible LRPT schedule if no search improved on INF.
    if (gOrders.empty()) {
        vector<vector<long long>> prio(J, vector<long long>(M));
        for (int j = 0; j < J; ++j) {
            long long acc = 0;
            for (int k = M - 1; k >= 0; --k) { acc += pmat[j][k]; prio[j][k] = acc; }
        }
        vector<vector<int>> orders;
        serial_sgs(prio, orders);
        gOrders = orders;
    }

    for (int m = 0; m < M; ++m) {
        const vector<int>& o = gOrders[m];
        for (int i = 0; i < J; ++i) {
            if (i + 1 < J) printf("%d ", o[i]);
            else printf("%d\n", o[i]);
        }
    }
    return 0;
}
