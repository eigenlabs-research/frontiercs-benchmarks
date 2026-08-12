#include <cstdio>
#include <ctime>
#include <vector>
#include <algorithm>
#include <climits>
using namespace std;
using ll = long long;

static int J, M, N;
static vector<int> machineOf, positionOnJob, orderOnMachine, slotOf, machineNext;
static vector<int> indegree, topo;
static vector<ll> duration, headTime, tailTime;

static void relink(int m, int lo, int hi) {
    lo = max(lo, 0); hi = min(hi, J - 1);
    for (int i = lo; i <= hi; ++i) {
        int u = orderOnMachine[m * J + i];
        machineNext[u] = i + 1 < J ? orderOnMachine[m * J + i + 1] : -1;
    }
}
static void relinkAll() { for (int m = 0; m < M; ++m) relink(m, 0, J - 1); }
static inline void swapAdjacent(int m, int i) {
    int a = m * J + i, b = a + 1, u = orderOnMachine[a], v = orderOnMachine[b];
    swap(orderOnMachine[a], orderOnMachine[b]);
    slotOf[u] = b; slotOf[v] = a;
    relink(m, i - 1, i + 2);
}

// Longest path in the oriented disjunctive graph. A negative result means cycle.
static ll evaluate(bool computeTails) {
    fill(indegree.begin(), indegree.end(), 0);
    for (int u = 0; u < N; ++u) {
        if (u % M + 1 < M) ++indegree[u + 1];
        if (machineNext[u] >= 0) ++indegree[machineNext[u]];
    }
    topo.clear();
    for (int u = 0; u < N; ++u) if (!indegree[u]) topo.push_back(u);
    fill(headTime.begin(), headTime.end(), 0);
    for (size_t qi = 0; qi < topo.size(); ++qi) {
        int u = topo[qi]; ll finish = headTime[u] + duration[u];
        if (u % M + 1 < M) {
            int v = u + 1; headTime[v] = max(headTime[v], finish);
            if (!--indegree[v]) topo.push_back(v);
        }
        int v = machineNext[u];
        if (v >= 0) {
            headTime[v] = max(headTime[v], finish);
            if (!--indegree[v]) topo.push_back(v);
        }
    }
    if ((int)topo.size() != N) return -1;
    ll makespan = 0;
    for (int u = 0; u < N; ++u) makespan = max(makespan, headTime[u] + duration[u]);
    if (computeTails) {
        fill(tailTime.begin(), tailTime.end(), 0);
        for (int qi = N - 1; qi >= 0; --qi) {
            int u = topo[qi]; ll t = 0;
            if (u % M + 1 < M) t = max(t, tailTime[u + 1] + duration[u + 1]);
            if (machineNext[u] >= 0) {
                int v = machineNext[u]; t = max(t, tailTime[v] + duration[v]);
            }
            tailTime[u] = t;
        }
    }
    return makespan;
}

// N5: reverse the first or last arc of each machine block on one critical path.
static void criticalMoves(ll makespan, vector<int>& moves) {
    moves.clear();
    int cur = -1;
    for (int u = 0; u < N; ++u)
        if (!headTime[u] && headTime[u] + duration[u] + tailTime[u] == makespan) { cur = u; break; }
    if (cur < 0) return;
    vector<int> path;
    while (cur >= 0) {
        path.push_back(cur); int next = -1;
        if (cur % M + 1 < M && headTime[cur] + duration[cur] == headTime[cur + 1]
            && headTime[cur + 1] + duration[cur + 1] + tailTime[cur + 1] == makespan) next = cur + 1;
        else if (machineNext[cur] >= 0) {
            int v = machineNext[cur];
            if (headTime[cur] + duration[cur] == headTime[v]
                && headTime[v] + duration[v] + tailTime[v] == makespan) next = v;
        }
        cur = next;
    }
    for (size_t i = 0; i < path.size();) {
        size_t e = i;
        while (e + 1 < path.size() && machineNext[path[e]] == path[e + 1]) ++e;
        if (e > i) {
            moves.push_back(slotOf[path[i]]);
            if (e - i >= 2) moves.push_back(slotOf[path[e - 1]]);
        }
        i = e + 1;
    }
}

static unsigned rngState = 0x9e3779b9u;
static inline unsigned random32() {
    rngState ^= rngState << 13; rngState ^= rngState >> 17; rngState ^= rngState << 5;
    return rngState;
}

// Giffler-Thompson active-schedule generator. Different modes only change the
// priority inside the conflict set, hence every returned schedule is feasible.
static vector<int> dispatch(int mode) {
    vector<int> nextOp(J), result(N), count(M);
    vector<ll> jobReady(J), machineReady(M), remaining(J);
    for (int j = 0; j < J; ++j) for (int k = 0; k < M; ++k) remaining[j] += duration[j * M + k];
    for (int done = 0; done < N; ++done) {
        ll earliestFinish = LLONG_MAX; int conflictMachine = -1;
        for (int j = 0; j < J; ++j) if (nextOp[j] < M) {
            int u = j * M + nextOp[j], m = machineOf[u];
            ll f = max(jobReady[j], machineReady[m]) + duration[u];
            if (f < earliestFinish) earliestFinish = f, conflictMachine = m;
        }
        int chosen = -1; ll bestPriority = LLONG_MIN;
        for (int j = 0; j < J; ++j) if (nextOp[j] < M) {
            int u = j * M + nextOp[j];
            if (machineOf[u] != conflictMachine || max(jobReady[j], machineReady[conflictMachine]) >= earliestFinish) continue;
            ll priority;
            if (mode == 0) priority = remaining[j];
            else if (mode == 1) priority = duration[u];
            else if (mode == 2) priority = -duration[u];
            else priority = (ll)random32();
            if (chosen < 0 || priority > bestPriority) chosen = j, bestPriority = priority;
        }
        int u = chosen * M + nextOp[chosen], m = machineOf[u];
        ll finish = max(jobReady[chosen], machineReady[m]) + duration[u];
        result[m * J + count[m]++] = u;
        jobReady[chosen] = machineReady[m] = finish;
        remaining[chosen] -= duration[u]; ++nextOp[chosen];
    }
    return result;
}

static void install(const vector<int>& candidate) {
    orderOnMachine = candidate;
    for (int i = 0; i < N; ++i) slotOf[orderOnMachine[i]] = i;
    relinkAll();
}

int main() {
    if (scanf("%d %d", &J, &M) != 2 || J <= 0 || M <= 0) return 0;
    N = J * M;
    machineOf.resize(N); duration.resize(N); positionOnJob.resize(N);
    for (int j = 0; j < J; ++j) for (int k = 0; k < M; ++k) {
        int m; ll p; if (scanf("%d %lld", &m, &p) != 2) return 0;
        machineOf[j * M + k] = m; duration[j * M + k] = p; positionOnJob[j * M + m] = k;
    }
    orderOnMachine.resize(N); slotOf.resize(N); machineNext.assign(N, -1);
    indegree.resize(N); headTime.resize(N); tailTime.resize(N); topo.reserve(N);

    // Four deterministic general seeds plus two randomized conflict priorities.
    vector<int> best;
    ll bestMakespan = LLONG_MAX;
    for (int mode = 0; mode < 6; ++mode) {
        vector<int> candidate = dispatch(mode < 3 ? mode : 3);
        install(candidate); ll c = evaluate(false);
        if (c >= 0 && c < bestMakespan) bestMakespan = c, best = candidate;
    }
    install(best);
    ll current = evaluate(true);

    if (J >= 2) {
        vector<int> tabu((size_t)N * N), moves;
        int iteration = 0, stagnant = 0, tenure = 8 + J / 4;
        clock_t start = clock();
        const clock_t limit = (clock_t)(0.78 * CLOCKS_PER_SEC);
        while (clock() - start < limit) {
            ++iteration; criticalMoves(current, moves);
            if (moves.empty()) break;
            ll chosenValue = LLONG_MAX; int chosenMove = -1;
            for (int encoded : moves) {
                int m = encoded / J, i = encoded % J;
                if (i + 1 >= J) continue;
                int u = orderOnMachine[encoded], v = orderOnMachine[encoded + 1];
                swapAdjacent(m, i); ll value = evaluate(false); swapAdjacent(m, i);
                if (value < 0) continue;
                bool forbidden = tabu[(size_t)u * N + v] > iteration;
                if (forbidden && value >= bestMakespan) continue;
                if (value < chosenValue) chosenValue = value, chosenMove = encoded;
            }
            if (chosenMove < 0) break;
            int m = chosenMove / J, i = chosenMove % J;
            int u = orderOnMachine[chosenMove], v = orderOnMachine[chosenMove + 1];
            swapAdjacent(m, i);
            tabu[(size_t)v * N + u] = iteration + tenure + (int)(random32() % 5);
            current = evaluate(true);
            if (current < bestMakespan) bestMakespan = current, best = orderOnMachine, stagnant = 0;
            else ++stagnant;

            // Restart from the incumbent and make a small feasible random walk.
            if (stagnant > 700) {
                install(best); current = evaluate(true);
                for (int tries = 0, accepted = 0; tries < 30 && accepted < 4; ++tries) {
                    int km = random32() % M, ki = random32() % (J - 1);
                    swapAdjacent(km, ki); ll value = evaluate(true);
                    if (value < 0) swapAdjacent(km, ki), current = evaluate(true);
                    else current = value, ++accepted;
                }
                stagnant = 0; tenure = 7 + (int)(random32() % 12);
                fill(tabu.begin(), tabu.end(), 0);
            }
        }
    }

    install(best);
    for (int m = 0; m < M; ++m) for (int i = 0; i < J; ++i) {
        int u = orderOnMachine[m * J + i];
        printf(i + 1 < J ? "%d " : "%d\n", u / M);
    }
    return 0;
}
