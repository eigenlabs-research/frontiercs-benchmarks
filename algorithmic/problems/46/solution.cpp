#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3,unroll-loops")
#endif
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <tuple>
#include <vector>
using namespace std;
using i64 = long long;
const i64 INF = numeric_limits<i64>::max();
struct Instance {
    int jobs = 0;
    int machines = 0;
    int operations = 0;
    vector<int> machineOf;
    vector<int> jobPred;
    vector<int> jobSucc;
    vector<i64> processing;
    vector<vector<int>> operationOnMachine;
};
struct Graph {
    const Instance* in = nullptr;
    vector<int> machinePred;
    vector<int> machineSucc;
    vector<int> topo;
    vector<i64> head;
    vector<i64> tail;
    vector<char> resolved;
    i64 makespan = -1;
    explicit Graph(const Instance& instance)
        : in(&instance),
          machinePred(instance.operations, -1),
          machineSucc(instance.operations, -1),
          head(instance.operations, 0),
          tail(instance.operations, 0),
          resolved(instance.machines, 0) {
        topo.reserve(instance.operations);
    }
    void clearMachine(int machine) {
        for (int job = 0; job < in->jobs; ++job) {
            int op = in->operationOnMachine[machine][job];
            machinePred[op] = -1;
            machineSucc[op] = -1;
        }
        resolved[machine] = 0;
    }
    void installMachine(int machine, const vector<int>& jobOrder) {
        clearMachine(machine);
        if (int(jobOrder.size()) != in->jobs) {
            return;
        }
        for (int i = 1; i < in->jobs; ++i) {
            int left = in->operationOnMachine[machine][jobOrder[i - 1]];
            int right = in->operationOnMachine[machine][jobOrder[i]];
            machineSucc[left] = right;
            machinePred[right] = left;
        }
        resolved[machine] = 1;
    }
    bool isComplete() const {
        return all_of(resolved.begin(), resolved.end(), [](char value) {
            return value != 0;
        });
    }
    bool evaluate(bool needTail) {
        const int n = in->operations;
        vector<int> indegree(n, 0);
        topo.clear();
        topo.reserve(n);
        fill(head.begin(), head.end(), 0);
        for (int op = 0; op < n; ++op) {
            if (in->jobPred[op] >= 0) {
                ++indegree[op];
            }
            if (machinePred[op] >= 0) {
                ++indegree[op];
            }
        }
        for (int op = 0; op < n; ++op) {
            if (indegree[op] == 0) {
                head[op] = in->processing[op];
                topo.push_back(op);
            }
        }
        size_t cursor = 0;
        while (cursor < topo.size()) {
            int op = topo[cursor++];
            auto relax = [&](int successor) {
                if (successor < 0) {
                    return;
                }
                head[successor] = max(
                    head[successor],
                    head[op] + in->processing[successor]
                );
                if (--indegree[successor] == 0) {
                    topo.push_back(successor);
                }
            };
            relax(in->jobSucc[op]);
            relax(machineSucc[op]);
        }
        if (int(topo.size()) != n) {
            makespan = -1;
            return false;
        }
        makespan = 0;
        for (i64 value : head) {
            makespan = max(makespan, value);
        }
        if (needTail) {
            fill(tail.begin(), tail.end(), 0);
            for (int index = n - 1; index >= 0; --index) {
                int op = topo[index];
                i64 after = 0;
                if (in->jobSucc[op] >= 0) {
                    after = max(after, tail[in->jobSucc[op]]);
                }
                if (machineSucc[op] >= 0) {
                    after = max(after, tail[machineSucc[op]]);
                }
                tail[op] = in->processing[op] + after;
            }
        }
        return true;
    }
};
struct Deadline {
    chrono::steady_clock::time_point end;
    bool expired() const {
        return chrono::steady_clock::now() >= end;
    }
};
struct MachineResult {
    vector<int> order;
    i64 objective = INF;
    bool valid = false;
};
static uint64_t splitmix64(uint64_t& state) {
    uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static i64 sequenceObjective(
    const Instance& in,
    int machine,
    const vector<int>& order,
    const vector<i64>& release,
    const vector<i64>& tailAfter
) {
    i64 time = 0;
    i64 objective = 0;
    for (int job : order) {
        int op = in.operationOnMachine[machine][job];
        time = max(time, release[op]) + in.processing[op];
        objective = max(objective, time + tailAfter[op]);
    }
    return objective;
}
static MachineResult schrageTail(
    const Instance& in,
    int machine,
    const vector<i64>& release,
    const vector<i64>& tailAfter,
    const vector<int>& tieRank,
    const Deadline& deadline
) {
    MachineResult result;
    result.objective = 0;
    vector<char> done(in.jobs, 0);
    i64 time = 0;
    result.order.reserve(in.jobs);
    while (int(result.order.size()) < in.jobs) {
        if (deadline.expired()) {
            return result;
        }
        int chosen = -1;
        i64 nextRelease = INF;
        for (int job = 0; job < in.jobs; ++job) {
            if (done[job]) {
                continue;
            }
            int op = in.operationOnMachine[machine][job];
            nextRelease = min(nextRelease, release[op]);
            if (release[op] > time) {
                continue;
            }
            if (chosen < 0) {
                chosen = job;
                continue;
            }
            int current = in.operationOnMachine[machine][chosen];
            if (tailAfter[op] > tailAfter[current] ||
                (tailAfter[op] == tailAfter[current] &&
                 in.processing[op] > in.processing[current]) ||
                (tailAfter[op] == tailAfter[current] &&
                 in.processing[op] == in.processing[current] &&
                 tieRank[op] < tieRank[current])) {
                chosen = job;
            }
        }
        if (chosen < 0) {
            if (nextRelease == INF) {
                return result;
            }
            time = max(time, nextRelease);
            continue;
        }
        int op = in.operationOnMachine[machine][chosen];
        done[chosen] = 1;
        time = max(time, release[op]) + in.processing[op];
        result.order.push_back(chosen);
        result.objective = max(result.objective, time + tailAfter[op]);
    }
    result.valid = true;
    return result;
}
static void carlierSearch(
    const Instance& in,
    int machine,
    const vector<i64>& release,
    const vector<i64>& tailAfter,
    const vector<i64>& processing,
    const vector<int>& tieRank,
    int& nodes,
    int nodeLimit,
    const Deadline& deadline,
    MachineResult& incumbent
) {
    if (nodes >= nodeLimit || deadline.expired()) {
        return;
    }
    ++nodes;
    MachineResult current = schrageTail(
        in, machine, release, tailAfter, tieRank, deadline
    );
    if (!current.valid) {
        return;
    }
    if (current.objective < incumbent.objective) {
        incumbent = current;
    }
    vector<i64> completion(in.jobs, 0);
    i64 time = 0;
    for (int index = 0; index < in.jobs; ++index) {
        int job = current.order[index];
        int op = in.operationOnMachine[machine][job];
        time = max(time, release[op]) + processing[op];
        completion[index] = time;
    }
    i64 currentObjective = 0;
    for (int index = 0; index < in.jobs; ++index) {
        int op = in.operationOnMachine[machine][current.order[index]];
        currentObjective = max(currentObjective, completion[index] + tailAfter[op]);
    }
    int b = -1;
    for (int index = in.jobs - 1; index >= 0; --index) {
        int op = in.operationOnMachine[machine][current.order[index]];
        if (completion[index] + tailAfter[op] == currentObjective) {
            b = index;
            break;
        }
    }
    if (b < 0) {
        return;
    }
    int a = -1;
    i64 blockProcessing = 0;
    int tailOp = in.operationOnMachine[machine][current.order[b]];
    for (int index = b; index >= 0; --index) {
        int op = in.operationOnMachine[machine][current.order[index]];
        blockProcessing += processing[op];
        if (release[op] + blockProcessing + tailAfter[tailOp] == currentObjective) {
            a = index;
        }
    }
    if (a < 0) {
        return;
    }
    int c = -1;
    for (int index = a; index < b; ++index) {
        int op = in.operationOnMachine[machine][current.order[index]];
        if (tailAfter[op] < tailAfter[tailOp]) {
            c = index;
        }
    }
    if (c < 0) {
        return;
    }
    i64 blockP = 0;
    i64 blockR = INF;
    i64 blockQ = INF;
    for (int index = c + 1; index <= b; ++index) {
        int op = in.operationOnMachine[machine][current.order[index]];
        blockP += processing[op];
        blockR = min(blockR, release[op]);
        blockQ = min(blockQ, tailAfter[op]);
    }
    int cOp = in.operationOnMachine[machine][current.order[c]];
    i64 lowerBound = max(
        blockR + blockP + blockQ,
        min(release[cOp], blockR) + processing[cOp] + blockP +
            min(tailAfter[cOp], blockQ)
    );
    if (lowerBound >= incumbent.objective) {
        return;
    }
    i64 raisedTail = blockP + blockQ;
    if (raisedTail > tailAfter[cOp]) {
        vector<i64> branch = tailAfter;
        branch[cOp] = raisedTail;
        carlierSearch(
            in, machine, release, branch, processing, tieRank,
            nodes, nodeLimit, deadline, incumbent
        );
    }
    if (nodes >= nodeLimit || deadline.expired()) {
        return;
    }
    vector<i64> raisedRelease = release;
    raisedRelease[cOp] = max(raisedRelease[cOp], blockR + blockP);
    if (raisedRelease[cOp] > release[cOp]) {
        carlierSearch(
            in, machine, raisedRelease, tailAfter, processing, tieRank,
            nodes, nodeLimit, deadline, incumbent
        );
    }
}
static vector<int> repairMachineOrder(
    const Graph& graph,
    const Instance& in,
    int machine,
    const vector<int>& proposed
) {
    if (int(proposed.size()) != in.jobs) {
        return {};
    }
    vector<int> rank(in.operations, in.jobs * 2);
    for (int index = 0; index < in.jobs; ++index) {
        int job = proposed[index];
        if (job < 0 || job >= in.jobs) {
            return {};
        }
        rank[in.operationOnMachine[machine][job]] = index;
    }
    vector<vector<char>> reachable(in.jobs, vector<char>(in.jobs, 0));
    for (int source = 0; source < in.jobs; ++source) {
        int sourceOp = in.operationOnMachine[machine][source];
        vector<int> stack(1, sourceOp);
        vector<char> seen(in.operations, 0);
        seen[sourceOp] = 1;
        while (!stack.empty()) {
            int op = stack.back();
            stack.pop_back();
            auto visit = [&](int successor) {
                if (successor >= 0 && !seen[successor]) {
                    seen[successor] = 1;
                    stack.push_back(successor);
                }
            };
            visit(in.jobSucc[op]);
            visit(graph.machineSucc[op]);
        }
        for (int target = 0; target < in.jobs; ++target) {
            int targetOp = in.operationOnMachine[machine][target];
            reachable[source][target] = seen[targetOp] && source != target;
        }
    }
    vector<int> indegree(in.jobs, 0);
    for (int source = 0; source < in.jobs; ++source) {
        for (int target = 0; target < in.jobs; ++target) {
            if (reachable[source][target]) {
                ++indegree[target];
            }
        }
    }
    vector<int> result;
    vector<char> used(in.jobs, 0);
    result.reserve(in.jobs);
    for (int step = 0; step < in.jobs; ++step) {
        int chosen = -1;
        for (int job = 0; job < in.jobs; ++job) {
            if (used[job] || indegree[job] != 0) {
                continue;
            }
            if (chosen < 0 || rank[in.operationOnMachine[machine][job]] <
                                  rank[in.operationOnMachine[machine][chosen]]) {
                chosen = job;
            }
        }
        if (chosen < 0) {
            return {};
        }
        used[chosen] = 1;
        result.push_back(chosen);
        for (int target = 0; target < in.jobs; ++target) {
            if (reachable[chosen][target]) {
                --indegree[target];
            }
        }
    }
    return result;
}
static MachineResult solveSingleMachine(
    Graph& graph,
    const Instance& in,
    int machine,
    uint64_t seed,
    const Deadline& deadline
) {
    MachineResult result;
    graph.clearMachine(machine);
    if (!graph.evaluate(true)) {
        return result;
    }
    vector<i64> release(in.operations, 0);
    vector<i64> tailAfter(in.operations, 0);
    for (int job = 0; job < in.jobs; ++job) {
        int op = in.operationOnMachine[machine][job];
        release[op] = graph.head[op] - in.processing[op];
        tailAfter[op] = graph.tail[op] - in.processing[op];
    }
    vector<int> tieRank(in.operations, 0);
    uint64_t randomState = seed ^ (uint64_t(machine) << 32);
    for (int job = 0; job < in.jobs; ++job) {
        int op = in.operationOnMachine[machine][job];
        tieRank[op] = int(splitmix64(randomState) & 0x7fffffff);
    }
    result = schrageTail(in, machine, release, tailAfter, tieRank, deadline);
    if (!result.valid) {
        return result;
    }
    int nodeLimit = max(32, min(512, in.jobs * in.jobs * 2));
    int nodes = 0;
    carlierSearch(
        in, machine, release, tailAfter,
        in.processing, tieRank,
        nodes, nodeLimit, deadline, result
    );
    vector<int> repaired = repairMachineOrder(graph, in, machine, result.order);
    if (repaired.empty()) {
        return MachineResult{};
    }
    result.order = std::move(repaired);
    result.objective = sequenceObjective(in, machine, result.order, release, tailAfter);
    result.valid = true;
    return result;
}
struct CandidateMachine {
    int machine = -1;
    vector<int> order;
    i64 partialMakespan = INF;
    i64 criticalWeight = 0;
    i64 load = 0;
    bool valid = false;
};
struct Schedule {
    Graph graph;
    vector<vector<int>> order;
    explicit Schedule(const Instance& in) : graph(in), order(in.machines) {}
};
struct Elite {
    Schedule schedule;
    uint64_t hash = 0;
    explicit Elite(const Instance& in) : schedule(in) {}
};
static i64 machineLoad(const Instance& in, int machine) {
    i64 load = 0;
    for (int job = 0; job < in.jobs; ++job) {
        load += in.processing[in.operationOnMachine[machine][job]];
    }
    return load;
}
static vector<int> identityOrder(int jobs) {
    vector<int> order(jobs);
    iota(order.begin(), order.end(), 0);
    return order;
}
static CandidateMachine solveUnresolvedMachine(
    Graph& partial,
    const Instance& in,
    int machine,
    uint64_t seed,
    const Deadline& deadline
) {
    CandidateMachine candidate;
    candidate.machine = machine;
    candidate.load = machineLoad(in, machine);
    vector<int> oldPred;
    vector<int> oldSucc;
    oldPred.reserve(in.jobs);
    oldSucc.reserve(in.jobs);
    for (int job = 0; job < in.jobs; ++job) {
        int op = in.operationOnMachine[machine][job];
        oldPred.push_back(partial.machinePred[op]);
        oldSucc.push_back(partial.machineSucc[op]);
    }
    char oldResolved = partial.resolved[machine];
    MachineResult oracle = solveSingleMachine(partial, in, machine, seed, deadline);
    if (oracle.valid) {
        partial.installMachine(machine, oracle.order);
        if (partial.evaluate(true)) {
            candidate.order = oracle.order;
            candidate.partialMakespan = partial.makespan;
            for (int job = 0; job < in.jobs; ++job) {
                int op = in.operationOnMachine[machine][job];
                if (partial.head[op] + partial.tail[op] -
                        in.processing[op] == partial.makespan) {
                    candidate.criticalWeight += in.processing[op];
                }
            }
            candidate.valid = true;
        }
    }
    for (int job = 0; job < in.jobs; ++job) {
        int op = in.operationOnMachine[machine][job];
        partial.machinePred[op] = oldPred[job];
        partial.machineSucc[op] = oldSucc[job];
    }
    partial.resolved[machine] = oldResolved;
    partial.evaluate(true);
    return candidate;
}
static vector<int> fallbackOrder(Graph& partial, const Instance& in, int machine) {
    partial.clearMachine(machine);
    if (!partial.evaluate(true)) {
        return {};
    }
    vector<int> proposed = identityOrder(in.jobs);
    return repairMachineOrder(partial, in, machine, proposed);
}
static vector<vector<int>> shiftingBottleneck(
    const Instance& in,
    uint64_t seed,
    const Deadline& deadline
) {
    Graph partial(in);
    vector<vector<int>> order(in.machines);
    vector<char> fixed(in.machines, 0);
    if (!partial.evaluate(true)) {
        return {};
    }
    for (int step = 0; step < in.machines; ++step) {
        if (deadline.expired()) {
            return {};
        }
        CandidateMachine best;
        for (int machine = 0; machine < in.machines; ++machine) {
            if (fixed[machine]) {
                continue;
            }
            CandidateMachine candidate = solveUnresolvedMachine(
                partial, in, machine,
                seed ^ (uint64_t(step + 1) << 40) ^
                    uint64_t(machine * 0x9e3779b9U),
                deadline
            );
            if (!candidate.valid) {
                continue;
            }
            if (!best.valid ||
                candidate.partialMakespan > best.partialMakespan ||
                (candidate.partialMakespan == best.partialMakespan &&
                 candidate.criticalWeight > best.criticalWeight) ||
                (candidate.partialMakespan == best.partialMakespan &&
                 candidate.criticalWeight == best.criticalWeight &&
                 candidate.load > best.load)) {
                best = std::move(candidate);
            }
        }
        if (!best.valid) {
            for (int machine = 0; machine < in.machines; ++machine) {
                if (fixed[machine]) {
                    continue;
                }
                vector<int> safe = fallbackOrder(partial, in, machine);
                if (!safe.empty()) {
                    best.machine = machine;
                    best.order = std::move(safe);
                    best.valid = true;
                    break;
                }
            }
        }
        if (!best.valid) {
            return {};
        }
        fixed[best.machine] = 1;
        order[best.machine] = best.order;
        partial.installMachine(best.machine, best.order);
        if (!partial.evaluate(true)) {
            return {};
        }
    }
    return partial.isComplete() ? order : vector<vector<int>>{};
}
static vector<int> orderFromGraph(
    const Graph& graph,
    const Instance& in,
    int machine
) {
    vector<int> result;
    result.reserve(in.jobs);
    int first = -1;
    for (int job = 0; job < in.jobs; ++job) {
        int op = in.operationOnMachine[machine][job];
        if (graph.machinePred[op] < 0) {
            first = op;
            break;
        }
    }
    vector<char> seen(in.jobs, 0);
    int op = first;
    while (op >= 0) {
        int job = op / in.machines;
        if (job < 0 || job >= in.jobs || seen[job]) {
            return {};
        }
        seen[job] = 1;
        result.push_back(job);
        op = graph.machineSucc[op];
    }
    return int(result.size()) == in.jobs ? result : vector<int>{};
}
static bool reoptimizeMachine(
    Graph& graph,
    const Instance& in,
    int machine,
    i64 incumbentMakespan,
    const Deadline& deadline
) {
    vector<int> oldOrder = orderFromGraph(graph, in, machine);
    if (oldOrder.empty()) {
        return false;
    }
    MachineResult candidate = solveSingleMachine(
        graph, in, machine,
        0x9e3779b97f4a7c15ULL ^ uint64_t(machine),
        deadline
    );
    if (!candidate.valid) {
        graph.installMachine(machine, oldOrder);
        graph.evaluate(true);
        return false;
    }
    graph.installMachine(machine, candidate.order);
    if (graph.evaluate(true) && graph.makespan < incumbentMakespan) {
        return true;
    }
    graph.installMachine(machine, oldOrder);
    graph.evaluate(true);
    return false;
}
static vector<int> chooseDestroySet(
    const Graph& graph,
    const Instance& in,
    uint64_t& rng,
    int round
) {
    struct Score {
        i64 critical = 0;
        i64 load = 0;
        int machine = -1;
    };
    vector<Score> scores;
    scores.reserve(in.machines);
    for (int machine = 0; machine < in.machines; ++machine) {
        i64 critical = 0;
        i64 load = machineLoad(in, machine);
        for (int job = 0; job < in.jobs; ++job) {
            int op = in.operationOnMachine[machine][job];
            if (graph.head[op] + graph.tail[op] -
                    in.processing[op] == graph.makespan) {
                critical += in.processing[op];
            }
        }
        scores.push_back({critical, load, machine});
    }
    sort(scores.begin(), scores.end(), [&](const Score& left, const Score& right) {
        if (left.critical != right.critical) {
            return left.critical > right.critical;
        }
        if (left.load != right.load) {
            return left.load > right.load;
        }
        return left.machine < right.machine;
    });
    int count = 1 + min(4, round / 3);
    count = min(count, in.machines);
    vector<int> result;
    result.reserve(count);
    for (int index = 0; index < count; ++index) {
        result.push_back(scores[index].machine);
    }
    if (count > 1 && (splitmix64(rng) & 3ULL) == 0) {
        swap(result.back(), result[(splitmix64(rng) % result.size())]);
    }
    return result;
}
static uint64_t orderHash(
    const Instance& in,
    const Graph& graph
) {
    uint64_t hash = 1469598103934665603ULL;
    for (int machine = 0; machine < in.machines; ++machine) {
        for (int job : orderFromGraph(graph, in, machine)) {
            hash ^= uint64_t(job + 1);
            hash *= 1099511628211ULL;
        }
    }
    return hash;
}
static void rememberElite(
    vector<Elite>& pool,
    const Instance& in,
    const Graph& graph
) {
    uint64_t hash = orderHash(in, graph);
    for (const Elite& elite : pool) {
        if (elite.hash == hash) {
            return;
        }
    }
    if (pool.size() >= 6) {
        pool.erase(pool.begin());
    }
    pool.emplace_back(in);
    Elite& elite = pool.back();
    elite.hash = hash;
    elite.schedule.graph = graph;
    elite.schedule.order.resize(in.machines);
    for (int machine = 0; machine < in.machines; ++machine) {
        elite.schedule.order[machine] = orderFromGraph(graph, in, machine);
    }
}
static bool topologicalRepair(
    Graph& graph,
    const Instance& in,
    vector<Elite>& pool,
    uint64_t& rng,
    const Deadline& deadline
) {
    if (!graph.isComplete() || !graph.evaluate(true) || deadline.expired()) {
        return false;
    }
    i64 oldMakespan = graph.makespan;
    vector<vector<int>> oldOrder(in.machines);
    for (int machine = 0; machine < in.machines; ++machine) {
        oldOrder[machine] = orderFromGraph(graph, in, machine);
    }
    static int roundCounter = 0;
    vector<int> removed = chooseDestroySet(graph, in, rng, roundCounter++);
    if (removed.empty()) {
        return false;
    }
    for (int machine : removed) {
        graph.clearMachine(machine);
    }
    if (!graph.evaluate(true)) {
        for (int machine = 0; machine < in.machines; ++machine) {
            graph.installMachine(machine, oldOrder[machine]);
        }
        graph.evaluate(true);
        return false;
    }
    vector<vector<int>> repaired(in.machines);
    bool foundRepair = false;
    i64 bestMakespan = oldMakespan;
    const int trials = 4;
    for (int trial = 0; trial < trials && !deadline.expired(); ++trial) {
        const int n = in.operations;
        vector<int> indegree(n, 0);
        vector<int> ready;
        ready.reserve(n);
        for (int op = 0; op < n; ++op) {
            indegree[op] = (in.jobPred[op] >= 0) +
                (graph.machinePred[op] >= 0);
            if (indegree[op] == 0) {
                ready.push_back(op);
            }
        }
        vector<i64> priority(n, 0);
        vector<i64> releaseTime(n, 0);
        vector<i64> machineReady(in.machines, 0);
        for (int op = 0; op < n; ++op) {
            i64 critical = graph.head[op] + graph.tail[op] -
                in.processing[op];
            if (trial == 0) {
                priority[op] = graph.tail[op];
            } else if (trial == 1) {
                priority[op] = critical;
            } else if (trial == 2) {
                priority[op] = -graph.head[op];
            } else {
                priority[op] = i64(splitmix64(rng));
            }
        }
        vector<vector<int>> order(in.machines);
        for (int machine : removed) {
            order[machine].reserve(in.jobs);
        }
        int popped = 0;
        while (!ready.empty() && !deadline.expired()) {
            auto finishAt = [&](int op) {
                return max(
                    releaseTime[op], machineReady[in.machineOf[op]]
                ) + in.processing[op];
            };
            i64 earliest = INF;
            int selectedMachine = -1;
            int firstSlot = -1;
            for (int slot = 0; slot < int(ready.size());
                 ++slot) {
                int op = ready[slot];
                i64 finish = finishAt(op);
                if (finish < earliest ||
                    (finish == earliest &&
                     (firstSlot < 0 || priority[op] >
                          priority[ready[firstSlot]]))) {
                    earliest = finish;
                    selectedMachine = in.machineOf[op];
                    firstSlot = slot;
                }
            }
            int conflictSlot = -1;
            for (int slot = 0; slot < int(ready.size());
                 ++slot) {
                int op = ready[slot];
                if (in.machineOf[op] != selectedMachine ||
                    releaseTime[op] >= earliest) {
                    continue;
                }
                if (conflictSlot < 0 || priority[op] >
                        priority[ready[conflictSlot]] ||
                    (priority[op] == priority[ready[conflictSlot]] &&
                     finishAt(op) < finishAt(ready[conflictSlot]))) {
                    conflictSlot = slot;
                }
            }
            int bestSlot = conflictSlot >= 0 ? conflictSlot : firstSlot;
            int op = ready[bestSlot];
            ready[bestSlot] = ready.back();
            ready.pop_back();
            ++popped;
            i64 finish = max(
                releaseTime[op], machineReady[in.machineOf[op]]
            ) + in.processing[op];
            machineReady[in.machineOf[op]] = finish;
            int machine = in.machineOf[op];
            if (find(removed.begin(), removed.end(), machine) !=
                    removed.end()) {
                order[machine].push_back(op / in.machines);
            }
            auto relax = [&](int successor) {
                if (successor >= 0) {
                    releaseTime[successor] = max(
                        releaseTime[successor], finish
                    );
                }
                if (successor >= 0 && --indegree[successor] == 0) {
                    ready.push_back(successor);
                }
            };
            relax(in.jobSucc[op]);
            relax(graph.machineSucc[op]);
        }
        if (popped != n) {
            continue;
        }
        Graph candidate = graph;
        bool complete = true;
        for (int machine : removed) {
            if (int(order[machine].size()) != in.jobs) {
                complete = false;
                break;
            }
            candidate.installMachine(machine, order[machine]);
        }
        if (complete && candidate.evaluate(true) &&
            candidate.makespan < bestMakespan) {
            bestMakespan = candidate.makespan;
            repaired = std::move(order);
            foundRepair = true;
        }
    }
    if (!foundRepair) {
        for (int machine = 0; machine < in.machines; ++machine) {
            graph.installMachine(machine, oldOrder[machine]);
        }
        graph.evaluate(true);
        return false;
    }
    for (int machine : removed) {
        graph.installMachine(machine, repaired[machine]);
    }
    graph.evaluate(true);
    if (graph.makespan < oldMakespan) {
        rememberElite(pool, in, graph);
        return true;
    }
    for (int machine = 0; machine < in.machines; ++machine) {
        graph.installMachine(machine, oldOrder[machine]);
    }
    graph.evaluate(true);
    return false;
}
static vector<vector<int>> criticalMoves(
    const Graph& graph,
    const Instance& in
) {
    vector<vector<int>> current(in.machines);
    for (int machine = 0; machine < in.machines; ++machine) {
        current[machine] = orderFromGraph(graph, in, machine);
    }
    return current;
}
static bool polish(
  Graph& graph,
  const Instance& in,
  const Deadline& deadline
) {
  if (!graph.isComplete() || !graph.evaluate(true)) {
    return false;
  }
  bool improvedAny = false;
  for (int pass = 0; pass < 24 && !deadline.expired(); ++pass) {
    vector<vector<int>> current = criticalMoves(graph, in);
    i64 base = graph.makespan;
    i64 best = base;
    int bestMachine = -1;
    vector<int> bestOrder;
    auto boundaryPolish = [&]() {
      i64 boundaryBest = base;
      int boundaryMachine = -1;
      vector<int> boundaryOrder;
      for (int machine = 0; machine < in.machines; ++machine) {
        const vector<int>& sequence = current[machine];
        int begin = 0;
        while (begin < in.jobs) {
          int beginOp = in.operationOnMachine[machine][sequence[begin]];
          if (graph.head[beginOp] + graph.tail[beginOp] -
              in.processing[beginOp] != graph.makespan) {
            ++begin;
            continue;
          }
          int end = begin;
          while (end + 1 < in.jobs) {
            int op = in.operationOnMachine[machine][sequence[end + 1]];
            if (graph.head[op] + graph.tail[op] -
                in.processing[op] != graph.makespan) {
              break;
            }
            ++end;
          }
          if (end > begin) {
            if (begin > 0) {
              vector<int> boundaryCandidate = sequence;
              int job = boundaryCandidate[begin];
              boundaryCandidate.erase(boundaryCandidate.begin() + begin);
              boundaryCandidate.insert(
                boundaryCandidate.begin() + begin - 1, job
              );
              graph.installMachine(machine, boundaryCandidate);
              if (graph.evaluate(false) &&
                  graph.makespan < boundaryBest) {
                boundaryBest = graph.makespan;
                boundaryMachine = machine;
                boundaryOrder = std::move(boundaryCandidate);
              }
              graph.installMachine(machine, sequence);
              graph.evaluate(true);
            }
            if (end + 1 < in.jobs) {
              vector<int> boundaryCandidate = sequence;
              int job = boundaryCandidate[end];
              boundaryCandidate.erase(boundaryCandidate.begin() + end);
              boundaryCandidate.insert(
                boundaryCandidate.begin() + end + 1, job
              );
              graph.installMachine(machine, boundaryCandidate);
              if (graph.evaluate(false) &&
                  graph.makespan < boundaryBest) {
                boundaryBest = graph.makespan;
                boundaryMachine = machine;
                boundaryOrder = std::move(boundaryCandidate);
              }
              graph.installMachine(machine, sequence);
              graph.evaluate(true);
            }
          }
          begin = end + 1;
        }
      }
      if (boundaryMachine < 0) {
        return false;
      }
      graph.installMachine(boundaryMachine, boundaryOrder);
      return graph.evaluate(true);
    };
    for (int machine = 0; machine < in.machines; ++machine) {
      const vector<int>& sequence = current[machine];
      int begin = 0;
      while (begin < in.jobs) {
        int beginOp = in.operationOnMachine[machine][sequence[begin]];
        bool beginCritical =
          graph.head[beginOp] + graph.tail[beginOp] -
            in.processing[beginOp] == graph.makespan;
        if (!beginCritical) {
          ++begin;
          continue;
        }
        int end = begin;
        while (end + 1 < in.jobs) {
          int op = in.operationOnMachine[machine][sequence[end + 1]];
          if (graph.head[op] + graph.tail[op] -
              in.processing[op] != graph.makespan) {
            break;
          }
          ++end;
        }
        if (end > begin) {
          vector<vector<int>> candidates;
          vector<int> swapped = sequence;
          swap(swapped[begin], swapped[begin + 1]);
          candidates.push_back(std::move(swapped));
          swapped = sequence;
          swap(swapped[end - 1], swapped[end]);
          candidates.push_back(std::move(swapped));
          for (int target = begin + 1; target <= end; ++target) {
            vector<int> moved = sequence;
            int job = moved[begin];
            moved.erase(moved.begin() + begin);
            moved.insert(moved.begin() + target, job);
            candidates.push_back(std::move(moved));
          }
          for (int target = begin; target < end; ++target) {
            vector<int> moved = sequence;
            int job = moved[end];
            moved.erase(moved.begin() + end);
            moved.insert(moved.begin() + target, job);
            candidates.push_back(std::move(moved));
          }
          for (const vector<int>& candidate : candidates) {
            graph.installMachine(machine, candidate);
            if (graph.evaluate(false) && graph.makespan < best) {
              best = graph.makespan;
              bestMachine = machine;
              bestOrder = candidate;
            }
          }
          graph.installMachine(machine, sequence);
          graph.evaluate(true);
        }
        begin = end + 1;
      }
    }
    if (bestMachine < 0) {
      if (!boundaryPolish()) {
        break;
      }
      improvedAny = true;
      continue;
    }
    graph.installMachine(bestMachine, bestOrder);
    if (!graph.evaluate(true)) {
      break;
    }
    improvedAny = true;
  }
  return improvedAny;
}
struct TabuMove {
    int machine = -1;
    int from = -1;
    int to = -1;
    int job = -1;
    int anchor = -1;
    i64 estimate = INF;
    int serial = 0;
};
static vector<int> movedOrder(
    const vector<int>& sequence,
    int from,
    int to
) {
    vector<int> result = sequence;
    int job = result[from];
    result.erase(result.begin() + from);
    result.insert(result.begin() + to, job);
    return result;
}
static i64 estimateMove(
    const Graph& graph,
    const Instance& in,
    int machine,
    const vector<int>& sequence,
    int from,
    int to
) {
    int lo = min(from, to);
    int hi = max(from, to);
    static vector<int> segment;
    static vector<i64> completion;
    segment.clear();
    segment.reserve(in.jobs);
    if (from < to) {
        segment.push_back(sequence[from]);
        for (int index = from + 1; index <= to; ++index) {
            segment.push_back(sequence[index]);
        }
        rotate(segment.begin(), segment.begin() + 1, segment.end());
    } else {
        segment.push_back(sequence[from]);
        for (int index = to; index < from; ++index) {
            segment.push_back(sequence[index]);
        }
    }
    completion.resize(segment.size());
    i64 previous = lo > 0
        ? graph.head[in.operationOnMachine[machine][sequence[lo - 1]]]
        : 0;
    for (size_t index = 0; index < segment.size(); ++index) {
        int op = in.operationOnMachine[machine][segment[index]];
        i64 release = in.jobPred[op] >= 0 ? graph.head[in.jobPred[op]] : 0;
        previous = max(previous, release) + in.processing[op];
        completion[index] = previous;
    }
    i64 nextTail = hi + 1 < in.jobs
        ? graph.tail[in.operationOnMachine[machine][sequence[hi + 1]]]
        : 0;
    i64 estimate = 0;
    for (int index = int(segment.size()) - 1; index >= 0; --index) {
        int op = in.operationOnMachine[machine][segment[index]];
        i64 successorTail = in.jobSucc[op] >= 0 ? graph.tail[in.jobSucc[op]] : 0;
        i64 after = max(nextTail, successorTail);
        i64 length = completion[index] - in.processing[op] +
                     in.processing[op] + after;
        estimate = max(estimate, length);
        nextTail = in.processing[op] + after;
    }
    return estimate;
}
static i64 incrementalAfterMove(
    Graph& graph,
    const Instance& in,
    const vector<vector<int>>& order,
    const TabuMove& move
) {
    int machine = move.machine;
    const vector<int>& sequence = order[machine];
    int lo = min(move.from, move.to);
    int hi = max(move.from, move.to);
    int first = max(0, lo - 1);
    int last = min(in.jobs - 1, hi + 1);
    for (int index = first; index <= last; ++index) {
        int op = in.operationOnMachine[machine][sequence[index]];
        graph.machinePred[op] = index > 0
            ? in.operationOnMachine[machine][sequence[index - 1]] : -1;
        graph.machineSucc[op] = index + 1 < in.jobs
            ? in.operationOnMachine[machine][sequence[index + 1]] : -1;
    }
    static vector<int> work;
    static vector<char> queued;
    if (int(queued.size()) != in.operations) {
        queued.assign(in.operations, 0);
    }
    work.clear();
    queued.assign(in.operations, 0);
    auto pushForward = [&](int op) {
        if (op >= 0 && !queued[op]) {
            queued[op] = 1;
            work.push_back(op);
        }
    };
    for (int index = lo; index <= min(in.jobs - 1, hi + 1); ++index) {
        pushForward(in.operationOnMachine[machine][sequence[index]]);
    }
    int pops = 0;
    size_t cursor = 0;
    const int limit = max(64, 16 * in.operations);
    while (cursor < work.size()) {
        int op = work[cursor++];
        queued[op] = 0;
        if (++pops > limit) {
            for (int pending : work) {
                queued[pending] = 0;
            }
            return -2;
        }
        i64 before = in.jobPred[op] >= 0 ? graph.head[in.jobPred[op]] : 0;
        if (graph.machinePred[op] >= 0) {
            before = max(before, graph.head[graph.machinePred[op]]);
        }
        i64 finish = before + in.processing[op];
        if (finish == graph.head[op]) {
            continue;
        }
        graph.head[op] = finish;
        pushForward(in.jobSucc[op]);
        pushForward(graph.machineSucc[op]);
    }
    work.clear();
    auto pushBackward = [&](int op) {
        if (op >= 0 && !queued[op]) {
            queued[op] = 1;
            work.push_back(op);
        }
    };
    for (int index = hi; index >= max(0, lo - 1); --index) {
        pushBackward(in.operationOnMachine[machine][sequence[index]]);
    }
    cursor = 0;
    pops = 0;
    while (cursor < work.size()) {
        int op = work[cursor++];
        queued[op] = 0;
        if (++pops > limit) {
            for (int pending : work) {
                queued[pending] = 0;
            }
            return -2;
        }
        i64 after = in.jobSucc[op] >= 0 ? graph.tail[in.jobSucc[op]] : 0;
        if (graph.machineSucc[op] >= 0) {
            after = max(after, graph.tail[graph.machineSucc[op]]);
        }
        i64 length = in.processing[op] + after;
        if (length == graph.tail[op]) {
            continue;
        }
        graph.tail[op] = length;
        pushBackward(in.jobPred[op]);
        pushBackward(graph.machinePred[op]);
    }
    graph.makespan = 0;
    for (i64 finish : graph.head) {
        graph.makespan = max(graph.makespan, finish);
    }
    return graph.makespan;
}
static bool tabuPolish(
  Graph& graph,
  const Instance& in,
  const Deadline& deadline
) {
  if (!graph.isComplete() || !graph.evaluate(true)) {
    return false;
  }
  vector<vector<int>> current(in.machines);
  for (int machine = 0; machine < in.machines; ++machine) {
    current[machine] = orderFromGraph(graph, in, machine);
  }
  vector<vector<int>> best = current;
  i64 bestMakespan = graph.makespan;
  vector<int> tabuJob(size_t(in.machines) * in.jobs, 0);
  int iteration = 0;
  auto lastImprove = chrono::steady_clock::now();
  auto moveIsTabu = [&](const TabuMove& move) {
    return tabuJob[size_t(move.machine) * in.jobs + move.job] > iteration;
  };
  auto markMove = [&](const TabuMove& move, int until) {
    tabuJob[size_t(move.machine) * in.jobs + move.job] = until;
  };
  vector<TabuMove> moves;
  moves.reserve(size_t(in.machines) * in.jobs * 2);
  uint64_t tabuRng = 0x9e3779b97f4a7c15ULL;
  auto nextTabuRandom = [&]() {
    tabuRng ^= tabuRng << 13;
    tabuRng ^= tabuRng >> 7;
    tabuRng ^= tabuRng << 17;
    return tabuRng;
  };
  auto dynamicTenure = [&]() {
    i64 stagnantMs = chrono::duration_cast<chrono::milliseconds>(
      chrono::steady_clock::now() - lastImprove
    ).count();
    int minimum = 8;
    int span = max(4, in.jobs / 3);
    if (stagnantMs >= 120) {
      minimum = 15 + int(min<i64>(6, (stagnantMs - 120) / 250));
      span = max(4, in.jobs / 2);
    } else if (stagnantMs > 30) {
      int fraction = int((stagnantMs - 30) * 100 / 90);
      minimum = 8 + (15 - 8) * fraction / 100;
      span = max(4, in.jobs / 3) +
        (max(4, in.jobs / 2) - max(4, in.jobs / 3)) * fraction / 100;
    }
    return minimum + int(nextTabuRandom() % uint64_t(span));
  };
  auto restart = [&]() {
    current = best;
    for (int kick = 0; kick < min(3, in.machines); ++kick) {
      int machine = (iteration + kick * 7) % in.machines;
      if (current[machine].size() < 2) {
        continue;
      }
      int from = (iteration + kick * 5) % current[machine].size();
      int to = (from + 1 + (iteration %
        max(1, int(current[machine].size()) - 1))) %
        current[machine].size();
      current[machine] = movedOrder(current[machine], from, to);
    }
    for (int machine = 0; machine < in.machines; ++machine) {
      graph.installMachine(machine, current[machine]);
    }
    if (!graph.evaluate(true)) {
      current = best;
      for (int machine = 0; machine < in.machines; ++machine) {
        graph.installMachine(machine, current[machine]);
      }
      graph.evaluate(true);
    }
    fill(tabuJob.begin(), tabuJob.end(), 0);
    lastImprove = chrono::steady_clock::now();
  };
  while (!deadline.expired() && iteration < 1000000) {
    ++iteration;
    moves.clear();
    for (int machine = 0; machine < in.machines; ++machine) {
      const vector<int>& sequence = current[machine];
      int begin = 0;
      while (begin < in.jobs) {
        int op = in.operationOnMachine[machine][sequence[begin]];
        if (graph.head[op] + graph.tail[op] -
            in.processing[op] != graph.makespan) {
          ++begin;
          continue;
        }
        int end = begin;
        while (end + 1 < in.jobs) {
          int nextOp = in.operationOnMachine[machine][sequence[end + 1]];
          if (graph.head[nextOp] + graph.tail[nextOp] -
              in.processing[nextOp] != graph.makespan) {
            break;
          }
          ++end;
        }
        if (end > begin) {
          for (int target = begin + 1; target <= end; ++target) {
            moves.push_back({
              machine, target, begin, sequence[target], sequence[begin],
              estimateMove(graph, in, machine, sequence, target, begin)
            });
          }
          for (int target = begin; target < end; ++target) {
            if (target == begin && end == begin + 1) {
              continue;
            }
            moves.push_back({
              machine, target, end, sequence[target], sequence[end],
              estimateMove(graph, in, machine, sequence, target, end)
            });
          }
        }
        begin = end + 1;
      }
    }
    if (moves.empty()) {
      restart();
      continue;
    }
    for (size_t index = 0; index < moves.size(); ++index) {
      moves[index].serial = int(index);
    }
    auto moveLess = [](const TabuMove& left, const TabuMove& right) {
      if (left.estimate != right.estimate) {
        return left.estimate < right.estimate;
      }
      return left.serial < right.serial;
    };
    size_t shortlist = min<size_t>(moves.size(), 24);
    if (shortlist == moves.size()) {
      sort(moves.begin(), moves.end(), moveLess);
    } else {
      partial_sort(moves.begin(), moves.begin() + shortlist, moves.end(), moveLess);
    }
    TabuMove selected;
    bool found = false;
    size_t limit = shortlist;
    for (int pass = 0; pass < 2 && !found; ++pass) {
      for (size_t index = 0; index < limit; ++index) {
        const TabuMove& move = moves[index];
        bool isTabu = moveIsTabu(move);
        if (pass == 0 && isTabu && move.estimate >= bestMakespan) {
          continue;
        }
        selected = move;
        found = true;
        break;
      }
    }
    if (!found) {
      restart();
      continue;
    }
    vector<int> oldOrder = current[selected.machine];
    current[selected.machine] = movedOrder(
      current[selected.machine], selected.from, selected.to
    );
    i64 candidateMakespan = incrementalAfterMove(graph, in, current, selected);
    if (candidateMakespan == -2) {
      candidateMakespan = graph.evaluate(true) ? graph.makespan : -1;
    }
    if (candidateMakespan < 0) {
      current[selected.machine] = oldOrder;
      graph.installMachine(selected.machine, oldOrder);
      graph.evaluate(true);
      continue;
    }
    markMove(selected, iteration + dynamicTenure());
    graph.makespan = candidateMakespan;
    if (candidateMakespan < bestMakespan) {
      bestMakespan = candidateMakespan;
      best = current;
      lastImprove = chrono::steady_clock::now();
    }
    i64 stagnantMs = chrono::duration_cast<chrono::milliseconds>(
      chrono::steady_clock::now() - lastImprove
    ).count();
    if (stagnantMs > 400) {
      restart();
    }
  }
  for (int machine = 0; machine < in.machines; ++machine) {
    graph.installMachine(machine, best[machine]);
  }
  graph.evaluate(true);
  return graph.makespan == bestMakespan;
}
static vector<vector<int>> gifflerThompson(
    const Instance& in,
    int mode,
    uint64_t& rng
) {
    vector<int> next(in.jobs, 0);
    vector<i64> jobReady(in.jobs, 0);
    vector<i64> machineReady(in.machines, 0);
    vector<i64> remaining(in.jobs, 0);
    vector<vector<int>> order(in.machines);
    for (int job = 0; job < in.jobs; ++job) {
        for (int k = 0; k < in.machines; ++k) {
            remaining[job] += in.processing[job * in.machines + k];
        }
    }
    for (int scheduled = 0; scheduled < in.operations; ++scheduled) {
        i64 earliestFinish = INF;
        int selectedMachine = -1;
        for (int job = 0; job < in.jobs; ++job) {
            if (next[job] >= in.machines) {
                continue;
            }
            int op = job * in.machines + next[job];
            int machine = in.machineOf[op];
            i64 start = max(jobReady[job], machineReady[machine]);
            i64 finish = start + in.processing[op];
            if (finish < earliestFinish) {
                earliestFinish = finish;
                selectedMachine = machine;
            }
        }
        int chosenJob = -1;
        i64 chosenPriority = -INF;
        for (int job = 0; job < in.jobs; ++job) {
            if (next[job] >= in.machines) {
                continue;
            }
            int op = job * in.machines + next[job];
            int machine = in.machineOf[op];
            if (machine != selectedMachine) {
                continue;
            }
            i64 start = max(jobReady[job], machineReady[machine]);
            if (start >= earliestFinish) {
                continue;
            }
            i64 priority = 0;
            if (mode == 0) {
                priority = remaining[job];
            } else if (mode == 1) {
                priority = in.processing[op];
            } else if (mode == 2) {
                priority = -in.processing[op];
            } else {
                priority = i64(splitmix64(rng) & 0x7fffffffffffffffULL);
            }
            if (chosenJob < 0 || priority > chosenPriority) {
                chosenJob = job;
                chosenPriority = priority;
            }
        }
        if (chosenJob < 0) {
            for (int job = 0; job < in.jobs; ++job) {
                if (next[job] < in.machines &&
                    in.machineOf[job * in.machines + next[job]] == selectedMachine) {
                    chosenJob = job;
                    break;
                }
            }
        }
        if (chosenJob < 0) {
            return {};
        }
        int op = chosenJob * in.machines + next[chosenJob];
        int machine = in.machineOf[op];
        i64 finish = max(jobReady[chosenJob], machineReady[machine]) +
                     in.processing[op];
        order[machine].push_back(chosenJob);
        jobReady[chosenJob] = finish;
        machineReady[machine] = finish;
        remaining[chosenJob] -= in.processing[op];
        ++next[chosenJob];
    }
    return order;
}
static i64 evaluateJobPermutation(
    const Instance& in,
    const vector<int>& permutation,
    int length
) {
    vector<i64> machineReady(in.machines, 0);
    i64 makespan = 0;
    for (int index = 0; index < length; ++index) {
        int job = permutation[index];
        i64 jobReady = 0;
        for (int k = 0; k < in.machines; ++k) {
            int op = job * in.machines + k;
            int machine = in.machineOf[op];
            i64 start = max(jobReady, machineReady[machine]);
            jobReady = start + in.processing[op];
            machineReady[machine] = jobReady;
        }
        makespan = max(makespan, jobReady);
    }
    return makespan;
}
static vector<int> nehBuild(
    const Instance& in,
    const vector<int>& inputOrder
) {
    vector<int> result;
    result.reserve(in.jobs);
    for (int job : inputOrder) {
        i64 bestMakespan = INF;
        int bestPosition = 0;
        for (int position = 0; position <= int(result.size()); ++position) {
            vector<int> candidate = result;
            candidate.insert(candidate.begin() + position, job);
            i64 makespan = evaluateJobPermutation(
                in, candidate, int(candidate.size())
            );
            if (makespan < bestMakespan) {
                bestMakespan = makespan;
                bestPosition = position;
            }
        }
        result.insert(result.begin() + bestPosition, job);
    }
    return result;
}
static void improveJobPermutation(
    const Instance& in,
    vector<int>& permutation,
    const Deadline& deadline
) {
    i64 current = evaluateJobPermutation(in, permutation, in.jobs);
    for (int pass = 0; pass < 8 && !deadline.expired(); ++pass) {
        bool improved = false;
        for (int from = 0; from < in.jobs && !deadline.expired(); ++from) {
            int job = permutation[from];
            vector<int> without = permutation;
            without.erase(without.begin() + from);
            i64 best = current;
            int bestPosition = from;
            for (int position = 0; position <= in.jobs - 1; ++position) {
                vector<int> candidate = without;
                candidate.insert(candidate.begin() + position, job);
                i64 makespan = evaluateJobPermutation(in, candidate, in.jobs);
                if (makespan < best) {
                    best = makespan;
                    bestPosition = position;
                }
            }
            if (best < current) {
                permutation = without;
                permutation.insert(permutation.begin() + bestPosition, job);
                current = best;
                improved = true;
            }
        }
        if (!improved) {
            break;
        }
    }
}
static Instance readInstance() {
    Instance in;
    if (!(cin >> in.jobs >> in.machines)) {
        return in;
    }
    in.operations = in.jobs * in.machines;
    in.machineOf.assign(in.operations, -1);
    in.jobPred.assign(in.operations, -1);
    in.jobSucc.assign(in.operations, -1);
    in.processing.assign(in.operations, 0);
    in.operationOnMachine.assign(in.machines, vector<int>(in.jobs, -1));
    for (int job = 0; job < in.jobs; ++job) {
        for (int k = 0; k < in.machines; ++k) {
            int machine;
            i64 processing;
            cin >> machine >> processing;
            int op = job * in.machines + k;
            in.machineOf[op] = machine;
            in.processing[op] = processing;
            in.operationOnMachine[machine][job] = op;
            if (k > 0) {
                in.jobPred[op] = op - 1;
            }
            if (k + 1 < in.machines) {
                in.jobSucc[op] = op + 1;
            }
        }
    }
    return in;
}
static vector<vector<int>> dispatchSeed(const Instance& in, int mode) {
    vector<vector<int>> order(in.machines);
    vector<int> next(in.jobs, 0);
    vector<i64> jobFree(in.jobs, 0);
    vector<i64> machineFree(in.machines, 0);
    vector<i64> remaining(in.jobs, 0);
    for (int job = 0; job < in.jobs; ++job) {
        for (int k = 0; k < in.machines; ++k) {
            remaining[job] += in.processing[job * in.machines + k];
        }
    }
    for (int scheduled = 0; scheduled < in.operations; ++scheduled) {
        int chosen = -1;
        i64 chosenStart = 0;
        i64 chosenFinish = 0;
        i64 chosenProcessing = 0;
        i64 chosenRemaining = 0;
        for (int job = 0; job < in.jobs; ++job) {
            if (next[job] >= in.machines) {
                continue;
            }
            int op = job * in.machines + next[job];
            int machine = in.machineOf[op];
            i64 start = max(jobFree[job], machineFree[machine]);
            i64 finish = start + in.processing[op];
            i64 processing = in.processing[op];
            bool better = chosen < 0 || finish < chosenFinish;
            if (!better && finish == chosenFinish) {
                if (mode == 1) {
                    better = processing > chosenProcessing;
                } else if (mode == 2) {
                    better = processing < chosenProcessing;
                } else if (mode == 3) {
                    better = remaining[job] > chosenRemaining;
                } else if (mode == 4) {
                    better = start < chosenStart;
                } else if (mode == 5) {
                    better = ((job * 1315423911u + machine * 2654435761u) &
                              1u) != 0;
                } else {
                    better = remaining[job] < chosenRemaining;
                }
            }
            if (better) {
                chosen = job;
                chosenStart = start;
                chosenFinish = finish;
                chosenProcessing = processing;
                chosenRemaining = remaining[job];
            }
        }
        if (chosen < 0) {
            return {};
        }
        int op = chosen * in.machines + next[chosen];
        int machine = in.machineOf[op];
        order[machine].push_back(chosen);
        jobFree[chosen] = chosenFinish;
        machineFree[machine] = chosenFinish;
        remaining[chosen] -= in.processing[op];
        ++next[chosen];
    }
    return order;
}
static void emit(const Instance& in, const vector<vector<int>>& order) {
    for (int machine = 0; machine < in.machines; ++machine) {
        for (int i = 0; i < in.jobs; ++i) {
            if (i) {
                cout << ' ';
            }
            cout << order[machine][i];
        }
        cout << '\n';
    }
}
int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);
  Instance in = readInstance();
  if (in.jobs <= 0 || in.machines <= 0) {
    return 0;
  }
  auto start = chrono::steady_clock::now();
  Deadline deadline{start + chrono::milliseconds(955)};
  Deadline seedDeadline{start + chrono::milliseconds(180)};
  Deadline constructionDeadline{start + chrono::milliseconds(300)};
  vector<vector<int>> best;
  i64 bestMakespan = INF;
  struct SeedCandidate {
    i64 makespan = INF;
    uint64_t hash = 0;
    vector<vector<int>> order;
  };
  vector<SeedCandidate> seedCandidates;
  auto consider = [&](const vector<vector<int>>& candidate) {
    if (int(candidate.size()) != in.machines) {
      return;
    }
    Graph graph(in);
    for (int machine = 0; machine < in.machines; ++machine) {
      graph.installMachine(machine, candidate[machine]);
    }
    if (graph.evaluate(false) && graph.makespan < bestMakespan) {
      bestMakespan = graph.makespan;
      best = candidate;
    }
    if (graph.makespan >= 0) {
      uint64_t hash = 1469598103934665603ULL;
      for (const vector<int>& sequence : candidate) {
        for (int job : sequence) {
          hash ^= uint64_t(job + 1);
          hash *= 1099511628211ULL;
        }
      }
      bool duplicate = false;
      for (const SeedCandidate& existing : seedCandidates) {
        if (existing.hash == hash) {
          duplicate = true;
          break;
        }
      }
      if (!duplicate) {
        seedCandidates.push_back({graph.makespan, hash, candidate});
        sort(seedCandidates.begin(), seedCandidates.end(),
          [](const SeedCandidate& left, const SeedCandidate& right) {
            return left.makespan < right.makespan;
          });
        if (seedCandidates.size() > 8) {
          seedCandidates.pop_back();
        }
      }
    }
  };
  for (int mode = 0; mode < 6; ++mode) {
    consider(dispatchSeed(in, mode));
  }
  uint64_t seedRng = 0x9e3779b97f4a7c15ULL;
  for (int mode = 0; mode < 4; ++mode) {
    consider(gifflerThompson(in, mode, seedRng));
  }
  vector<i64> totalWork(in.jobs, 0);
  vector<i64> frontWork(in.jobs, 0);
  vector<i64> maxOperation(in.jobs, 0);
  for (int job = 0; job < in.jobs; ++job) {
    for (int k = 0; k < in.machines; ++k) {
      i64 processing = in.processing[job * in.machines + k];
      totalWork[job] += processing;
      maxOperation[job] = max(maxOperation[job], processing);
      if (k < in.machines / 2) {
        frontWork[job] += processing;
      }
    }
  }
  vector<int> baseJobs(in.jobs);
  iota(baseJobs.begin(), baseJobs.end(), 0);
  for (int variant = 0; variant < 4; ++variant) {
    sort(baseJobs.begin(), baseJobs.end(), [&](int left, int right) {
      i64 leftKey = variant == 0 ? totalWork[left]
        : variant == 1 ? -totalWork[left]
        : variant == 2 ? frontWork[left] : maxOperation[left];
      i64 rightKey = variant == 0 ? totalWork[right]
        : variant == 1 ? -totalWork[right]
        : variant == 2 ? frontWork[right] : maxOperation[right];
      if (leftKey != rightKey) {
        return leftKey > rightKey;
      }
      return left < right;
    });
    vector<int> permutation = nehBuild(in, baseJobs);
    improveJobPermutation(in, permutation, seedDeadline);
    consider(vector<vector<int>>(in.machines, permutation));
    reverse(permutation.begin(), permutation.end());
    consider(vector<vector<int>>(in.machines, permutation));
  }
  const uint64_t seeds[] = {
    0x9e3779b97f4a7c15ULL,
    0xd1b54a32d192ed03ULL,
    0x94d049bb133111ebULL,
  };
  for (uint64_t seed : seeds) {
    if (constructionDeadline.expired()) {
      break;
    }
    consider(shiftingBottleneck(in, seed, constructionDeadline));
  }
  if (best.empty()) {
    best.assign(in.machines, vector<int>(in.jobs));
    for (int machine = 0; machine < in.machines; ++machine) {
      iota(best[machine].begin(), best[machine].end(), 0);
    }
    Graph fallback(in);
    for (int machine = 0; machine < in.machines; ++machine) {
      fallback.installMachine(machine, best[machine]);
    }
    if (!fallback.evaluate(false)) {
      return 1;
    }
  }
  Graph working(in);
  for (int machine = 0; machine < in.machines; ++machine) {
    working.installMachine(machine, best[machine]);
  }
  if (!working.evaluate(true)) {
    return 1;
  }
  Graph overallBest = working;
  auto updateOverallBest = [&]() {
    if (working.makespan < overallBest.makespan) {
      overallBest = working;
    }
  };
  vector<Elite> pool;
  rememberElite(pool, in, working);
  for (const SeedCandidate& seed : seedCandidates) {
    Graph candidateGraph(in);
    for (int machine = 0; machine < in.machines; ++machine) {
      candidateGraph.installMachine(machine, seed.order[machine]);
    }
    if (candidateGraph.evaluate(true)) {
      rememberElite(pool, in, candidateGraph);
    }
  }
  for (int pass = 0; pass < 2 && !deadline.expired(); ++pass) {
    if (pass == 0) {
      for (int machine = 0; machine < in.machines; ++machine) {
        if (deadline.expired()) {
          break;
        }
        reoptimizeMachine(
          working, in, machine, working.makespan, deadline
        );
        updateOverallBest();
      }
    } else {
      for (int machine = in.machines - 1; machine >= 0; --machine) {
        if (deadline.expired()) {
          break;
        }
        reoptimizeMachine(
          working, in, machine, working.makespan, deadline
        );
        updateOverallBest();
      }
    }
  }
  uint64_t searchRng = 0xd1b54a32d192ed03ULL;
  for (int round = 0; round < 24 && !deadline.expired(); ++round) {
    topologicalRepair(working, in, pool, searchRng, deadline);
    updateOverallBest();
  }
  polish(working, in, deadline);
  updateOverallBest();
  tabuPolish(working, in, deadline);
  updateOverallBest();
  rememberElite(pool, in, overallBest);
  auto pathRelink = [&]() {
    if(deadline.expired()||pool.size()<2)return;
    vector<vector<int>> cur(in.machines);
    for(int m=0;m<in.machines;++m)cur[m]=orderFromGraph(working,in,m);
    size_t target = 0;
    int far = -1;
    for(size_t index=0;index<pool.size();++index){
      int d = 0;
      for(int m=0;m<in.machines;++m)for(int i=0;i<in.jobs;++i){
        auto it=find(cur[m].begin(),cur[m].end(),pool[index].schedule.order[m][i]);
        d+=abs(int(it-cur[m].begin())-i);
      }
      if(d>far){
        far=d;
        target=index;
      }
    }
    vector<vector<int>> targetOrder=pool[target].schedule.order;
    Graph path=working;
    for(int step=0;step<min(in.operations,64)&&!deadline.expired();++step){
      Graph next=path;
      i64 nm=INF;
      int mm=-1,mp=-1;
      for(int m=0;m<in.machines;++m){
        for(int p=0;p+1<in.jobs;++p){
          auto a=find(targetOrder[m].begin(),targetOrder[m].end(),cur[m][p]);
          auto b=find(targetOrder[m].begin(),targetOrder[m].end(),cur[m][p+1]);
          if(a<=b)continue;
          vector<int> co=cur[m];
          swap(co[p],co[p+1]);
          Graph candidate=path;
          candidate.installMachine(m,co);
          if(candidate.evaluate(true)&&candidate.makespan<nm){
            next=std::move(candidate);
            nm=next.makespan;
            mm=m;
            mp=p;
          }
          break;
        }
      }
      if(mm<0)break;
      swap(cur[mm][mp],cur[mm][mp+1]);
      path=std::move(next);
      working=path;
      updateOverallBest();
    }
    working = path;
    if(!deadline.expired()){
      tabuPolish(working, in, deadline);
      updateOverallBest();
    }
  };
  pathRelink();
  uint64_t lateRng = 0x94d049bb133111ebULL;
  for (int round = 0; round < 8 && !deadline.expired(); ++round) {
    if (round > 0 && !pool.empty() && (round & 1) == 0) {
      size_t index = splitmix64(lateRng) % pool.size();
      working = pool[index].schedule.graph;
      working.evaluate(true);
    }
    bool changed = topologicalRepair(working, in, pool, lateRng, deadline);
    if (!changed && !deadline.expired()) {
      int machine = (round * 7 + 3) % in.machines;
      reoptimizeMachine(working, in, machine, working.makespan, deadline);
    }
    if (!deadline.expired()) {
      tabuPolish(working, in, deadline);
      updateOverallBest();
    }
  }
  working = overallBest;
  if (!working.evaluate(true)) {
    return 1;
  }
  best.resize(in.machines);
  for (int machine = 0; machine < in.machines; ++machine) {
    best[machine] = orderFromGraph(working, in, machine);
  }
  emit(in, best);
  return 0;
}
