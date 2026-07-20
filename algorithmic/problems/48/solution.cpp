#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

namespace {

using Clock = chrono::steady_clock;
using i64 = long long;

#ifndef SOLVER_TIME_MS
#define SOLVER_TIME_MS 700
#endif

struct RNG {
    uint64_t state = 0x46c0ffee12345678ULL;

    uint64_t next() {
        state += 0x9e3779b97f4a7c15ULL;
        uint64_t value = state;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }

    int integer(int bound) { return static_cast<int>(next() % static_cast<uint64_t>(bound)); }
};

struct Problem {
    int jobs = 0;
    int machines = 0;
    int operations = 0;
    vector<int> machine;
    vector<i64> duration;
    vector<vector<int>> node_on_machine;
    vector<i64> remaining_work;
    vector<i64> job_load;
    vector<i64> machine_load;
    i64 lower_bound = 0;

    int node(int job, int operation) const { return job * machines + operation; }
    int job(int node_id) const { return node_id / machines; }
    int operation(int node_id) const { return node_id % machines; }
};

struct Schedule {
    vector<vector<int>> order;
    vector<int> position;
    vector<int> critical_path;
    i64 cost = numeric_limits<i64>::max();
};

uint64_t mix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

void rebuild_positions(const Problem& problem, Schedule& schedule) {
    schedule.position.assign(problem.operations, -1);
    for (int machine = 0; machine < problem.machines; ++machine) {
        for (int index = 0; index < problem.jobs; ++index) {
            schedule.position[schedule.order[machine][index]] = index;
        }
    }
}

void swap_adjacent(Schedule& schedule, int machine, int index) {
    int first = schedule.order[machine][index];
    int second = schedule.order[machine][index + 1];
    swap(schedule.order[machine][index], schedule.order[machine][index + 1]);
    schedule.position[first] = index + 1;
    schedule.position[second] = index;
}

void relocate(Schedule& schedule, int machine, int from, int to) {
    if (from == to) return;
    vector<int>& sequence = schedule.order[machine];
    int node_id = sequence[from];
    sequence.erase(sequence.begin() + from);
    sequence.insert(sequence.begin() + to, node_id);
    int first = min(from, to);
    int last = max(from, to);
    for (int index = first; index <= last; ++index) {
        schedule.position[sequence[index]] = index;
    }
}

class Evaluator {
public:
    explicit Evaluator(const Problem& problem)
        : problem_(problem),
          indegree_(problem.operations),
          head_(problem.operations),
          parent_(problem.operations),
          queue_(problem.operations) {}

    bool evaluate(Schedule& schedule, bool trace_path, uint64_t path_salt = 0) {
        const int n = problem_.operations;
        fill(indegree_.begin(), indegree_.end(), 0);
        fill(head_.begin(), head_.end(), 0);
        if (trace_path) fill(parent_.begin(), parent_.end(), -1);

        for (int node_id = 0; node_id < n; ++node_id) {
            int operation = problem_.operation(node_id);
            if (operation > 0) ++indegree_[node_id];
            if (schedule.position[node_id] > 0) ++indegree_[node_id];
        }

        int front = 0;
        int back = 0;
        for (int node_id = 0; node_id < n; ++node_id) {
            if (indegree_[node_id] == 0) {
                head_[node_id] = problem_.duration[node_id];
                queue_[back++] = node_id;
            }
        }

        int visited = 0;
        while (front < back) {
            int current = queue_[front++];
            ++visited;
            int operation = problem_.operation(current);
            if (operation + 1 < problem_.machines) {
                relax(current, current + 1, trace_path, path_salt, back);
            }
            int machine = problem_.machine[current];
            int position = schedule.position[current];
            if (position + 1 < problem_.jobs) {
                int successor = schedule.order[machine][position + 1];
                relax(current, successor, trace_path, path_salt, back);
            }
        }
        if (visited != n) {
            schedule.cost = numeric_limits<i64>::max();
            if (trace_path) schedule.critical_path.clear();
            return false;
        }

        int finish = 0;
        for (int node_id = 1; node_id < n; ++node_id) {
            if (head_[node_id] > head_[finish] ||
                (trace_path && head_[node_id] == head_[finish] &&
                 tie_rank(node_id, path_salt) < tie_rank(finish, path_salt))) {
                finish = node_id;
            }
        }
        schedule.cost = head_[finish];

        if (trace_path) {
            schedule.critical_path.clear();
            for (int node_id = finish; node_id != -1; node_id = parent_[node_id]) {
                schedule.critical_path.push_back(node_id);
            }
            reverse(schedule.critical_path.begin(), schedule.critical_path.end());
        }
        return true;
    }

private:
    const Problem& problem_;
    vector<int> indegree_;
    vector<i64> head_;
    vector<int> parent_;
    vector<int> queue_;

    uint64_t tie_rank(int node_id, uint64_t salt) const {
        return mix64(static_cast<uint64_t>(node_id) ^ salt);
    }

    void relax(
        int current,
        int successor,
        bool trace_path,
        uint64_t path_salt,
        int& queue_back
    ) {
        i64 candidate = head_[current] + problem_.duration[successor];
        if (candidate > head_[successor] ||
            (trace_path && candidate == head_[successor] &&
             (parent_[successor] == -1 ||
              tie_rank(current, path_salt) < tie_rank(parent_[successor], path_salt)))) {
            head_[successor] = candidate;
            if (trace_path) parent_[successor] = current;
        }
        if (--indegree_[successor] == 0) queue_[queue_back++] = successor;
    }
};

Schedule serial_dispatch(const Problem& problem, bool longest_processing_time) {
    vector<int> next_operation(problem.jobs, 0);
    vector<i64> job_ready(problem.jobs, 0);
    vector<i64> machine_ready(problem.machines, 0);
    Schedule schedule;
    schedule.order.assign(problem.machines, {});
    for (auto& sequence : schedule.order) sequence.reserve(problem.jobs);

    for (int step = 0; step < problem.operations; ++step) {
        tuple<i64, i64, int, int> best{
            numeric_limits<i64>::max(), numeric_limits<i64>::max(), -1, -1
        };
        for (int job = 0; job < problem.jobs; ++job) {
            int operation = next_operation[job];
            if (operation == problem.machines) continue;
            int node_id = problem.node(job, operation);
            int machine = problem.machine[node_id];
            i64 start = max(job_ready[job], machine_ready[machine]);
            i64 priority = longest_processing_time ? -problem.duration[node_id]
                                                   : problem.duration[node_id];
            best = min(best, make_tuple(start, priority, job, node_id));
        }
        auto [start, ignored, job, node_id] = best;
        (void)ignored;
        int machine = problem.machine[node_id];
        i64 finish = start + problem.duration[node_id];
        schedule.order[machine].push_back(node_id);
        ++next_operation[job];
        job_ready[job] = finish;
        machine_ready[machine] = finish;
    }
    rebuild_positions(problem, schedule);
    return schedule;
}

Schedule giffler_thompson(const Problem& problem, int rule, uint64_t seed) {
    vector<int> next_operation(problem.jobs, 0);
    vector<i64> job_ready(problem.jobs, 0);
    vector<i64> machine_ready(problem.machines, 0);
    Schedule schedule;
    schedule.order.assign(problem.machines, {});
    for (auto& sequence : schedule.order) sequence.reserve(problem.jobs);

    i64 maximum_job_load = *max_element(problem.job_load.begin(), problem.job_load.end());
    i64 maximum_duration = *max_element(problem.duration.begin(), problem.duration.end());
    auto coefficient = [&](int shift, int bias) {
        return bias + static_cast<int>((seed >> shift) & 255ULL) - 128;
    };
    const int remaining_coefficient = coefficient(0, -64);
    const int duration_coefficient = coefficient(8, 0);
    const int ready_coefficient = coefficient(16, 32);
    const int operations_coefficient = coefficient(24, -32);
    const int load_coefficient = coefficient(32, -16);

    for (int step = 0; step < problem.operations; ++step) {
        int pivot = -1;
        i64 pivot_finish = numeric_limits<i64>::max();
        for (int job = 0; job < problem.jobs; ++job) {
            int operation = next_operation[job];
            if (operation == problem.machines) continue;
            int node_id = problem.node(job, operation);
            int machine = problem.machine[node_id];
            i64 start = max(job_ready[job], machine_ready[machine]);
            i64 finish = start + problem.duration[node_id];
            if (finish < pivot_finish ||
                (finish == pivot_finish && node_id < pivot)) {
                pivot_finish = finish;
                pivot = node_id;
            }
        }

        int pivot_machine = problem.machine[pivot];
        int chosen = -1;
        tuple<i64, uint64_t, int> chosen_key{
            numeric_limits<i64>::max(), numeric_limits<uint64_t>::max(), -1
        };
        for (int job = 0; job < problem.jobs; ++job) {
            int operation = next_operation[job];
            if (operation == problem.machines) continue;
            int node_id = problem.node(job, operation);
            if (problem.machine[node_id] != pivot_machine) continue;
            i64 start = max(job_ready[job], machine_ready[pivot_machine]);
            if (start >= pivot_finish) continue;

            i64 primary = 0;
            if (rule == 0) primary = -problem.remaining_work[node_id];
            else if (rule == 1) primary = problem.duration[node_id];
            else if (rule == 2) primary = -problem.duration[node_id];
            else if (rule == 3) {
                primary = -(problem.remaining_work[node_id] - problem.duration[node_id]);
            } else if (rule == 4) {
                primary = problem.remaining_work[node_id];
            } else if (rule == 5) {
                primary = job_ready[job];
            } else if (rule == 6) {
                i64 remaining = problem.remaining_work[node_id];
                i64 noise = static_cast<i64>(
                    mix64(seed ^ static_cast<uint64_t>(node_id)) %
                    static_cast<uint64_t>(1 + max<i64>(1, remaining / 2))
                );
                primary = -remaining + noise;
            } else if (rule == 7) {
                primary = problem.duration[node_id] + static_cast<i64>(
                    mix64(seed ^ static_cast<uint64_t>(node_id)) %
                    static_cast<uint64_t>(1 + problem.duration[node_id])
                );
            } else if (rule == 8) {
                i64 tail = problem.remaining_work[node_id] - problem.duration[node_id];
                i64 noise = static_cast<i64>(
                    mix64(seed ^ static_cast<uint64_t>(node_id)) %
                    static_cast<uint64_t>(1 + max<i64>(1, tail / 2))
                );
                primary = -tail + noise;
            } else {
                i64 remaining_feature =
                    4096 * problem.remaining_work[node_id] / maximum_job_load;
                i64 duration_feature =
                    4096 * problem.duration[node_id] / maximum_duration;
                i64 ready_feature = 4096 * job_ready[job] / maximum_job_load;
                i64 operations_feature =
                    4096 * (problem.machines - operation) / problem.machines;
                i64 load_feature = 4096 * problem.job_load[job] / maximum_job_load;
                primary = remaining_coefficient * remaining_feature
                    + duration_coefficient * duration_feature
                    + ready_coefficient * ready_feature
                    + operations_coefficient * operations_feature
                    + load_coefficient * load_feature;
            }
            uint64_t secondary = mix64(
                seed + static_cast<uint64_t>(node_id) * 0x9e3779b97f4a7c15ULL
            );
            auto key = make_tuple(primary, secondary, node_id);
            if (key < chosen_key) {
                chosen_key = key;
                chosen = node_id;
            }
        }

        int job = problem.job(chosen);
        int machine = problem.machine[chosen];
        i64 start = max(job_ready[job], machine_ready[machine]);
        i64 finish = start + problem.duration[chosen];
        schedule.order[machine].push_back(chosen);
        ++next_operation[job];
        job_ready[job] = finish;
        machine_ready[machine] = finish;
    }
    rebuild_positions(problem, schedule);
    return schedule;
}

Problem reversed_problem(const Problem& problem) {
    Problem reversed;
    reversed.jobs = problem.jobs;
    reversed.machines = problem.machines;
    reversed.operations = problem.operations;
    reversed.machine.resize(problem.operations);
    reversed.duration.resize(problem.operations);
    reversed.node_on_machine.assign(
        problem.jobs, vector<int>(problem.machines, -1)
    );
    reversed.job_load = problem.job_load;
    reversed.machine_load = problem.machine_load;
    reversed.lower_bound = problem.lower_bound;
    for (int job = 0; job < problem.jobs; ++job) {
        for (int operation = 0; operation < problem.machines; ++operation) {
            int reversed_node = reversed.node(job, operation);
            int original_node = problem.node(
                job, problem.machines - 1 - operation
            );
            reversed.machine[reversed_node] = problem.machine[original_node];
            reversed.duration[reversed_node] = problem.duration[original_node];
            reversed.node_on_machine[job][reversed.machine[reversed_node]] =
                reversed_node;
        }
    }
    reversed.remaining_work.assign(problem.operations, 0);
    for (int job = 0; job < problem.jobs; ++job) {
        i64 remaining = 0;
        for (int operation = problem.machines - 1; operation >= 0; --operation) {
            int node_id = reversed.node(job, operation);
            remaining += reversed.duration[node_id];
            reversed.remaining_work[node_id] = remaining;
        }
    }
    return reversed;
}

Schedule mirror_reversed_schedule(
    const Problem& problem,
    const Schedule& reversed_schedule
) {
    Schedule schedule;
    schedule.order.assign(problem.machines, {});
    for (int machine = 0; machine < problem.machines; ++machine) {
        schedule.order[machine].reserve(problem.jobs);
        for (auto iterator = reversed_schedule.order[machine].rbegin();
             iterator != reversed_schedule.order[machine].rend(); ++iterator) {
            int reversed_node = *iterator;
            int job = reversed_node / problem.machines;
            int reversed_operation = reversed_node % problem.machines;
            schedule.order[machine].push_back(problem.node(
                job, problem.machines - 1 - reversed_operation
            ));
        }
    }
    rebuild_positions(problem, schedule);
    return schedule;
}

vector<int> schrage_order(const Problem& problem, int fixed_machine) {
    struct Task {
        int node = -1;
        i64 release = 0;
        i64 tail = 0;
        bool done = false;
    };
    vector<Task> tasks(problem.jobs);
    for (int job = 0; job < problem.jobs; ++job) {
        int node_id = problem.node_on_machine[job][fixed_machine];
        int operation = problem.operation(node_id);
        i64 release = 0;
        for (int k = 0; k < operation; ++k) {
            release += problem.duration[problem.node(job, k)];
        }
        i64 tail = 0;
        for (int k = operation + 1; k < problem.machines; ++k) {
            tail += problem.duration[problem.node(job, k)];
        }
        tasks[job] = {node_id, release, tail, false};
    }

    vector<int> result;
    result.reserve(problem.jobs);
    i64 now = numeric_limits<i64>::max();
    for (const Task& task : tasks) now = min(now, task.release);
    while (static_cast<int>(result.size()) < problem.jobs) {
        int chosen = -1;
        for (int job = 0; job < problem.jobs; ++job) {
            const Task& task = tasks[job];
            if (task.done || task.release > now) continue;
            if (chosen == -1 || task.tail > tasks[chosen].tail ||
                (task.tail == tasks[chosen].tail &&
                 problem.duration[task.node] > problem.duration[tasks[chosen].node])) {
                chosen = job;
            }
        }
        if (chosen == -1) {
            now = numeric_limits<i64>::max();
            for (const Task& task : tasks) {
                if (!task.done) now = min(now, task.release);
            }
            continue;
        }
        tasks[chosen].done = true;
        result.push_back(tasks[chosen].node);
        now = max(now, tasks[chosen].release) + problem.duration[tasks[chosen].node];
    }
    return result;
}

Schedule complete_fixed_machine(
    const Problem& problem,
    int,
    const vector<int>& fixed_order
) {
    vector<int> fixed_predecessor(problem.operations, -1);
    vector<int> fixed_successor(problem.operations, -1);
    for (int index = 0; index + 1 < problem.jobs; ++index) {
        fixed_successor[fixed_order[index]] = fixed_order[index + 1];
        fixed_predecessor[fixed_order[index + 1]] = fixed_order[index];
    }

    vector<int> indegree(problem.operations, 0);
    vector<char> scheduled(problem.operations, false);
    vector<i64> job_ready(problem.jobs, 0);
    vector<i64> machine_ready(problem.machines, 0);
    for (int node_id = 0; node_id < problem.operations; ++node_id) {
        if (problem.operation(node_id) > 0) ++indegree[node_id];
        if (fixed_predecessor[node_id] != -1) ++indegree[node_id];
    }

    Schedule schedule;
    schedule.order.assign(problem.machines, {});
    for (auto& sequence : schedule.order) sequence.reserve(problem.jobs);
    for (int step = 0; step < problem.operations; ++step) {
        int chosen = -1;
        tuple<i64, i64, int> best{
            numeric_limits<i64>::max(), numeric_limits<i64>::max(), -1
        };
        for (int node_id = 0; node_id < problem.operations; ++node_id) {
            if (scheduled[node_id] || indegree[node_id] != 0) continue;
            int job = problem.job(node_id);
            int machine = problem.machine[node_id];
            i64 start = max(job_ready[job], machine_ready[machine]);
            auto key = make_tuple(start, -problem.remaining_work[node_id], node_id);
            if (key < best) {
                best = key;
                chosen = node_id;
            }
        }

        int job = problem.job(chosen);
        int machine = problem.machine[chosen];
        i64 start = max(job_ready[job], machine_ready[machine]);
        i64 finish = start + problem.duration[chosen];
        scheduled[chosen] = true;
        schedule.order[machine].push_back(chosen);
        job_ready[job] = finish;
        machine_ready[machine] = finish;
        if (problem.operation(chosen) + 1 < problem.machines) --indegree[chosen + 1];
        if (fixed_successor[chosen] != -1) --indegree[fixed_successor[chosen]];
    }
    rebuild_positions(problem, schedule);
    return schedule;
}

class InsertionEvaluator {
public:
    explicit InsertionEvaluator(const Problem& problem)
        : problem_(problem),
          machine_successor_(problem.operations),
          indegree_(problem.operations),
          head_(problem.operations),
          tail_(problem.operations),
          queue_(problem.operations),
          topological_(problem.operations) {}

    bool evaluate(const vector<vector<int>>& partial_order, int focus) {
        fill(machine_successor_.begin(), machine_successor_.end(), -1);
        fill(indegree_.begin(), indegree_.end(), 0);
        fill(head_.begin(), head_.end(), 0);
        for (int node_id = 0; node_id < problem_.operations; ++node_id) {
            if (problem_.operation(node_id) > 0) ++indegree_[node_id];
        }
        for (const vector<int>& sequence : partial_order) {
            for (int index = 0; index + 1 < static_cast<int>(sequence.size()); ++index) {
                machine_successor_[sequence[index]] = sequence[index + 1];
                ++indegree_[sequence[index + 1]];
            }
        }

        int front = 0;
        int back = 0;
        for (int node_id = 0; node_id < problem_.operations; ++node_id) {
            if (indegree_[node_id] == 0) {
                head_[node_id] = problem_.duration[node_id];
                queue_[back++] = node_id;
            }
        }
        int visited = 0;
        while (front < back) {
            int current = queue_[front++];
            topological_[visited++] = current;
            if (problem_.operation(current) + 1 < problem_.machines) {
                relax_head(current, current + 1, back);
            }
            if (machine_successor_[current] != -1) {
                relax_head(current, machine_successor_[current], back);
            }
        }
        if (visited != problem_.operations) return false;

        fill(tail_.begin(), tail_.end(), 0);
        for (int index = problem_.operations - 1; index >= 0; --index) {
            int current = topological_[index];
            i64 successor_tail = 0;
            if (problem_.operation(current) + 1 < problem_.machines) {
                successor_tail = max(successor_tail, tail_[current + 1]);
            }
            if (machine_successor_[current] != -1) {
                successor_tail = max(successor_tail, tail_[machine_successor_[current]]);
            }
            tail_[current] = problem_.duration[current] + successor_tail;
        }
        through_focus_ = head_[focus] + tail_[focus] - problem_.duration[focus];
        makespan_ = *max_element(head_.begin(), head_.end());
        return true;
    }

    i64 through_focus() const { return through_focus_; }
    i64 makespan() const { return makespan_; }

private:
    const Problem& problem_;
    vector<int> machine_successor_;
    vector<int> indegree_;
    vector<i64> head_;
    vector<i64> tail_;
    vector<int> queue_;
    vector<int> topological_;
    i64 through_focus_ = 0;
    i64 makespan_ = 0;

    void relax_head(int current, int successor, int& queue_back) {
        head_[successor] = max(
            head_[successor], head_[current] + problem_.duration[successor]
        );
        if (--indegree_[successor] == 0) queue_[queue_back++] = successor;
    }
};

Schedule greedy_insertion_seed(
    const Problem& problem,
    Evaluator& full_evaluator,
    int first_job,
    int beam_width,
    const Clock::time_point& deadline
) {
    vector<vector<int>> partial_order(problem.machines);
    vector<char> inserted(problem.operations, false);
    for (int operation = 0; operation < problem.machines; ++operation) {
        int node_id = problem.node(first_job, operation);
        partial_order[problem.machine[node_id]].push_back(node_id);
        inserted[node_id] = true;
    }

    vector<int> remaining;
    remaining.reserve(problem.operations - problem.machines);
    for (int node_id = 0; node_id < problem.operations; ++node_id) {
        if (!inserted[node_id]) remaining.push_back(node_id);
    }
    stable_sort(remaining.begin(), remaining.end(), [&](int left, int right) {
        if (problem.duration[left] != problem.duration[right]) {
            return problem.duration[left] > problem.duration[right];
        }
        return problem.remaining_work[left] > problem.remaining_work[right];
    });

    struct BeamState {
        vector<vector<int>> order;
        tuple<i64, i64, int> key;
    };
    vector<BeamState> beam;
    beam.push_back({std::move(partial_order), {0, 0, 0}});
    InsertionEvaluator evaluator(problem);
    for (int node_id : remaining) {
        if (Clock::now() >= deadline) return serial_dispatch(problem, false);
        int machine = problem.machine[node_id];
        vector<BeamState> next_beam;
        for (BeamState& state : beam) {
            vector<int>& sequence = state.order[machine];
            for (int position = 0; position <= static_cast<int>(sequence.size()); ++position) {
                sequence.insert(sequence.begin() + position, node_id);
                if (evaluator.evaluate(state.order, node_id)) {
                    auto key = make_tuple(
                        evaluator.through_focus(), evaluator.makespan(), position
                    );
                    if (static_cast<int>(next_beam.size()) < beam_width ||
                        key < next_beam.back().key) {
                        BeamState candidate{state.order, key};
                        auto insertion = lower_bound(
                            next_beam.begin(),
                            next_beam.end(),
                            key,
                            [](const BeamState& existing, const tuple<i64, i64, int>& value) {
                                return existing.key < value;
                            }
                        );
                        next_beam.insert(insertion, std::move(candidate));
                        if (static_cast<int>(next_beam.size()) > beam_width) {
                            next_beam.pop_back();
                        }
                    }
                }
                sequence.erase(sequence.begin() + position);
            }
        }
        if (next_beam.empty()) return serial_dispatch(problem, false);
        beam = std::move(next_beam);
    }

    Schedule best_schedule;
    for (BeamState& state : beam) {
        Schedule schedule;
        schedule.order = std::move(state.order);
        rebuild_positions(problem, schedule);
        full_evaluator.evaluate(schedule, false);
        if (schedule.cost < best_schedule.cost) best_schedule = std::move(schedule);
    }
    return best_schedule;
}

Schedule common_permutation(const Problem& problem, const vector<int>& permutation) {
    Schedule schedule;
    schedule.order.assign(problem.machines, vector<int>(problem.jobs));
    for (int machine = 0; machine < problem.machines; ++machine) {
        for (int index = 0; index < problem.jobs; ++index) {
            int job = permutation[index];
            schedule.order[machine][index] = problem.node_on_machine[job][machine];
        }
    }
    rebuild_positions(problem, schedule);
    return schedule;
}

i64 common_subset_makespan(
    const Problem& problem,
    const vector<int>& permutation
) {
    vector<int> indegree(problem.operations, 0);
    vector<int> machine_successor(problem.operations, -1);
    vector<char> active(problem.operations, false);
    for (int job : permutation) {
        for (int operation = 0; operation < problem.machines; ++operation) {
            int node_id = problem.node(job, operation);
            active[node_id] = true;
            if (operation > 0) ++indegree[node_id];
        }
    }
    for (int machine = 0; machine < problem.machines; ++machine) {
        for (int index = 0; index + 1 < static_cast<int>(permutation.size()); ++index) {
            int first = problem.node_on_machine[permutation[index]][machine];
            int second = problem.node_on_machine[permutation[index + 1]][machine];
            machine_successor[first] = second;
            ++indegree[second];
        }
    }

    vector<int> queue;
    queue.reserve(permutation.size() * problem.machines);
    vector<i64> completion(problem.operations, 0);
    for (int node_id = 0; node_id < problem.operations; ++node_id) {
        if (active[node_id] && indegree[node_id] == 0) {
            completion[node_id] = problem.duration[node_id];
            queue.push_back(node_id);
        }
    }
    for (int front = 0; front < static_cast<int>(queue.size()); ++front) {
        int current = queue[front];
        auto relax = [&](int successor) {
            completion[successor] = max(
                completion[successor],
                completion[current] + problem.duration[successor]
            );
            if (--indegree[successor] == 0) queue.push_back(successor);
        };
        if (problem.operation(current) + 1 < problem.machines) {
            relax(current + 1);
        }
        if (machine_successor[current] != -1) {
            relax(machine_successor[current]);
        }
    }
    return *max_element(completion.begin(), completion.end());
}

vector<int> neh_permutation(
    const Problem& problem,
    const vector<int>& job_order,
    const Clock::time_point& deadline
) {
    vector<int> permutation;
    permutation.reserve(problem.jobs);
    for (int job : job_order) {
        if (Clock::now() >= deadline) {
            permutation.push_back(job);
            continue;
        }
        i64 best_cost = numeric_limits<i64>::max();
        int best_position = 0;
        for (int position = 0; position <= static_cast<int>(permutation.size()); ++position) {
            vector<int> candidate = permutation;
            candidate.insert(candidate.begin() + position, job);
            i64 cost = common_subset_makespan(problem, candidate);
            if (cost < best_cost) {
                best_cost = cost;
                best_position = position;
            }
        }
        permutation.insert(permutation.begin() + best_position, job);
    }
    return permutation;
}

Schedule improve_common_permutation(
    const Problem& problem,
    Evaluator& evaluator,
    vector<int> permutation,
    const Clock::time_point& construction_deadline
) {
    Schedule best = common_permutation(problem, permutation);
    evaluator.evaluate(best, false);
    bool improved = true;
    for (int pass = 0; pass < 3 && improved; ++pass) {
        improved = false;
        for (int old_position = 0; old_position < problem.jobs; ++old_position) {
            if (Clock::now() >= construction_deadline) return best;
            vector<int> base = permutation;
            int job = base[old_position];
            base.erase(base.begin() + old_position);
            int best_position = old_position;
            i64 best_cost = best.cost;
            for (int new_position = 0; new_position < problem.jobs; ++new_position) {
                vector<int> candidate_permutation = base;
                candidate_permutation.insert(candidate_permutation.begin() + new_position, job);
                Schedule candidate = common_permutation(problem, candidate_permutation);
                evaluator.evaluate(candidate, false);
                if (candidate.cost < best_cost) {
                    best_cost = candidate.cost;
                    best_position = new_position;
                }
            }
            if (best_position != old_position) {
                permutation = base;
                permutation.insert(permutation.begin() + best_position, job);
                best = common_permutation(problem, permutation);
                evaluator.evaluate(best, false);
                improved = true;
            }
        }
    }
    return best;
}

uint64_t schedule_hash(const Schedule& schedule) {
    uint64_t hash = 1469598103934665603ULL;
    for (const auto& sequence : schedule.order) {
        for (int node_id : sequence) {
            hash = (hash ^ static_cast<uint64_t>(node_id + 1)) * 1099511628211ULL;
        }
    }
    return hash;
}

void add_seed(
    const Problem&,
    Evaluator& evaluator,
    vector<Schedule>& seeds,
    vector<uint64_t>& hashes,
    Schedule schedule
) {
    if (!evaluator.evaluate(schedule, false)) return;
    uint64_t hash = schedule_hash(schedule);
    if (find(hashes.begin(), hashes.end(), hash) != hashes.end()) return;
    hashes.push_back(hash);
    seeds.push_back(std::move(schedule));
}

struct Move {
    int machine = -1;
    int index = -1;
};

struct RelocateMove {
    int machine = -1;
    int from = -1;
    int to = -1;
};

struct Block {
    int first = 0;
    int last = 0;
};

vector<Block> critical_blocks(const Problem& problem, const Schedule& schedule) {
    vector<Block> blocks;
    const vector<int>& path = schedule.critical_path;
    for (int first = 0; first < static_cast<int>(path.size());) {
        int last = first;
        int machine = problem.machine[path[first]];
        while (last + 1 < static_cast<int>(path.size()) &&
               problem.machine[path[last + 1]] == machine &&
               schedule.position[path[last + 1]] == schedule.position[path[last]] + 1) {
            ++last;
        }
        blocks.push_back({first, last});
        first = last + 1;
    }
    return blocks;
}

vector<Move> n5_moves(const Problem& problem, const Schedule& schedule) {
    vector<Block> blocks = critical_blocks(problem, schedule);
    vector<Move> moves;
    for (int block_index = 0; block_index < static_cast<int>(blocks.size()); ++block_index) {
        const Block& block = blocks[block_index];
        if (block.first == block.last) continue;
        int machine = problem.machine[schedule.critical_path[block.first]];
        int first_index = schedule.position[schedule.critical_path[block.first]];
        int last_index = schedule.position[schedule.critical_path[block.last - 1]];
        if (blocks.size() == 1 || block_index > 0) {
            moves.push_back({machine, first_index});
        }
        if ((blocks.size() == 1 || block_index + 1 < static_cast<int>(blocks.size())) &&
            last_index != first_index) {
            moves.push_back({machine, last_index});
        }
    }
    return moves;
}

vector<Move> n1_moves(const Problem& problem, const Schedule& schedule) {
    vector<Move> moves;
    for (const Block& block : critical_blocks(problem, schedule)) {
        int machine = problem.machine[schedule.critical_path[block.first]];
        for (int path_index = block.first; path_index < block.last; ++path_index) {
            moves.push_back({machine, schedule.position[schedule.critical_path[path_index]]});
        }
    }
    return moves;
}

vector<RelocateMove> n7_moves(const Problem& problem, const Schedule& schedule) {
    vector<RelocateMove> moves;
    for (const Block& block : critical_blocks(problem, schedule)) {
        int machine = problem.machine[schedule.critical_path[block.first]];
        int first = schedule.position[schedule.critical_path[block.first]];
        int last = schedule.position[schedule.critical_path[block.last]];
        for (int index = first + 1; index < last; ++index) {
            moves.push_back({machine, index, first});
            moves.push_back({machine, index, last});
            moves.push_back({machine, first, index});
            moves.push_back({machine, last, index});
        }
    }
    sort(moves.begin(), moves.end(), [](const RelocateMove& left, const RelocateMove& right) {
        return tie(left.machine, left.from, left.to) < tie(right.machine, right.from, right.to);
    });
    moves.erase(
        unique(moves.begin(), moves.end(), [](const RelocateMove& left, const RelocateMove& right) {
            return left.machine == right.machine && left.from == right.from && left.to == right.to;
        }),
        moves.end()
    );
    return moves;
}

bool perturb(
    const Problem& problem,
    Evaluator& evaluator,
    Schedule& schedule,
    RNG& rng,
    int steps,
    const Clock::time_point& deadline
) {
    for (int step = 0; step < steps; ++step) {
        if (Clock::now() >= deadline) return false;
        evaluator.evaluate(schedule, true, rng.next());
        vector<Move> moves = n1_moves(problem, schedule);
        if (moves.empty()) return false;
        int offset = rng.integer(static_cast<int>(moves.size()));
        bool applied = false;
        for (int attempt = 0; attempt < static_cast<int>(moves.size()); ++attempt) {
            Move move = moves[(offset + attempt) % moves.size()];
            swap_adjacent(schedule, move.machine, move.index);
            if (evaluator.evaluate(schedule, true, rng.next())) {
                applied = true;
                break;
            }
            swap_adjacent(schedule, move.machine, move.index);
        }
        if (!applied) {
            evaluator.evaluate(schedule, true, rng.next());
            return false;
        }
    }
    return true;
}

Schedule annealing_search(
    const Problem& problem,
    Evaluator& evaluator,
    vector<Schedule> seeds,
    RNG& rng,
    const Clock::time_point& deadline
) {
    sort(seeds.begin(), seeds.end(), [](const Schedule& left, const Schedule& right) {
        return left.cost < right.cost;
    });
    Schedule best = seeds.front();
    if (best.cost == problem.lower_bound) return best;
    Schedule current = best;
    double average_duration = accumulate(
        problem.duration.begin(), problem.duration.end(), 0.0
    ) / max(1, problem.operations);
    double temperature = average_duration;
    int stage_iteration = 0;

    while (Clock::now() < deadline) {
        evaluator.evaluate(current, true, rng.next());
        vector<Move> moves = n5_moves(problem, current);
        if (moves.empty()) {
            current = best;
            perturb(problem, evaluator, current, rng, 5, deadline);
            continue;
        }

        const Move& selected = moves[rng.integer(static_cast<int>(moves.size()))];
        i64 original_cost = current.cost;
        swap_adjacent(current, selected.machine, selected.index);
        bool valid = evaluator.evaluate(current, false);
        i64 candidate_cost = current.cost;
        bool accept = valid && candidate_cost <= original_cost;
        if (valid && !accept) {
            double probability = exp(
                -static_cast<double>(candidate_cost - original_cost) /
                max(1.0, temperature)
            );
            double draw = static_cast<double>(rng.next() >> 11) /
                static_cast<double>(uint64_t{1} << 53);
            accept = draw < probability;
        }
        if (!accept) {
            swap_adjacent(current, selected.machine, selected.index);
            current.cost = original_cost;
        } else if (current.cost < best.cost) {
            best = current;
            if (best.cost == problem.lower_bound) return best;
        }

        temperature *= 0.9999;
        if (++stage_iteration >= 16000 ||
            temperature < average_duration * 0.05) {
            stage_iteration = 0;
            current = best;
            perturb(
                problem,
                evaluator,
                current,
                rng,
                2 + rng.integer(6),
                deadline
            );
            temperature = average_duration;
        }
    }
    return best;
}

Schedule iterated_jump_search(
    const Problem& problem,
    Evaluator& evaluator,
    vector<Schedule> seeds,
    RNG& rng,
    const Clock::time_point& deadline
) {
    sort(seeds.begin(), seeds.end(), [](const Schedule& left, const Schedule& right) {
        return left.cost < right.cost;
    });
    Schedule best = seeds.front();
    if (best.cost == problem.lower_bound) return best;
    Schedule current = best;
    int local_optima = 0;
    int unchanged_optima = 0;
    i64 last_optimum = numeric_limits<i64>::max();

    while (Clock::now() < deadline) {
        bool descended = false;
        while (Clock::now() < deadline) {
            evaluator.evaluate(current, true, rng.next());
            vector<Move> moves = n5_moves(problem, current);
            if (moves.empty()) break;
            i64 original_cost = current.cost;
            int best_move = -1;
            i64 best_cost = original_cost;
            int offset = rng.integer(static_cast<int>(moves.size()));
            for (int attempt = 0; attempt < static_cast<int>(moves.size()); ++attempt) {
                int index = (offset + attempt) % static_cast<int>(moves.size());
                const Move& move = moves[index];
                swap_adjacent(current, move.machine, move.index);
                bool valid = evaluator.evaluate(current, false);
                i64 candidate_cost = current.cost;
                swap_adjacent(current, move.machine, move.index);
                if (valid && candidate_cost < best_cost) {
                    best_cost = candidate_cost;
                    best_move = index;
                }
            }
            current.cost = original_cost;
            if (best_move == -1) break;
            Move move = moves[best_move];
            swap_adjacent(current, move.machine, move.index);
            evaluator.evaluate(current, true, rng.next());
            descended = true;
            if (current.cost < best.cost) {
                best = current;
                if (best.cost == problem.lower_bound) return best;
            }
        }
        if (Clock::now() >= deadline) break;

        ++local_optima;
        if (current.cost == last_optimum) ++unchanged_optima;
        else unchanged_optima = 0;
        last_optimum = current.cost;

        if (unchanged_optima >= 4 || local_optima % 30 == 0) {
            perturb(problem, evaluator, current, rng, 5 + min(15, unchanged_optima), deadline);
            unchanged_optima = 0;
            continue;
        }

        int jumps = 1 + rng.integer(2);
        for (int jump = 0; jump < jumps && Clock::now() < deadline; ++jump) {
            evaluator.evaluate(current, true, rng.next());
            vector<Move> moves = n5_moves(problem, current);
            if (moves.empty()) moves = n1_moves(problem, current);
            if (moves.empty()) break;
            i64 original_cost = current.cost;
            int offset = rng.integer(static_cast<int>(moves.size()));
            int selected = -1;
            i64 mildest_worsening = numeric_limits<i64>::max();
            for (int attempt = 0; attempt < static_cast<int>(moves.size()); ++attempt) {
                int index = (offset + attempt) % static_cast<int>(moves.size());
                const Move& move = moves[index];
                swap_adjacent(current, move.machine, move.index);
                bool valid = evaluator.evaluate(current, false);
                i64 candidate_cost = current.cost;
                swap_adjacent(current, move.machine, move.index);
                if (valid && candidate_cost >= original_cost &&
                    candidate_cost < mildest_worsening) {
                    mildest_worsening = candidate_cost;
                    selected = index;
                }
            }
            current.cost = original_cost;
            if (selected == -1) selected = offset;
            Move move = moves[selected];
            swap_adjacent(current, move.machine, move.index);
            if (!evaluator.evaluate(current, true, rng.next())) {
                swap_adjacent(current, move.machine, move.index);
                evaluator.evaluate(current, true, rng.next());
                break;
            }
        }

        if (!descended && local_optima % 12 == 0) {
            current = best;
            perturb(problem, evaluator, current, rng, 3 + rng.integer(5), deadline);
        }
    }
    return best;
}

Schedule n7_tabu_search(
    const Problem& problem,
    Evaluator& evaluator,
    Schedule current,
    RNG& rng,
    const Clock::time_point& deadline
) {
    Schedule best = current;
    evaluator.evaluate(current, true, rng.next());
    vector<int> tabu_until(problem.operations, 0);
    int iteration = 0;
    int since_best = 0;
    while (Clock::now() < deadline) {
        ++iteration;
        vector<RelocateMove> moves = n7_moves(problem, current);
        for (const Move& move : n5_moves(problem, current)) {
            moves.push_back({move.machine, move.index, move.index + 1});
        }
        for (const Block& block : critical_blocks(problem, current)) {
            int first_node = current.critical_path[block.first];
            int last_node = current.critical_path[block.last];
            int machine = problem.machine[first_node];
            int first = current.position[first_node];
            int last = current.position[last_node];
            if (last + 1 < problem.jobs) {
                int to = last + 1 + rng.integer(problem.jobs - last - 1);
                moves.push_back({machine, first, to});
            }
            if (first > 0) {
                int to = rng.integer(first);
                moves.push_back({machine, last, to});
            }
            if (block.last - block.first >= 2) {
                int path_index = block.first + 1 +
                    rng.integer(block.last - block.first - 1);
                int from = current.position[current.critical_path[path_index]];
                int to = rng.integer(problem.jobs);
                if (from != to) moves.push_back({machine, from, to});
            }
        }
        if (moves.empty()) break;

        i64 original_cost = current.cost;
        i64 chosen_cost = numeric_limits<i64>::max();
        int chosen_index = -1;
        int ties = 0;
        for (int index = 0; index < static_cast<int>(moves.size()); ++index) {
            const RelocateMove& move = moves[index];
            int node_id = current.order[move.machine][move.from];
            relocate(current, move.machine, move.from, move.to);
            bool valid = evaluator.evaluate(current, false);
            i64 candidate_cost = current.cost;
            relocate(current, move.machine, move.to, move.from);
            if (!valid) continue;
            bool admissible = tabu_until[node_id] <= iteration || candidate_cost < best.cost;
            if (!admissible) continue;
            if (candidate_cost < chosen_cost) {
                chosen_cost = candidate_cost;
                chosen_index = index;
                ties = 1;
            } else if (candidate_cost == chosen_cost && rng.integer(++ties) == 0) {
                chosen_index = index;
            }
            if ((index & 15) == 15 && Clock::now() >= deadline) break;
        }
        current.cost = original_cost;
        if (chosen_index == -1 || Clock::now() >= deadline) break;

        const RelocateMove& chosen = moves[chosen_index];
        int node_id = current.order[chosen.machine][chosen.from];
        relocate(current, chosen.machine, chosen.from, chosen.to);
        if (!evaluator.evaluate(current, true, rng.next())) break;
        tabu_until[node_id] = iteration + 10 + rng.integer(6);

        if (current.cost < best.cost) {
            best = current;
            since_best = 0;
            if (best.cost == problem.lower_bound) break;
        } else {
            ++since_best;
        }
        if (since_best >= 120) {
            current = best;
            evaluator.evaluate(current, true, rng.next());
            perturb(problem, evaluator, current, rng, 5, deadline);
            fill(tabu_until.begin(), tabu_until.end(), 0);
            since_best = 0;
        }
    }
    return best;
}

Problem read_problem() {
    Problem problem;
    cin >> problem.jobs >> problem.machines;
    problem.operations = problem.jobs * problem.machines;
    problem.machine.resize(problem.operations);
    problem.duration.resize(problem.operations);
    problem.node_on_machine.assign(
        problem.jobs, vector<int>(problem.machines, -1)
    );
    problem.job_load.assign(problem.jobs, 0);
    problem.machine_load.assign(problem.machines, 0);
    for (int job = 0; job < problem.jobs; ++job) {
        for (int operation = 0; operation < problem.machines; ++operation) {
            int machine;
            i64 duration;
            cin >> machine >> duration;
            int node_id = problem.node(job, operation);
            problem.machine[node_id] = machine;
            problem.duration[node_id] = duration;
            problem.node_on_machine[job][machine] = node_id;
            problem.job_load[job] += duration;
            problem.machine_load[machine] += duration;
        }
    }
    problem.remaining_work.assign(problem.operations, 0);
    for (int job = 0; job < problem.jobs; ++job) {
        i64 remaining = 0;
        for (int operation = problem.machines - 1; operation >= 0; --operation) {
            int node_id = problem.node(job, operation);
            remaining += problem.duration[node_id];
            problem.remaining_work[node_id] = remaining;
        }
    }
    for (i64 load : problem.job_load) problem.lower_bound = max(problem.lower_bound, load);
    for (i64 load : problem.machine_load) problem.lower_bound = max(problem.lower_bound, load);
    return problem;
}

void print_schedule(const Problem& problem, const Schedule& schedule) {
    for (int machine = 0; machine < problem.machines; ++machine) {
        for (int index = 0; index < problem.jobs; ++index) {
            if (index) cout << ' ';
            cout << problem.job(schedule.order[machine][index]);
        }
        cout << '\n';
    }
}

}  // namespace

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    const Clock::time_point started = Clock::now();
    Problem problem = read_problem();
    int route_matches = 0;
    for (int job = 1; job < problem.jobs; ++job) {
        for (int operation = 0; operation < problem.machines; ++operation) {
            route_matches +=
                problem.machine[problem.node(job, operation)] ==
                problem.machine[problem.node(0, operation)];
        }
    }
    double flow_similarity = problem.jobs == 1
        ? 1.0
        : static_cast<double>(route_matches) /
          static_cast<double>((problem.jobs - 1) * problem.machines);
    Evaluator evaluator(problem);
    RNG rng;
    vector<Schedule> seeds;
    vector<uint64_t> hashes;

    add_seed(problem, evaluator, seeds, hashes, serial_dispatch(problem, false));
    add_seed(problem, evaluator, seeds, hashes, serial_dispatch(problem, true));
    for (int rule = 0; rule < 6; ++rule) {
        add_seed(
            problem,
            evaluator,
            seeds,
            hashes,
            giffler_thompson(problem, rule, rng.next())
        );
    }
    for (int attempt = 0; attempt < 512; ++attempt) {
        add_seed(
            problem,
            evaluator,
            seeds,
            hashes,
            giffler_thompson(problem, 6 + attempt % 4, rng.next())
        );
    }
    if (flow_similarity >= 0.15) {
        Problem backwards = reversed_problem(problem);
        add_seed(
            problem,
            evaluator,
            seeds,
            hashes,
            mirror_reversed_schedule(
                problem, serial_dispatch(backwards, false)
            )
        );
        add_seed(
            problem,
            evaluator,
            seeds,
            hashes,
            mirror_reversed_schedule(
                problem, serial_dispatch(backwards, true)
            )
        );
        for (int rule = 0; rule < 6; ++rule) {
            add_seed(
                problem,
                evaluator,
                seeds,
                hashes,
                mirror_reversed_schedule(
                    problem, giffler_thompson(backwards, rule, rng.next())
                )
            );
        }
        for (int attempt = 0; attempt < 128; ++attempt) {
            add_seed(
                problem,
                evaluator,
                seeds,
                hashes,
                mirror_reversed_schedule(
                    problem,
                    giffler_thompson(backwards, 6 + attempt % 4, rng.next())
                )
            );
        }
    }
#ifdef LOCAL_DEBUG
    cerr << "after_dispatch_ms="
         << chrono::duration_cast<chrono::milliseconds>(Clock::now() - started).count()
         << '\n';
#endif

    vector<int> machine_rank(problem.machines);
    iota(machine_rank.begin(), machine_rank.end(), 0);
    sort(machine_rank.begin(), machine_rank.end(), [&](int left, int right) {
        return problem.machine_load[left] > problem.machine_load[right];
    });
    for (int index = 0; index < min(2, problem.machines); ++index) {
        int machine = machine_rank[index];
        vector<int> order = schrage_order(problem, machine);
        add_seed(
            problem,
            evaluator,
            seeds,
            hashes,
            complete_fixed_machine(problem, machine, order)
        );
    }

    vector<int> insertion_jobs(problem.jobs);
    iota(insertion_jobs.begin(), insertion_jobs.end(), 0);
    sort(insertion_jobs.begin(), insertion_jobs.end(), [&](int left, int right) {
        return problem.job_load[left] > problem.job_load[right];
    });
    add_seed(
        problem,
        evaluator,
        seeds,
        hashes,
        greedy_insertion_seed(
            problem,
            evaluator,
            insertion_jobs.front(),
            1,
            started + chrono::milliseconds(200)
        )
    );
#ifdef LOCAL_DEBUG
    cerr << "after_insertion_ms="
         << chrono::duration_cast<chrono::milliseconds>(Clock::now() - started).count()
         << '\n';
#endif
    if (flow_similarity >= 0.15) {
        const Clock::time_point construction_deadline =
            started + chrono::milliseconds(280);
        vector<int> descending(problem.jobs);
        iota(descending.begin(), descending.end(), 0);
        sort(descending.begin(), descending.end(), [&](int left, int right) {
            return problem.job_load[left] > problem.job_load[right];
        });
        vector<int> neh = neh_permutation(
            problem,
            descending,
            started + chrono::milliseconds(220)
        );
        add_seed(
            problem,
            evaluator,
            seeds,
            hashes,
            common_permutation(problem, neh)
        );
        add_seed(
            problem,
            evaluator,
            seeds,
            hashes,
            improve_common_permutation(
                problem, evaluator, neh, construction_deadline
            )
        );
        add_seed(
            problem,
            evaluator,
            seeds,
            hashes,
            improve_common_permutation(
                problem, evaluator, descending, construction_deadline
            )
        );
        reverse(descending.begin(), descending.end());
        add_seed(
            problem,
            evaluator,
            seeds,
            hashes,
            improve_common_permutation(
                problem, evaluator, descending, construction_deadline
            )
        );
    }

#ifdef LOCAL_DEBUG
    {
        vector<i64> costs;
        for (const Schedule& seed : seeds) costs.push_back(seed.cost);
        sort(costs.begin(), costs.end());
        cerr << "seeds=" << costs.size() << " costs:";
        for (i64 cost : costs) cerr << ' ' << cost;
        cerr << " lb=" << problem.lower_bound << '\n';
    }
#endif
    const Clock::time_point deadline = started + chrono::milliseconds(SOLVER_TIME_MS);
    RNG search_rng;
    search_rng.state = 0xdeadbeefcafebabeULL;
    int annealing_ms = 500;
    if (flow_similarity < 0.15 && 4 * problem.machines > 3 * problem.jobs) {
        annealing_ms = 600;
    }
    Schedule best = annealing_search(
        problem,
        evaluator,
        std::move(seeds),
        search_rng,
        started + chrono::milliseconds(annealing_ms)
    );
    if (best.cost != problem.lower_bound && Clock::now() < deadline) {
        if (flow_similarity >= 0.15) {
            best = n7_tabu_search(
                problem, evaluator, std::move(best), search_rng, deadline
            );
        } else {
            vector<Schedule> final_seeds{best};
            best = iterated_jump_search(
                problem, evaluator, std::move(final_seeds), search_rng, deadline
            );
        }
    }
#ifdef LOCAL_DEBUG
    cerr << "final=" << best.cost << '\n';
#endif
    print_schedule(problem, best);
    return 0;
}
