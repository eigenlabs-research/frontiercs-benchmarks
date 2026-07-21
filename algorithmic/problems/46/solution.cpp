#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>
using namespace std;
namespace {
using Clock=chrono::steady_clock;
using i64=long long;
#ifndef SOLVER_TIME_MS
#define SOLVER_TIME_MS 820
#endif
struct RNG {
uint64_t state=0x46c0ffee12345678ULL;
uint64_t next(){
state += 0x9e3779b97f4a7c15ULL;
uint64_t value=state;
value=(value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
value=(value ^ (value >> 27)) * 0x94d049bb133111ebULL;
return value ^ (value >> 31);
}
int integer(int bound){ return static_cast<int>(next() % static_cast<uint64_t>(bound));}
};
struct Problem {
int jobs=0;
int machines=0;
int operations=0;
vector<int> machine;
vector<i64> duration;
vector<vector<int>> node_on_machine;
vector<i64> remaining_work;
vector<i64> job_load;
vector<i64> machine_load;
i64 lower_bound=0;
int node(int job,int operation) const { return job * machines + operation;}
int job(int node_id) const { return node_id / machines;}
int operation(int node_id) const { return node_id % machines;}
};
struct Schedule {
vector<vector<int>> order;
vector<int> position;
vector<int> critical_path;
i64 cost=numeric_limits<i64>::max();
};
uint64_t mix64(uint64_t value){
value += 0x9e3779b97f4a7c15ULL;
value=(value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
value=(value ^ (value >> 27)) * 0x94d049bb133111ebULL;
return value ^ (value >> 31);
}
void rebuild_positions(const Problem& problem,Schedule& schedule){
schedule.position.assign(problem.operations,-1);
for(int machine=0;machine<problem.machines;++machine){
for(int index=0;index<problem.jobs;++index){
schedule.position[schedule.order[machine][index]]=index;
}
}
}
void swap_adjacent(Schedule& schedule,int machine,int index){
int first=schedule.order[machine][index];
int second=schedule.order[machine][index + 1];
swap(schedule.order[machine][index],schedule.order[machine][index + 1]);
schedule.position[first]=index + 1;
schedule.position[second]=index;
}
class Evaluator {
public:
explicit Evaluator(const Problem& problem)
: problem_(problem),indegree_(problem.operations),head_(problem.operations),tail_(problem.operations),parent_(problem.operations),queue_(problem.operations){}
bool evaluate(Schedule& schedule,bool trace_path,uint64_t path_salt=0){
const int n=problem_.operations;
fill(indegree_.begin(),indegree_.end(),0);
fill(head_.begin(),head_.end(),0);
if(trace_path) fill(parent_.begin(),parent_.end(),-1);
for(int node_id=0;node_id<n;++node_id){
int operation=problem_.operation(node_id);
if(operation>0) ++indegree_[node_id];
if(schedule.position[node_id]>0) ++indegree_[node_id];
}
int front=0;
int back=0;
for(int node_id=0;node_id<n;++node_id){
if(indegree_[node_id]==0){
head_[node_id]=problem_.duration[node_id];
queue_[back++]=node_id;
}
}
int visited=0;
while(front<back){
int current=queue_[front++];
++visited;
int operation=problem_.operation(current);
if(operation + 1<problem_.machines){
relax(current,current + 1,trace_path,path_salt,back);
}
int machine=problem_.machine[current];
int position=schedule.position[current];
if(position + 1<problem_.jobs){
int successor=schedule.order[machine][position + 1];
relax(current,successor,trace_path,path_salt,back);
}
}
if(visited!=n){
schedule.cost=numeric_limits<i64>::max();
if(trace_path) schedule.critical_path.clear();
return false;
}
int finish=0;
for(int node_id=1;node_id<n;++node_id){
if(head_[node_id]>head_[finish] ||
(trace_path&&head_[node_id]==head_[finish] &&
tie_rank(node_id,path_salt)<tie_rank(finish,path_salt))){
finish=node_id;
}
}
schedule.cost=head_[finish];
if(trace_path){
schedule.critical_path.clear();
for(int node_id=finish;node_id!=-1;node_id=parent_[node_id]){
schedule.critical_path.push_back(node_id);
}
reverse(schedule.critical_path.begin(),schedule.critical_path.end());
}
return true;
}
void compute_tails(const Schedule& schedule){
fill(tail_.begin(),tail_.end(),0);
for(int index=problem_.operations - 1;index>=0;--index){
int current=queue_[index];
i64 successor_tail=0;
if(problem_.operation(current) + 1<problem_.machines){
successor_tail=max(successor_tail,tail_[current + 1]);
}
int machine=problem_.machine[current];
int position=schedule.position[current];
if(position + 1<problem_.jobs){
successor_tail=max(
successor_tail,tail_[schedule.order[machine][position + 1]]
);
}
tail_[current]=problem_.duration[current] + successor_tail;
}
}
i64 head(int node_id) const { return head_[node_id];}
i64 tail(int node_id) const { return tail_[node_id];}
private:
const Problem& problem_;
vector<int> indegree_;
vector<i64> head_;
vector<i64> tail_;
vector<int> parent_;
vector<int> queue_;
uint64_t tie_rank(int node_id,uint64_t salt) const {
return mix64(static_cast<uint64_t>(node_id) ^ salt);
}
void relax(
int current,int successor,bool trace_path,uint64_t path_salt,int& queue_back
){
i64 candidate=head_[current] + problem_.duration[successor];
if(candidate>head_[successor] ||
(trace_path&&candidate==head_[successor] &&
(parent_[successor]==-1 ||
tie_rank(current,path_salt)<tie_rank(parent_[successor],path_salt)))){
head_[successor]=candidate;
if(trace_path) parent_[successor]=current;
}
if(--indegree_[successor]==0) queue_[queue_back++]=successor;
}
};
class IncrementalEvaluator {
public:
IncrementalEvaluator(const Problem& problem,const Evaluator& base)
: problem_(problem),base_(base),mark_(problem.operations,0),indegree_(problem.operations,0),value_(problem.operations,0),queue_(problem.operations){
affected_.reserve(problem.operations);
}
pair<bool,i64> evaluate(
const Schedule& candidate,const vector<int>& roots
){
if(++generation_==0){
fill(mark_.begin(),mark_.end(),0);
generation_=1;
}
affected_.clear();
int front=0;
int back=0;
for(int root : roots){
if(marked(root)) continue;
mark_[root]=generation_;
queue_[back++]=root;
affected_.push_back(root);
}
while(front<back){
int current=queue_[front++];
int operation=problem_.operation(current);
if(operation + 1<problem_.machines){
add_descendant(current + 1,back);
}
int machine=problem_.machine[current];
int position=candidate.position[current];
if(position + 1<problem_.jobs){
add_descendant(candidate.order[machine][position + 1],back);
}
}
front=0;
back=0;
for(int node_id : affected_){
int degree=0;
int operation=problem_.operation(node_id);
if(operation>0&&marked(node_id - 1)) ++degree;
int machine=problem_.machine[node_id];
int position=candidate.position[node_id];
if(position>0 &&
marked(candidate.order[machine][position - 1])){
++degree;
}
indegree_[node_id]=degree;
if(degree==0) queue_[back++]=node_id;
}
i64 cost=0;
int visited=0;
while(front<back){
int node_id=queue_[front++];
++visited;
i64 predecessor=0;
int operation=problem_.operation(node_id);
if(operation>0){
int previous=node_id - 1;
predecessor=max(
predecessor,marked(previous) ? value_[previous] : base_.head(previous)
);
}
int machine=problem_.machine[node_id];
int position=candidate.position[node_id];
if(position>0){
int previous=candidate.order[machine][position - 1];
predecessor=max(
predecessor,marked(previous) ? value_[previous] : base_.head(previous)
);
}
value_[node_id]=predecessor + problem_.duration[node_id];
cost=max(cost,value_[node_id]);
if(operation + 1<problem_.machines){
release(node_id + 1,back);
}
if(position + 1<problem_.jobs){
release(candidate.order[machine][position + 1],back);
}
}
if(visited!=static_cast<int>(affected_.size())){
return {false,numeric_limits<i64>::max()};
}
for(int node_id=0;node_id<problem_.operations;++node_id){
if(!marked(node_id)) cost=max(cost,base_.head(node_id));
}
return {true,cost};
}
private:
const Problem& problem_;
const Evaluator& base_;
vector<uint32_t> mark_;
vector<int> indegree_;
vector<i64> value_;
vector<int> queue_;
vector<int> affected_;
uint32_t generation_=0;
bool marked(int node_id) const { return mark_[node_id]==generation_;}
void add_descendant(int node_id,int& back){
if(marked(node_id)) return;
mark_[node_id]=generation_;
queue_[back++]=node_id;
affected_.push_back(node_id);
}
void release(int node_id,int& back){
if(marked(node_id)&&--indegree_[node_id]==0){
queue_[back++]=node_id;
}
}
};
Schedule serial_dispatch(const Problem& problem,bool longest_processing_time){
vector<int> next_operation(problem.jobs,0);
vector<i64> job_ready(problem.jobs,0);
vector<i64> machine_ready(problem.machines,0);
Schedule schedule;
schedule.order.assign(problem.machines,{});
for(auto& sequence : schedule.order) sequence.reserve(problem.jobs);
for(int step=0;step<problem.operations;++step){
tuple<i64,i64,int,int> best{
numeric_limits<i64>::max(),numeric_limits<i64>::max(),-1,-1
};
for(int job=0;job<problem.jobs;++job){
int operation=next_operation[job];
if(operation==problem.machines) continue;
int node_id=problem.node(job,operation);
int machine=problem.machine[node_id];
i64 start=max(job_ready[job],machine_ready[machine]);
i64 priority=longest_processing_time ? -problem.duration[node_id]
: problem.duration[node_id];
best=min(best,make_tuple(start,priority,job,node_id));
}
auto [start,ignored,job,node_id]=best;
(void)ignored;
int machine=problem.machine[node_id];
i64 finish=start + problem.duration[node_id];
schedule.order[machine].push_back(node_id);
++next_operation[job];
job_ready[job]=finish;
machine_ready[machine]=finish;
}
rebuild_positions(problem,schedule);
return schedule;
}
Schedule giffler_thompson(const Problem& problem,int rule,uint64_t seed){
vector<int> next_operation(problem.jobs,0);
vector<i64> job_ready(problem.jobs,0);
vector<i64> machine_ready(problem.machines,0);
Schedule schedule;
schedule.order.assign(problem.machines,{});
for(auto& sequence : schedule.order) sequence.reserve(problem.jobs);
i64 maximum_job_load=*max_element(problem.job_load.begin(),problem.job_load.end());
i64 maximum_duration=*max_element(problem.duration.begin(),problem.duration.end());
auto coefficient=[&](int shift,int bias){
return bias + static_cast<int>((seed >> shift) & 255ULL) - 128;
};
const int remaining_coefficient=coefficient(0,-64);
const int duration_coefficient=coefficient(8,0);
const int ready_coefficient=coefficient(16,32);
const int operations_coefficient=coefficient(24,-32);
const int load_coefficient=coefficient(32,-16);
for(int step=0;step<problem.operations;++step){
int pivot=-1;
i64 pivot_finish=numeric_limits<i64>::max();
for(int job=0;job<problem.jobs;++job){
int operation=next_operation[job];
if(operation==problem.machines) continue;
int node_id=problem.node(job,operation);
int machine=problem.machine[node_id];
i64 start=max(job_ready[job],machine_ready[machine]);
i64 finish=start + problem.duration[node_id];
if(finish<pivot_finish ||
(finish==pivot_finish&&node_id<pivot)){
pivot_finish=finish;
pivot=node_id;
}
}
int pivot_machine=problem.machine[pivot];
int chosen=-1;
tuple<i64,uint64_t,int> chosen_key{
numeric_limits<i64>::max(),numeric_limits<uint64_t>::max(),-1
};
for(int job=0;job<problem.jobs;++job){
int operation=next_operation[job];
if(operation==problem.machines) continue;
int node_id=problem.node(job,operation);
if(problem.machine[node_id]!=pivot_machine) continue;
i64 start=max(job_ready[job],machine_ready[pivot_machine]);
if(start>=pivot_finish) continue;
i64 primary=0;
if(rule==0) primary=-problem.remaining_work[node_id];
else if(rule==1) primary=problem.duration[node_id];
else if(rule==2) primary=-problem.duration[node_id];
else if(rule==3){
primary=-(problem.remaining_work[node_id] - problem.duration[node_id]);
} else if(rule==4){
primary=problem.remaining_work[node_id];
} else if(rule==5){
primary=job_ready[job];
} else if(rule==6){
i64 remaining=problem.remaining_work[node_id];
i64 noise=static_cast<i64>(
mix64(seed ^ static_cast<uint64_t>(node_id)) %
static_cast<uint64_t>(1 + max<i64>(1,remaining / 2))
);
primary=-remaining + noise;
} else if(rule==7){
primary=problem.duration[node_id] + static_cast<i64>(
mix64(seed ^ static_cast<uint64_t>(node_id)) %
static_cast<uint64_t>(1 + problem.duration[node_id])
);
} else if(rule==8){
i64 tail=problem.remaining_work[node_id] - problem.duration[node_id];
i64 noise=static_cast<i64>(
mix64(seed ^ static_cast<uint64_t>(node_id)) %
static_cast<uint64_t>(1 + max<i64>(1,tail / 2))
);
primary=-tail + noise;
}else{
i64 remaining_feature =
4096.L * problem.remaining_work[node_id] / maximum_job_load;
i64 duration_feature =
4096.L * problem.duration[node_id] / maximum_duration;
i64 ready_feature=4096.L * job_ready[job] / maximum_job_load;
i64 operations_feature =
4096 * (problem.machines - operation) / problem.machines;
i64 load_feature=4096.L * problem.job_load[job] / maximum_job_load;
primary=remaining_coefficient * remaining_feature
+ duration_coefficient * duration_feature
+ ready_coefficient * ready_feature
+ operations_coefficient * operations_feature
+ load_coefficient * load_feature;
}
uint64_t secondary=mix64(
seed + static_cast<uint64_t>(node_id) * 0x9e3779b97f4a7c15ULL
);
auto key=make_tuple(primary,secondary,node_id);
if(key<chosen_key){
chosen_key=key;
chosen=node_id;
}
}
int job=problem.job(chosen);
int machine=problem.machine[chosen];
i64 start=max(job_ready[job],machine_ready[machine]);
i64 finish=start + problem.duration[chosen];
schedule.order[machine].push_back(chosen);
++next_operation[job];
job_ready[job]=finish;
machine_ready[machine]=finish;
}
rebuild_positions(problem,schedule);
return schedule;
}
class PriorityScheduleEvaluator {
public:
explicit PriorityScheduleEvaluator(const Problem& problem)
: problem_(problem),next_operation_(problem.jobs),job_ready_(problem.jobs),machine_ready_(problem.machines),bucket_head_(problem.machines),next_job_(problem.jobs),previous_job_(problem.jobs){}
i64 evaluate(const vector<int>& rank,Schedule* output=nullptr){
fill(next_operation_.begin(),next_operation_.end(),0);
fill(job_ready_.begin(),job_ready_.end(),0);
fill(machine_ready_.begin(),machine_ready_.end(),0);
fill(bucket_head_.begin(),bucket_head_.end(),-1);
fill(next_job_.begin(),next_job_.end(),-1);
fill(previous_job_.begin(),previous_job_.end(),-1);
for(int job=problem_.jobs - 1;job>=0;--job){
int machine=problem_.machine[problem_.node(job,0)];
next_job_[job]=bucket_head_[machine];
if(bucket_head_[machine]!=-1){
previous_job_[bucket_head_[machine]]=job;
}
bucket_head_[machine]=job;
}
if(output!=nullptr){
output->order.assign(problem_.machines,{});
for(auto& sequence : output->order){
sequence.reserve(problem_.jobs);
}
}
for(int step=0;step<problem_.operations;++step){
int pivot=-1;
i64 pivot_finish=numeric_limits<i64>::max();
for(int job=0;job<problem_.jobs;++job){
int operation=next_operation_[job];
if(operation==problem_.machines) continue;
int node_id=problem_.node(job,operation);
int machine=problem_.machine[node_id];
i64 finish=max(job_ready_[job],machine_ready_[machine]) +
problem_.duration[node_id];
if(finish<pivot_finish ||
(finish==pivot_finish&&node_id<pivot)){
pivot_finish=finish;
pivot=node_id;
}
}
int pivot_machine=problem_.machine[pivot];
int chosen=-1;
pair<int,int> best{numeric_limits<int>::max(),-1};
for(int job=bucket_head_[pivot_machine];job!=-1;
job=next_job_[job]){
int node_id=problem_.node(job,next_operation_[job]);
i64 start=max(job_ready_[job],machine_ready_[pivot_machine]);
if(start>=pivot_finish) continue;
pair<int,int> key{rank[job],node_id};
if(key<best){
best=key;
chosen=node_id;
}
}
int job=problem_.job(chosen);
int machine=problem_.machine[chosen];
i64 finish=max(job_ready_[job],machine_ready_[machine]) +
problem_.duration[chosen];
if(output!=nullptr) output->order[machine].push_back(chosen);
int previous=previous_job_[job];
int next=next_job_[job];
if(previous==-1) bucket_head_[machine]=next;
else next_job_[previous]=next;
if(next!=-1) previous_job_[next]=previous;
++next_operation_[job];
job_ready_[job]=finish;
machine_ready_[machine]=finish;
next_job_[job]=-1;
previous_job_[job]=-1;
if(next_operation_[job]<problem_.machines){
int next_machine=problem_.machine[
problem_.node(job,next_operation_[job])
];
next_job_[job]=bucket_head_[next_machine];
if(bucket_head_[next_machine]!=-1){
previous_job_[bucket_head_[next_machine]]=job;
}
bucket_head_[next_machine]=job;
}
}
i64 cost=*max_element(job_ready_.begin(),job_ready_.end());
if(output!=nullptr){
rebuild_positions(problem_,*output);
output->cost=cost;
}
return cost;
}
private:
const Problem& problem_;
vector<int> next_operation_;
vector<i64> job_ready_;
vector<i64> machine_ready_;
vector<int> bucket_head_;
vector<int> next_job_;
vector<int> previous_job_;
};
Schedule priority_permutation_search(
const Problem& problem,const vector<vector<int>>& starts,RNG& rng,const Clock::time_point& deadline
){
PriorityScheduleEvaluator evaluator(problem);
vector<int> rank(problem.jobs);
vector<int> best_permutation;
i64 best_cost=numeric_limits<i64>::max();
for(const vector<int>& permutation : starts){
for(int position=0;position<problem.jobs;++position){
rank[permutation[position]]=position;
}
i64 cost=evaluator.evaluate(rank);
if(cost<best_cost){
best_cost=cost;
best_permutation=permutation;
}
}
vector<int> permutation=best_permutation;
i64 current_cost=best_cost;
double average_duration=accumulate(
problem.duration.begin(),problem.duration.end(),0.0
) / max(1,problem.operations);
const double temperature=max(1.0,0.30 * average_duration);
int since_best=0;
while(Clock::now()<deadline){
vector<int> candidate=permutation;
if((rng.next() & 7)==0){
swap(
candidate[rng.integer(problem.jobs)],candidate[rng.integer(problem.jobs)]
);
}else{
int old_position=rng.integer(problem.jobs);
int job=candidate[old_position];
candidate.erase(candidate.begin() + old_position);
candidate.insert(candidate.begin() + rng.integer(problem.jobs),job);
}
for(int position=0;position<problem.jobs;++position){
rank[candidate[position]]=position;
}
i64 candidate_cost=evaluator.evaluate(rank);
double draw=static_cast<double>(rng.next() >> 11) /
static_cast<double>(uint64_t{1} << 53);
bool accept=candidate_cost<=current_cost ||
draw<exp(
-static_cast<double>(candidate_cost - current_cost) /
temperature
);
if(candidate_cost<best_cost){
best_cost=candidate_cost;
best_permutation=candidate;
since_best=0;
}else{
++since_best;
}
if(accept){
permutation=std::move(candidate);
current_cost=candidate_cost;
}
if(since_best>400){
permutation=best_permutation;
for(int count=0;count<min(6,problem.jobs);++count){
swap(
permutation[rng.integer(problem.jobs)],permutation[rng.integer(problem.jobs)]
);
}
for(int position=0;position<problem.jobs;++position){
rank[permutation[position]]=position;
}
current_cost=evaluator.evaluate(rank);
since_best=0;
}
}
for(int position=0;position<problem.jobs;++position){
rank[best_permutation[position]]=position;
}
Schedule result;
evaluator.evaluate(rank,&result);
return result;
}
vector<vector<int>> flow_priority_starts(
const Problem& problem,const vector<int>& neh,const vector<int>& descending
){
vector<int> consensus(problem.machines);
iota(consensus.begin(),consensus.end(),0);
vector<int> stage_sum(problem.machines,0);
for(int job=0;job<problem.jobs;++job){
for(int operation=0;operation<problem.machines;++operation){
stage_sum[problem.machine[problem.node(job,operation)]] += operation;
}
}
sort(consensus.begin(),consensus.end(),[&](int left,int right){
return tie(stage_sum[left],left)<tie(stage_sum[right],right);
});
vector<vector<int>> starts{neh,descending};
vector<int> reversed=descending;
reverse(reversed.begin(),reversed.end());
starts.push_back(std::move(reversed));
vector<int> palmer(problem.jobs);
vector<i64> slope(problem.jobs,0);
iota(palmer.begin(),palmer.end(),0);
for(int job=0;job<problem.jobs;++job){
for(int stage=0;stage<problem.machines;++stage){
slope[job] += static_cast<i64>(
problem.machines - 1 - 2 * stage
) * problem.duration[
problem.node_on_machine[job][consensus[stage]]
];
}
}
sort(palmer.begin(),palmer.end(),[&](int left,int right){
return tie(slope[left],left)>tie(slope[right],right);
});
starts.push_back(std::move(palmer));
for(int split=1;split<problem.machines;++split){
vector<pair<i64,int>> early;
vector<pair<i64,int>> late;
for(int job=0;job<problem.jobs;++job){
i64 first_sum=0;
i64 second_sum=0;
for(int stage=0;stage<split;++stage){
first_sum += problem.duration[
problem.node_on_machine[job][consensus[stage]]
];
}
for(int stage=split;stage<problem.machines;++stage){
second_sum += problem.duration[
problem.node_on_machine[job][consensus[stage]]
];
}
if(first_sum<second_sum) early.push_back({first_sum,job});
else late.push_back({second_sum,job});
}
sort(early.begin(),early.end());
sort(late.begin(),late.end(),greater<pair<i64,int>>());
vector<int> candidate;
candidate.reserve(problem.jobs);
for(auto [ignored,job] : early) candidate.push_back(job);
for(auto [ignored,job] : late) candidate.push_back(job);
starts.push_back(std::move(candidate));
}
return starts;
}
Problem reversed_problem(const Problem& problem){
Problem reversed;
reversed.jobs=problem.jobs;
reversed.machines=problem.machines;
reversed.operations=problem.operations;
reversed.machine.resize(problem.operations);
reversed.duration.resize(problem.operations);
reversed.node_on_machine.assign(
problem.jobs,vector<int>(problem.machines,-1)
);
reversed.job_load=problem.job_load;
reversed.machine_load=problem.machine_load;
reversed.lower_bound=problem.lower_bound;
for(int job=0;job<problem.jobs;++job){
for(int operation=0;operation<problem.machines;++operation){
int reversed_node=reversed.node(job,operation);
int original_node=problem.node(
job,problem.machines - 1 - operation
);
reversed.machine[reversed_node]=problem.machine[original_node];
reversed.duration[reversed_node]=problem.duration[original_node];
reversed.node_on_machine[job][reversed.machine[reversed_node]] =
reversed_node;
}
}
reversed.remaining_work.assign(problem.operations,0);
for(int job=0;job<problem.jobs;++job){
i64 remaining=0;
for(int operation=problem.machines - 1;operation>=0;--operation){
int node_id=reversed.node(job,operation);
remaining += reversed.duration[node_id];
reversed.remaining_work[node_id]=remaining;
}
}
return reversed;
}
Schedule mirror_reversed_schedule(
const Problem& problem,const Schedule& reversed_schedule
){
Schedule schedule;
schedule.order.assign(problem.machines,{});
for(int machine=0;machine<problem.machines;++machine){
schedule.order[machine].reserve(problem.jobs);
for(auto iterator=reversed_schedule.order[machine].rbegin();
iterator!=reversed_schedule.order[machine].rend();++iterator){
int reversed_node=*iterator;
int job=reversed_node / problem.machines;
int reversed_operation=reversed_node % problem.machines;
schedule.order[machine].push_back(problem.node(
job,problem.machines - 1 - reversed_operation
));
}
}
rebuild_positions(problem,schedule);
return schedule;
}
vector<int> schrage_order(const Problem& problem,int fixed_machine){
struct Task {
int node=-1;
i64 release=0;
i64 tail=0;
bool done=false;
};
vector<Task> tasks(problem.jobs);
for(int job=0;job<problem.jobs;++job){
int node_id=problem.node_on_machine[job][fixed_machine];
int operation=problem.operation(node_id);
i64 release=0;
for(int k=0;k<operation;++k){
release += problem.duration[problem.node(job,k)];
}
i64 tail=0;
for(int k=operation + 1;k<problem.machines;++k){
tail += problem.duration[problem.node(job,k)];
}
tasks[job]={node_id,release,tail,false};
}
vector<int> result;
result.reserve(problem.jobs);
i64 now=numeric_limits<i64>::max();
for(const Task& task : tasks) now=min(now,task.release);
while(static_cast<int>(result.size())<problem.jobs){
int chosen=-1;
for(int job=0;job<problem.jobs;++job){
const Task& task=tasks[job];
if(task.done||task.release>now) continue;
if(chosen==-1||task.tail>tasks[chosen].tail ||
(task.tail==tasks[chosen].tail &&
problem.duration[task.node]>problem.duration[tasks[chosen].node])){
chosen=job;
}
}
if(chosen==-1){
now=numeric_limits<i64>::max();
for(const Task& task : tasks){
if(!task.done) now=min(now,task.release);
}
continue;
}
tasks[chosen].done=true;
result.push_back(tasks[chosen].node);
now=max(now,tasks[chosen].release) + problem.duration[tasks[chosen].node];
}
return result;
}
Schedule complete_fixed_machine(
const Problem& problem,int,const vector<int>& fixed_order
){
vector<int> fixed_predecessor(problem.operations,-1);
vector<int> fixed_successor(problem.operations,-1);
for(int index=0;index + 1<problem.jobs;++index){
fixed_successor[fixed_order[index]]=fixed_order[index + 1];
fixed_predecessor[fixed_order[index + 1]]=fixed_order[index];
}
vector<int> indegree(problem.operations,0);
vector<char> scheduled(problem.operations,false);
vector<i64> job_ready(problem.jobs,0);
vector<i64> machine_ready(problem.machines,0);
for(int node_id=0;node_id<problem.operations;++node_id){
if(problem.operation(node_id)>0) ++indegree[node_id];
if(fixed_predecessor[node_id]!=-1) ++indegree[node_id];
}
Schedule schedule;
schedule.order.assign(problem.machines,{});
for(auto& sequence : schedule.order) sequence.reserve(problem.jobs);
for(int step=0;step<problem.operations;++step){
int chosen=-1;
tuple<i64,i64,int> best{
numeric_limits<i64>::max(),numeric_limits<i64>::max(),-1
};
for(int node_id=0;node_id<problem.operations;++node_id){
if(scheduled[node_id]||indegree[node_id]!=0) continue;
int job=problem.job(node_id);
int machine=problem.machine[node_id];
i64 start=max(job_ready[job],machine_ready[machine]);
auto key=make_tuple(start,-problem.remaining_work[node_id],node_id);
if(key<best){
best=key;
chosen=node_id;
}
}
int job=problem.job(chosen);
int machine=problem.machine[chosen];
i64 start=max(job_ready[job],machine_ready[machine]);
i64 finish=start + problem.duration[chosen];
scheduled[chosen]=true;
schedule.order[machine].push_back(chosen);
job_ready[job]=finish;
machine_ready[machine]=finish;
if(problem.operation(chosen) + 1<problem.machines) --indegree[chosen + 1];
if(fixed_successor[chosen]!=-1) --indegree[fixed_successor[chosen]];
}
rebuild_positions(problem,schedule);
return schedule;
}
class InsertionEvaluator {
public:
explicit InsertionEvaluator(const Problem& problem)
: problem_(problem),machine_successor_(problem.operations),indegree_(problem.operations),head_(problem.operations),tail_(problem.operations),queue_(problem.operations),topological_(problem.operations){}
bool evaluate(const vector<vector<int>>& partial_order,int focus){
fill(machine_successor_.begin(),machine_successor_.end(),-1);
fill(indegree_.begin(),indegree_.end(),0);
fill(head_.begin(),head_.end(),0);
for(int node_id=0;node_id<problem_.operations;++node_id){
if(problem_.operation(node_id)>0) ++indegree_[node_id];
}
for(const vector<int>& sequence : partial_order){
for(int index=0;index + 1<static_cast<int>(sequence.size());++index){
machine_successor_[sequence[index]]=sequence[index + 1];
++indegree_[sequence[index + 1]];
}
}
int front=0;
int back=0;
for(int node_id=0;node_id<problem_.operations;++node_id){
if(indegree_[node_id]==0){
head_[node_id]=problem_.duration[node_id];
queue_[back++]=node_id;
}
}
int visited=0;
while(front<back){
int current=queue_[front++];
topological_[visited++]=current;
if(problem_.operation(current) + 1<problem_.machines){
relax_head(current,current + 1,back);
}
if(machine_successor_[current]!=-1){
relax_head(current,machine_successor_[current],back);
}
}
if(visited!=problem_.operations) return false;
fill(tail_.begin(),tail_.end(),0);
for(int index=problem_.operations - 1;index>=0;--index){
int current=topological_[index];
i64 successor_tail=0;
if(problem_.operation(current) + 1<problem_.machines){
successor_tail=max(successor_tail,tail_[current + 1]);
}
if(machine_successor_[current]!=-1){
successor_tail=max(successor_tail,tail_[machine_successor_[current]]);
}
tail_[current]=problem_.duration[current] + successor_tail;
}
through_focus_=head_[focus] + tail_[focus] - problem_.duration[focus];
makespan_=*max_element(head_.begin(),head_.end());
return true;
}
i64 through_focus() const { return through_focus_;}
i64 makespan() const { return makespan_;}
private:
const Problem& problem_;
vector<int> machine_successor_;
vector<int> indegree_;
vector<i64> head_;
vector<i64> tail_;
vector<int> queue_;
vector<int> topological_;
i64 through_focus_=0;
i64 makespan_=0;
void relax_head(int current,int successor,int& queue_back){
head_[successor]=max(
head_[successor],head_[current] + problem_.duration[successor]
);
if(--indegree_[successor]==0) queue_[queue_back++]=successor;
}
};
Schedule greedy_insertion_seed(
const Problem& problem,Evaluator& full_evaluator,int first_job,int beam_width,const Clock::time_point& deadline
){
vector<vector<int>> partial_order(problem.machines);
vector<char> inserted(problem.operations,false);
for(int operation=0;operation<problem.machines;++operation){
int node_id=problem.node(first_job,operation);
partial_order[problem.machine[node_id]].push_back(node_id);
inserted[node_id]=true;
}
vector<int> remaining;
remaining.reserve(problem.operations - problem.machines);
for(int node_id=0;node_id<problem.operations;++node_id){
if(!inserted[node_id]) remaining.push_back(node_id);
}
stable_sort(remaining.begin(),remaining.end(),[&](int left,int right){
if(problem.duration[left]!=problem.duration[right]){
return problem.duration[left]>problem.duration[right];
}
return problem.remaining_work[left]>problem.remaining_work[right];
});
struct BeamState {
vector<vector<int>> order;
tuple<i64,i64,int> key;
};
vector<BeamState> beam;
beam.push_back({std::move(partial_order),{0,0,0}});
InsertionEvaluator evaluator(problem);
for(int node_id : remaining){
if(Clock::now()>=deadline) return serial_dispatch(problem,false);
int machine=problem.machine[node_id];
vector<BeamState> next_beam;
for(BeamState& state : beam){
vector<int>& sequence=state.order[machine];
for(int position=0;position<=static_cast<int>(sequence.size());++position){
sequence.insert(sequence.begin() + position,node_id);
if(evaluator.evaluate(state.order,node_id)){
auto key=make_tuple(
evaluator.through_focus(),evaluator.makespan(),position
);
if(static_cast<int>(next_beam.size())<beam_width ||
key<next_beam.back().key){
BeamState candidate{state.order,key};
auto insertion=lower_bound(
next_beam.begin(),next_beam.end(),key,[](const BeamState& existing,const tuple<i64,i64,int>& value){
return existing.key<value;
}
);
next_beam.insert(insertion,std::move(candidate));
if(static_cast<int>(next_beam.size())>beam_width){
next_beam.pop_back();
}
}
}
sequence.erase(sequence.begin() + position);
}
}
if(next_beam.empty()) return serial_dispatch(problem,false);
beam=std::move(next_beam);
}
Schedule best_schedule;
for(BeamState& state : beam){
Schedule schedule;
schedule.order=std::move(state.order);
rebuild_positions(problem,schedule);
full_evaluator.evaluate(schedule,false);
if(schedule.cost<best_schedule.cost) best_schedule=std::move(schedule);
}
return best_schedule;
}
Schedule common_permutation(const Problem& problem,const vector<int>& permutation){
Schedule schedule;
schedule.order.assign(problem.machines,vector<int>(problem.jobs));
for(int machine=0;machine<problem.machines;++machine){
for(int index=0;index<problem.jobs;++index){
int job=permutation[index];
schedule.order[machine][index]=problem.node_on_machine[job][machine];
}
}
rebuild_positions(problem,schedule);
return schedule;
}
i64 common_cost(
const Problem& problem,const vector<int>& permutation
){
vector<i64> machine_ready(problem.machines,0);
i64 makespan=0;
for(int job : permutation){
i64 job_ready=0;
for(int operation=0;operation<problem.machines;++operation){
int node_id=problem.node(job,operation);
int machine=problem.machine[node_id];
job_ready=max(job_ready,machine_ready[machine]) + problem.duration[node_id];
machine_ready[machine]=job_ready;
makespan=max(makespan,job_ready);
}
}
return makespan;
}
i64 priority_cost(
const Problem& problem,const vector<int>& permutation
){
vector<int> next_operation(problem.jobs,0);
vector<i64> job_ready(problem.jobs,0),machine_ready(problem.machines,0);
for(int step=0;step<static_cast<int>(permutation.size()) * problem.machines;++step){
int pivot_machine=-1;
i64 pivot_finish=numeric_limits<i64>::max();
for(int job : permutation){
int operation=next_operation[job];
if(operation==problem.machines) continue;
int node_id=problem.node(job,operation),machine=problem.machine[node_id];
i64 finish=max(job_ready[job],machine_ready[machine]) + problem.duration[node_id];
if(finish<pivot_finish){
pivot_finish=finish;
pivot_machine=machine;
}
}
int chosen=-1;
for(int job : permutation){
int operation=next_operation[job];
if(operation==problem.machines) continue;
int node_id=problem.node(job,operation);
if(problem.machine[node_id]==pivot_machine &&
max(job_ready[job],machine_ready[pivot_machine])<pivot_finish){
chosen=job;
break;
}
}
int node_id=problem.node(chosen,next_operation[chosen]++);
job_ready[chosen]=max(job_ready[chosen],machine_ready[pivot_machine]) +
problem.duration[node_id];
machine_ready[pivot_machine]=job_ready[chosen];
}
return *max_element(job_ready.begin(),job_ready.end());
}
vector<int> neh_permutation(
const Problem& problem,const vector<int>& job_order,const Clock::time_point& deadline,bool priority_decode=false
){
vector<int> permutation;
permutation.reserve(problem.jobs);
for(int job : job_order){
if(Clock::now()>=deadline){
permutation.push_back(job);
continue;
}
i64 best_cost=numeric_limits<i64>::max();
int best_position=0;
for(int position=0;position<=static_cast<int>(permutation.size());++position){
if(Clock::now()>=deadline) break;
vector<int> candidate=permutation;
candidate.insert(candidate.begin() + position,job);
i64 cost=priority_decode
? priority_cost(problem,candidate)
: common_cost(problem,candidate);
if(cost<best_cost){
best_cost=cost;
best_position=position;
}
}
permutation.insert(permutation.begin() + best_position,job);
}
return permutation;
}
vector<int> improve_perm(
const Problem& problem,vector<int> permutation,const Clock::time_point& search_end,bool priority_decode
){
vector<int> best_permutation=permutation;
auto cost=[&](const vector<int>& order){
return priority_decode ? priority_cost(problem,order)
: common_cost(problem,order);
};
i64 best_cost=cost(permutation);
i64 current_cost=best_cost;
RNG rng;
rng.state=0x3141592653589793ULL;
while(Clock::now()<search_end){
vector<int> candidate=permutation,removed;
int destruction=min(5,problem.jobs);
for(int count=0;count<destruction;++count){
int position=rng.integer(candidate.size());
removed.push_back(candidate[position]);
candidate.erase(candidate.begin() + position);
}
for(int job : removed){
i64 insertion_cost=numeric_limits<i64>::max();
int insertion_position=0;
for(int position=0;position<=static_cast<int>(candidate.size());++position){
if(Clock::now()>=search_end) break;
candidate.insert(candidate.begin() + position,job);
i64 candidate_insertion_cost=cost(candidate);
candidate.erase(candidate.begin() + position);
if(candidate_insertion_cost<insertion_cost){
insertion_cost=candidate_insertion_cost;
insertion_position=position;
}
}
candidate.insert(candidate.begin() + insertion_position,job);
}
i64 candidate_cost=cost(candidate);
if(candidate_cost<best_cost){
best_cost=candidate_cost;
best_permutation=candidate;
}
if(candidate_cost<=current_cost||rng.integer(16)==0){
permutation=std::move(candidate);
current_cost=candidate_cost;
}else if((rng.next() & 31)==0){
permutation=best_permutation;
current_cost=best_cost;
}
}
return best_permutation;
}
uint64_t schedule_hash(const Schedule& schedule){
uint64_t hash=1469598103934665603ULL;
for(const auto& sequence : schedule.order){
for(int node_id : sequence){
hash=(hash ^ static_cast<uint64_t>(node_id + 1)) * 1099511628211ULL;
}
}
return hash;
}
void add_seed(
const Problem&,Evaluator& evaluator,vector<Schedule>& seeds,vector<uint64_t>& hashes,Schedule schedule
){
if(!evaluator.evaluate(schedule,false)) return;
uint64_t hash=schedule_hash(schedule);
if(find(hashes.begin(),hashes.end(),hash)!=hashes.end()) return;
hashes.push_back(hash);
seeds.push_back(std::move(schedule));
}
struct Move {
int machine=-1;
int index=-1;
};
struct RelocateMove {
int machine=-1;
int from=-1;
int to=-1;
};
struct Block {
int first=0;
int last=0;
};
vector<Block> critical_blocks(const Problem& problem,const Schedule& schedule){
vector<Block> blocks;
const vector<int>& path=schedule.critical_path;
for(int first=0;first<static_cast<int>(path.size());){
int last=first;
int machine=problem.machine[path[first]];
while(last + 1<static_cast<int>(path.size()) &&
problem.machine[path[last + 1]]==machine &&
schedule.position[path[last + 1]]==schedule.position[path[last]] + 1){
++last;
}
blocks.push_back({first,last});
first=last + 1;
}
return blocks;
}
vector<Move> n5_moves(const Problem& problem,const Schedule& schedule){
vector<Block> blocks=critical_blocks(problem,schedule);
vector<Move> moves;
for(int block_index=0;block_index<static_cast<int>(blocks.size());++block_index){
const Block& block=blocks[block_index];
if(block.first==block.last) continue;
int machine=problem.machine[schedule.critical_path[block.first]];
int first_index=schedule.position[schedule.critical_path[block.first]];
int last_index=schedule.position[schedule.critical_path[block.last - 1]];
if(blocks.size()==1||block_index>0){
moves.push_back({machine,first_index});
}
if((blocks.size()==1||block_index + 1<static_cast<int>(blocks.size())) &&
last_index!=first_index){
moves.push_back({machine,last_index});
}
}
return moves;
}
vector<Move> n1_moves(const Problem& problem,const Schedule& schedule){
vector<Move> moves;
for(const Block& block : critical_blocks(problem,schedule)){
int machine=problem.machine[schedule.critical_path[block.first]];
for(int path_index=block.first;path_index<block.last;++path_index){
moves.push_back({machine,schedule.position[schedule.critical_path[path_index]]});
}
}
return moves;
}
void relocate_operation(Schedule& schedule,const RelocateMove& move){
vector<int>& sequence=schedule.order[move.machine];
int moved=sequence[move.from];
if(move.from<move.to){
for(int index=move.from;index<move.to;++index){
sequence[index]=sequence[index + 1];
schedule.position[sequence[index]]=index;
}
}else{
for(int index=move.from;index>move.to;--index){
sequence[index]=sequence[index - 1];
schedule.position[sequence[index]]=index;
}
}
sequence[move.to]=moved;
schedule.position[moved]=move.to;
}
vector<RelocateMove> n7_moves(
const Problem& problem,const Schedule& schedule
){
const vector<Block> blocks=critical_blocks(problem,schedule);
vector<RelocateMove> moves;
vector<uint64_t> seen;
auto add=[&](int machine,int from,int to){
if(from==to) return;
uint64_t key=(static_cast<uint64_t>(machine) << 32) |
(static_cast<uint64_t>(from) << 16) |
static_cast<uint64_t>(to);
if(find(seen.begin(),seen.end(),key)!=seen.end()) return;
seen.push_back(key);
moves.push_back({machine,from,to});
};
for(int block_index=0;
block_index<static_cast<int>(blocks.size());++block_index){
const Block& block=blocks[block_index];
if(block.first==block.last) continue;
int machine=problem.machine[schedule.critical_path[block.first]];
int first=schedule.position[schedule.critical_path[block.first]];
int last=schedule.position[schedule.critical_path[block.last]];
if(block_index>0) add(machine,first,first + 1);
if(block_index + 1<static_cast<int>(blocks.size())){
add(machine,last,last - 1);
}
for(int index=first + 1;index<last;++index){
if(block_index>0) add(machine,index,first);
if(block_index + 1<static_cast<int>(blocks.size())){
add(machine,index,last);
}
}
for(int index=first + 2;index<=last;++index){
add(machine,first,index);
}
for(int index=first;index + 1<last;++index){
add(machine,last,index);
}
}
return moves;
}
bool perturb(
const Problem& problem,Evaluator& evaluator,Schedule& schedule,RNG& rng,int steps,const Clock::time_point& deadline
){
for(int step=0;step<steps;++step){
if(Clock::now()>=deadline) return false;
evaluator.evaluate(schedule,true,rng.next());
vector<Move> moves=n1_moves(problem,schedule);
if(moves.empty()) return false;
int offset=rng.integer(static_cast<int>(moves.size()));
bool applied=false;
for(int attempt=0;attempt<static_cast<int>(moves.size());++attempt){
if((attempt&7)==0&&Clock::now()>=deadline) return false;
Move move=moves[(offset + attempt) % moves.size()];
swap_adjacent(schedule,move.machine,move.index);
if(evaluator.evaluate(schedule,true,rng.next())){
applied=true;
break;
}
swap_adjacent(schedule,move.machine,move.index);
}
if(!applied){
evaluator.evaluate(schedule,true,rng.next());
return false;
}
}
return true;
}
Schedule annealing_search(
const Problem& problem,Evaluator& evaluator,vector<Schedule> seeds,RNG& rng,const Clock::time_point& deadline
){
sort(seeds.begin(),seeds.end(),[](const Schedule& left,const Schedule& right){
return left.cost<right.cost;
});
Schedule best=seeds.front();
if(best.cost==problem.lower_bound) return best;
Schedule current=best;
double average_duration=accumulate(
problem.duration.begin(),problem.duration.end(),0.0
) / max(1,problem.operations);
double temperature=average_duration;
int stage_iteration=0;
while(Clock::now()<deadline){
evaluator.evaluate(current,true,rng.next());
vector<Move> moves=n5_moves(problem,current);
if(moves.empty()){
current=best;
perturb(problem,evaluator,current,rng,5,deadline);
continue;
}
const Move& selected=moves[rng.integer(static_cast<int>(moves.size()))];
i64 original_cost=current.cost;
swap_adjacent(current,selected.machine,selected.index);
bool valid=evaluator.evaluate(current,false);
i64 candidate_cost=current.cost;
bool accept=valid&&candidate_cost<=original_cost;
if(valid&&!accept){
double probability=exp(
-static_cast<double>(candidate_cost - original_cost) /
max(1.0,temperature)
);
double draw=static_cast<double>(rng.next() >> 11) /
static_cast<double>(uint64_t{1} << 53);
accept=draw<probability;
}
if(!accept){
swap_adjacent(current,selected.machine,selected.index);
current.cost=original_cost;
} else if(current.cost<best.cost){
best=current;
if(best.cost==problem.lower_bound) return best;
}
temperature *= 0.9999;
if(++stage_iteration>=16000 ||
temperature<average_duration * 0.05){
stage_iteration=0;
current=best;
perturb(
problem,evaluator,current,rng,2 + rng.integer(6),deadline
);
temperature=average_duration;
}
}
return best;
}
i64 estimate_adjacent_swap(
const Problem& problem,const Schedule& schedule,const Evaluator& evaluator,const Move& move
){
int first=schedule.order[move.machine][move.index];
int second=schedule.order[move.machine][move.index + 1];
i64 machine_predecessor_head=move.index==0
? 0
: evaluator.head(schedule.order[move.machine][move.index - 1]);
i64 first_job_predecessor_head=problem.operation(first)==0
? 0
: evaluator.head(first - 1);
i64 second_job_predecessor_head=problem.operation(second)==0
? 0
: evaluator.head(second - 1);
i64 second_start=max(machine_predecessor_head,second_job_predecessor_head);
i64 second_completion=second_start + problem.duration[second];
i64 first_start=max(first_job_predecessor_head,second_completion);
i64 machine_successor_tail=move.index + 2>=problem.jobs
? 0
: evaluator.tail(schedule.order[move.machine][move.index + 2]);
i64 first_job_successor_tail=problem.operation(first) + 1==problem.machines
? 0
: evaluator.tail(first + 1);
i64 second_job_successor_tail=problem.operation(second) + 1==problem.machines
? 0
: evaluator.tail(second + 1);
i64 first_tail_after=max(machine_successor_tail,first_job_successor_tail);
i64 second_tail_after=max(
second_job_successor_tail,problem.duration[first] + first_tail_after
);
return max(
second_start + problem.duration[second] + second_tail_after,first_start + problem.duration[first] + first_tail_after
);
}
i64 estimate_relocation(
const Problem& problem,const Schedule& schedule,const Evaluator& evaluator,const RelocateMove& move
){
int first=min(move.from,move.to);
int last=max(move.from,move.to);
vector<int> changed(
schedule.order[move.machine].begin() + first,schedule.order[move.machine].begin() + last + 1
);
int local_from=move.from - first;
int local_to=move.to - first;
int moved=changed[local_from];
changed.erase(changed.begin() + local_from);
changed.insert(changed.begin() + local_to,moved);
vector<i64> changed_head(changed.size());
i64 machine_head=first==0
? 0
: evaluator.head(schedule.order[move.machine][first - 1]);
for(int index=0;index<static_cast<int>(changed.size());++index){
int node_id=changed[index];
i64 job_head=problem.operation(node_id)==0
? 0
: evaluator.head(node_id - 1);
changed_head[index]=max(machine_head,job_head) +
problem.duration[node_id];
machine_head=changed_head[index];
}
i64 estimate=0;
i64 machine_tail=last + 1==problem.jobs
? 0
: evaluator.tail(schedule.order[move.machine][last + 1]);
for(int index=static_cast<int>(changed.size()) - 1;index>=0;--index){
int node_id=changed[index];
i64 job_tail=problem.operation(node_id) + 1==problem.machines
? 0
: evaluator.tail(node_id + 1);
i64 changed_tail=problem.duration[node_id] +
max(machine_tail,job_tail);
estimate=max(
estimate,changed_head[index] + changed_tail - problem.duration[node_id]
);
machine_tail=changed_tail;
}
return estimate;
}
bool relocation_is_tabu(
const Problem& problem,const Schedule& schedule,const RelocateMove& move,const vector<int>& tabu_until,int iteration
){
const vector<int>& sequence=schedule.order[move.machine];
int moved=sequence[move.from];
if(move.from<move.to){
for(int index=move.from + 1;index<=move.to;++index){
int crossed=sequence[index];
if(tabu_until[crossed * problem.operations + moved]>iteration){
return true;
}
}
}else{
for(int index=move.to;index<move.from;++index){
int crossed=sequence[index];
if(tabu_until[moved * problem.operations + crossed]>iteration){
return true;
}
}
}
return false;
}
void forbid_relocation_inverse(
const Problem& problem,const Schedule& schedule,const RelocateMove& move,vector<int>& tabu_until,int expires
){
const vector<int>& sequence=schedule.order[move.machine];
int moved=sequence[move.from];
if(move.from<move.to){
for(int index=move.from + 1;index<=move.to;++index){
int crossed=sequence[index];
tabu_until[moved * problem.operations + crossed]=expires;
}
}else{
for(int index=move.to;index<move.from;++index){
int crossed=sequence[index];
tabu_until[crossed * problem.operations + moved]=expires;
}
}
}
Schedule n7_tabu_search(
const Problem& problem,Evaluator& evaluator,vector<Schedule> seeds,RNG& rng,const Clock::time_point& deadline
){
sort(seeds.begin(),seeds.end(),[](const Schedule& left,const Schedule& right){
return left.cost<right.cost;
});
Schedule best=seeds.front();
Schedule current=best;
if(best.cost==problem.lower_bound) return best;
const int n=problem.operations;
const int tenure_base=max(9,(problem.jobs + problem.machines) / 5);
const int tenure_spread=max(7,problem.jobs / 3);
const int restart_after=max(600,min(2400,32 * problem.jobs));
const int elite_count=min<int>(16,seeds.size());
vector<int> tabu_until(n * n,0);
int iteration=0;
int since_best=0;
int next_elite=1;
while(Clock::now()<deadline){
++iteration;
if(!evaluator.evaluate(current,true,rng.next())){
current=best;
rebuild_positions(problem,current);
continue;
}
evaluator.compute_tails(current);
vector<RelocateMove> moves=n7_moves(problem,current);
if(moves.empty()){
current=best;
perturb(problem,evaluator,current,rng,3 + rng.integer(6),deadline);
fill(tabu_until.begin(),tabu_until.end(),0);
continue;
}
struct RankedMove {
i64 estimate;
uint64_t tie;
bool tabu;
RelocateMove move;
};
vector<RankedMove> ranked;
ranked.reserve(moves.size());
for(int move_index=0;move_index<static_cast<int>(moves.size());++move_index){
if((move_index & 15)==0&&Clock::now()>=deadline) return best;
const RelocateMove& move=moves[move_index];
ranked.push_back({
estimate_relocation(problem,current,evaluator,move),rng.next(),relocation_is_tabu(
problem,current,move,tabu_until,iteration
),move
});
}
sort(ranked.begin(),ranked.end(),[](const RankedMove& left,const RankedMove& right){
return tie(left.estimate,left.tabu,left.tie) <
tie(right.estimate,right.tabu,right.tie);
});
bool found=false;
i64 chosen_cost=numeric_limits<i64>::max();
uint64_t chosen_tie=numeric_limits<uint64_t>::max();
RelocateMove chosen;
int exact_non_tabu=0;
int exact_tabu=0;
const i64 original_cost=current.cost;
for(const RankedMove& candidate : ranked){
if(Clock::now()>=deadline) break;
if(candidate.tabu){
if((problem.operations>=160&&problem.operations<=250&&candidate.estimate>=best.cost)||exact_tabu>=1) continue;
++exact_tabu;
}else{
if(exact_non_tabu>=(problem.operations>=160&&problem.operations<=250 ? 4 : 1)) continue;
++exact_non_tabu;
}
relocate_operation(current,candidate.move);
bool feasible=evaluator.evaluate(current,false);
i64 exact_cost=current.cost;
RelocateMove inverse{
candidate.move.machine,candidate.move.to,candidate.move.from
};
relocate_operation(current,inverse);
current.cost=original_cost;
bool admissible=feasible &&
(!candidate.tabu||exact_cost<best.cost);
if(admissible &&
(exact_cost<chosen_cost ||
(exact_cost==chosen_cost&&candidate.tie<chosen_tie))){
found=true;
chosen_cost=exact_cost;
chosen_tie=candidate.tie;
chosen=candidate.move;
}
}
if(!found){
current=best;
perturb(problem,evaluator,current,rng,4 + rng.integer(8),deadline);
fill(tabu_until.begin(),tabu_until.end(),0);
since_best=0;
continue;
}
int expires=iteration + tenure_base + rng.integer(tenure_spread);
forbid_relocation_inverse(
problem,current,chosen,tabu_until,expires
);
relocate_operation(current,chosen);
if(!evaluator.evaluate(current,false)){
current=best;
fill(tabu_until.begin(),tabu_until.end(),0);
continue;
}
if(current.cost<best.cost){
best=current;
since_best=0;
if(best.cost==problem.lower_bound) break;
}else{
++since_best;
}
if(since_best>=restart_after){
if(next_elite<elite_count){
current=seeds[next_elite++];
}else{
current=best;
perturb(
problem,evaluator,current,rng,5 + rng.integer(11),deadline
);
next_elite=1;
}
fill(tabu_until.begin(),tabu_until.end(),0);
since_best=0;
}
}
return best;
}
Schedule estimated_n5_tabu(
const Problem& problem,Evaluator& evaluator,vector<Schedule> seeds,RNG& rng,const Clock::time_point& deadline,bool long_trajectory
){
sort(seeds.begin(),seeds.end(),[](const Schedule& left,const Schedule& right){
return left.cost<right.cost;
});
Schedule best=seeds.front();
if(best.cost==problem.lower_bound) return best;
Schedule current=best;
int iteration=0;
int since_best=0;
const int n=problem.operations;
const int tabu_base=long_trajectory ? 10 : 5;
const int tabu_spread=long_trajectory ? 8 : 6;
const int restart_after=long_trajectory ? 1600 : 400;
vector<int> tabu_until(n * n,0);
while(Clock::now()<deadline){
++iteration;
if(!evaluator.evaluate(current,true,rng.next())){
current=best;
continue;
}
evaluator.compute_tails(current);
vector<Move> moves=n5_moves(problem,current);
if(moves.empty()){
current=best;
perturb(problem,evaluator,current,rng,3 + rng.integer(6),deadline);
fill(tabu_until.begin(),tabu_until.end(),0);
continue;
}
struct RankedMove {
i64 estimate;
uint64_t tie;
Move move;
};
vector<RankedMove> ranked;
ranked.reserve(moves.size());
for(int move_index=0;move_index<static_cast<int>(moves.size());++move_index){
if((move_index & 15)==0&&Clock::now()>=deadline) return best;
const Move& move=moves[move_index];
int first=current.order[move.machine][move.index];
int second=current.order[move.machine][move.index + 1];
int resulting_arc=second * n + first;
if(tabu_until[resulting_arc]>iteration) continue;
ranked.push_back({
estimate_adjacent_swap(problem,current,evaluator,move),rng.next(),move
});
}
sort(ranked.begin(),ranked.end(),[](const RankedMove& left,const RankedMove& right){
return tie(left.estimate,left.tie)<tie(right.estimate,right.tie);
});
bool applied=false;
for(const RankedMove& candidate : ranked){
if(Clock::now()>=deadline) return best;
const Move move=candidate.move;
int first=current.order[move.machine][move.index];
int second=current.order[move.machine][move.index + 1];
swap_adjacent(current,move.machine,move.index);
if(evaluator.evaluate(current,true,rng.next())){
tabu_until[first * n + second] =
iteration + tabu_base + rng.integer(tabu_spread);
applied=true;
break;
}
swap_adjacent(current,move.machine,move.index);
}
if(!applied){
current=best;
perturb(problem,evaluator,current,rng,4 + rng.integer(8),deadline);
fill(tabu_until.begin(),tabu_until.end(),0);
continue;
}
if(current.cost<best.cost){
best=current;
since_best=0;
if(best.cost==problem.lower_bound) break;
}else{
++since_best;
}
if(since_best>=restart_after){
current=best;
perturb(problem,evaluator,current,rng,4 + rng.integer(10),deadline);
fill(tabu_until.begin(),tabu_until.end(),0);
since_best=0;
}
}
return best;
}
Schedule incremental_n5_tabu(
const Problem& problem,Evaluator& evaluator,vector<Schedule> seeds,RNG& rng,const Clock::time_point& deadline,bool long_trajectory
){
sort(seeds.begin(),seeds.end(),[](const Schedule& left,const Schedule& right){
return left.cost<right.cost;
});
Schedule best=seeds.front();
if(best.cost==problem.lower_bound) return best;
Schedule current=best;
const int n=problem.operations;
vector<int> tabu_until(n * n,0);
int iteration=0;
int since_best=0;
int restart=0;
const int restart_after=long_trajectory ? 900 : 500;
const int top_seeds=min<int>(12,seeds.size());
IncrementalEvaluator incremental(problem,evaluator);
vector<int> roots;
roots.reserve(2);
while(Clock::now()<deadline){
++iteration;
if(!evaluator.evaluate(current,true,rng.next())){
current=best;
continue;
}
evaluator.compute_tails(current);
vector<Move> moves=n5_moves(problem,current);
if(moves.empty()){
current=best;
perturb(problem,evaluator,current,rng,3 + rng.integer(6),deadline);
fill(tabu_until.begin(),tabu_until.end(),0);
continue;
}
struct Candidate {
long double metric=numeric_limits<long double>::infinity();
i64 cost=numeric_limits<i64>::max();
uint64_t tie=numeric_limits<uint64_t>::max();
Move move;
int first=-1;
int second=-1;
} chosen;
for(int move_index=0;move_index<static_cast<int>(moves.size());++move_index){
if((move_index & 7)==0&&Clock::now()>=deadline) return best;
const Move& move=moves[move_index];
int first=current.order[move.machine][move.index];
int second=current.order[move.machine][move.index + 1];
i64 estimate=estimate_adjacent_swap(
problem,current,evaluator,move
);
swap_adjacent(current,move.machine,move.index);
roots.clear();
roots.push_back(second);
roots.push_back(first);
auto [valid,candidate_cost]=incremental.evaluate(current,roots);
swap_adjacent(current,move.machine,move.index);
if(!valid) continue;
bool tabu=tabu_until[second * n + first]>iteration;
if(tabu&&candidate_cost>=best.cost) continue;
uint64_t tie_value=rng.next();
long double metric=100.L*estimate+50.L*candidate_cost;
if(make_pair(metric,tie_value) <
make_pair(chosen.metric,chosen.tie)){
chosen={
metric,candidate_cost,tie_value,move,first,second
};
}
}
if(chosen.first==-1){
current=best;
perturb(problem,evaluator,current,rng,4 + rng.integer(8),deadline);
fill(tabu_until.begin(),tabu_until.end(),0);
since_best=0;
continue;
}
swap_adjacent(current,chosen.move.machine,chosen.move.index);
if(!evaluator.evaluate(current,true,rng.next())){
swap_adjacent(current,chosen.move.machine,chosen.move.index);
evaluator.evaluate(current,true,rng.next());
continue;
}
int base_tenure=long_trajectory ? 18 : 11;
int spread=long_trajectory ? 15 : 12;
if(since_best>restart_after / 2) base_tenure += 5;
tabu_until[chosen.first * n + chosen.second] =
iteration + base_tenure + rng.integer(spread);
if(current.cost<best.cost){
best=current;
since_best=0;
if(best.cost==problem.lower_bound) break;
}else{
++since_best;
}
if(since_best>=restart_after){
++restart;
if(restart % 3==0&&top_seeds>1){
current=seeds[1 + (restart / 3 - 1) % (top_seeds - 1)];
}else{
current=best;
perturb(
problem,evaluator,current,rng,5 + rng.integer(long_trajectory ? 15 : 10),deadline
);
}
fill(tabu_until.begin(),tabu_until.end(),0);
since_best=0;
}
}
return best;
}
Schedule single_continuation(
const Problem& problem,const Schedule& incumbent,const Clock::time_point& deadline
){
if(Clock::now()>=deadline||incumbent.cost==problem.lower_bound){
return incumbent;
}
Evaluator local_evaluator(problem);
RNG local_rng;
local_rng.state=mix64(schedule_hash(incumbent) ^ 0x9e3779b97f4a7c15ULL);
vector<Schedule> local_seeds{incumbent};
Schedule perturbed=incumbent;
RNG perturb_rng=local_rng;
perturb_rng.state^=0x53a9bf128ce476d0ULL;
bool ok=perturb(problem,local_evaluator,perturbed,perturb_rng,5,deadline);
if(ok) local_seeds.push_back(std::move(perturbed));
Schedule candidate=estimated_n5_tabu(
problem,local_evaluator,std::move(local_seeds),local_rng,deadline,false
);
return candidate.cost<incumbent.cost ? candidate : incumbent;
}
Problem read_problem(){
Problem problem;
cin >> problem.jobs >> problem.machines;
problem.operations=problem.jobs * problem.machines;
problem.machine.resize(problem.operations);
problem.duration.resize(problem.operations);
problem.node_on_machine.assign(
problem.jobs,vector<int>(problem.machines,-1)
);
problem.job_load.assign(problem.jobs,0);
problem.machine_load.assign(problem.machines,0);
for(int job=0;job<problem.jobs;++job){
for(int operation=0;operation<problem.machines;++operation){
int machine;
i64 duration;
cin >> machine >> duration;
int node_id=problem.node(job,operation);
problem.machine[node_id]=machine;
problem.duration[node_id]=duration;
problem.node_on_machine[job][machine]=node_id;
problem.job_load[job] += duration;
problem.machine_load[machine] += duration;
}
}
problem.remaining_work.assign(problem.operations,0);
for(int job=0;job<problem.jobs;++job){
i64 remaining=0;
for(int operation=problem.machines - 1;operation>=0;--operation){
int node_id=problem.node(job,operation);
remaining += problem.duration[node_id];
problem.remaining_work[node_id]=remaining;
}
}
for(i64 load : problem.job_load) problem.lower_bound=max(problem.lower_bound,load);
for(i64 load : problem.machine_load) problem.lower_bound=max(problem.lower_bound,load);
return problem;
}
void print_schedule(const Problem& problem,const Schedule& schedule){
for(int machine=0;machine<problem.machines;++machine){
for(int index=0;index<problem.jobs;++index){
if(index) cout << ' ';
cout << problem.job(schedule.order[machine][index]);
}
cout << '\n';
}
}
}
int main(){
ios::sync_with_stdio(false);
cin.tie(nullptr);
const Clock::time_point started=Clock::now();
Problem problem=read_problem();
vector<int> consensus(problem.machines),stage_sum(problem.machines,0);
iota(consensus.begin(),consensus.end(),0);
for(int job=0;job<problem.jobs;++job){
for(int machine=0;machine<problem.machines;++machine){
stage_sum[machine] += problem.operation(problem.node_on_machine[job][machine]);
}
}
sort(consensus.begin(),consensus.end(),[&](int left,int right){
return tie(stage_sum[left],left)<tie(stage_sum[right],right);
});
i64 route_agree=0,route_pairs=0;
int stage_agree=0;
for(int job=0;job<problem.jobs;++job){
for(int first=0;first<problem.machines;++first){
stage_agree += problem.machine[problem.node(job,first)]==consensus[first];
for(int second=first + 1;second<problem.machines;++second){
route_agree += problem.operation(problem.node_on_machine[job][consensus[first]]) <
problem.operation(problem.node_on_machine[job][consensus[second]]);
++route_pairs;
}
}
}
bool flow_like=route_pairs==0||5 * route_agree>=4 * route_pairs;
bool same_route=stage_agree==problem.operations;
bool strong_flow=flow_like&&(same_route||2 * stage_agree>=problem.operations);
i64 total_work=accumulate(problem.machine_load.begin(),problem.machine_load.end(),i64{0});
i64 max_machine_load=*max_element(problem.machine_load.begin(),problem.machine_load.end());
bool skewed=max_machine_load * problem.machines>=2 * total_work;
Evaluator evaluator(problem);
RNG rng;
vector<Schedule> seeds;
vector<uint64_t> hashes;
add_seed(problem,evaluator,seeds,hashes,serial_dispatch(problem,false));
add_seed(problem,evaluator,seeds,hashes,serial_dispatch(problem,true));
for(int rule=0;rule<6;++rule){
add_seed(
problem,evaluator,seeds,hashes,giffler_thompson(problem,rule,rng.next())
);
}
int dispatches=skewed ? 96 : flow_like ? 192 : 320;
for(int attempt=0;attempt<dispatches;++attempt){
add_seed(
problem,evaluator,seeds,hashes,giffler_thompson(problem,6 + attempt % 4,rng.next())
);
}
if(flow_like){
Problem backwards=reversed_problem(problem);
add_seed(
problem,evaluator,seeds,hashes,mirror_reversed_schedule(
problem,serial_dispatch(backwards,false)
)
);
add_seed(
problem,evaluator,seeds,hashes,mirror_reversed_schedule(
problem,serial_dispatch(backwards,true)
)
);
for(int rule=0;rule<6;++rule){
add_seed(
problem,evaluator,seeds,hashes,mirror_reversed_schedule(
problem,giffler_thompson(backwards,rule,rng.next())
)
);
}
for(int attempt=0;attempt<64;++attempt){
add_seed(
problem,evaluator,seeds,hashes,mirror_reversed_schedule(
problem,giffler_thompson(backwards,6 + attempt % 4,rng.next())
)
);
}
}
vector<int> machine_rank(problem.machines);
iota(machine_rank.begin(),machine_rank.end(),0);
sort(machine_rank.begin(),machine_rank.end(),[&](int left,int right){
return problem.machine_load[left]>problem.machine_load[right];
});
for(int index=0;index<min(2,problem.machines);++index){
int machine=machine_rank[index];
vector<int> order=schrage_order(problem,machine);
add_seed(
problem,evaluator,seeds,hashes,complete_fixed_machine(problem,machine,order)
);
}
vector<int> insertion_jobs(problem.jobs);
iota(insertion_jobs.begin(),insertion_jobs.end(),0);
sort(insertion_jobs.begin(),insertion_jobs.end(),[&](int left,int right){
return problem.job_load[left]>problem.job_load[right];
});
add_seed(
problem,evaluator,seeds,hashes,greedy_insertion_seed(
problem,evaluator,insertion_jobs.front(),1,started + chrono::milliseconds(180)
)
);
if(flow_like){
vector<int> descending(problem.jobs);
iota(descending.begin(),descending.end(),0);
sort(descending.begin(),descending.end(),[&](int left,int right){
return problem.job_load[left]>problem.job_load[right];
});
vector<int> neh=neh_permutation(
problem,descending,started + chrono::milliseconds(300),!same_route
);
if(same_route){
add_seed(problem,evaluator,seeds,hashes,common_permutation(problem,neh));
neh=improve_perm(problem,neh,started + chrono::milliseconds(400),false);
add_seed(problem,evaluator,seeds,hashes,common_permutation(problem,neh));
}else if(strong_flow){
neh=improve_perm(problem,neh,started + chrono::milliseconds(380),true);
}
vector<vector<int>> priority_starts=flow_priority_starts(
problem,neh,descending
);
RNG priority_rng;
priority_rng.state=0x369dea0f31a53f85ULL;
add_seed(
problem,evaluator,seeds,hashes,priority_permutation_search(
problem,priority_starts,priority_rng,started + chrono::milliseconds(270)
)
);
}
auto seed_best=min_element(seeds.begin(),seeds.end(),[](const Schedule& left,const Schedule& right){
return left.cost<right.cost;
});
i64 seed_gap=seed_best->cost - problem.lower_bound;
if(seed_gap<=seed_best->cost / (skewed ? 200 : 500)){
print_schedule(problem,*seed_best);
return 0;
}
int budget_ms=flow_like ? SOLVER_TIME_MS : skewed ? min(SOLVER_TIME_MS,400)
: min(SOLVER_TIME_MS,650);
const Clock::time_point deadline=started + chrono::milliseconds(budget_ms);
const Clock::time_point portfolio_start =
started + chrono::milliseconds(problem.operations<=800 ? 650 : 600);
RNG search_rng;
search_rng.state=strong_flow&&problem.operations>=1000
?0x123456789abcdef0ULL:0xdeadbeefcafebabeULL;
Schedule best;
if(flow_like&&!strong_flow){
best=annealing_search(
problem,evaluator,std::move(seeds),search_rng,started + chrono::milliseconds(
problem.operations>=1000
? 450
: problem.operations>=500 ? 350 : 580
)
);
if(best.cost!=problem.lower_bound&&Clock::now()<deadline){
vector<Schedule> final_seeds{best};
best=n7_tabu_search(
problem,evaluator,std::move(final_seeds),search_rng,deadline
);
}
}else{
bool long_trajectory=!flow_like &&
4 * problem.machines>3 * problem.jobs;
bool incremental_flow=strong_flow &&
problem.operations>=160&&problem.operations<=800;
if(incremental_flow){
best=incremental_n5_tabu(
problem,evaluator,std::move(seeds),search_rng,min(deadline,portfolio_start),long_trajectory
);
}else{
best=estimated_n5_tabu(
problem,evaluator,std::move(seeds),search_rng,strong_flow
? min(deadline,portfolio_start)
: deadline,long_trajectory
);
}
if(strong_flow&&best.cost!=problem.lower_bound &&
Clock::now()<deadline){
if(incremental_flow){
Problem backwards=reversed_problem(problem);
Schedule reverse_seed=mirror_reversed_schedule(backwards,best);
Evaluator reverse_evaluator(backwards);
reverse_evaluator.evaluate(reverse_seed,false);
vector<Schedule> reverse_seeds{reverse_seed};
RNG reverse_rng;
reverse_rng.state=schedule_hash(reverse_seed) ^ 0x2718281828459045ULL;
Schedule reverse_best=n7_tabu_search(backwards,reverse_evaluator,std::move(reverse_seeds),reverse_rng,deadline);
Schedule candidate=mirror_reversed_schedule(problem,reverse_best);
if(evaluator.evaluate(candidate,false)&&candidate.cost<best.cost) best=std::move(candidate);
}else{
best=single_continuation(problem,best,deadline);
}
}
}
print_schedule(problem,best);
return 0;
}
