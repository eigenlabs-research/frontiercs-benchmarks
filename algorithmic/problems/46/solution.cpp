#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
using namespace std;
using namespace std::chrono;

struct Op { int m, p; };

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int J, M;
    if (!(cin >> J >> M)) return 0;
    
    vector<vector<Op>> jobs(J, vector<Op>(M));
    for(int j=0; j<J; ++j) {
        for(int m=0; m<M; ++m) {
            cin >> jobs[j][m].m >> jobs[j][m].p;
        }
    }
    
    auto start_time = high_resolution_clock::now();
    mt19937 rng(1337);
    
    int best_makespan = 2e9;
    vector<vector<int>> best_machine_orders(M);
    
    int iterations = 0;
    while(true) {
        if ((iterations & 63) == 0) {
            auto now = high_resolution_clock::now();
            if (duration_cast<milliseconds>(now - start_time).count() > 1850) break;
        }
        iterations++;
        
        vector<int> next_op(J, 0);
        vector<int> job_time(J, 0);
        vector<int> mach_time(M, 0);
        
        vector<vector<int>> mach_orders(M);
        
        // Active schedule generation with SPT / random bias
        int completed = 0;
        while(completed < J * M) {
            // Find valid operations
            vector<int> available_jobs;
            for(int j=0; j<J; ++j) {
                if(next_op[j] < M) {
                    available_jobs.push_back(j);
                }
            }
            
            // Heuristic: pick the one that finishes earliest, with some randomness
            int best_j = -1;
            int best_finish = 2e9;
            
            // Randomly select a job to schedule
            // or use SPT:
            int num_cand = available_jobs.size();
            
            int cand_idx = rng() % num_cand;
            // with 80% chance, pick the absolute best finish time
            if (rng() % 100 < 80) {
                for(int i=0; i<num_cand; ++i) {
                    int j = available_jobs[i];
                    int m = jobs[j][next_op[j]].m;
                    int p = jobs[j][next_op[j]].p;
                    int start = max(job_time[j], mach_time[m]);
                    int finish = start + p;
                    if (finish < best_finish) {
                        best_finish = finish;
                        best_j = j;
                    }
                }
            } else {
                best_j = available_jobs[cand_idx];
            }
            
            int j = best_j;
            int m = jobs[j][next_op[j]].m;
            int p = jobs[j][next_op[j]].p;
            
            int start = max(job_time[j], mach_time[m]);
            int finish = start + p;
            
            job_time[j] = finish;
            mach_time[m] = finish;
            
            mach_orders[m].push_back(j);
            next_op[j]++;
            completed++;
        }
        
        int current_makespan = 0;
        for(int m=0; m<M; ++m) {
            current_makespan = max(current_makespan, mach_time[m]);
        }
        
        if (current_makespan < best_makespan) {
            best_makespan = current_makespan;
            best_machine_orders = mach_orders;
        }
    }
    
    for(int m=0; m<M; ++m) {
        for(int i=0; i<J; ++i) {
            cout << best_machine_orders[m][i] << (i==J-1 ? "" : " ");
        }
        cout << "
";
    }
    
    return 0;
}
