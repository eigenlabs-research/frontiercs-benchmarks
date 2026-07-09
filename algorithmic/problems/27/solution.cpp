#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <chrono>
using namespace std;
using namespace std::chrono;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int n, m;
    if (!(cin >> n >> m)) return 0;
    
    auto start_time = high_resolution_clock::now();
    mt19937 rng(1337);
    
    vector<int> all_cells(n * m);
    iota(all_cells.begin(), all_cells.end(), 0);
    
    int best_k = -1;
    vector<int> best_solution;
    
    int iterations = 0;
    while (true) {
        if ((iterations & 7) == 0) {
            auto now = high_resolution_clock::now();
            if (duration_cast<milliseconds>(now - start_time).count() > 1850) break;
        }
        iterations++;
        
        shuffle(all_cells.begin(), all_cells.end(), rng);
        
        vector<bool> grid(n * m, false);
        vector<vector<int>> rows(n);
        vector<vector<int>> cols(m);
        
        vector<int> current_solution;
        int k = 0;
        
        for (int cell : all_cells) {
            int r = cell / m;
            int c = cell % m;
            
            bool forms_rect = false;
            // To form a rectangle with (r, c), there must exist r' != r and c' != c
            // such that (r', c), (r, c'), and (r', c') are all placed.
            // We iterate over all c' placed in row r.
            // For each c', we check if there's an r' placed in col c, such that (r', c') is also placed.
            // If the number of elements is large, we should iterate over the smaller set.
            
            if (rows[r].size() < cols[c].size()) {
                for (int c_prime : rows[r]) {
                    for (int r_prime : cols[c]) {
                        if (grid[r_prime * m + c_prime]) {
                            forms_rect = true;
                            break;
                        }
                    }
                    if (forms_rect) break;
                }
            } else {
                for (int r_prime : cols[c]) {
                    for (int c_prime : rows[r]) {
                        if (grid[r_prime * m + c_prime]) {
                            forms_rect = true;
                            break;
                        }
                    }
                    if (forms_rect) break;
                }
            }
            
            if (!forms_rect) {
                grid[cell] = true;
                rows[r].push_back(c);
                cols[c].push_back(r);
                current_solution.push_back(cell);
                k++;
            }
        }
        
        if (k > best_k) {
            best_k = k;
            best_solution = current_solution;
        }
    }
    
    cout << best_k << "
";
    for (int cell : best_solution) {
        cout << (cell / m) + 1 << " " << (cell % m) + 1 << "
";
    }
    return 0;
}
