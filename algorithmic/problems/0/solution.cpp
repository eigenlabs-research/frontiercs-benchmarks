
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <cmath>
using namespace std;
using namespace std::chrono;

struct Cell { int x, y; };

struct Poly {
    int id;
    int k;
    vector<Cell> cells;
};

// Transform Cell based on Rotation (0..3) and Flip (0..1)
Cell transform_cell(Cell c, int R, int F) {
    if (F == 1) c.x = -c.x;
    if (R == 1) return {c.y, -c.x};
    if (R == 2) return {-c.x, -c.y};
    if (R == 3) return {-c.y, c.x};
    return c;
}

const int MAX_GRID = 2000;
bool grid2d[MAX_GRID][MAX_GRID];

struct Result {
    int X, Y, R, F;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int n;
    if (!(cin >> n)) return 0;
    
    vector<Poly> polys(n);
    int total_cells = 0;
    for(int i=0; i<n; ++i) {
        polys[i].id = i;
        cin >> polys[i].k;
        total_cells += polys[i].k;
        polys[i].cells.resize(polys[i].k);
        for(int j=0; j<polys[i].k; ++j) {
            cin >> polys[i].cells[j].x >> polys[i].cells[j].y;
        }
    }
    
    // Base State
    vector<int> order(n);
    vector<int> flips(n);
    vector<int> rots(n);
    for(int i=0; i<n; ++i) {
        order[i] = i;
        flips[i] = 0;
        rots[i] = 0;
    }
    
    // Initial sort by area descending
    sort(order.begin(), order.end(), [&](int a, int b) {
        return polys[a].k > polys[b].k;
    });

    auto start_time = high_resolution_clock::now();
    mt19937 rng(1337);
    
    auto evaluate = [&](const vector<int>& curr_order, const vector<int>& curr_flips, const vector<int>& curr_rots, vector<Result>& out_res) {
        for(int x=0; x<MAX_GRID; ++x) {
            for(int y=0; y<MAX_GRID; ++y) {
                grid2d[x][y] = false;
            }
        }
        int max_x = 0;
        int max_y = 0;
        out_res.assign(n, {-1, -1, 0, 0});
        
        for(int idx : curr_order) {
            const Poly& p = polys[idx];
            int R = curr_rots[idx];
            int F = curr_flips[idx];
            
            // find bounding box
            int min_tx = 1e9, min_ty = 1e9;
            vector<Cell> tcells(p.k);
            for(int j=0; j<p.k; ++j) {
                Cell tc = transform_cell(p.cells[j], R, F);
                tcells[j] = tc;
                min_tx = min(min_tx, tc.x);
                min_ty = min(min_ty, tc.y);
            }
            
            bool placed = false;
            // Bottom-Left heuristic: prioritize smaller Y + X
            // We search in a diamond shape or just standard nested loops
            for(int sum = 0; sum < MAX_GRID * 2 && !placed; ++sum) {
                for(int startY = 0; startY <= sum && !placed; ++startY) {
                    int startX = sum - startY;
                    if(startX >= MAX_GRID || startY >= MAX_GRID) continue;
                    
                    int X_trans = startX - min_tx;
                    int Y_trans = startY - min_ty;
                    
                    bool ok = true;
                    for(int j=0; j<p.k; ++j) {
                        int final_x = tcells[j].x + X_trans;
                        int final_y = tcells[j].y + Y_trans;
                        if (final_x >= MAX_GRID || final_y >= MAX_GRID || final_x < 0 || final_y < 0) {
                            ok = false; break;
                        }
                        if (grid2d[final_x][final_y]) {
                            ok = false; break;
                        }
                    }
                    if(ok) {
                        placed = true;
                        out_res[idx] = {X_trans, Y_trans, R, F};
                        for(int j=0; j<p.k; ++j) {
                            int final_x = tcells[j].x + X_trans;
                            int final_y = tcells[j].y + Y_trans;
                            grid2d[final_x][final_y] = true;
                            max_x = max(max_x, final_x);
                            max_y = max(max_y, final_y);
                        }
                    }
                }
            }
        }
        return (max_x + 1) * (max_y + 1);
    };
    
    vector<Result> best_results;
    int best_area = evaluate(order, flips, rots, best_results);
    vector<int> best_order = order;
    vector<int> best_flips = flips;
    vector<int> best_rots = rots;
    
    double T = 1000.0;
    double alpha = 0.99;
    
    vector<int> curr_order = order;
    vector<int> curr_flips = flips;
    vector<int> curr_rots = rots;
    int curr_area = best_area;
    
    int iterations = 0;
    while(true) {
        iterations++;
        if ((iterations & 7) == 0) {
            auto now = high_resolution_clock::now();
            double elapsed = duration_cast<milliseconds>(now - start_time).count() / 1000.0;
            if (elapsed > 1.85) break;
        }
        
        int type = rng() % 3;
        vector<int> next_order = curr_order;
        vector<int> next_flips = curr_flips;
        vector<int> next_rots = curr_rots;
        
        if (type == 0 && n > 1) { // swap two pieces
            int i = rng() % n;
            int j = rng() % n;
            swap(next_order[i], next_order[j]);
        } else if (type == 1) { // change rotation
            int i = rng() % n;
            next_rots[i] = rng() % 4;
        } else { // change flip
            int i = rng() % n;
            next_flips[i] = rng() % 2;
        }
        
        vector<Result> next_res;
        int next_area = evaluate(next_order, next_flips, next_rots, next_res);
        
        int delta = next_area - curr_area;
        if (delta < 0 || exp(-delta / T) > (double)rng() / rng.max()) {
            curr_area = next_area;
            curr_order = next_order;
            curr_flips = next_flips;
            curr_rots = next_rots;
            
            if (curr_area < best_area) {
                best_area = curr_area;
                best_results = next_res;
                best_order = curr_order;
                best_flips = curr_flips;
                best_rots = curr_rots;
            }
        }
        
        T *= alpha;
        if(T < 0.001) T = 0.001;
    }
    
    int max_x = 0, max_y = 0;
    for(int i=0; i<n; ++i) {
        const Poly& p = polys[i];
        int R = best_results[i].R;
        int F = best_results[i].F;
        int X = best_results[i].X;
        int Y = best_results[i].Y;
        for(int j=0; j<p.k; ++j) {
            Cell tc = transform_cell(p.cells[j], R, F);
            max_x = max(max_x, tc.x + X);
            max_y = max(max_y, tc.y + Y);
        }
    }
    
    cout << max_x + 1 << " " << max_y + 1 << "\n";
    for(int i=0; i<n; ++i) {
        cout << best_results[i].X << " " << best_results[i].Y << " " << best_results[i].R << " " << best_results[i].F << "\n";
    }
    
    return 0;
}

