#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

struct Cell { int x, y; };

struct Poly {
    int id;
    int k;
    vector<Cell> cells;
};

// 0: none, 1: rot90, 2: rot180, 3: rot270
// F: 0 no ref, 1 ref across y axis (x = -x)
Cell transform_cell(Cell c, int R, int F) {
    if (F == 1) c.x = -c.x;
    if (R == 1) return {c.y, -c.x};
    if (R == 2) return {-c.x, -c.y};
    if (R == 3) return {-c.y, c.x};
    return c;
}

const int MAX_GRID = 3000;
bool grid2d[MAX_GRID][MAX_GRID];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int n;
    if (!(cin >> n)) return 0;
    
    vector<Poly> polys(n);
    for(int i=0; i<n; ++i) {
        polys[i].id = i;
        cin >> polys[i].k;
        polys[i].cells.resize(polys[i].k);
        for(int j=0; j<polys[i].k; ++j) {
            cin >> polys[i].cells[j].x >> polys[i].cells[j].y;
        }
    }
    
    // Sort by area (number of cells) descending
    vector<Poly> sorted_polys = polys;
    sort(sorted_polys.begin(), sorted_polys.end(), [](const Poly& a, const Poly& b){
        return a.k > b.k;
    });
    
    struct Result {
        int X, Y, R, F;
    };
    vector<Result> results(n);
    
    int max_x_used = 0;
    int max_y_used = 0;
    
    for(int i=0; i<n; ++i) {
        const Poly& p = sorted_polys[i];
        
        bool placed = false;
        int best_X = -1, best_Y = -1, best_R = -1, best_F = -1;
        vector<Cell> best_placed_cells;
        
        // Try to place as bottom-left as possible
        // We will scan Y from 0 up, then X from 0 up
        for(int startY = 0; startY < MAX_GRID && !placed; ++startY) {
            for(int startX = 0; startX < MAX_GRID && !placed; ++startX) {
                // Try all 8 orientations
                for(int F = 0; F < 2 && !placed; ++F) {
                    for(int R = 0; R < 4 && !placed; ++R) {
                        
                        // find bounding box of transformed cells
                        int min_tx = 1e9, min_ty = 1e9;
                        for(auto c : p.cells) {
                            Cell tc = transform_cell(c, R, F);
                            min_tx = min(min_tx, tc.x);
                            min_ty = min(min_ty, tc.y);
                        }
                        
                        // the translation needed to put the minimums at (startX, startY)
                        int X_trans = startX - min_tx;
                        int Y_trans = startY - min_ty;
                        
                        bool ok = true;
                        vector<Cell> test_cells;
                        for(auto c : p.cells) {
                            Cell tc = transform_cell(c, R, F);
                            int final_x = tc.x + X_trans;
                            int final_y = tc.y + Y_trans;
                            if (final_x >= MAX_GRID || final_y >= MAX_GRID || final_x < 0 || final_y < 0) {
                                ok = false; break;
                            }
                            if (grid2d[final_x][final_y]) {
                                ok = false; break;
                            }
                            test_cells.push_back({final_x, final_y});
                        }
                        
                        if (ok) {
                            placed = true;
                            best_X = X_trans;
                            best_Y = Y_trans;
                            best_R = R;
                            best_F = F;
                            best_placed_cells = test_cells;
                        }
                    }
                }
            }
        }
        
        for(auto c : best_placed_cells) {
            grid2d[c.x][c.y] = true;
            max_x_used = max(max_x_used, c.x);
            max_y_used = max(max_y_used, c.y);
        }
        
        results[p.id] = {best_X, best_Y, best_R, best_F};
    }
    
    cout << max_x_used + 1 << " " << max_y_used + 1 << "
";
    for(int i=0; i<n; ++i) {
        cout << results[i].X << " " << results[i].Y << " " << results[i].R << " " << results[i].F << "
";
    }
    
    return 0;
}
