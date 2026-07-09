#include <bits/stdc++.h>
using namespace std;

struct Point {
    int x, y;
    bool operator<(const Point& p) const { return tie(x, y) < tie(p.x, p.y); }
    bool operator==(const Point& p) const { return x == p.x && y == p.y; }
};

struct Polyomino {
    int k;
    vector<Point> cells;
    vector<vector<Point>> transforms;
    
    void generate_transforms() {
        transforms.clear();
        vector<Point> original = cells;
        
        for (int flip = 0; flip <= 1; flip++) {
            for (int rot = 0; rot < 4; rot++) {
                vector<Point> current = original;
                
                if (flip) {
                    for (auto& p : current) p.x = -p.x;
                }
                
                for (int r = 0; r < rot; r++) {
                    for (auto& p : current) {
                        int nx = -p.y;
                        int ny = p.x;
                        p.x = nx;
                        p.y = ny;
                    }
                }
                
                // Normalize
                int min_x = INT_MAX, min_y = INT_MAX;
                for (auto& p : current) {
                    min_x = min(min_x, p.x);
                    min_y = min(min_y, p.y);
                }
                for (auto& p : current) {
                    p.x -= min_x;
                    p.y -= min_y;
                }
                
                sort(current.begin(), current.end());
                current.erase(unique(current.begin(), current.end()), current.end());
                
                bool duplicate = false;
                for (const auto& t : transforms) {
                    if (t == current) {
                        duplicate = true;
                        break;
                    }
                }
                
                if (!duplicate) {
                    transforms.push_back(current);
                }
            }
        }
    }
};

vector<Polyomino> polys;
int total_cells = 0;

struct Placement {
    int x, y, rot, flip;
    int transform_idx;
    vector<Point> cells;
};

vector<Placement> best_placements;
int best_W = 0, best_H = 0;

bool is_valid(const vector<Placement>& placements, int W, int H) {
    set<Point> occupied;
    for (const auto& pl : placements) {
        for (const auto& p : pl.cells) {
            int x = p.x + pl.x;
            int y = p.y + pl.y;
            if (x < 0 || x >= W || y < 0 || y >= H) return false;
            if (occupied.count({x, y})) return false;
            occupied.insert({x, y});
        }
    }
    return true;
}

void solve() {
    // Initialize with a simple packing
    int W = ceil(sqrt(total_cells));
    int H = ceil(total_cells * 1.0 / W);
    
    best_W = W;
    best_H = H;
    
    vector<Placement> placements;
    set<Point> occupied;
    
    // Sort polyominoes by size descending
    vector<int> indices(polys.size());
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(), [&](int a, int b) {
        return polys[a].k > polys[b].k;
    });
    
    for (int idx : indices) {
        const auto& poly = polys[idx];
        bool placed = false;
        
        // Try all possible transforms in random order
        vector<int> transform_order(poly.transforms.size());
        iota(transform_order.begin(), transform_order.end(), 0);
        random_shuffle(transform_order.begin(), transform_order.end());
        
        for (int ti : transform_order) {
            const auto& cells = poly.transforms[ti];
            
            // Find flip and rot from transform index
            int flip = ti / 4;
            int rot = ti % 4;
            
            // Try to place it in the first available position
            for (int y = 0; y < H && !placed; y++) {
                for (int x = 0; x < W && !placed; x++) {
                    bool can_place = true;
                    for (const auto& p : cells) {
                        int nx = x + p.x;
                        int ny = y + p.y;
                        if (nx >= W || ny >= H || occupied.count({nx, ny})) {
                            can_place = false;
                            break;
                        }
                    }
                    
                    if (can_place) {
                        Placement pl;
                        pl.x = x;
                        pl.y = y;
                        pl.rot = rot;
                        pl.flip = flip;
                        pl.transform_idx = ti;
                        pl.cells = cells;
                        placements.push_back(pl);
                        
                        for (const auto& p : cells) {
                            occupied.insert({x + p.x, y + p.y});
                        }
                        placed = true;
                    }
                }
            }
        }
        
        if (!placed) {
            // Need to expand the grid
            if (W <= H) {
                W++;
            } else {
                H++;
            }
            
            // Reset and try again
            placements.clear();
            occupied.clear();
            idx = 0; // Restart from beginning
        }
    }
    
    best_placements = placements;
    best_W = W;
    best_H = H;
    
    // Try to optimize the packing
    for (int iter = 0; iter < 100; iter++) {
        int new_W = best_W;
        int new_H = best_H;
        
        if (iter % 2 == 0 && new_W * new_H > total_cells) {
            if (new_H > 1 && (new_W) * (new_H - 1) >= total_cells) {
                new_H--;
            } else if (new_W > 1 && (new_W - 1) * (new_H) >= total_cells) {
                new_W--;
            }
        }
        
        if (new_W * new_H >= best_W * best_H) continue;
        
        vector<Placement> new_placements;
        set<Point> new_occupied;
        bool valid = true;
        
        for (const auto& pl : best_placements) {
            bool placed = false;
            for (const auto& p : pl.cells) {
                int x = pl.x + p.x;
                int y = pl.y + p.y;
                if (x >= new_W || y >= new_H) {
                    valid = false;
                    break;
                }
            }
            if (!valid) break;
        }
        
        if (valid && is_valid(best_placements, new_W, new_H)) {
            best_W = new_W;
            best_H = new_H;
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int n;
    cin >> n;
    
    polys.resize(n);
    for (int i = 0; i < n; i++) {
        int k;
        cin >> k;
        polys[i].k = k;
        total_cells += k;
        
        for (int j = 0; j < k; j++) {
            int x, y;
            cin >> x >> y;
            polys[i].cells.push_back({x, y});
        }
        
        polys[i].generate_transforms();
    }
    
    solve();
    
    cout << best_W << " " << best_H << "\n";
    for (const auto& pl : best_placements) {
        cout << pl.x << " " << pl.y << " " << pl.rot << " " << pl.flip << "\n";
    }
    
    return 0;
}