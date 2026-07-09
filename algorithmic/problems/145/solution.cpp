#include <bits/stdc++.h>
using namespace std;

struct Point {
    int x, y;
    bool operator<(const Point& p) const {
        return tie(x, y) < tie(p.x, p.y);
    }
    bool operator==(const Point& p) const {
        return x == p.x && y == p.y;
    }
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
                for (auto p : current) {
                    min_x = min(min_x, p.x);
                    min_y = min(min_y, p.y);
                }
                for (auto& p : current) {
                    p.x -= min_x;
                    p.y -= min_y;
                }
                sort(current.begin(), current.end());
                
                // Check if unique
                bool unique = true;
                for (const auto& t : transforms) {
                    if (t == current) {
                        unique = false;
                        break;
                    }
                }
                if (unique) {
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

vector<Placement> placements;
set<Point> occupied;
int W, H;

void solve() {
    // Initialize with a square-ish rectangle
    int area = total_cells;
    W = sqrt(area);
    while (area % W != 0) W++;
    H = area / W;
    
    // Sort polyominoes by size (largest first)
    vector<int> order(polys.size());
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
        return polys[a].k > polys[b].k;
    });
    
    placements.resize(polys.size());
    
    // Try placing each polyomino
    for (int idx : order) {
        auto& poly = polys[idx];
        bool placed = false;
        
        // Try all possible transforms in random order
        vector<int> torder(poly.transforms.size());
        iota(torder.begin(), torder.end(), 0);
        random_shuffle(torder.begin(), torder.end());
        
        for (int ti : torder) {
            auto& t = poly.transforms[ti];
            
            // Determine flip and rotation from transform index
            int flip = (ti / 4) % 2;
            int rot = ti % 4;
            
            // Find possible positions
            for (int y = 0; y <= H; y++) {
                for (int x = 0; x <= W; x++) {
                    bool can_place = true;
                    vector<Point> new_cells;
                    
                    for (auto p : t) {
                        Point np = {x + p.x, y + p.y};
                        if (np.x >= W || np.y >= H || occupied.count(np)) {
                            can_place = false;
                            break;
                        }
                        new_cells.push_back(np);
                    }
                    
                    if (can_place) {
                        placements[idx] = {x, y, rot, flip, ti, new_cells};
                        for (auto p : new_cells) {
                            occupied.insert(p);
                        }
                        placed = true;
                        goto next_poly;
                    }
                }
            }
        }
        
        // If we get here, couldn't place - need to expand
        if (!placed) {
            // Simple expansion strategy: increase the smaller dimension
            if (W <= H) {
                W++;
            } else {
                H++;
            }
            // Reset and try again
            occupied.clear();
            placements.clear();
            placements.resize(polys.size());
            idx = order[0]; // Start from beginning
        }
        next_poly:;
    }
    
    // Try to optimize the rectangle
    while (true) {
        bool improved = false;
        
        // Try reducing height
        if (H > 1) {
            set<Point> test_occ;
            bool possible = true;
            for (auto& pl : placements) {
                for (auto p : pl.cells) {
                    if (p.y >= H-1) {
                        possible = false;
                        break;
                    }
                }
                if (!possible) break;
            }
            if (possible) {
                H--;
                improved = true;
            }
        }
        
        // Try reducing width
        if (W > 1) {
            set<Point> test_occ;
            bool possible = true;
            for (auto& pl : placements) {
                for (auto p : pl.cells) {
                    if (p.x >= W-1) {
                        possible = false;
                        break;
                    }
                }
                if (!possible) break;
            }
            if (possible) {
                W--;
                improved = true;
            }
        }
        
        if (!improved) break;
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
        
        // Generate all unique transformations
        polys[i].generate_transforms();
    }
    
    solve();
    
    // Output the solution
    cout << W << " " << H << "\n";
    for (int i = 0; i < n; i++) {
        auto& pl = placements[i];
        cout << pl.x << " " << pl.y << " " << pl.rot << " " << pl.flip << "\n";
    }
    
    return 0;
}