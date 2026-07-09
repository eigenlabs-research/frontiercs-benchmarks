#include <bits/stdc++.h>
using namespace std;

struct Point {
    int x, y;
    bool operator<(const Point& other) const {
        return x < other.x || (x == other.x && y < other.y);
    }
    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }
};

struct Polyomino {
    int k;
    vector<Point> cells;
    vector<vector<Point>> transformations;

    void generate_transformations() {
        transformations.clear();
        vector<Point> current = cells;

        for (int flip = 0; flip < 2; ++flip) {
            for (int rot = 0; rot < 4; ++rot) {
                vector<Point> transformed = current;
                
                // Normalize by translating to minimum (0,0)
                int min_x = INT_MAX, min_y = INT_MAX;
                for (auto& p : transformed) {
                    min_x = min(min_x, p.x);
                    min_y = min(min_y, p.y);
                }
                for (auto& p : transformed) {
                    p.x -= min_x;
                    p.y -= min_y;
                }
                sort(transformed.begin(), transformed.end());
                
                // Check if this transformation is unique
                bool unique = true;
                for (const auto& t : transformations) {
                    if (t == transformed) {
                        unique = false;
                        break;
                    }
                }
                if (unique) {
                    transformations.push_back(transformed);
                }

                // Rotate for next iteration
                for (auto& p : current) {
                    int new_x = -p.y;
                    int new_y = p.x;
                    p.x = new_x;
                    p.y = new_y;
                }
            }
            // Flip for next iteration
            for (auto& p : current) {
                p.x = -p.x;
            }
        }
    }
};

vector<Polyomino> polyominoes;
int total_cells = 0;

struct Placement {
    int x, y, rot, flip;
    int variant_idx;
};

vector<Placement> best_placements;
int best_W = 0, best_H = 0;
int best_area = INT_MAX;

void solve() {
    // Sort polyominoes by size (largest first)
    sort(polyominoes.begin(), polyominoes.end(), [](const Polyomino& a, const Polyomino& b) {
        return a.k > b.k;
    });

    // Initialize grid
    int W = ceil(sqrt(total_cells));
    int H = ceil(sqrt(total_cells));
    vector<vector<bool>> grid(W, vector<bool>(H, false));
    vector<Placement> placements;

    for (const auto& poly : polyominoes) {
        bool placed = false;
        
        // Try all possible positions and transformations
        for (int y = 0; y < H && !placed; ++y) {
            for (int x = 0; x < W && !placed; ++x) {
                for (size_t t = 0; t < poly.transformations.size() && !placed; ++t) {
                    const auto& variant = poly.transformations[t];
                    
                    // Check if this placement is valid
                    bool valid = true;
                    for (const auto& p : variant) {
                        int nx = x + p.x;
                        int ny = y + p.y;
                        if (nx >= W || ny >= H) {
                            valid = false;
                            break;
                        }
                        if (grid[nx][ny]) {
                            valid = false;
                            break;
                        }
                    }
                    
                    if (valid) {
                        // Place the polyomino
                        for (const auto& p : variant) {
                            int nx = x + p.x;
                            int ny = y + p.y;
                            grid[nx][ny] = true;
                        }
                        
                        // Determine flip and rotation from variant index
                        int flip = (t / 4) % 2;
                        int rot = t % 4;
                        
                        placements.push_back({x, y, rot, flip, (int)t});
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
            grid.resize(W, vector<bool>(H, false));
            // Retry placement
            x--;
            continue;
        }
    }
    
    // Update best solution
    int area = W * H;
    if (area < best_area || 
        (area == best_area && H < best_H) || 
        (area == best_area && H == best_H && W < best_W)) {
        best_area = area;
        best_W = W;
        best_H = H;
        best_placements = placements;
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int n;
    cin >> n;
    polyominoes.resize(n);
    
    for (int i = 0; i < n; ++i) {
        int k;
        cin >> k;
        polyominoes[i].k = k;
        total_cells += k;
        
        for (int j = 0; j < k; ++j) {
            int x, y;
            cin >> x >> y;
            polyominoes[i].cells.push_back({x, y});
        }
        
        // Generate all unique transformations
        polyominoes[i].generate_transformations();
    }
    
    // Try multiple times with different random seeds
    const int ATTEMPTS = 5;
    for (int attempt = 0; attempt < ATTEMPTS; ++attempt) {
        // Shuffle the order of polyominoes
        random_shuffle(polyominoes.begin(), polyominoes.end());
        solve();
    }
    
    // Output the best solution found
    cout << best_W << " " << best_H << "\n";
    
    // Reconstruct the placements in original order
    vector<Placement> original_order_placements(n);
    for (size_t i = 0; i < polyominoes.size(); ++i) {
        original_order_placements[i] = best_placements[i];
    }
    
    for (const auto& p : original_order_placements) {
        cout << p.x << " " << p.y << " " << p.rot << " " << p.flip << "\n";
    }
    
    return 0;
}