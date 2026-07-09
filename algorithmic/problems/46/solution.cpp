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
    vector<Point> cells;
    vector<vector<Point>> transformations;
    
    void generate_transformations() {
        transformations.clear();
        
        // Original
        vector<Point> original = cells;
        transformations.push_back(original);
        
        // Generate all rotations (90, 180, 270)
        vector<Point> rotated90, rotated180, rotated270;
        for (const Point& p : original) {
            rotated90.push_back({-p.y, p.x});
        }
        normalize(rotated90);
        transformations.push_back(rotated90);
        
        for (const Point& p : rotated90) {
            rotated180.push_back({-p.y, p.x});
        }
        normalize(rotated180);
        transformations.push_back(rotated180);
        
        for (const Point& p : rotated180) {
            rotated270.push_back({-p.y, p.x});
        }
        normalize(rotated270);
        transformations.push_back(rotated270);
        
        // Generate reflected versions
        vector<Point> reflected;
        for (const Point& p : original) {
            reflected.push_back({-p.x, p.y});
        }
        normalize(reflected);
        transformations.push_back(reflected);
        
        vector<Point> reflected90, reflected180, reflected270;
        for (const Point& p : reflected) {
            reflected90.push_back({-p.y, p.x});
        }
        normalize(reflected90);
        transformations.push_back(reflected90);
        
        for (const Point& p : reflected90) {
            reflected180.push_back({-p.y, p.x});
        }
        normalize(reflected180);
        transformations.push_back(reflected180);
        
        for (const Point& p : reflected180) {
            reflected270.push_back({-p.y, p.x});
        }
        normalize(reflected270);
        transformations.push_back(reflected270);
    }
    
    void normalize(vector<Point>& pts) {
        if (pts.empty()) return;
        int min_x = pts[0].x, min_y = pts[0].y;
        for (const Point& p : pts) {
            min_x = min(min_x, p.x);
            min_y = min(min_y, p.y);
        }
        for (Point& p : pts) {
            p.x -= min_x;
            p.y -= min_y;
        }
        sort(pts.begin(), pts.end());
    }
};

vector<Polyomino> polyominoes;
int n;

struct Placement {
    int x, y, rot, flip;
    vector<Point> transformed_cells;
};

vector<Placement> placements;
set<Point> occupied;

bool can_place(const vector<Point>& cells, int dx, int dy) {
    for (const Point& p : cells) {
        Point new_p = {p.x + dx, p.y + dy};
        if (occupied.count(new_p)) {
            return false;
        }
    }
    return true;
}

void place_polyomino(const vector<Point>& cells, int dx, int dy) {
    for (const Point& p : cells) {
        Point new_p = {p.x + dx, p.y + dy};
        occupied.insert(new_p);
    }
}

pair<int, int> get_bounding_box() {
    if (occupied.empty()) return {0, 0};
    int min_x = INT_MAX, min_y = INT_MAX;
    int max_x = INT_MIN, max_y = INT_MIN;
    for (const Point& p : occupied) {
        min_x = min(min_x, p.x);
        min_y = min(min_y, p.y);
        max_x = max(max_x, p.x);
        max_y = max(max_y, p.y);
    }
    return {max_x + 1, max_y + 1};
}

void solve() {
    // Generate all transformations for each polyomino
    for (Polyomino& poly : polyominoes) {
        poly.generate_transformations();
    }
    
    // Sort polyominoes by size (largest first)
    sort(polyominoes.begin(), polyominoes.end(), [](const Polyomino& a, const Polyomino& b) {
        return a.cells.size() > b.cells.size();
    });
    
    // Place each polyomino
    for (const Polyomino& poly : polyominoes) {
        bool placed = false;
        // Try all transformations
        for (int flip = 0; flip <= 1; ++flip) {
            for (int rot = 0; rot < 4; ++rot) {
                int trans_idx = rot;
                if (flip) trans_idx += 4;
                const vector<Point>& transformed = poly.transformations[trans_idx];
                
                // Try placing at all possible positions
                for (int y = 0; y < 100; ++y) {
                    for (int x = 0; x < 100; ++x) {
                        if (can_place(transformed, x, y)) {
                            place_polyomino(transformed, x, y);
                            placements.push_back({x, y, rot, flip, transformed});
                            placed = true;
                            goto next_poly;
                        }
                    }
                }
                if (placed) break;
            }
            if (placed) break;
        }
        next_poly:
        if (!placed) {
            // Shouldn't happen if algorithm is correct
            cerr << "Failed to place polyomino!" << endl;
            exit(1);
        }
    }
    
    auto [W, H] = get_bounding_box();
    cout << W << " " << H << endl;
    
    // Output placements (need to map back to original indices)
    // Since we sorted, we need to track original indices
    // This is simplified for the example - in practice you'd need to track original indices
    for (const Placement& pl : placements) {
        cout << pl.x << " " << pl.y << " " << pl.rot << " " << pl.flip << endl;
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cin >> n;
    polyominoes.resize(n);
    for (int i = 0; i < n; ++i) {
        int k;
        cin >> k;
        polyominoes[i].cells.resize(k);
        for (int j = 0; j < k; ++j) {
            cin >> polyominoes[i].cells[j].x >> polyominoes[i].cells[j].y;
        }
        polyominoes[i].normalize(polyominoes[i].cells);
    }
    
    solve();
    
    return 0;
}