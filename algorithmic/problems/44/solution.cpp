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
        
        // Generate all 8 possible transformations (4 rotations × 2 reflections)
        for (int flip = 0; flip <= 1; ++flip) {
            for (int rot = 0; rot < 4; ++rot) {
                vector<Point> transformed;
                for (const auto& p : cells) {
                    Point np = p;
                    // Flip if needed
                    if (flip) np.x = -np.x;
                    // Rotate
                    for (int r = 0; r < rot; ++r) {
                        int nx = -np.y;
                        int ny = np.x;
                        np.x = nx;
                        np.y = ny;
                    }
                    transformed.push_back(np);
                }
                // Normalize by translating to minimum (0,0)
                int min_x = transformed[0].x, min_y = transformed[0].y;
                for (const auto& p : transformed) {
                    min_x = min(min_x, p.x);
                    min_y = min(min_y, p.y);
                }
                for (auto& p : transformed) {
                    p.x -= min_x;
                    p.y -= min_y;
                }
                sort(transformed.begin(), transformed.end());
                transformations.push_back(transformed);
            }
        }
    }
};

vector<Polyomino> polys;
int total_cells = 0;

struct Placement {
    int x, y, rot, flip;
    vector<Point> cells;
};

vector<Placement> best_placements;
int best_w = 0, best_h = 0;

void solve() {
    // Sort polyominoes by size descending
    sort(polys.begin(), polys.end(), [](const Polyomino& a, const Polyomino& b) {
        return a.k > b.k;
    });

    // Generate all transformations for each polyomino
    for (auto& poly : polys) {
        poly.generate_transformations();
    }

    // Simple greedy packing
    int current_w = 1, current_h = 1;
    vector<Placement> placements;
    set<Point> occupied;

    for (const auto& poly : polys) {
        bool placed = false;
        Placement best_placement;
        int best_area = INT_MAX;
        int best_h_for_area = INT_MAX;

        // Try all transformations
        for (int flip = 0; flip <= 1; ++flip) {
            for (int rot = 0; rot < 4; ++rot) {
                int tidx = flip * 4 + rot;
                const auto& transformed = poly.transformations[tidx];

                // Find bounding box of transformed polyomino
                int pw = 0, ph = 0;
                for (const auto& p : transformed) {
                    pw = max(pw, p.x + 1);
                    ph = max(ph, p.y + 1);
                }

                // Try to place it in current grid
                for (int y = 0; y <= current_h + ph; ++y) {
                    for (int x = 0; x <= current_w + pw; ++x) {
                        bool can_place = true;
                        vector<Point> abs_cells;
                        for (const auto& p : transformed) {
                            Point np = {x + p.x, y + p.y};
                            if (occupied.count(np)) {
                                can_place = false;
                                break;
                            }
                            abs_cells.push_back(np);
                        }
                        if (can_place) {
                            // Calculate new bounding box
                            int new_w = current_w;
                            int new_h = current_h;
                            for (const auto& p : abs_cells) {
                                new_w = max(new_w, p.x + 1);
                                new_h = max(new_h, p.y + 1);
                            }
                            int new_area = new_w * new_h;
                            if (new_area < best_area || 
                                (new_area == best_area && new_h < best_h_for_area)) {
                                best_area = new_area;
                                best_h_for_area = new_h;
                                best_placement = {x, y, rot, flip, abs_cells};
                                placed = true;
                            }
                        }
                    }
                }
            }
        }

        if (placed) {
            placements.push_back(best_placement);
            for (const auto& p : best_placement.cells) {
                occupied.insert(p);
                current_w = max(current_w, p.x + 1);
                current_h = max(current_h, p.y + 1);
            }
        } else {
            // Shouldn't happen as we expand the grid
            cerr << "Failed to place polyomino" << endl;
            exit(1);
        }
    }

    best_placements = placements;
    best_w = current_w;
    best_h = current_h;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    polys.resize(n);
    for (int i = 0; i < n; ++i) {
        int k;
        cin >> k;
        polys[i].k = k;
        total_cells += k;
        for (int j = 0; j < k; ++j) {
            int x, y;
            cin >> x >> y;
            polys[i].cells.push_back({x, y});
        }
    }

    solve();

    cout << best_w << " " << best_h << "\n";
    for (const auto& placement : best_placements) {
        cout << placement.x << " " << placement.y << " " 
             << placement.rot << " " << placement.flip << "\n";
    }

    return 0;
}