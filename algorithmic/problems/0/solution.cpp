#include <bits/stdc++.h>
using namespace std;

struct Polyomino {
    int k;
    vector<pair<int, int>> cells;
    vector<vector<pair<int, int>>> transforms;

    void generate_transforms() {
        transforms.clear();
        
        // Generate all 8 possible transformations (4 rotations × 2 reflections)
        for (int reflect = 0; reflect <= 1; reflect++) {
            for (int rotate = 0; rotate < 4; rotate++) {
                vector<pair<int, int>> current = cells;
                
                // Apply reflection (across y-axis)
                if (reflect) {
                    for (auto& p : current) p.first = -p.first;
                }
                
                // Apply rotation
                for (int r = 0; r < rotate; r++) {
                    for (auto& p : current) {
                        int x = p.first, y = p.second;
                        p.first = -y;
                        p.second = x;
                    }
                }
                
                // Normalize by translating to minimum (0,0)
                int min_x = INT_MAX, min_y = INT_MAX;
                for (auto& p : current) {
                    min_x = min(min_x, p.first);
                    min_y = min(min_y, p.second);
                }
                for (auto& p : current) {
                    p.first -= min_x;
                    p.second -= min_y;
                }
                
                transforms.push_back(current);
            }
        }
    }
};

struct Placement {
    int x, y, transform_idx;
    vector<pair<int, int>> occupied;
};

vector<Polyomino> polys;
int n;
int best_W = 0, best_H = 0;
vector<Placement> best_placements;

void solve() {
    // Generate all transformations for each polyomino
    for (auto& poly : polys) {
        poly.generate_transforms();
    }

    // Sort polyominoes by size (largest first)
    sort(polys.begin(), polys.end(), [](const Polyomino& a, const Polyomino& b) {
        return a.k > b.k;
    });

    // Initialize grid and placements
    vector<vector<bool>> grid(100, vector<bool>(100, false));
    vector<Placement> placements;

    // Place each polyomino greedily
    for (const auto& poly : polys) {
        Placement best_placement;
        int best_score = INT_MAX;
        int best_w = 0, best_h = 0;

        // Try all transformations
        for (int t = 0; t < poly.transforms.size(); t++) {
            const auto& transformed = poly.transforms[t];
            
            // Find bounding box of transformed polyomino
            int max_x = 0, max_y = 0;
            for (const auto& p : transformed) {
                max_x = max(max_x, p.first);
                max_y = max(max_y, p.second);
            }

            // Try placing at all possible positions
            for (int y = 0; y + max_y < 100; y++) {
                for (int x = 0; x + max_x < 100; x++) {
                    // Check if placement is valid
                    bool valid = true;
                    for (const auto& p : transformed) {
                        int px = x + p.first;
                        int py = y + p.second;
                        if (px >= 100 || py >= 100 || grid[py][px]) {
                            valid = false;
                            break;
                        }
                    }
                    if (!valid) continue;

                    // Calculate score (current bounding box)
                    int min_x = x, max_x_placed = x + max_x;
                    int min_y = y, max_y_placed = y + max_y;
                    for (const auto& pl : placements) {
                        for (const auto& p : pl.occupied) {
                            min_x = min(min_x, p.first);
                            max_x_placed = max(max_x_placed, p.first);
                            min_y = min(min_y, p.second);
                            max_y_placed = max(max_y_placed, p.second);
                        }
                    }
                    int w = max_x_placed - min_x + 1;
                    int h = max_y_placed - min_y + 1;
                    int score = w * h * 1000 + h * 10 + w;

                    if (score < best_score) {
                        best_score = score;
                        best_placement.x = x;
                        best_placement.y = y;
                        best_placement.transform_idx = t;
                        best_w = w;
                        best_h = h;
                    }
                }
            }
        }

        // Place the polyomino
        const auto& transformed = poly.transforms[best_placement.transform_idx];
        for (const auto& p : transformed) {
            int px = best_placement.x + p.first;
            int py = best_placement.y + p.second;
            grid[py][px] = true;
            best_placement.occupied.emplace_back(px, py);
        }
        placements.push_back(best_placement);
    }

    // Calculate final bounding box
    int min_x = INT_MAX, max_x = INT_MIN;
    int min_y = INT_MAX, max_y = INT_MIN;
    for (const auto& pl : placements) {
        for (const auto& p : pl.occupied) {
            min_x = min(min_x, p.first);
            max_x = max(max_x, p.first);
            min_y = min(min_y, p.second);
            max_y = max(max_y, p.second);
        }
    }
    best_W = max_x - min_x + 1;
    best_H = max_y - min_y + 1;

    // Adjust coordinates to start at (0,0)
    for (auto& pl : placements) {
        for (auto& p : pl.occupied) {
            p.first -= min_x;
            p.second -= min_y;
        }
        pl.x -= min_x;
        pl.y -= min_y;
    }

    best_placements = placements;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> n;
    polys.resize(n);
    for (int i = 0; i < n; i++) {
        cin >> polys[i].k;
        polys[i].cells.resize(polys[i].k);
        for (int j = 0; j < polys[i].k; j++) {
            cin >> polys[i].cells[j].first >> polys[i].cells[j].second;
        }
    }

    solve();

    // Output solution
    cout << best_W << " " << best_H << "\n";
    for (int i = 0; i < n; i++) {
        int orig_idx = -1;
        // Find original index (since we sorted)
        for (int j = 0; j < n; j++) {
            if (polys[j].k == polys[i].k && polys[j].cells == polys[i].cells) {
                orig_idx = j;
                break;
            }
        }
        const auto& pl = best_placements[orig_idx];
        int t = pl.transform_idx;
        int reflect = t / 4;
        int rotate = t % 4;
        cout << pl.x << " " << pl.y << " " << rotate << " " << reflect << "\n";
    }

    return 0;
}