#include <bits/stdc++.h>
using namespace std;

struct Point {
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
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
        vector<Point> original = cells;
        for (int reflect = 0; reflect <= 1; ++reflect) {
            vector<Point> reflected;
            if (reflect) {
                for (auto& p : original) {
                    reflected.emplace_back(-p.x, p.y);
                }
            } else {
                reflected = original;
            }
            for (int rotate = 0; rotate < 4; ++rotate) {
                vector<Point> rotated = reflected;
                for (int r = 0; r < rotate; ++r) {
                    for (auto& p : rotated) {
                        int nx = -p.y;
                        int ny = p.x;
                        p.x = nx;
                        p.y = ny;
                    }
                }
                // Normalize
                int min_x = rotated[0].x, min_y = rotated[0].y;
                for (auto& p : rotated) {
                    min_x = min(min_x, p.x);
                    min_y = min(min_y, p.y);
                }
                for (auto& p : rotated) {
                    p.x -= min_x;
                    p.y -= min_y;
                }
                sort(rotated.begin(), rotated.end());
                transformations.push_back(rotated);
            }
        }
        // Remove duplicates
        sort(transformations.begin(), transformations.end());
        transformations.erase(unique(transformations.begin(), transformations.end()), transformations.end());
    }
};

vector<Polyomino> polyominoes;
int n;

struct Placement {
    int x, y, r, f;
    vector<Point> cells;
};

vector<Placement> best_placements;
int best_W = 0, best_H = 0;
int best_area = INT_MAX;

void solve() {
    // Sort polyominoes by size in descending order
    sort(polyominoes.begin(), polyominoes.end(), [](const Polyomino& a, const Polyomino& b) {
        return a.k > b.k;
    });

    // Initialize with a simple greedy packing
    int current_W = 0, current_H = 0;
    vector<Placement> placements;
    vector<Point> occupied;

    for (auto& poly : polyominoes) {
        poly.generate_transformations();
        bool placed = false;
        for (auto& trans : poly.transformations) {
            int trans_idx = &trans - &poly.transformations[0];
            int r = trans_idx % 4;
            int f = trans_idx / 4;
            vector<Point> transformed = trans;
            // Try to place at the first possible position
            for (int y = 0; y <= current_H + 10; ++y) {
                for (int x = 0; x <= current_W + 10; ++x) {
                    bool can_place = true;
                    for (auto& p : transformed) {
                        Point np(p.x + x, p.y + y);
                        if (find(occupied.begin(), occupied.end(), np) != occupied.end()) {
                            can_place = false;
                            break;
                        }
                    }
                    if (can_place) {
                        // Place it here
                        Placement pl;
                        pl.x = x;
                        pl.y = y;
                        pl.r = r;
                        pl.f = f;
                        for (auto& p : transformed) {
                            pl.cells.emplace_back(p.x + x, p.y + y);
                            occupied.push_back(pl.cells.back());
                        }
                        placements.push_back(pl);
                        // Update bounding box
                        for (auto& p : pl.cells) {
                            current_W = max(current_W, p.x + 1);
                            current_H = max(current_H, p.y + 1);
                        }
                        placed = true;
                        goto next_poly;
                    }
                }
            }
            next_poly:
            if (placed) break;
        }
        if (!placed) {
            // Shouldn't happen if algorithm is correct
            cerr << "Failed to place polyomino" << endl;
            exit(1);
        }
    }

    // Try to optimize the bounding box
    int area = current_W * current_H;
    if (area < best_area || (area == best_area && current_H < best_H) || 
        (area == best_area && current_H == best_H && current_W < best_W)) {
        best_area = area;
        best_W = current_W;
        best_H = current_H;
        best_placements = placements;
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> n;
    polyominoes.resize(n);
    for (int i = 0; i < n; ++i) {
        int k;
        cin >> k;
        polyominoes[i].k = k;
        polyominoes[i].cells.resize(k);
        for (int j = 0; j < k; ++j) {
            cin >> polyominoes[i].cells[j].x >> polyominoes[i].cells[j].y;
        }
    }

    // Run multiple times with different random seeds to get better results
    const int ITERATIONS = 5;
    for (int iter = 0; iter < ITERATIONS; ++iter) {
        if (iter > 0) {
            // Shuffle the order of polyominoes
            random_shuffle(polyominoes.begin(), polyominoes.end());
        }
        solve();
    }

    // Output the best solution found
    cout << best_W << " " << best_H << "\n";
    for (auto& pl : best_placements) {
        cout << pl.x << " " << pl.y << " " << pl.r << " " << pl.f << "\n";
    }

    return 0;
}