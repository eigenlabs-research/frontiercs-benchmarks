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

struct Placement {
    int x, y, rot, flip;
    vector<Point> cells;
};

vector<Polyomino> polys;
int n;
int best_area = INT_MAX;
int best_w = 0, best_h = 0;
vector<Placement> best_placements;

void solve() {
    // Sort polyominoes by size in descending order
    sort(polys.begin(), polys.end(), [](const Polyomino& a, const Polyomino& b) {
        return a.k > b.k;
    });
    
    // Initialize with greedy placement
    int current_w = 0, current_h = 0;
    vector<Placement> placements;
    vector<Point> occupied;
    
    for (auto& poly : polys) {
        int best_x = 0, best_y = 0, best_rot = 0, best_flip = 0;
        int best_fit = INT_MAX;
        vector<Point> best_cells;
        
        for (int t = 0; t < poly.transforms.size(); t++) {
            auto& cells = poly.transforms[t];
            int flip = (t >= 4) ? 1 : 0;
            int rot = t % 4;
            
            // Try to place at (0,0) first
            int min_x = 0, min_y = 0;
            for (auto& p : cells) {
                min_x = min(min_x, p.x);
                min_y = min(min_y, p.y);
            }
            
            // Try to find the first available position
            for (int y = 0; y <= current_h + 10; y++) {
                for (int x = 0; x <= current_w + 10; x++) {
                    bool valid = true;
                    vector<Point> transformed;
                    for (auto& p : cells) {
                        Point np = {x + p.x, y + p.y};
                        transformed.push_back(np);
                        
                        if (find(occupied.begin(), occupied.end(), np) != occupied.end()) {
                            valid = false;
                            break;
                        }
                    }
                    
                    if (valid) {
                        int new_w = max(current_w, x + cells.back().x + 1);
                        int new_h = max(current_h, y + cells.back().y + 1);
                        int area = new_w * new_h;
                        
                        if (area < best_fit || 
                            (area == best_fit && new_h < best_h) ||
                            (area == best_fit && new_h == best_h && new_w < best_w)) {
                            best_fit = area;
                            best_x = x;
                            best_y = y;
                            best_rot = rot;
                            best_flip = flip;
                            best_cells = transformed;
                        }
                    }
                }
            }
        }
        
        if (best_fit == INT_MAX) {
            // Shouldn't happen if algorithm is correct
            cerr << "Failed to place polyomino" << endl;
            return;
        }
        
        placements.push_back({best_x, best_y, best_rot, best_flip, best_cells});
        occupied.insert(occupied.end(), best_cells.begin(), best_cells.end());
        current_w = max(current_w, best_x + best_cells.back().x + 1);
        current_h = max(current_h, best_y + best_cells.back().y + 1);
    }
    
    best_area = current_w * current_h;
    best_w = current_w;
    best_h = current_h;
    best_placements = placements;
    
    // Try to improve by swapping placements
    random_device rd;
    mt19937 gen(rd());
    
    auto compute_bbox = [&]() {
        int min_x = INT_MAX, min_y = INT_MAX;
        int max_x = INT_MIN, max_y = INT_MIN;
        for (auto& p : occupied) {
            min_x = min(min_x, p.x);
            min_y = min(min_y, p.y);
            max_x = max(max_x, p.x);
            max_y = max(max_y, p.y);
        }
        return make_tuple(min_x, min_y, max_x, max_y);
    };
    
    for (int iter = 0; iter < 1000; iter++) {
        int i = uniform_int_distribution<>(0, n-1)(gen);
        int j = uniform_int_distribution<>(0, n-1)(gen);
        if (i == j) continue;
        
        auto old_i = placements[i];
        auto old_j = placements[j];
        
        // Remove their cells
        for (auto& p : old_i.cells) {
            occupied.erase(find(occupied.begin(), occupied.end(), p));
        }
        for (auto& p : old_j.cells) {
            occupied.erase(find(occupied.begin(), occupied.end(), p));
        }
        
        // Try to place j in i's position and vice versa
        bool can_place_j = true;
        vector<Point> new_j_cells;
        for (auto& p : old_j.cells) {
            Point np = {old_i.x + (p.x - old_j.x), old_i.y + (p.y - old_j.y)};
            if (find(occupied.begin(), occupied.end(), np) != occupied.end()) {
                can_place_j = false;
                break;
            }
            new_j_cells.push_back(np);
        }
        
        bool can_place_i = true;
        vector<Point> new_i_cells;
        for (auto& p : old_i.cells) {
            Point np = {old_j.x + (p.x - old_i.x), old_j.y + (p.y - old_i.y)};
            if (find(occupied.begin(), occupied.end(), np) != occupied.end()) {
                can_place_i = false;
                break;
            }
            new_i_cells.push_back(np);
        }
        
        if (can_place_j && can_place_i) {
            placements[i] = {old_j.x, old_j.y, old_j.rot, old_j.flip, new_i_cells};
            placements[j] = {old_i.x, old_i.y, old_i.rot, old_i.flip, new_j_cells};
            occupied.insert(occupied.end(), new_i_cells.begin(), new_i_cells.end());
            occupied.insert(occupied.end(), new_j_cells.begin(), new_j_cells.end());
            
            auto [min_x, min_y, max_x, max_y] = compute_bbox();
            int new_w = max_x - min_x + 1;
            int new_h = max_y - min_y + 1;
            int new_area = new_w * new_h;
            
            if (new_area < best_area || 
                (new_area == best_area && new_h < best_h) ||
                (new_area == best_area && new_h == best_h && new_w < best_w)) {
                best_area = new_area;
                best_w = new_w;
                best_h = new_h;
            }
        } else {
            // Revert changes
            occupied.insert(occupied.end(), old_i.cells.begin(), old_i.cells.end());
            occupied.insert(occupied.end(), old_j.cells.begin(), old_j.cells.end());
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cin >> n;
    polys.resize(n);
    for (int i = 0; i < n; i++) {
        int k;
        cin >> k;
        polys[i].k = k;
        polys[i].cells.resize(k);
        for (int j = 0; j < k; j++) {
            cin >> polys[i].cells[j].x >> polys[i].cells[j].y;
        }
        polys[i].generate_transforms();
    }
    
    solve();
    
    cout << best_w << " " << best_h << "\n";
    for (auto& p : best_placements) {
        cout << p.x << " " << p.y << " " << p.rot << " " << p.flip << "\n";
    }
    
    return 0;
}