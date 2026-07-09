#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>
using namespace std;
using namespace std::chrono;

struct Point { double x, y; };

bool is_prime(int n) {
    if(n < 2) return false;
    for(int i=2; i*i<=n; ++i) if(n%i==0) return false;
    return true;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> pts(n);
    vector<bool> prime(n, false);
    for(int i=0; i<n; ++i) {
        cin >> pts[i].x >> pts[i].y;
        prime[i] = is_prime(i);
    }
    
    vector<int> path(n);
    for(int i=0; i<n; ++i) path[i] = i; // Input order baseline
    
    auto get_dist = [&](int a, int b) {
        double dx = pts[a].x - pts[b].x;
        double dy = pts[a].y - pts[b].y;
        return sqrt(dx*dx + dy*dy);
    };
    
    auto eval_full = [&]() {
        double cost = 0;
        for(int t=1; t<=n; ++t) {
            int u = path[t-1];
            int v = (t == n) ? 0 : path[t];
            double m = (t % 10 == 0 && !prime[u]) ? 1.1 : 1.0;
            cost += m * get_dist(u, v);
        }
        return cost;
    };

    double best_cost = eval_full();
    auto start_time = high_resolution_clock::now();

    while (true) {
        auto now = high_resolution_clock::now();
        if (duration_cast<milliseconds>(now - start_time).count() > 1850) break;
        
        bool improved = false;
        for (int i = 1; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                // Approximate 2-opt check without 10% penalty for speed
                int a = path[i-1], b = path[i];
                int c = path[j], d = (j+1 == n) ? 0 : path[j+1];
                double d1 = get_dist(a, b) + get_dist(c, d);
                double d2 = get_dist(a, c) + get_dist(b, d);
                if (d2 + 1e-9 < d1) {
                    reverse(path.begin() + i, path.begin() + j + 1);
                    improved = true;
                }
            }
        }
        if (!improved) break;
    }
    
    path.push_back(0);
    cout << path.size() << "\n";
    for(int i=0; i<path.size(); ++i) cout << path[i] << "\n";
    return 0;
}
