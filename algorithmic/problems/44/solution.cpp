#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <algorithm>
using namespace std;
using namespace std::chrono;

struct Point { double x, y; };

bool is_prime(int n) {
    if(n < 2) return false;
    if(n == 2 || n == 3) return true;
    if(n % 2 == 0 || n % 3 == 0) return false;
    for(int i=5; i*i<=n; i+=6) {
        if(n%i==0 || n%(i+2)==0) return false;
    }
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
    for(int i=0; i<n; ++i) path[i] = i; 
    
    auto get_dist = [&](int a, int b) {
        double dx = pts[a].x - pts[b].x;
        double dy = pts[a].y - pts[b].y;
        return sqrt(dx*dx + dy*dy);
    };
    
    auto eval_full = [&](const vector<int>& p) {
        double cost = 0;
        for(int t=1; t<=n; ++t) {
            int u = p[t-1];
            int v = (t == n) ? 0 : p[t];
            double m = (t % 10 == 0 && !prime[u]) ? 1.1 : 1.0;
            cost += m * get_dist(u, v);
        }
        return cost;
    };

    double best_cost = eval_full(path);
    vector<int> best_path = path;
    double current_cost = best_cost;
    
    auto start_time = high_resolution_clock::now();
    mt19937 rng(1337);
    
    double T = 10000.0;
    double T_min = 0.001;
    double alpha = 0.99999;
    
    int iterations = 0;
    while (true) {
        if ((iterations & 255) == 0) {
            auto now = high_resolution_clock::now();
            if (duration_cast<milliseconds>(now - start_time).count() > 1850) break;
        }
        iterations++;
        
        int i = 1 + rng() % (n - 2);
        int j = 1 + rng() % (n - 2);
        if (i == j) continue;
        if (i > j) swap(i, j);
        
        int a = path[i-1], b = path[i];
        int c = path[j], d = (j+1 == n) ? 0 : path[j+1];
        
        double delta = 0;
        
        for(int t = i; t <= j+1; ++t) {
            int u = path[t-1];
            int v = (t == n) ? 0 : path[t];
            double m = (t % 10 == 0 && !prime[u]) ? 1.1 : 1.0;
            delta -= m * get_dist(u, v);
        }
        
        reverse(path.begin() + i, path.begin() + j + 1);
        
        for(int t = i; t <= j+1; ++t) {
            int u = path[t-1];
            int v = (t == n) ? 0 : path[t];
            double m = (t % 10 == 0 && !prime[u]) ? 1.1 : 1.0;
            delta += m * get_dist(u, v);
        }
        
        if (delta < 0.0 || exp(-delta / T) > (double)rng() / rng.max()) {
            current_cost += delta;
            if (current_cost < best_cost) {
                best_cost = current_cost;
                best_path = path;
            }
        } else {
            reverse(path.begin() + i, path.begin() + j + 1);
        }
        
        T *= alpha;
        if (T < T_min) T = T_min;
    }
    
    best_path.push_back(0);
    cout << best_path.size() << "\n";
    for(size_t i=0; i<best_path.size(); ++i) {
        cout << best_path[i] << "\n";
    }
    return 0;
}
