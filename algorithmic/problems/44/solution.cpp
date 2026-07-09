#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
using namespace std;

// TSP baseline O(N^2) Nearest Neighbor
struct Point { double x, y; };

bool is_prime(int n) {
    if(n < 2) return false;
    for(int i=2; i*i<=n; ++i) if(n%i==0) return false;
    return true;
}

int main() {
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> pts(n);
    for(int i=0; i<n; ++i) {
        double x, y; cin >> x >> y;
        pts[i] = {x, y};
    }
    
    vector<int> path;
    vector<bool> vis(n, false);
    path.push_back(0);
    vis[0] = true;
    
    for(int step=1; step<n; ++step) {
        int curr = path.back();
        int best_next = -1;
        double best_d = 1e18;
        
        for(int i=0; i<n; ++i) {
            if(!vis[i]) {
                double dx = pts[curr].x - pts[i].x;
                double dy = pts[curr].y - pts[i].y;
                double d = dx*dx + dy*dy;
                if(d < best_d) { best_d = d; best_next = i; }
            }
        }
        path.push_back(best_next);
        vis[best_next] = true;
    }
    
    for(int i=0; i<n; ++i) cout << path[i] << (i==n-1 ? "" : "\n");
    cout << endl;
    return 0;
}
