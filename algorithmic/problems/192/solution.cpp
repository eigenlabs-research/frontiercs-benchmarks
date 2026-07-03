#include <cstdio>
#include <vector>
#include <random>
#include <chrono>
using namespace std;

// Max-Cut: multi-restart randomized local search under a ~0.9s budget.
// Each restart: random partition, then first-improvement flip sweeps until a
// local optimum; plateau kicks (random flips) escape zero-gain traps.
int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    vector<vector<int>> adj(n + 1);
    vector<pair<int,int>> edges(m);
    for (int i = 0; i < m; i++){
        int u, v; scanf("%d %d", &u, &v);
        adj[u].push_back(v); adj[v].push_back(u);
        edges[i] = {u, v};
    }

    auto t0 = chrono::steady_clock::now();
    auto elapsed = [&]{ return chrono::duration<double>(chrono::steady_clock::now() - t0).count(); };
    const double BUDGET = 0.90;

    mt19937 rng(12345);
    vector<int> side(n + 1), best(n + 1, 0);
    int bestCut = -1;

    auto cutOf = [&](const vector<int>& s){
        int c = 0;
        for (auto& e : edges) c += s[e.first] != s[e.second];
        return c;
    };

    // cross[v] = # neighbors on the other side; flipping v changes cut by cross_after - cross_before
    vector<int> crossCnt(n + 1);
    while (elapsed() < BUDGET){
        for (int v = 1; v <= n; v++) side[v] = rng() & 1;
        for (int v = 1; v <= n; v++){
            int c = 0;
            for (int w : adj[v]) c += side[w] != side[v];
            crossCnt[v] = c;
        }
        int kicks = 3;
        while (elapsed() < BUDGET){
            bool improved = false;
            for (int v = 1; v <= n; v++){
                int deg = (int)adj[v].size();
                // gain of flipping v = (deg - cross) - cross
                if (deg - 2 * crossCnt[v] > 0){
                    side[v] ^= 1;
                    crossCnt[v] = deg - crossCnt[v];
                    for (int w : adj[v]){
                        if (side[w] != side[v]) crossCnt[w]++;
                        else                    crossCnt[w]--;
                    }
                    improved = true;
                }
            }
            if (!improved){
                if (kicks-- <= 0) break;
                // plateau kick: flip a few random vertices, then keep searching
                for (int k = 0; k < 1 + n / 50; k++){
                    int v = 1 + (int)(rng() % n);
                    side[v] ^= 1;
                    crossCnt[v] = (int)adj[v].size() - crossCnt[v];
                    for (int w : adj[v]){
                        if (side[w] != side[v]) crossCnt[w]++;
                        else                    crossCnt[w]--;
                    }
                }
            }
        }
        int c = cutOf(side);
        if (c > bestCut){ bestCut = c; best = side; }
    }

    if (bestCut < 0) best.assign(n + 1, 0);
    for (int v = 1; v <= n; v++) printf("%d%c", best[v], v == n ? '\n' : ' ');
    return 0;
}
