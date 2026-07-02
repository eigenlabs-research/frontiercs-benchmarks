#include <cstdio>
#include <vector>
using namespace std;
int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    vector<vector<int>> adj(n + 1);
    for (int i = 0; i < m; i++){ int u, v; scanf("%d %d", &u, &v); adj[u].push_back(v); adj[v].push_back(u); }
    vector<int> side(n + 1, 0);
    for (int v = 1; v <= n; v++) side[v] = v & 1;          // seed: alternate
    // Local search: flip any vertex that has more same-side than cross-side neighbors.
    bool improved = true; int passes = 0;
    while (improved && passes++ < 50){
        improved = false;
        for (int v = 1; v <= n; v++){
            int same = 0, cross = 0;
            for (int w : adj[v]) (side[w] == side[v] ? same : cross)++;
            if (same > cross){ side[v] ^= 1; improved = true; }
        }
    }
    for (int v = 1; v <= n; v++) printf("%d%c", side[v], v == n ? '\n' : ' ');
    return 0;
}
