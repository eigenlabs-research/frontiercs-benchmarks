#include <cstdio>
#include <vector>
using namespace std;
int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    vector<int> eu(m), ev(m); vector<int> deg(n + 1, 0);
    for (int i = 0; i < m; i++){ scanf("%d %d", &eu[i], &ev[i]); deg[eu[i]]++; deg[ev[i]]++; }
    vector<char> in(n + 1, 0);
    for (int i = 0; i < m; i++){
        int u = eu[i], v = ev[i];
        if (in[u] || in[v]) continue;                 // already covered
        if (deg[u] >= deg[v]) in[u] = 1; else in[v] = 1;
    }
    for (int v = 1; v <= n; v++) printf("%d\n", (int)in[v]);
    return 0;
}
