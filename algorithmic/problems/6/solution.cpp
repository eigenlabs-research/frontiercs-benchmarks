#include <bits/stdc++.h>
using namespace std;

static vector<vector<int>> build_map(int N, int M, const vector<int>& A, const vector<int>& B) {
    if (N == 1) return vector<vector<int>>(1, vector<int>(1, 1));

    vector<vector<int>> adj(N + 1);
    vector<vector<char>> has(N + 1, vector<char>(N + 1, 0));
    for (int i = 0; i < M; ++i) {
        int u = A[i], v = B[i];
        if (!has[u][v]) {
            has[u][v] = has[v][u] = 1;
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
    }
    for (int i = 1; i <= N; ++i) sort(adj[i].begin(), adj[i].end());

    auto prune_walk = [&](vector<int> w) {
        vector<int> cnt(N + 1, 0);
        for (int x : w) cnt[x]++;
        bool changed = true;
        while (changed) {
            changed = false;
            for (int i = 0; i < (int)w.size(); ++i) {
                int x = w[i];
                if (cnt[x] <= 1) continue;
                bool ok = false;
                if (i == 0 || i + 1 == (int)w.size()) ok = true;
                else {
                    int a = w[i - 1], b = w[i + 1];
                    ok = (a == b) || has[a][b];
                }
                if (ok) {
                    w.erase(w.begin() + i);
                    cnt[x]--;
                    changed = true;
                    break;
                }
            }
        }
        return w;
    };

    // Guaranteed safe walk: Euler tour of a BFS spanning tree (length at most 2N-1).
    vector<int> par(N + 1, -1), order;
    queue<int> q;
    par[1] = 0; q.push(1);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        order.push_back(u);
        for (int v : adj[u]) if (par[v] == -1) par[v] = u, q.push(v);
    }
    // The statement promises a drawable instance; for N>1 this should be connected.
    for (int i = 1; i <= N; ++i) if (par[i] == -1) {
        // Defensive fallback for malformed disconnected inputs.
        return vector<vector<int>>(1, vector<int>(1, 1));
    }
    vector<vector<int>> child(N + 1);
    for (int v = 2; v <= N; ++v) child[par[v]].push_back(v);
    vector<int> safe;
    function<void(int)> dfs = [&](int u) {
        safe.push_back(u);
        for (int v : child[u]) { dfs(v); safe.push_back(u); }
    };
    dfs(1);
    safe = prune_walk(safe);
    vector<int> best = safe;

    // Try shorter graph walks: repeatedly go to the nearest unvisited country.
    for (int st = 1; st <= N; ++st) {
        vector<int> w(1, st), seen(N + 1, 0); seen[st] = 1;
        int left = N - 1, cur = st;
        while (left) {
            vector<int> dist(N + 1, -1), pre(N + 1, -1);
            queue<int> qq; dist[cur] = 0; qq.push(cur);
            while (!qq.empty()) {
                int u = qq.front(); qq.pop();
                for (int v : adj[u]) if (dist[v] == -1) {
                    dist[v] = dist[u] + 1; pre[v] = u; qq.push(v);
                }
            }
            int target = -1;
            for (int v = 1; v <= N; ++v) if (!seen[v] && (target == -1 || dist[v] < dist[target] || (dist[v] == dist[target] && v < target))) target = v;
            vector<int> path;
            for (int x = target; x != cur; x = pre[x]) path.push_back(x);
            reverse(path.begin(), path.end());
            for (int x : path) { w.push_back(x); if (!seen[x]) seen[x] = 1, --left; }
            cur = target;
        }
        w = prune_walk(w);
        if (w.size() < best.size()) best = w;
    }

    vector<vector<char>> already(N + 1, vector<char>(N + 1, 0));
    for (int i = 1; i < (int)best.size(); ++i) {
        int u = best[i - 1], v = best[i];
        if (u != v) already[u][v] = already[v][u] = 1;
    }

    vector<pair<int,int>> need;
    for (int u = 1; u <= N; ++u)
        for (int v : adj[u]) if (u < v && !already[u][v]) need.push_back({u, v});

    struct Dinic {
        struct E { int to, rev, cap; };
        vector<vector<E>> g;
        vector<int> level, it;
        Dinic(int n=0): g(n), level(n), it(n) {}
        void addEdge(int v, int to, int cap) {
            E a{to, (int)g[to].size(), cap}, b{v, (int)g[v].size(), 0};
            g[v].push_back(a); g[to].push_back(b);
        }
        bool bfs(int s, int t) {
            fill(level.begin(), level.end(), -1);
            queue<int> q; level[s] = 0; q.push(s);
            while (!q.empty()) {
                int v = q.front(); q.pop();
                for (auto &e : g[v]) if (e.cap && level[e.to] < 0) level[e.to] = level[v] + 1, q.push(e.to);
            }
            return level[t] >= 0;
        }
        int dfs(int v, int t, int f) {
            if (v == t) return f;
            for (int &i = it[v]; i < (int)g[v].size(); ++i) {
                E &e = g[v][i];
                if (!e.cap || level[e.to] != level[v] + 1) continue;
                int ret = dfs(e.to, t, min(f, e.cap));
                if (ret) { e.cap -= ret; g[e.to][e.rev].cap += ret; return ret; }
            }
            return 0;
        }
        int flow(int s, int t) {
            int res = 0, f;
            while (bfs(s, t)) { fill(it.begin(), it.end(), 0); while ((f = dfs(s, t, 1e9))) res += f; }
            return res;
        }
    };

    vector<vector<int>> gadgets(N + 1);
    int bestB = 0;
    for (int Bcap = 0; Bcap <= (int)need.size(); ++Bcap) {
        int Ecnt = need.size(), S = 0, edgeBase = 1, vertBase = edgeBase + Ecnt, T = vertBase + N;
        Dinic din(T + 1);
        for (int i = 0; i < Ecnt; ++i) {
            din.addEdge(S, edgeBase + i, 1);
            din.addEdge(edgeBase + i, vertBase + need[i].first - 1, 1);
            din.addEdge(edgeBase + i, vertBase + need[i].second - 1, 1);
        }
        for (int u = 1; u <= N; ++u) din.addEdge(vertBase + u - 1, T, Bcap);
        if (din.flow(S, T) == Ecnt) {
            bestB = Bcap;
            for (int i = 0; i < Ecnt; ++i) {
                int node = edgeBase + i;
                for (auto &e : din.g[node]) if (e.to >= vertBase && e.to < vertBase + N && e.cap == 0) {
                    int u = e.to - vertBase + 1;
                    int a = need[i].first, b = need[i].second;
                    gadgets[u].push_back(a ^ b ^ u);
                    break;
                }
            }
            break;
        }
    }

    // The min-max flow above minimizes width, but it may spread gadgets over many
    // countries.  Since each used country costs two extra rows, a more concentrated
    // assignment can be better when height is the bottleneck.  Try all per-country
    // capacities with a greedy capacitated vertex-cover assignment and keep it if
    // it improves the actual square side estimate.
    auto side_estimate = [&](const vector<vector<int>>& gg, int widthLoad) {
        int active = 0;
        for (int u = 1; u <= N; ++u) if (!gg[u].empty()) ++active;
        return max((int)best.size() + 2 * active, 2 * widthLoad + 1);
    };
    int chosenLoad = bestB;
    int chosenSide = side_estimate(gadgets, chosenLoad);
    int Eneed = (int)need.size();
    for (int cap = 1; cap <= Eneed && 2 * cap + 1 <= chosenSide; ++cap) {
        vector<char> done(Eneed, 0);
        vector<vector<int>> cand(N + 1);
        int left = Eneed, maxLoad = 0;
        while (left > 0) {
            int bu = -1, bc = -1;
            for (int u = 1; u <= N; ++u) {
                if ((int)cand[u].size() >= cap) continue;
                int c = 0;
                for (int i = 0; i < Eneed; ++i) if (!done[i] && (need[i].first == u || need[i].second == u)) ++c;
                c = min(c, cap - (int)cand[u].size());
                if (c > bc) bc = c, bu = u;
            }
            if (bc <= 0) break; // should not happen for valid cap>=1
            for (int i = 0; i < Eneed && (int)cand[bu].size() < cap; ++i) if (!done[i] && (need[i].first == bu || need[i].second == bu)) {
                done[i] = 1; --left;
                cand[bu].push_back(need[i].first ^ need[i].second ^ bu);
            }
            maxLoad = max(maxLoad, (int)cand[bu].size());
        }
        if (left == 0) {
            int side = side_estimate(cand, maxLoad);
            if (side < chosenSide || (side == chosenSide && maxLoad < chosenLoad)) {
                chosenSide = side;
                chosenLoad = maxLoad;
                gadgets.swap(cand);
            }
        }
    }

    vector<int> first_pos(N + 1, -1);
    for (int g = 0; g < (int)best.size(); ++g) if (first_pos[best[g]] == -1) first_pos[best[g]] = g;
    vector<int> h(best.size(), 1);
    for (int u = 1; u <= N; ++u) if (!gadgets[u].empty()) h[first_pos[u]] = 3;
    int rows = 0;
    for (int x : h) rows += x;
    int K = max(rows, 2 * bestB + 1);
    K = max(K, 1);
    // Safe-walk fallback gives rows <= (2N-1)+2N = 4N-1 <= 159, and width <= 2M/N+1 <= 79.

    vector<vector<int>> C(K, vector<int>(K, best.back()));
    vector<int> start(best.size(), 0);
    int rcur = 0;
    for (int g = 0; g < (int)best.size(); ++g) {
        start[g] = rcur;
        for (int r = rcur; r < rcur + h[g]; ++r) fill(C[r].begin(), C[r].end(), best[g]);
        rcur += h[g];
    }
    for (int u = 1; u <= N; ++u) if (!gadgets[u].empty()) {
        int r = start[first_pos[u]] + 1;
        for (int idx = 0; idx < (int)gadgets[u].size(); ++idx) C[r][2 * idx + 1] = gadgets[u][idx];
    }
    return C;
}

vector<vector<int>> create_map(int N, int M, vector<int> A, vector<int> B) {
    return build_map(N, M, A, B);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int T;
    if (!(cin >> T)) return 0;
    while (T--) {
        int N, M; cin >> N >> M;
        vector<int> A(M), B(M);
        for (int i = 0; i < M; ++i) cin >> A[i] >> B[i];
        auto C = create_map(N, M, A, B);
        int K = (int)C.size();
        cout << K << '\n';
        for (int i = 0; i < K; ++i) cout << K << (i + 1 == K ? '\n' : ' ');
        for (int i = 0; i < K; ++i) {
            for (int j = 0; j < K; ++j) cout << C[i][j] << (j + 1 == K ? '\n' : ' ');
        }
    }
    return 0;
}
