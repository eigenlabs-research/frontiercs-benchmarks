#include <bits/stdc++.h>
using namespace std;

// Maximum clique in a graph on <= 40 vertices, used on the complement graph
// to obtain a large (usually maximum) independent set of the original graph.
struct MaxClique {
    int n;
    vector<unsigned long long> adj;
    vector<int> best, cur;

    void expand(vector<int> P, vector<int> col) {
        while (!P.empty()) {
            if ((int)cur.size() + col.back() <= (int)best.size()) return;
            int v = P.back();
            P.pop_back();
            col.pop_back();
            cur.push_back(v);
            vector<int> P2;
            unsigned long long av = adj[v];
            for (int u : P) if (av & (1ULL << u)) P2.push_back(u);
            if (P2.empty()) {
                if (cur.size() > best.size()) best = cur;
            } else {
                vector<int> col2;
                color_sort(P2, col2);
                expand(P2, col2);
            }
            cur.pop_back();
        }
    }

    void color_sort(vector<int>& P, vector<int>& col) {
        int m = P.size();
        vector<int> rem = P, ordered, colors;
        int c = 0;
        while (!rem.empty()) {
            ++c;
            unsigned long long forbidden = 0;
            vector<int> next;
            for (int v : rem) {
                if (!(forbidden & (1ULL << v))) {
                    ordered.push_back(v);
                    colors.push_back(c);
                    forbidden |= adj[v];
                    forbidden |= (1ULL << v);
                } else next.push_back(v);
            }
            rem.swap(next);
        }
        P.swap(ordered);
        col.swap(colors);
        (void)m;
    }

    vector<int> solve(const vector<unsigned long long>& g) {
        adj = g; n = (int)g.size(); best.clear(); cur.clear();
        vector<int> P(n), col;
        iota(P.begin(), P.end(), 0);
        color_sort(P, col);
        expand(P, col);
        return best;
    }
};

static vector<int> dfs_walk_from(const vector<vector<int>>& adj, int start) {
    int n = (int)adj.size() - 1;
    vector<int> seen(n + 1, 0), walk;
    function<void(int)> dfs = [&](int v) {
        seen[v] = 1;
        walk.push_back(v);
        vector<int> ns = adj[v];
        sort(ns.begin(), ns.end(), [&](int a, int b){
            if (adj[a].size() != adj[b].size()) return adj[a].size() < adj[b].size();
            return a < b;
        });
        for (int to : ns) if (!seen[to]) {
            dfs(to);
            walk.push_back(v);
        }
    };
    dfs(start);
    return walk;
}

static vector<int> greedy_cover_walk(const vector<vector<int>>& adj, int start, int mode) {
    int n = (int)adj.size() - 1;
    vector<int> used(n + 1, 0), walk;
    int cur = start, cnt = 1;
    used[cur] = 1;
    walk.push_back(cur);

    while (cnt < n) {
        vector<int> cand;
        for (int to : adj[cur]) if (!used[to]) cand.push_back(to);
        if (!cand.empty()) {
            sort(cand.begin(), cand.end(), [&](int a, int b) {
                int da = (int)adj[a].size(), db = (int)adj[b].size();
                if (mode == 0) return da < db || (da == db && a < b);
                if (mode == 1) return da > db || (da == db && a < b);
                return a < b;
            });
            cur = cand[0];
            used[cur] = 1; cnt++;
            walk.push_back(cur);
            continue;
        }

        vector<int> par(n + 1, -1), dist(n + 1, -1);
        queue<int> q;
        q.push(cur); par[cur] = cur; dist[cur] = 0;
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (int to : adj[v]) if (dist[to] == -1) {
                dist[to] = dist[v] + 1;
                par[to] = v;
                q.push(to);
            }
        }
        int best = -1;
        for (int v = 1; v <= n; ++v) if (!used[v] && dist[v] != -1) {
            if (best == -1 || dist[v] < dist[best] ||
                (dist[v] == dist[best] && adj[v].size() > adj[best].size())) best = v;
        }
        if (best == -1) break;
        vector<int> p;
        for (int v = best; v != cur; v = par[v]) p.push_back(v);
        reverse(p.begin(), p.end());
        for (int v : p) {
            cur = v;
            walk.push_back(cur);
            if (!used[cur]) { used[cur] = 1; cnt++; }
        }
    }
    return walk;
}

static int row_cost(const vector<int>& walk, const vector<int>& inCover) {
    int ncover = 0;
    for (int i = 1; i < (int)inCover.size(); ++i) ncover += inCover[i];
    int savings = 0;
    if (!walk.empty() && inCover[walk.front()]) savings++;
    if (walk.size() > 1 && inCover[walk.back()] && walk.back() != walk.front()) savings++;
    return (int)walk.size() + 2 * ncover - savings;
}

vector<vector<int>> create_map(int N, int M, vector<int> A, vector<int> B) {
    if (N == 1) return {{1}};

    vector<vector<int>> adj(N + 1);
    vector<vector<int>> has(N + 1, vector<int>(N + 1, 0));
    for (int i = 0; i < M; ++i) {
        int a = A[i], b = B[i];
        if (!has[a][b]) {
            has[a][b] = has[b][a] = 1;
            adj[a].push_back(b);
            adj[b].push_back(a);
        }
    }
    for (int i = 1; i <= N; ++i) sort(adj[i].begin(), adj[i].end());

    // Choose a vertex cover: patterns are needed only for these countries, yet
    // every edge is represented because at least one endpoint is in the cover.
    vector<unsigned long long> comp(N, 0);
    for (int i = 1; i <= N; ++i) for (int j = 1; j <= N; ++j) {
        if (i != j && !has[i][j]) comp[i-1] |= 1ULL << (j-1);
    }
    MaxClique mc;
    vector<int> indep0 = mc.solve(comp);
    vector<int> inInd(N + 1, 0), inCover(N + 1, 1);
    for (int x : indep0) inInd[x + 1] = 1;
    for (int i = 1; i <= N; ++i) inCover[i] = !inInd[i];
    // Defensive repair if the heuristic/exact routine ever misses an edge.
    for (int i = 0; i < M; ++i) if (!inCover[A[i]] && !inCover[B[i]]) {
        if (adj[A[i]].size() <= adj[B[i]].size()) inCover[A[i]] = 1;
        else inCover[B[i]] = 1;
    }

    vector<int> best;
    int bestCost = INT_MAX;
    for (int s = 1; s <= N; ++s) {
        vector<vector<int>> cands;
        cands.push_back(dfs_walk_from(adj, s));
        for (int mode = 0; mode < 3; ++mode) cands.push_back(greedy_cover_walk(adj, s, mode));
        for (auto &w : cands) {
            vector<int> seen(N + 1, 0);
            for (int v : w) seen[v] = 1;
            bool ok = true;
            for (int v = 1; v <= N; ++v) if (!seen[v]) ok = false;
            if (!ok) continue;
            int c = row_cost(w, inCover);
            if (c < bestCost || (c == bestCost && w.size() < best.size())) {
                bestCost = c;
                best = w;
            }
        }
    }
    if (best.empty()) best = dfs_walk_from(adj, 1);

    int maxPattern = 1;
    for (int v = 1; v <= N; ++v) if (inCover[v]) maxPattern = max(maxPattern, 2 * (int)adj[v].size() + 1);

    vector<int> patternAt(best.size(), 0), done(N + 1, 0);
    for (int pass = 0; pass < 2; ++pass) {
        for (int i = 0; i < (int)best.size(); ++i) {
            int v = best[i];
            if (!inCover[v] || done[v]) continue;
            bool endpoint = (i == 0 || i + 1 == (int)best.size());
            if ((pass == 0) != endpoint) continue;
            patternAt[i] = 1; done[v] = 1;
        }
    }
    int rows = 0;
    for (int i = 0; i < (int)best.size(); ++i) {
        if (!patternAt[i]) rows += 1;
        else rows += (i == 0 || i + 1 == (int)best.size()) ? 2 : 3;
    }
    int K = max(rows, maxPattern);
    if (K > 240) K = 240;

    vector<vector<int>> C;
    C.reserve(K);
    auto solid = [&](int v) { return vector<int>(K, v); };
    auto pattern = [&](int v) {
        vector<int> r(K, v);
        int c = 0;
        r[c++] = v;
        for (int to : adj[v]) {
            if (c + 1 >= K) break;
            r[c++] = to;
            r[c++] = v;
        }
        return r;
    };

    for (int i = 0; i < (int)best.size(); ++i) {
        int v = best[i];
        if (!patternAt[i]) {
            C.push_back(solid(v));
        } else if (i == 0) {
            C.push_back(pattern(v));
            C.push_back(solid(v));
        } else if (i + 1 == (int)best.size()) {
            C.push_back(solid(v));
            C.push_back(pattern(v));
        } else {
            C.push_back(solid(v));
            C.push_back(pattern(v));
            C.push_back(solid(v));
        }
    }
    while ((int)C.size() < K) C.push_back(solid(best.back()));
    return C;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int T;
    if (!(cin >> T)) return 0;
    while (T--) {
        int N, M;
        cin >> N >> M;
        vector<int> A(M), B(M);
        for (int i = 0; i < M; ++i) cin >> A[i] >> B[i];
        vector<vector<int>> C = create_map(N, M, A, B);
        int K = (int)C.size();
        cout << K << '\n';
        for (int i = 0; i < K; ++i) cout << K << (i + 1 == K ? '\n' : ' ');
        for (int i = 0; i < K; ++i) {
            for (int j = 0; j < K; ++j) {
                if (j) cout << ' ';
                cout << C[i][j];
            }
            cout << '\n';
        }
    }
    return 0;
}
