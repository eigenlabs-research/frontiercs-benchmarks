#include <bits/stdc++.h>
using namespace std;

static vector<vector<int>> make_walk_map(int N, const vector<pair<int,int>>& edges, const vector<vector<int>>& adj) {
    int M = (int)edges.size();
    vector<vector<pair<int,int>>> g(N+1);
    for (int i=0;i<M;i++) {
        auto [u,v]=edges[i];
        g[u].push_back({v,i}); g[v].push_back({u,i});
    }
    vector<int> seq;
    vector<char> used(M,0);
    function<void(int)> dfs = [&](int u){
        seq.push_back(u);
        for (auto [v,id]: g[u]) if(!used[id]) {
            used[id]=1;
            dfs(v);
            seq.push_back(u);
        }
    };
    int st=1;
    for(int i=1;i<=N;i++) if(!g[i].empty()) { st=i; break; }
    dfs(st);
    if (N==1) seq={1};

    // Put the edge-covering walk on anti-diagonals instead of repeating it in
    // every row.  Adjacent cells then differ by one walk step both horizontally
    // and vertically, so all contacts are legal, while the needed side length is
    // roughly half of the walk length.
    int L=(int)seq.size();
    int K=max(N, (L+1)/2); // need 2*K-1 >= L
    vector<vector<int>> C(K, vector<int>(K));
    for(int i=0;i<K;i++) for(int j=0;j<K;j++) C[i][j]=seq[(i+j)%L];
    return C;
}

static bool verify(int N, const vector<pair<int,int>>& edges, const vector<vector<int>>& C,
                   const vector<vector<char>>& ok) {
    int K=C.size();
    vector<vector<char>> seen(N+1, vector<char>(N+1,0));
    vector<char> has(N+1,0);
    for(int i=0;i<K;i++) for(int j=0;j<K;j++) {
        int a=C[i][j];
        if(a<1||a>N) return false;
        has[a]=1;
        if(i+1<K){ int b=C[i+1][j]; if(a!=b && !ok[a][b]) return false; seen[a][b]=seen[b][a]=1; }
        if(j+1<K){ int b=C[i][j+1]; if(a!=b && !ok[a][b]) return false; seen[a][b]=seen[b][a]=1; }
    }
    for(int i=1;i<=N;i++) if(!has[i]) return false;
    for(auto [u,v]: edges) if(!seen[u][v]) return false;
    return true;
}

vector<vector<int>> create_map(int N, int M, vector<int> A, vector<int> B) {
    vector<pair<int,int>> edges;
    vector<vector<int>> adj(N+1);
    vector<vector<char>> ok(N+1, vector<char>(N+1,0));
    for(int i=1;i<=N;i++) ok[i][i]=1;
    for(int i=0;i<M;i++){
        int u=A[i], v=B[i];
        edges.push_back({u,v});
        adj[u].push_back(v); adj[v].push_back(u);
        ok[u][v]=ok[v][u]=1;
    }

    // For sparse graphs a repeated edge-covering DFS walk is simple and exact.
    if (M==0) return vector<vector<int>>(1, vector<int>(1,1));
    if (2*M+1 <= 479) {
        auto C = make_walk_map(N, edges, adj);
        if (verify(N, edges, C, ok)) return C;
    }

    // Dense/medium fallback: build a graph homomorphic grid greedily.  Every new cell is
    // chosen from the common closed neighbourhood of its already-written top/left cells;
    // among legal choices we bias toward still-unrepresented borders.
    int K = min(240, max(N, min(6*N-1, 5*N)));
    if (K < N) K = N;
    mt19937 rng(71236721u + 1009u*N + 9176u*M);
    vector<vector<int>> best;
    int bestMiss = INT_MAX;

    auto eval_missing = [&](const vector<vector<int>>& C){
        vector<vector<char>> seen(N+1, vector<char>(N+1,0));
        vector<char> has(N+1,0);
        for(int i=0;i<K;i++) for(int j=0;j<K;j++){
            int a=C[i][j]; has[a]=1;
            if(i+1<K){int b=C[i+1][j]; seen[a][b]=seen[b][a]=1;}
            if(j+1<K){int b=C[i][j+1]; seen[a][b]=seen[b][a]=1;}
        }
        int miss=0;
        for(int i=1;i<=N;i++) if(!has[i]) miss += 1000;
        for(auto [u,v]: edges) if(!seen[u][v]) miss++;
        return miss;
    };

    int attempts = 180;
    for(int at=0; at<attempts; at++){
        vector<vector<int>> C(K, vector<int>(K,1));
        vector<vector<char>> seen(N+1, vector<char>(N+1,0));
        vector<int> forced;
        for(int i=1;i<=N;i++) forced.push_back(i);
        shuffle(forced.begin(), forced.end(), rng);
        int fpos=0;
        for(int i=0;i<K;i++) for(int j=0;j<K;j++){
            vector<int> cand;
            for(int c=1;c<=N;c++){
                if(i && !ok[c][C[i-1][j]]) continue;
                if(j && !ok[c][C[i][j-1]]) continue;
                cand.push_back(c);
            }
            if(cand.empty()) cand.push_back(1); // should not happen for connected feasible inputs
            int pick=cand[uniform_int_distribution<int>(0,(int)cand.size()-1)(rng)];
            int bestw=-1;
            // Occasionally force a not-yet-used colour if it is legal.
            if(fpos<N){
                for(int t=fpos;t<N;t++){
                    int c=forced[t];
                    if(find(cand.begin(), cand.end(), c)!=cand.end()) { pick=c; swap(forced[fpos], forced[t]); fpos++; bestw=INT_MAX/2; break; }
                }
            }
            if(bestw<INT_MAX/2){
                for(int c: cand){
                    int w=uniform_int_distribution<int>(0,99)(rng);
                    if(i){ int b=C[i-1][j]; if(c!=b && ok[c][b] && !seen[c][b]) w += 10000; }
                    if(j){ int b=C[i][j-1]; if(c!=b && ok[c][b] && !seen[c][b]) w += 10000; }
                    // prefer colours with larger degree to keep future intersections nonempty
                    w += (int)adj[c].size()*3;
                    if(w>bestw){ bestw=w; pick=c; }
                }
            }
            C[i][j]=pick;
            if(i){ int b=C[i-1][j]; seen[pick][b]=seen[b][pick]=1; }
            if(j){ int b=C[i][j-1]; seen[pick][b]=seen[b][pick]=1; }
        }
        int miss=eval_missing(C);
        if(miss < bestMiss){ bestMiss=miss; best=C; if(miss==0 && verify(N,edges,C,ok)) return C; }
    }

    // Last-resort: if random did not cover everything, return the best valid homomorphic grid.
    // (On the intended dense cases above usually succeeds; for very sparse cases the DFS path did.)
    return best.empty() ? vector<vector<int>>(N, vector<int>(N,1)) : best;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int T; if(!(cin>>T)) return 0;
    // The benchmark checker supplies one scenario per file.  If more are present, emit maps in order.
    for(int tc=0; tc<T; tc++){
        int N,M; cin>>N>>M;
        vector<int>A(M),B(M);
        for(int i=0;i<M;i++) cin>>A[i]>>B[i];
        auto C=create_map(N,M,A,B);
        int P=C.size();
        cout << P << '\n';
        for(int i=0;i<P;i++){ if(i) cout << ' '; cout << (int)C[i].size(); } cout << '\n';
        for(int i=0;i<P;i++){
            for(int j=0;j<(int)C[i].size();j++){ if(j) cout << ' '; cout << C[i][j]; }
            cout << '\n';
        }
    }
    return 0;
}
