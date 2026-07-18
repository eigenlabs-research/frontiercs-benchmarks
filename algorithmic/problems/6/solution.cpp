#include <bits/stdc++.h>
using namespace std;

static vector<vector<int>> create_map(int N, int M, vector<int> A, vector<int> B) {
    vector<vector<int>> adj(N+1, vector<int>(N+1, 0));
    vector<vector<int>> g(N+1);
    for (int i=0;i<M;i++) {
        int u=A[i], v=B[i];
        adj[u][v]=adj[v][u]=1;
        g[u].push_back(v); g[v].push_back(u);
    }
    if (N==1) return {{1}};

    mt19937 rng(1234567);

    auto make_sequence = [&](const vector<int>& S)->vector<int>{
        int SZ = (int)S.size();
        vector<int> inS(N+1,0); for(int v:S) inS[v]=1;
        vector<int> best;
        auto try_path = [&](vector<int> order)->vector<int>{
            vector<int> path;
            vector<int> used(N+1,0);
            path.push_back(order[0]); used[order[0]]=1;
            int rem=SZ-1;
            while(rem>0){
                struct Move{int v,pos;};
                vector<Move> moves;
                for(int v: order) if(!used[v]){
                    int L=path.size();
                    if(adj[v][path[0]]) moves.push_back({v,0});
                    if(adj[v][path[L-1]]) moves.push_back({v,L});
                    for(int p=1;p<L;p++) if(adj[v][path[p-1]] && adj[v][path[p]]) moves.push_back({v,p});
                }
                if(moves.empty()) break;
                Move mv = moves[rng()%moves.size()];
                path.insert(path.begin()+mv.pos, mv.v);
                used[mv.v]=1; rem--;
            }
            return path;
        };
        vector<int> ord=S;
        for(int s=0;s<SZ;s++){
            vector<int> o=ord;
            rotate(o.begin(), o.begin()+s, o.end());
            auto p=try_path(o);
            if(p.size()>best.size()) best=p;
            if((int)best.size()==SZ) return best;
        }
        for(int it=0; it<1200 && (int)best.size()<SZ; ++it){
            shuffle(ord.begin(), ord.end(), rng);
            auto p=try_path(ord);
            if(p.size()>best.size()) best=p;
        }
        if((int)best.size()==SZ) return best;

        // Connected-subgraph fallback: DFS tour of a spanning tree of S.
        vector<int> seq, seen(N+1,0);
        function<void(int)> dfs = [&](int u){
            seen[u]=1; seq.push_back(u);
            vector<int> nb;
            for(int v:g[u]) if(inS[v]) nb.push_back(v);
            sort(nb.begin(), nb.end(), [&](int x,int y){return g[x].size()>g[y].size();});
            for(int v:nb) if(!seen[v]) { dfs(v); seq.push_back(u); }
        };
        dfs(S[0]);
        for(int v:S) if(!seen[v]) seq.push_back(v); // should not happen for our candidates
        return seq;
    };

    auto build_for_cover = [&](const vector<int>& S)->vector<vector<int>>{
        vector<int> inS(N+1,0); for(int v:S) inS[v]=1;
        vector<int> seq = make_sequence(S);
        int L = (int)seq.size();
        vector<int> firstOcc(N+1, -1);
        for(int i=0;i<L;i++) if(firstOcc[seq[i]]<0) firstOcc[seq[i]]=i;

        vector<vector<int>> place(N+1);
        bool ok=true;
        for(int i=0;i<M;i++){
            int u=A[i], v=B[i];
            if(inS[v]) place[v].push_back(u);
            else if(inS[u]) place[u].push_back(v);
            else ok=false;
        }
        for(int v=1; v<=N; ++v) if(!inS[v]){
            bool appears=false;
            for(int h:S) for(int x:place[h]) if(x==v) appears=true;
            if(!appears) ok=false;
        }
        if(!ok) return {};
        int maxNeed=0;
        for(int v=1; v<=N; ++v) maxNeed=max(maxNeed, (int)place[v].size());
        int K = max(3*L, 2*maxNeed + 3);
        if(K>240) return {};
        vector<vector<int>> C(K, vector<int>(K, seq[0]));
        for(int b=0;b<L;b++){
            for(int r=3*b; r<min(K,3*b+3); ++r)
                fill(C[r].begin(), C[r].end(), seq[b]);
        }
        for(int r=3*L; r<K; ++r) fill(C[r].begin(), C[r].end(), seq.back());

        for(int v:S){
            int b=firstOcc[v];
            if(b<0) continue;
            int r=3*b+1;
            int idx=0;
            for(int u:place[v]){
                int c=1+2*idx;
                if(c>=K-1) break; // K was sized to prevent this
                C[r][c]=u;
                idx++;
            }
        }
        return C;
    };

    vector<int> all(N); iota(all.begin(), all.end(), 1);
    vector<vector<int>> bestC = build_for_cover(all);

    auto connectedize_cover = [&](vector<int> indep)->vector<int>{
        vector<int> inI(N+1,0); for(int v:indep) inI[v]=1;
        vector<int> inS(N+1,0); for(int v=1; v<=N; ++v) if(!inI[v]) inS[v]=1;
        auto count_comp = [&](){
            vector<int> comp(N+1,-1); int cc=0;
            for(int s=1;s<=N;s++) if(inS[s] && comp[s]<0){
                queue<int> q; q.push(s); comp[s]=cc;
                while(!q.empty()){
                    int u=q.front(); q.pop();
                    for(int v:g[u]) if(inS[v] && comp[v]<0){comp[v]=cc; q.push(v);} 
                }
                cc++;
            }
            return pair<int,vector<int>>(cc,comp);
        };
        for(int iter=0; iter<N; ++iter){
            auto [cc,comp]=count_comp();
            if(cc<=1) break;
            int best=-1, bestcnt=-1;
            for(int v=1; v<=N; ++v) if(inI[v]){
                set<int> touched;
                for(int u:g[v]) if(inS[u] && comp[u]>=0) touched.insert(comp[u]);
                if((int)touched.size()>bestcnt){ bestcnt=touched.size(); best=v; }
            }
            if(best<0 || bestcnt<=0) { for(int v=1; v<=N; ++v) inS[v]=1; break; }
            inI[best]=0; inS[best]=1;
        }
        vector<int> S; for(int v=1; v<=N; ++v) if(inS[v]) S.push_back(v);
        if(S.empty()) S.push_back(1);
        return S;
    };

    // Try connected vertex covers obtained as complements of large independent sets.
    for(int trial=0; trial<220; ++trial){
        vector<int> ord=all;
        if(trial==0) sort(ord.begin(), ord.end(), [&](int x,int y){return g[x].size()<g[y].size();});
        else if(trial==1) sort(ord.begin(), ord.end(), [&](int x,int y){return g[x].size()>g[y].size();});
        else shuffle(ord.begin(), ord.end(), rng);
        vector<int> indep, inI(N+1,0);
        for(int v:ord){
            bool can=true;
            for(int u:g[v]) if(inI[u]) { can=false; break; }
            if(can){ indep.push_back(v); inI[v]=1; }
        }
        vector<int> S = connectedize_cover(indep);
        if((int)S.size()==N) continue;
        auto C = build_for_cover(S);
        if(!C.empty() && (bestC.empty() || C.size()<bestC.size())) bestC.swap(C);
    }
    return bestC;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int T; if(!(cin>>T)) return 0;
    for(int tc=0; tc<T; ++tc){
        int N,M; cin>>N>>M;
        vector<int>A(M),B(M);
        for(int i=0;i<M;i++) cin>>A[i]>>B[i];
        auto C=create_map(N,M,A,B);
        int P=C.size();
        cout << P << '\n';
        for(int i=0;i<P;i++) cout << P << (i+1==P?'\n':' ');
        for(int i=0;i<P;i++){
            for(int j=0;j<P;j++) cout << C[i][j] << (j+1==P?'\n':' ');
        }
    }
    return 0;
}
