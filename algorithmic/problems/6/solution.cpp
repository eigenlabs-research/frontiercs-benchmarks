#include <bits/stdc++.h>
using namespace std;

static const int INF = 1e9;

vector<vector<int>> build_complete(int N, const vector<pair<int,int>>& edges){
    vector<int> s;
    s.reserve(edges.size()*2 + 5);
    for(auto [a,b]: edges){ s.push_back(a); s.push_back(b); }
    if(s.empty()) s.push_back(1);
    int K = max(N, (int)ceil(sqrt((double)s.size())));
    K = max(K, 1);
    vector<vector<int>> C(K, vector<int>(K, 1));
    int p=0;
    for(int i=0;i<K;i++){
        if(i%2==0){
            for(int j=0;j<K;j++) if(p<(int)s.size()) C[i][j]=s[p++];
        }else{
            for(int j=K-1;j>=0;j--) if(p<(int)s.size()) C[i][j]=s[p++];
        }
    }
    return C;
}

vector<int> euler_walk(int N, const vector<pair<int,int>>& orig, const vector<vector<int>>& extraPaths){
    vector<vector<pair<int,int>>> adj(N+1);
    int id=0;
    auto addEdge=[&](int u,int v){
        adj[u].push_back({v,id}); adj[v].push_back({u,id}); id++;
    };
    for(auto [u,v]: orig) addEdge(u,v);
    for(auto &p: extraPaths) for(int i=0;i+1<(int)p.size();i++) addEdge(p[i],p[i+1]);
    int start=1;
    for(int i=1;i<=N;i++) if(adj[i].size()%2){ start=i; break; }
    if(start==1) for(int i=1;i<=N;i++) if(!adj[i].empty()){ start=i; break; }
    vector<int> it(N+1,0), used(id,0), st, ans;
    st.push_back(start);
    while(!st.empty()){
        int v=st.back();
        while(it[v]<(int)adj[v].size() && used[adj[v][it[v]].second]) it[v]++;
        if(it[v]==(int)adj[v].size()){
            ans.push_back(v); st.pop_back();
        }else{
            auto [to, eid]=adj[v][it[v]++];
            if(!used[eid]){ used[eid]=1; st.push_back(to); }
        }
    }
    reverse(ans.begin(), ans.end());
    return ans;
}

vector<vector<int>> create_map(int N, int M, vector<int> A, vector<int> B){
    if(N==1) return {{1}};
    vector<pair<int,int>> edges;
    vector<vector<int>> g(N+1);
    vector<vector<int>> has(N+1, vector<int>(N+1,0));
    for(int i=0;i<M;i++){
        int a=A[i], b=B[i];
        edges.push_back({a,b});
        g[a].push_back(b); g[b].push_back(a);
        has[a][b]=has[b][a]=1;
    }
    if(M == N*(N-1)/2) return build_complete(N, edges);

    // All-pairs shortest paths in the input graph (unweighted).
    vector<vector<int>> dist(N+1, vector<int>(N+1, INF)), par(N+1, vector<int>(N+1,-1));
    for(int s=1;s<=N;s++){
        queue<int> q; dist[s][s]=0; par[s][s]=s; q.push(s);
        while(!q.empty()){
            int v=q.front(); q.pop();
            for(int to:g[v]) if(dist[s][to]==INF){
                dist[s][to]=dist[s][v]+1; par[s][to]=v; q.push(to);
            }
        }
    }
    auto path_between = [&](int s,int t){
        vector<int> p;
        if(dist[s][t]>=INF) return p;
        int cur=t; p.push_back(cur);
        while(cur!=s){ cur=par[s][cur]; p.push_back(cur); }
        reverse(p.begin(), p.end());
        return p;
    };

    vector<int> odd;
    for(int i=1;i<=N;i++) if(g[i].size()%2) odd.push_back(i);
    int O=odd.size();

    vector<pair<int,int>> pairs;
    if(O <= 22){
        int FULL = (1<<O) - 1;
        vector<int> dp(1<<O, INF), choice(1<<O, -1), leave(1<<O, 0);
        dp[0]=0;
        for(int mask=1; mask<=FULL; mask++){
            int bits=__builtin_popcount((unsigned)mask);
            if(bits==1 || bits==2){ dp[mask]=0; leave[mask]=1; continue; } // endpoints of open trail
            int i=__builtin_ctz((unsigned)mask);
            int m2 = mask ^ (1<<i);
            for(int j=i+1;j<O;j++) if(m2&(1<<j)){
                int nm=m2^(1<<j), val=dp[nm]+dist[odd[i]][odd[j]];
                if(val<dp[mask]){ dp[mask]=val; choice[mask]=(i<<8)|j; }
            }
        }
        int mask=FULL;
        while(mask && !leave[mask]){
            int ch=choice[mask];
            if(ch<0) break;
            int i=ch>>8, j=ch&255;
            pairs.push_back({odd[i], odd[j]});
            mask ^= (1<<i); mask ^= (1<<j);
        }
    }else{
        // Large odd set: leave a far-apart pair as trail endpoints, greedily pair the rest.
        vector<int> alive(O,1); int bi=0,bj=1,bd=-1;
        for(int i=0;i<O;i++) for(int j=i+1;j<O;j++) if(dist[odd[i]][odd[j]]>bd){bd=dist[odd[i]][odd[j]];bi=i;bj=j;}
        alive[bi]=alive[bj]=0;
        while(true){
            int ai=-1,aj=-1,ad=INF;
            for(int i=0;i<O;i++) if(alive[i]) for(int j=i+1;j<O;j++) if(alive[j] && dist[odd[i]][odd[j]]<ad){ad=dist[odd[i]][odd[j]];ai=i;aj=j;}
            if(ai<0) break;
            alive[ai]=alive[aj]=0; pairs.push_back({odd[ai], odd[aj]});
        }
    }

    vector<vector<int>> extra;
    for(auto [u,v]: pairs) extra.push_back(path_between(u,v));
    vector<int> walk = euler_walk(N, edges, extra);

    int L = (int)walk.size();
    int K = max(N, (L+2)/2); // diagonals 0..2K-2 contain the whole walk
    K = min(K, 240);
    vector<vector<int>> C(K, vector<int>(K, walk.empty()?1:walk.back()));
    for(int i=0;i<K;i++) for(int j=0;j<K;j++){
        int d=i+j;
        if(d<L) C[i][j]=walk[d];
        else C[i][j]=walk.back();
    }
    return C;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int T; if(!(cin>>T)) return 0;
    for(int tc=0; tc<T; tc++){
        int N,M; cin>>N>>M;
        vector<int>A(M),B(M);
        for(int i=0;i<M;i++) cin>>A[i]>>B[i];
        auto C = create_map(N,M,A,B);
        int K=C.size();
        cout << K << "\n";
        for(int i=0;i<K;i++) cout << K << (i+1==K?'\n':' ');
        for(int i=0;i<K;i++){
            for(int j=0;j<K;j++) cout << C[i][j] << (j+1==K?'\n':' ');
        }
    }
    return 0;
}
