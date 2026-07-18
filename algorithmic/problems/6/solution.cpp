#include <bits/stdc++.h>
using namespace std;

static vector<vector<int>> adjm;

static vector<int> bfs_path(int s, int t, const vector<vector<int>>& adj){
    int n=adj.size()-1;
    vector<int> par(n+1,-1); queue<int> q; q.push(s); par[s]=0;
    while(!q.empty()){
        int u=q.front(); q.pop(); if(u==t) break;
        for(int v: adj[u]) if(par[v]==-1){ par[v]=u; q.push(v); }
    }
    vector<int> p;
    if(par[t]==-1) return p;
    for(int x=t;x;x=par[x]){ p.push_back(x); if(x==s) break; }
    reverse(p.begin(),p.end()); return p;
}

static vector<int> dfs_tour_order(int n, const vector<vector<int>>& adj){
    vector<int> vis(n+1,0), ord;
    function<void(int)> dfs = [&](int u){
        vis[u]=1; ord.push_back(u);
        vector<int> nb=adj[u];
        sort(nb.begin(), nb.end(), [&](int a,int b){return adj[a].size()>adj[b].size();});
        for(int v: nb) if(!vis[v]){ dfs(v); ord.push_back(u); }
    };
    dfs(1);
    for(int i=1;i<=n;i++) if(!vis[i]) { // should not happen for valid instances with n>1
        if(!ord.empty()){
            auto p=bfs_path(ord.back(), i, adj);
            for(size_t k=1;k<p.size();k++) ord.push_back(p[k]);
        } else dfs(i);
        if(!vis[i]) dfs(i);
    }
    return ord;
}

static bool valid_order(const vector<int>& o){
    for(size_t i=1;i<o.size();++i) if(o[i]!=o[i-1] && !adjm[o[i]][o[i-1]]) return false;
    return true;
}

static vector<int> hamiltonian_like_order(int n, const vector<vector<int>>& adj){
    vector<int> best;
    vector<int> verts(n); iota(verts.begin(), verts.end(), 1);
    sort(verts.begin(), verts.end(), [&](int a,int b){return adj[a].size()>adj[b].size();});

    auto improve_insert = [&](vector<int> path){
        vector<int> used(n+1,0); for(int x:path) used[x]=1;
        bool changed=true;
        while(changed){
            changed=false;
            vector<int> unused;
            for(int v=1;v<=n;v++) if(!used[v]) unused.push_back(v);
            sort(unused.begin(), unused.end(), [&](int a,int b){return adj[a].size()>adj[b].size();});
            for(int v: unused){
                int pos=-1;
                if(adjm[v][path.front()]) pos=0;
                else if(adjm[v][path.back()]) pos=(int)path.size();
                else{
                    for(int i=0;i+1<(int)path.size();i++) if(adjm[path[i]][v] && adjm[v][path[i+1]]){ pos=i+1; break; }
                }
                if(pos!=-1){ path.insert(path.begin()+pos, v); used[v]=1; changed=true; }
            }
        }
        return path;
    };

    for(int st: verts){
        vector<int> path(1, st), used(n+1,0); used[st]=1;
        while(true){
            int L=path.front(), R=path.back();
            int bv=-1, side=0, bs=1e9;
            for(int v=1;v<=n;v++) if(!used[v]){
                if(adjm[v][L]){
                    int cnt=0; for(int w:adj[v]) if(!used[w]) cnt++;
                    if(cnt<bs){bs=cnt; bv=v; side=-1;}
                }
                if(adjm[R][v]){
                    int cnt=0; for(int w:adj[v]) if(!used[w]) cnt++;
                    if(cnt<bs){bs=cnt; bv=v; side=1;}
                }
            }
            if(bv==-1) break;
            used[bv]=1;
            if(side<0) path.insert(path.begin(), bv); else path.push_back(bv);
            path=improve_insert(path);
            fill(used.begin(), used.end(), 0); for(int x:path) used[x]=1;
        }
        path=improve_insert(path);
        if(path.size()>best.size()) best=path;
        if((int)best.size()==n) break;
    }
    if((int)best.size()==n && valid_order(best)) return best;

    // Connect any remaining vertices by shortest paths, allowing repeats; then fall back to a DFS tour if needed.
    vector<int> used(n+1,0); for(int x:best) used[x]=1;
    if(best.empty()) best.push_back(1), used[1]=1;
    while(true){
        int target=0; for(int v=1;v<=n;v++) if(!used[v]){target=v; break;}
        if(!target) break;
        auto p=bfs_path(best.back(), target, adj);
        if(p.empty()) return dfs_tour_order(n, adj);
        for(size_t i=1;i<p.size();i++){ best.push_back(p[i]); used[p[i]]=1; }
    }
    if(valid_order(best)) return best;
    return dfs_tour_order(n, adj);
}

struct Candidate { vector<vector<int>> rows; int K; };

static vector<int> make_row(int bg, const vector<int>& vals, int K){
    vector<int> r(K, bg);
    int col=1;
    for(int v: vals){ if(col>=K) break; r[col]=v; col+=2; }
    return r;
}

static bool check_grid(const vector<vector<int>>& g, int n, const vector<pair<int,int>>& edges){
    int K=g.size(); if(K<1 || K>240) return false;
    vector<vector<int>> seen(n+1, vector<int>(n+1,0));
    vector<int> present(n+1,0);
    for(int i=0;i<K;i++) for(int j=0;j<K;j++){
        int c=g[i][j]; if(c<1||c>n) return false; present[c]=1;
        if(i+1<K && c!=g[i+1][j]){ if(!adjm[c][g[i+1][j]]) return false; seen[c][g[i+1][j]]=seen[g[i+1][j]][c]=1; }
        if(j+1<K && c!=g[i][j+1]){ if(!adjm[c][g[i][j+1]]) return false; seen[c][g[i][j+1]]=seen[g[i][j+1]][c]=1; }
    }
    for(int c=1;c<=n;c++) if(!present[c]) return false;
    for(auto [a,b]: edges) if(!seen[a][b]) return false;
    return true;
}

static Candidate build_three_band(int n, const vector<vector<int>>& adj, const vector<int>& order){
    int maxdeg=0; for(int i=1;i<=n;i++) maxdeg=max(maxdeg,(int)adj[i].size());
    int K=max((int)order.size()*3, 2*maxdeg+1);
    Candidate cand; cand.K=K;
    for(int a: order){
        cand.rows.push_back(vector<int>(K,a));
        cand.rows.push_back(make_row(a, adj[a], K));
        cand.rows.push_back(vector<int>(K,a));
    }
    while((int)cand.rows.size()<K) cand.rows.push_back(vector<int>(K, order.back()));
    return cand;
}

static Candidate build_compressed(int n, const vector<vector<int>>& adj, const vector<pair<int,int>>& edges, const vector<int>& order){
    int maxdeg=0; for(int i=1;i<=n;i++) maxdeg=max(maxdeg,(int)adj[i].size());
    int K=240; // temporary width for row construction; shrink afterwards
    set<pair<int,int>> unc;
    for(auto e: edges) unc.insert(e);
    vector<int> lastpos(n+1,-1);
    for(int i=0;i<(int)order.size();i++) lastpos[order[i]]=i;
    vector<vector<int>> rows;
    auto erase_edge=[&](int a,int b){ if(a>b) swap(a,b); unc.erase({a,b}); };
    auto is_unc=[&](int a,int b){ if(a>b) swap(a,b); return unc.count({a,b}); };
    for(int i=0;i<(int)order.size();i++){
        int a=order[i];
        int prev = (i?order[i-1]:0), next=(i+1<(int)order.size()?order[i+1]:0);
        vector<int> after;
        (void)prev;
        rows.push_back(vector<int>(K,a));

        if(i==lastpos[a]){
            vector<int> rem;
            for(int b: adj[a]) if(is_unc(a,b)) rem.push_back(b);
            if(!rem.empty()){
                rows.push_back(make_row(a, rem, K));
                rows.push_back(vector<int>(K,a));
                for(int b: rem) erase_edge(a,b);
            }
        }
        if(next){
            for(int b: adj[a]) if(is_unc(a,b) && (b==next || adjm[b][next])) after.push_back(b);
        }
        for(int b: after) erase_edge(a,b);
        if(!after.empty()) rows.push_back(make_row(a, after, K)); // between U_a and U_next
    }
    int needW=1;
    for(auto &r: rows){ int last=0; for(int j=0;j<(int)r.size();j++) if(r[j]!=r[0]) last=j; needW=max(needW,last+2); }
    int finalK=max((int)rows.size(), needW);
    Candidate cand; cand.K=finalK;
    for(auto &r: rows){ r.resize(finalK); cand.rows.push_back(r); }
    while((int)cand.rows.size()<finalK) cand.rows.push_back(vector<int>(finalK, order.back()));
    return cand;
}

vector<vector<int>> create_map(int N, int M, vector<int> A, vector<int> B){
    vector<vector<int>> adj(N+1);
    adjm.assign(N+1, vector<int>(N+1,0));
    vector<pair<int,int>> edges;
    for(int i=0;i<M;i++){
        int a=A[i], b=B[i]; adj[a].push_back(b); adj[b].push_back(a); adjm[a][b]=adjm[b][a]=1; edges.push_back({a,b});
    }
    if(N==1) return vector<vector<int>>(1, vector<int>(1,1));
    for(int i=1;i<=N;i++) sort(adj[i].begin(), adj[i].end(), [&](int x,int y){return adj[x].size()>adj[y].size();});
    vector<int> order=hamiltonian_like_order(N, adj);
    if(!valid_order(order)) order=dfs_tour_order(N, adj);

    vector<Candidate> cands;
    cands.push_back(build_three_band(N, adj, order));
    cands.push_back(build_compressed(N, adj, edges, order));
    vector<vector<int>> best;
    int bestK=1000000;
    for(auto &c: cands){
        if(c.K<=240 && check_grid(c.rows, N, edges) && c.K<bestK){ bestK=c.K; best=c.rows; }
    }
    if(best.empty()) best=cands[0].rows;
    return best;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int T; if(!(cin>>T)) return 0;
    while(T--){
        int N,M; cin>>N>>M;
        vector<int>A(M),B(M);
        for(int i=0;i<M;i++) cin>>A[i]>>B[i];
        auto C=create_map(N,M,A,B);
        int K=C.size();
        cout<<K<<"\n";
        for(int i=0;i<K;i++) cout<<K<<(i+1==K?'\n':' ');
        for(int i=0;i<K;i++){
            for(int j=0;j<K;j++) cout<<C[i][j]<<(j+1==K?'\n':' ');
        }
    }
    return 0;
}
