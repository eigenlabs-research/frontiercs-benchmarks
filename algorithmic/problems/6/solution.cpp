#include <cstdio>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <unordered_set>
#include <random>
#include <chrono>
#include <cstring>
#include <array>
using namespace std;
typedef long long ll;
typedef pair<int,int> pii;

static int N, M;
static bool ADJ[45][45];
static vector<pii> EDGES;
static mt19937 rng(987654321u);

// ---------------- Hamiltonian path (Warnsdorff + restarts) ----------------
static vector<int> bestHP, cur;
static vector<char> inpath;
static bool hpFound;
static chrono::steady_clock::time_point deadline;
static inline bool timeUp(){ return chrono::steady_clock::now() > deadline; }

static bool dfs_hp(int v, int d){
    if(d==N){ hpFound=true; bestHP=cur; return true; }
    if((d & 7)==0 && timeUp()) return false;
    int nb[45], c=0;
    for(int u=1;u<=N;u++) if(ADJ[v][u] && !inpath[u]) nb[c++]=u;
    sort(nb, nb+c, [&](int a,int b){
        int ca=0,cb=0;
        for(int w=1;w<=N;w++){ if(ADJ[a][w]&&!inpath[w])ca++; if(ADJ[b][w]&&!inpath[w])cb++; }
        return ca<cb;
    });
    for(int i=0;i<c;i++){
        int u=nb[i];
        inpath[u]=1; cur.push_back(u);
        if(dfs_hp(u,d+1)) return true;
        cur.pop_back(); inpath[u]=0;
        if(hpFound) return true;
        if(timeUp()) return false;
    }
    return false;
}

static bool findHP(double tl){
    hpFound=false; bestHP.clear();
    deadline = chrono::steady_clock::now() + chrono::milliseconds((int)(tl*1000));
    inpath.assign(N+1,0);
    vector<int> starts(N); for(int i=0;i<N;i++) starts[i]=i+1;
    shuffle(starts.begin(), starts.end(), rng);
    for(int s: starts){
        if(hpFound || timeUp()) break;
        cur.clear(); fill(inpath.begin(), inpath.end(), 0);
        inpath[s]=1; cur.push_back(s);
        dfs_hp(s,1);
    }
    return hpFound;
}

// ---------------- DFS walk of spanning tree (fallback) ----------------
static vector<int> buildWalk(){
    vector<int> par(N+1,0); vector<char> vis(N+1,0);
    vector<int> q; q.push_back(1); vis[1]=1; int head=0;
    while(head<(int)q.size()){
        int v=q[head++];
        for(int u=1;u<=N;u++) if(ADJ[v][u] && !vis[u]){ vis[u]=1; par[u]=v; q.push_back(u); }
    }
    vector<vector<int>> ch(N+1);
    for(int v=1;v<=N;v++) if(par[v]) ch[par[v]].push_back(v);
    vector<int> w; w.push_back(1);
    vector<pair<int,int>> st; st.push_back({1,0});
    while(!st.empty()){
        auto& bk = st.back();
        int v=bk.first, &ci=bk.second;
        if(ci < (int)ch[v].size()){
            int c=ch[v][ci++];
            w.push_back(c);
            st.push_back({c,0});
        } else {
            st.pop_back();
            if(!st.empty()) w.push_back(st.back().first);
        }
    }
    return w;
}

// ---------------- Construct from a given path ----------------
// returns K and fills grid (K x K). Also returns K for selection.
static int constructFrom(const vector<int>& path, vector<vector<int>>& grid){
    set<pii> pathEdges;
    for(int i=0;i+1<(int)path.size();i++){
        int a=path[i], b=path[i+1]; if(a>b) swap(a,b);
        pathEdges.insert({a,b});
    }
    vector<pii> chords;
    for(auto& e: EDGES) if(!pathEdges.count(e)) chords.push_back(e);

    vector<vector<int>> chordNbr(N+1);
    for(auto& c: chords){ chordNbr[c.first].push_back(c.second); chordNbr[c.second].push_back(c.first); }

    // vertex cover of chord-graph: multiple randomized greedy passes, keep smallest
    int C = (int)chords.size();
    vector<char> bestS(N+1,0); int bestSize = 1000000;
    int passes = 40;
    for(int pass=0; pass<passes; pass++){
        vector<char> inS(N+1,0);
        vector<char> covered(C,0);
        int remaining = C;
        while(remaining>0){
            // find max uncovered chord-degree, random tie-break
            vector<int> cd(N+1,0);
            for(int i=0;i<C;i++) if(!covered[i]){ cd[chords[i].first]++; cd[chords[i].second]++; }
            int bd=0; for(int v=1;v<=N;v++) bd=max(bd,cd[v]);
            if(bd==0) break;
            vector<int> cand; for(int v=1;v<=N;v++) if(cd[v]==bd) cand.push_back(v);
            int bv = cand[rng()%cand.size()];
            inS[bv]=1;
            for(int i=0;i<C;i++) if(!covered[i] && (chords[i].first==bv||chords[i].second==bv)){ covered[i]=1; remaining--; }
        }
        // prune redundant
        bool changed=true;
        while(changed){
            changed=false;
            for(int v=1;v<=N;v++) if(inS[v]){
                bool allInS=true;
                for(int u: chordNbr[v]) if(!inS[u]){ allInS=false; break; }
                if(allInS){ inS[v]=0; changed=true; }
            }
        }
        int sz=0; for(int v=1;v<=N;v++) sz+=inS[v];
        if(sz<bestSize){ bestSize=sz; bestS=inS; }
    }
    vector<char> inS = bestS;
    // assign chords to an endpoint in S (balance blob counts)
    vector<int> blobCount(N+1,0);
    vector<vector<pii>> blockBlobs(N+1);
    for(auto& c: chords){
        int a=c.first, b=c.second, chosen, other;
        if(inS[a] && inS[b]) chosen = (blobCount[a] <= blobCount[b]) ? a : b;
        else if(inS[a]) chosen = a;
        else chosen = b;
        other = (chosen==a)? b : a;
        int col = 2*blobCount[chosen];
        blockBlobs[chosen].push_back({col, other});
        blobCount[chosen]++;
    }
    // row sequence with 3-row expansions for S vertices (at first occurrence)
    vector<int> rows;
    vector<char> expanded(N+1,0);
    map<int,int> blockMiddle;
    for(int v: path){
        if(inS[v] && !expanded[v]){
            rows.push_back(v); rows.push_back(v); rows.push_back(v);
            expanded[v]=1;
            blockMiddle[v] = (int)rows.size()-2;
        } else rows.push_back(v);
    }
    int R = (int)rows.size();
    int maxBlobs = 0;
    for(int v=1;v<=N;v++) maxBlobs = max(maxBlobs, blobCount[v]);
    int W = max(1, 2*maxBlobs - 1);
    int K = max(R, W);
    if(K > 240) K = 240;
    grid.assign(K, vector<int>(K, 0));
    for(int r=0;r<K;r++){
        int cc = (r<R)? rows[r] : rows[R-1];
        for(int c=0;c<K;c++) grid[r][c]=cc;
    }
    for(int v=1;v<=N;v++) if(inS[v] && blockMiddle.count(v)){
        int mid = blockMiddle[v];
        for(auto& bp: blockBlobs[v]){
            int col = bp.first, other = bp.second;
            if(mid>=0 && mid<K && col>=0 && col<K) grid[mid][col]=other;
        }
    }
    return K;
}

// ---------------- min-conflicts local search for a fixed K ----------------
// seed: K x K grid (colors 1..N). Reach cost 0 within deadline. Returns true if valid.
static bool localSearch(vector<vector<int>>& seedGrid, int K,
                        chrono::steady_clock::time_point dl,
                        vector<vector<int>>& result){
    if((int)seedGrid.size() < K || (int)seedGrid[0].size() < K) return false;
    vector<int> g(K*K);
    for(int r=0;r<K;r++) for(int c=0;c<K;c++) g[r*K+c]=seedGrid[r][c];
    static int cnt[45][45];
    for(int a=0;a<45;a++) for(int b=0;b<45;b++) cnt[a][b]=0;
    int forbidden=0, missing=0, missingColors=0;
    vector<int> colorCount(N+2,0);
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
        int col=g[r*K+c]; colorCount[col]++;
        if(c+1<K){ int col2=g[r*K+c+1]; if(col!=col2){ int a=min(col,col2),b=max(col,col2); cnt[a][b]++; if(!ADJ[a][b])forbidden++; } }
        if(r+1<K){ int col2=g[(r+1)*K+c]; if(col!=col2){ int a=min(col,col2),b=max(col,col2); cnt[a][b]++; if(!ADJ[a][b])forbidden++; } }
    }
    for(auto& e: EDGES){ if(cnt[e.first][e.second]==0) missing++; }
    for(int c=1;c<=N;c++) if(colorCount[c]==0) missingColors++;
    if(forbidden==0 && missing==0 && missingColors==0){
        result.assign(K, vector<int>(K));
        for(int r=0;r<K;r++) for(int c=0;c<K;c++) result[r][c]=g[r*K+c];
        return true;
    }
    vector<char> inBad(K*K,0); vector<int> badVec;
    const int DR[4]={1,-1,0,0}, DC[4]={0,0,1,-1};
    auto cellInForbidden=[&](int p)->bool{
        int r=p/K, c=p%K, col=g[p];
        for(int d=0;d<4;d++){ int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int nc2=g[nr*K+nc]; if(col!=nc2 && !ADJ[min(col,nc2)][max(col,nc2)]) return true; }
        return false;
    };
    auto updateBad=[&](int p){ bool b=cellInForbidden(p); if(b && !inBad[p]){ inBad[p]=1; badVec.push_back(p);} else if(!b && inBad[p]) inBad[p]=0; };
    for(int p=0;p<K*K;p++) updateBad(p);
    vector<vector<int>> cellsOf(N+2);
    auto rebuildCellsOf=[&](){ for(int c=1;c<=N;c++) cellsOf[c].clear(); for(int p=0;p<K*K;p++) cellsOf[g[p]].push_back(p); };
    rebuildCellsOf();
    std::uniform_real_distribution<double> prob(0.0,1.0);
    long long iters=0;
    auto now=[&](){ return chrono::steady_clock::now(); };
    while(now() < dl){
        iters++;
        if((iters & 511)==0 && forbidden==0 && missing==0 && missingColors==0){
            result.assign(K, vector<int>(K));
            for(int r=0;r<K;r++) for(int c=0;c<K;c++) result[r][c]=g[r*K+c];
            return true;
        }
        int p=-1;
        if(prob(rng) < 0.02){
            p = rng() % (K*K);
        } else if(forbidden>0){
            for(int tries=0; tries<12 && !badVec.empty(); tries++){ int idx=rng()%badVec.size(); int cand=badVec[idx]; if(inBad[cand]){ p=cand; break; } else { badVec[idx]=badVec.back(); badVec.pop_back(); } }
            if(p<0) p = rng()%(K*K);
        } else if(missing>0){
            int eidx=-1;
            for(int tries=0; tries<20; tries++){ int e=rng()%EDGES.size(); if(cnt[EDGES[e].first][EDGES[e].second]==0){ eidx=e; break; } }
            if(eidx>=0){
                int a=EDGES[eidx].first,b=EDGES[eidx].second;
                int ep = (prob(rng)<0.5)?a:b; int other=(ep==a)?b:a;
                if(!cellsOf[ep].empty()){
                    int baseCell = cellsOf[ep][rng()%cellsOf[ep].size()];
                    if(g[baseCell]!=ep){ /* stale */ }
                    int r=baseCell/K, c=baseCell%K;
                    int d=rng()%4; int nr=r+DR[d],nc=c+DC[d];
                    if(nr>=0&&nr<K&&nc>=0&&nc<K) p=nr*K+nc;
                    (void)other;
                }
            }
            if(p<0) p=rng()%(K*K);
        } else { p=rng()%(K*K); }
        if(p<0) continue;
        int r=p/K, c=p%K, o=g[p];
        int nbr[4], nnbr=0;
        for(int d=0;d<4;d++){ int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; nbr[nnbr++]=g[nr*K+nc]; }
        long long bestDelta=0; int bestX=o, ties=1;
        bool first=true;
        for(int x=1;x<=N;x++){
            if(x==o) continue;
            int df=0, dmc=0;
            int touched[8], de[8], nt=0;
            for(int i=0;i<nnbr;i++){
                int nc=nbr[i];
                if(o!=nc){
                    int a=min(o,nc),b=max(o,nc);
                    if(ADJ[a][b]){ int e=a*45+b; int j; for(j=0;j<nt;j++) if(touched[j]==e) break; if(j==nt){touched[nt]=e;de[nt]=0;nt++;} de[j]-=1; }
                    else df--;
                }
                if(x!=nc){
                    int a=min(x,nc),b=max(x,nc);
                    if(ADJ[a][b]){ int e=a*45+b; int j; for(j=0;j<nt;j++) if(touched[j]==e) break; if(j==nt){touched[nt]=e;de[nt]=0;nt++;} de[j]+=1; }
                    else df++;
                }
            }
            int dm=0;
            for(int j=0;j<nt;j++){ int e=touched[j]; int a=e/45,b=e%45; int oldc=cnt[a][b]; int newc=oldc+de[j]; if(oldc>0&&newc==0)dm++; else if(oldc==0&&newc>0)dm--; }
            if(colorCount[o]==1) dmc++; if(colorCount[x]==0) dmc--;
            long long delta = 1000LL*df + 50LL*dm + 50LL*dmc;
            if(first){ bestDelta=delta; bestX=x; ties=1; first=false; }
            else if(delta < bestDelta){ bestDelta=delta; bestX=x; ties=1; }
            else if(delta==bestDelta){ ties++; if((int)(rng()%ties)==0) bestX=x; }
        }
        int x=bestX;
        if(!first && x!=o){
            for(int i=0;i<nnbr;i++){
                int nc=nbr[i];
                if(o!=nc){ int a=min(o,nc),b=max(o,nc); int oldc=cnt[a][b]; cnt[a][b]=oldc-1; if(ADJ[a][b]){ if(oldc==1) missing++; } else forbidden--; }
                if(x!=nc){ int a=min(x,nc),b=max(x,nc); int oldc=cnt[a][b]; cnt[a][b]=oldc+1; if(ADJ[a][b]){ if(oldc==0) missing--; } else forbidden++; }
            }
            colorCount[o]--; if(colorCount[o]==0) missingColors++;
            colorCount[x]++; if(colorCount[x]==1) missingColors--;
            g[p]=x;
            updateBad(p);
            for(int d=0;d<4;d++){ int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; updateBad(nr*K+nc); }
            if((iters & 127)==0) rebuildCellsOf();
        }
    }
    if(forbidden==0 && missing==0 && missingColors==0){
        result.assign(K, vector<int>(K));
        for(int r=0;r<K;r++) for(int c=0;c<K;c++) result[r][c]=g[r*K+c];
        return true;
    }
    return false;
}

// produce a (K-1)x(K-1) seed from a KxK grid by removing the least-valuable row and column
static vector<vector<int>> shrinkByOne(const vector<vector<int>>& G){
    int K=(int)G.size();
    if(K<=1){ return G; }
    // per-color total counts and per-row / per-col counts
    vector<int> totalC(N+2,0);
    vector<vector<int>> rowC(K, vector<int>(N+2,0)), colC(K, vector<int>(N+2,0));
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){ int v=G[r][c]; totalC[v]++; rowC[r][v]++; colC[c][v]++; }
    auto rowLoss=[&](int r)->long long{
        long long loss=0;
        if(r>0) for(int c=0;c<K;c++) if(G[r-1][c]!=G[r][c]) loss++;
        if(r+1<K) for(int c=0;c<K;c++) if(G[r][c]!=G[r+1][c]) loss++;
        for(int v=1;v<=N;v++) if(rowC[r][v]==totalC[v] && totalC[v]>0) loss+=100000; // would lose a color
        return loss;
    };
    auto colLoss=[&](int c)->long long{
        long long loss=0;
        if(c>0) for(int r=0;r<K;r++) if(G[r][c-1]!=G[r][c]) loss++;
        if(c+1<K) for(int r=0;r<K;r++) if(G[r][c]!=G[r][c+1]) loss++;
        for(int v=1;v<=N;v++) if(colC[c][v]==totalC[v] && totalC[v]>0) loss+=100000;
        return loss;
    };
    int bestR=0; long long lr=rowLoss(0);
    for(int r=1;r<K;r++){ long long x=rowLoss(r); if(x<lr){ lr=x; bestR=r; } }
    int bestC=0; long long lc=colLoss(0);
    for(int c=1;c<K;c++){ long long x=colLoss(c); if(x<lc){ lc=x; bestC=c; } }
    // build (K-1)x(K-1) excluding bestR and bestC
    vector<vector<int>> out(K-1, vector<int>(K-1));
    int rr=0;
    for(int r=0;r<K;r++){ if(r==bestR) continue; int cc=0; for(int c=0;c<K;c++){ if(c==bestC) continue; out[rr][cc]=G[r][c]; cc++; } rr++; }
    return out;
}

// verify a grid satisfies all checker conditions (defensive)
static bool verifyGrid(const vector<vector<int>>& grid){
    int K=(int)grid.size();
    if(K<1 || K>240) return false;
    for(auto& row: grid){
        if((int)row.size()!=K) return false;
        for(int v: row) if(v<1 || v>N) return false;
    }
    vector<char> present(N+1,0);
    for(auto& row: grid) for(int v: row) present[v]=1;
    for(int c=1;c<=N;c++) if(!present[c]) return false;
    vector<char> realized(M,0);
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
        int v=grid[r][c];
        if(c+1<K){ int u=grid[r][c+1]; if(v!=u){ if(!ADJ[min(v,u)][max(v,u)]) return false; } }
        if(r+1<K){ int u=grid[r+1][c]; if(v!=u){ if(!ADJ[min(v,u)][max(v,u)]) return false; } }
    }
    // mark realized edges
    static map<pii,int> edgeId;
    edgeId.clear();
    for(int i=0;i<M;i++) edgeId[EDGES[i]]=i;
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
        int v=grid[r][c];
        if(c+1<K){ int u=grid[r][c+1]; if(v!=u){ auto it=edgeId.find({min(v,u),max(v,u)}); if(it!=edgeId.end()) realized[it->second]=1; } }
        if(r+1<K){ int u=grid[r+1][c]; if(v!=u){ auto it=edgeId.find({min(v,u),max(v,u)}); if(it!=edgeId.end()) realized[it->second]=1; } }
    }
    for(int i=0;i<M;i++) if(!realized[i]) return false;
    return true;
}

int main(){
    int T; if(scanf("%d",&T)!=1) return 0;
    while(T--){
        scanf("%d %d",&N,&M);
        memset(ADJ,0,sizeof(ADJ));
        EDGES.clear();
        for(int i=0;i<M;i++){
            int a,b; scanf("%d %d",&a,&b);
            ADJ[a][b]=ADJ[b][a]=true;
            EDGES.push_back({min(a,b),max(a,b)});
        }
        // gather candidate paths (variety helps local search escape minima)
        vector<vector<int>> cands;
        auto hpEnd = chrono::steady_clock::now() + chrono::milliseconds(420);
        if(findHP(0.20)) cands.push_back(bestHP);
        for(int attempt=0; attempt<4 && chrono::steady_clock::now()<hpEnd; attempt++){
            if(findHP(0.05)){ cands.push_back(bestHP); }
        }
        cands.push_back(buildWalk()); // fallback always available

        int bestK = 1000000; vector<vector<int>> bestGrid, constructiveGrid;
        for(auto& p : cands){
            if((int)p.size() < N) continue; // walk covers all; safety
            vector<vector<int>> g;
            int K = constructFrom(p, g);
            if(K < bestK){ bestK = K; bestGrid = g; constructiveGrid = g; }
        }
        if(!verifyGrid(bestGrid)){
            // constructive should always be valid; if not, re-derive from walk
            vector<vector<int>> g; constructFrom(buildWalk(), g); bestGrid=g; constructiveGrid=g; bestK=(int)g.size();
        }
        // ---- local search: progressively shrink K ----
        auto tStart = chrono::steady_clock::now();
        auto tEnd = tStart + chrono::milliseconds(870);
        int failStreak=0;
        for(int K=bestK-1; K>=2 && chrono::steady_clock::now()<tEnd; K--){
            auto rem = chrono::duration_cast<chrono::milliseconds>(tEnd - chrono::steady_clock::now()).count();
            long long slice = max(10, (int)(rem/4));
            bool success=false; vector<vector<int>> res;
            // attempt 1: top-left crop
            {
                vector<vector<int>> seed(K, vector<int>(K));
                for(int r=0;r<K;r++) for(int c=0;c<K;c++) seed[r][c]=bestGrid[r][c];
                auto sl = chrono::steady_clock::now() + chrono::milliseconds((int)(slice*3/5));
                if(localSearch(seed, K, sl, res) && verifyGrid(res)) success=true;
            }
            // attempt 2: shrinkByOne (if crop failed and time remains)
            if(!success && chrono::steady_clock::now()<tEnd){
                vector<vector<int>> seed = shrinkByOne(bestGrid);
                if((int)seed.size()==K){
                    auto sl = chrono::steady_clock::now() + chrono::milliseconds((int)(slice*2/5));
                    if(localSearch(seed, K, sl, res) && verifyGrid(res)) success=true;
                }
            }
            if(success){ bestK=K; bestGrid=res; failStreak=0; }
            else { if(++failStreak>=4) break; }
        }
        if(!verifyGrid(bestGrid)) bestGrid = constructiveGrid; // ultimate fallback
        // output (buffered for speed)
        {
            int K = (int)bestGrid.size();
            string out;
            out.reserve((size_t)K*K*4 + K*8);
            out += to_string(K); out += '\n';
            for(int j=0;j<K;j++){ out += to_string(K); if(j+1<K) out += ' '; } out += '\n';
            for(int i=0;i<K;i++){
                for(int j=0;j<K;j++){ out += to_string(bestGrid[i][j]); if(j+1<K) out += ' '; }
                out += '\n';
            }
            fwrite(out.data(), 1, out.size(), stdout);
        }
    }
    return 0;
}
