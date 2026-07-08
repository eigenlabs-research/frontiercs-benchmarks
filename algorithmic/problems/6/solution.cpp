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
static mt19937 rng(246813579u);
static chrono::steady_clock::time_point HARD_DL; // global hard deadline

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
    if(deadline > HARD_DL) deadline = HARD_DL;
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

// ---------------- vertex cover of chord set (randomized greedy) ----------------
static vector<char> computeVC(const vector<pii>& chords){
    int C = (int)chords.size();
    vector<vector<int>> chordNbr(N+1);
    for(auto& c: chords){ chordNbr[c.first].push_back(c.second); chordNbr[c.second].push_back(c.first); }
    vector<char> bestS(N+1,0); int bestSize = 1000000;
    int passes = 30;
    for(int pass=0; pass<passes; pass++){
        vector<char> inS(N+1,0);
        vector<char> covered(C,0);
        int remaining = C;
        while(remaining>0){
            vector<int> cd(N+1,0);
            for(int i=0;i<C;i++) if(!covered[i]){ cd[chords[i].first]++; cd[chords[i].second]++; }
            int bd=0; for(int v=1;v<=N;v++) bd=max(bd,cd[v]);
            if(bd==0) break;
            vector<int> cand; for(int v=1;v<=N;v++) if(cd[v]==bd) cand.push_back(v);
            int bv = cand[rng()%cand.size()];
            inS[bv]=1;
            for(int i=0;i<C;i++) if(!covered[i] && (chords[i].first==bv||chords[i].second==bv)){ covered[i]=1; remaining--; }
        }
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
    return bestS;
}

// ---------------- stripes construction (rows; safe fallback) ----------------
static int constructFrom(const vector<int>& path, vector<vector<int>>& grid){
    set<pii> pathEdges;
    for(int i=0;i+1<(int)path.size();i++){
        int a=path[i], b=path[i+1]; if(a>b) swap(a,b);
        pathEdges.insert({a,b});
    }
    vector<pii> chords;
    for(auto& e: EDGES) if(!pathEdges.count(e)) chords.push_back(e);
    vector<char> inS(N+1,0);
    if(!chords.empty()) inS = computeVC(chords);
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
    vector<int> rows;
    vector<char> expanded(N+1,0);
    map<int,int> blockMiddle;
    for(int v: path){
        if(inS[v] && blobCount[v]>0 && !expanded[v]){
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

// ---------------- diagonal construction ----------------
// Lay a walk along anti-diagonals: cell (r,c) = seq[min(r+c, L-1)].
// Cover vertices of the chord set are tripled; chords realized as isolated
// blob cells on the (pure) middle diagonal of the triple.
static int constructDiag(const vector<int>& walk, vector<vector<int>>& grid){
    int W = (int)walk.size();
    if(W==0) return 100000;
    for(int i=0;i+1<W;i++) if(walk[i]==walk[i+1]) return 100000;
    set<pii> wE;
    for(int i=0;i+1<W;i++){ int a=walk[i], b=walk[i+1]; wE.insert({min(a,b),max(a,b)}); }
    vector<pii> chords;
    for(auto& e: EDGES) if(!wE.count(e)) chords.push_back(e);
    vector<char> inS(N+1,0);
    if(!chords.empty()) inS = computeVC(chords);
    vector<int> blobCount(N+1,0);
    vector<vector<int>> blobs(N+1);
    for(auto& c: chords){
        int a=c.first, b=c.second, chosen;
        if(inS[a] && inS[b]) chosen = (blobCount[a] <= blobCount[b]) ? a : b;
        else if(inS[a]) chosen = a;
        else chosen = b;
        int other = (chosen==a)? b : a;
        blobs[chosen].push_back(other); blobCount[chosen]++;
    }
    // build tripled sequence for an orientation; return needed front padding
    auto tryBuild=[&](bool rev, vector<int>& seq, vector<int>& midPos)->int{
        vector<int> w = walk; if(rev) reverse(w.begin(), w.end());
        seq.clear(); midPos.assign(N+1,-1);
        vector<char> done(N+1,0);
        for(int v: w){
            if(blobCount[v]>0 && !done[v]){
                done[v]=1;
                seq.push_back(v); seq.push_back(v); seq.push_back(v);
                midPos[v]=(int)seq.size()-2;
            } else seq.push_back(v);
        }
        int pad=0;
        for(int v=1;v<=N;v++) if(blobCount[v]>0 && midPos[v]>=0)
            pad=max(pad, blobCount[v]-1-midPos[v]);
        return pad;
    };
    vector<int> seqA, midA, seqB, midB;
    int padA=tryBuild(false, seqA, midA), padB=tryBuild(true, seqB, midB);
    vector<int>& seq = (padB<padA)? seqB : seqA;
    vector<int>& midPos = (padB<padA)? midB : midA;
    int pad = min(padA, padB);
    if(pad>0){
        vector<int> ns((size_t)pad, seq[0]);
        ns.insert(ns.end(), seq.begin(), seq.end());
        seq = ns;
        for(int v=1;v<=N;v++) if(midPos[v]>=0) midPos[v]+=pad;
    }
    int L=(int)seq.size();
    int K=(L+2)/2; if(K<1) K=1;
    auto capOK=[&](int k)->bool{
        for(int v=1;v<=N;v++) if(blobCount[v]>0){
            int t=midPos[v];
            int len = min(min(t+1,k), 2*k-1-t);
            if(len < blobCount[v]) return false;
        }
        return true;
    };
    while(K<=240 && !capOK(K)) K++;
    if(K>240) return 100000;
    grid.assign(K, vector<int>(K));
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
        int t=r+c; grid[r][c]=seq[min(t,L-1)];
    }
    for(int v=1;v<=N;v++) if(blobCount[v]>0){
        int t=midPos[v];
        int i0=max(0,t-K+1), i1=min(t,K-1);
        int idx=0;
        for(int i=i0;i<=i1 && idx<(int)blobs[v].size();i++)
            grid[i][t-i]=blobs[v][idx++];
        if(idx<(int)blobs[v].size()) return 100000;
    }
    return K;
}

// ---------------- generalized diagonal construction (v2) ----------------
// Blobs may sit on any diagonal t (1<=t<=L-2) whose side colors seq[t-1],
// seq[t+1] are compatible (== blob or adjacent) with one side the chord
// partner. Even-diagonal blobs pack from the top end of the diagonal, odd
// ones from the bottom end; margin constraints keep blobs non-adjacent and
// keep a pure contact pair across every diagonal boundary. Vertices are
// tripled only when a chord has no anchor at all.
static int constructDiagV2(const vector<int>& walk, vector<vector<int>>& grid){
    int W = (int)walk.size();
    if(W==0) return 100000;
    for(int i=0;i+1<W;i++) if(walk[i]==walk[i+1]) return 100000;
    set<pii> wE;
    for(int i=0;i+1<W;i++){ int a=walk[i], b=walk[i+1]; wE.insert({min(a,b),max(a,b)}); }
    vector<pii> chords;
    for(auto& e: EDGES) if(!wE.count(e)) chords.push_back(e);
    int C=(int)chords.size();
    vector<char> tripled(N+1,0);
    int pad=0;
    for(int round=0; round<12; round++){
        // ---- build sequence ----
        vector<int> seq((size_t)pad, walk[0]);
        vector<char> done(N+1,0);
        for(int v: walk){
            if(tripled[v] && !done[v]){
                done[v]=1;
                seq.push_back(v); seq.push_back(v); seq.push_back(v);
            } else seq.push_back(v);
        }
        int L=(int)seq.size();
        if(L>460) return 100000;
        // ---- candidate diagonals per chord (independent of K) ----
        vector<vector<int>> posOf(N+1);
        for(int i=0;i<L;i++) posOf[seq[i]].push_back(i);
        vector<vector<pii>> cand(C); // (t, blobColor)
        for(int ci=0; ci<C; ci++){
            set<pii> cs;
            for(int dir=0; dir<2; dir++){
                int anchor = dir? chords[ci].second : chords[ci].first;
                int blob   = dir? chords[ci].first  : chords[ci].second;
                for(int idx: posOf[anchor]){
                    for(int dt=-1; dt<=1; dt+=2){
                        int t=idx+dt;
                        if(t<1 || t>L-2) continue;
                        int s1=seq[t-1], s2=seq[t+1];
                        bool ok1 = (s1==blob) || ADJ[min(s1,blob)][max(s1,blob)];
                        bool ok2 = (s2==blob) || ADJ[min(s2,blob)][max(s2,blob)];
                        if(ok1 && ok2) cs.insert({t, blob});
                    }
                }
            }
            cand[ci].assign(cs.begin(), cs.end());
        }
        // ---- try assignment for increasing K ----
        int K0=(L+2)/2;
        vector<int> order(C); for(int i=0;i<C;i++) order[i]=i;
        sort(order.begin(), order.end(), [&](int a,int b){ return cand[a].size()<cand[b].size(); });
        bool anyAnchorless=false;
        for(int i=0;i<C;i++) if(cand[i].empty()) anyAnchorless=true;
        if(!anyAnchorless){
            for(int K=K0; K<=min(240,K0+8); K++){
                auto I0=[&](int t){ return max(0,t-K+1); };
                auto I1=[&](int t){ return min(t,K-1); };
                auto diaglen=[&](int t){ return I1(t)-I0(t)+1; };
                bool useAssignFallback = (N>=22 && N<=26 && 4LL*M > 1LL*N*(N-1));
                if(!useAssignFallback){
                    vector<int> load(L,0);
                    auto marginAfter=[&](int t)->int{
                        int ld=load[t]+1;
                        int m = (diaglen(t)-1) - ld;
                        for(int dt=-1; dt<=1; dt+=2){
                            int u=t+dt; if(u<0 || u>=L) continue;
                            int e = (t%2==0)? t : u;
                            int o = (t%2==0)? u : t;
                            int le = (e==t)? ld : load[e];
                            int lo = (o==t)? ld : load[o];
                            int lim;
                            if(e < o) lim = I1(o)-I0(e)-1;
                            else      lim = I1(o)-I0(e);
                            m = min(m, lim - (le+lo));
                        }
                        return m;
                    };
                    vector<pii> assign(C,{-1,-1});
                    bool ok=true;
                    for(int oi=0; oi<C && ok; oi++){
                        int ci=order[oi];
                        int bt=-1, bblob=-1, bslack=-1;
                        for(auto& tc: cand[ci]){
                            int slack = marginAfter(tc.first);
                            if(slack>bslack){ bslack=slack; bt=tc.first; bblob=tc.second; }
                        }
                        if(bt<0 || bslack<0){ ok=false; break; }
                        load[bt]++; assign[ci]={bt,bblob};
                    }
                    if(!ok) continue;
                    grid.assign(K, vector<int>(K));
                    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
                        int t=r+c; grid[r][c]=seq[min(t,L-1)];
                    }
                    vector<int> used(L,0);
                    for(int ci=0; ci<C; ci++){
                        int t=assign[ci].first, blob=assign[ci].second;
                        int i;
                        if(t%2==0) i = I0(t)+used[t];   // even: pack from top
                        else       i = I1(t)-used[t];   // odd: pack from bottom
                        used[t]++;
                        grid[i][t-i]=blob;
                    }
                    return K;
                }
                int attempts = (C>120 ? 3 : 5);
                for(int attempt=0; attempt<attempts; attempt++){
                    vector<int> load(L,0);
                    auto marginAfter=[&](int t)->int{
                        int ld=load[t]+1;
                        int m = (diaglen(t)-1) - ld;
                        for(int dt=-1; dt<=1; dt+=2){
                            int u=t+dt; if(u<0 || u>=L) continue;
                            int e = (t%2==0)? t : u;
                            int o = (t%2==0)? u : t;
                            int le = (e==t)? ld : load[e];
                            int lo = (o==t)? ld : load[o];
                            int lim;
                            if(e < o) lim = I1(o)-I0(e)-1;
                            else      lim = I1(o)-I0(e);
                            m = min(m, lim - (le+lo));
                        }
                        return m;
                    };
                    vector<int> ord=order;
                    if(attempt>0){
                        shuffle(ord.begin(), ord.end(), rng);
                        stable_sort(ord.begin(), ord.end(), [&](int a,int b){ return cand[a].size()<cand[b].size(); });
                    }
                    vector<pii> assign(C,{-1,-1});
                    bool ok=true;
                    for(int oi=0; oi<C && ok; oi++){
                        int ci=ord[oi];
                        int bt=-1, bblob=-1, bestScore=-1000000000;
                        for(auto& tc: cand[ci]){
                            int slack = marginAfter(tc.first);
                            if(slack<0) continue;
                            int score = slack*1000;
                            if(attempt>0){
                                score -= load[tc.first]*37;
                                score += (int)(rng()%97);
                            }
                            if(score>bestScore){
                                bestScore=score; bt=tc.first; bblob=tc.second;
                            }
                        }
                        if(bt<0){ ok=false; break; }
                        load[bt]++; assign[ci]={bt,bblob};
                    }
                    if(!ok) continue;
                    grid.assign(K, vector<int>(K));
                    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
                        int t=r+c; grid[r][c]=seq[min(t,L-1)];
                    }
                    vector<int> used(L,0);
                    for(int ci=0; ci<C; ci++){
                        int t=assign[ci].first, blob=assign[ci].second;
                        int i;
                        if(t%2==0) i = I0(t)+used[t];
                        else       i = I1(t)-used[t];
                        used[t]++;
                        grid[i][t-i]=blob;
                    }
                    return K;
                }
            }
        }
        // ---- failure: fix and retry ----
        if(anyAnchorless){
            // greedy vertex cover of anchorless chords among untripled vertices
            vector<pii> bad;
            for(int i=0;i<C;i++) if(cand[i].empty()) bad.push_back(chords[i]);
            vector<int> cd(N+1,0);
            for(auto& b: bad){ cd[b.first]++; cd[b.second]++; }
            vector<char> chosen(N+1,0);
            while(true){
                int bv=-1, bd=0;
                for(int v=1;v<=N;v++) if(!tripled[v] && !chosen[v] && cd[v]>bd){ bd=cd[v]; bv=v; }
                if(bv<0) break;
                chosen[bv]=1;
                bool any=false;
                vector<int> cd2(N+1,0);
                for(auto& b: bad){
                    if(chosen[b.first]||chosen[b.second]) continue;
                    cd2[b.first]++; cd2[b.second]++; any=true;
                }
                cd=cd2;
                if(!any) break;
            }
            bool progress=false;
            for(int v=1;v<=N;v++) if(chosen[v]){ tripled[v]=1; progress=true; }
            if(!progress) return 100000;
        } else {
            pad += 2; // capacity trouble on early diagonals: shift right
        }
    }
    return 100000;
}

// ---------------- greedy short covering walk over all vertices ----------------
static vector<int> greedyCoverWalk(int start){
    vector<char> vis(N+1,0);
    vector<int> walk; walk.push_back(start); vis[start]=1;
    int cur=start, cnt=1;
    while(cnt<N){
        vector<int> par(N+1,0), dist(N+1,-1);
        vector<int> q; q.push_back(cur); dist[cur]=0; int head=0;
        int foundDist=-1;
        vector<int> found;
        while(head<(int)q.size()){
            int v=q[head++];
            if(!vis[v]){
                if(foundDist<0) foundDist=dist[v];
                if(dist[v]>foundDist) break;
                found.push_back(v);
                continue; // don't expand beyond the frontier of interest
            }
            if(foundDist>=0 && dist[v]>=foundDist) continue;
            for(int u=1;u<=N;u++) if(ADJ[v][u] && dist[u]<0){ dist[u]=dist[v]+1; par[u]=v; q.push_back(u); }
        }
        if(found.empty()) return {};
        // prefer targets with the fewest unvisited neighbors (clear dead-ends first)
        int bestT=-1, bestScore=1<<29, ties=0;
        for(int v: found){
            int s=0; for(int u=1;u<=N;u++) if(ADJ[v][u] && !vis[u]) s++;
            if(s<bestScore){ bestScore=s; bestT=v; ties=1; }
            else if(s==bestScore){ ties++; if((int)(rng()%ties)==0) bestT=v; }
        }
        vector<int> path;
        for(int x=bestT; x!=cur; x=par[x]) path.push_back(x);
        reverse(path.begin(), path.end());
        for(int v: path){ walk.push_back(v); if(!vis[v]){ vis[v]=1; cnt++; } }
        cur=bestT;
    }
    return walk;
}

// ---------------- postman covering walk (short walk covering ALL edges) ----------------
static vector<int> postmanWalk(){
    if(M==0) return {};
    static int dist[41][41], par[41][41];
    for(int s=1;s<=N;s++){
        for(int v=1;v<=N;v++){ dist[s][v]=-1; par[s][v]=0; }
        vector<int> q; q.push_back(s); dist[s][s]=0; int head=0;
        while(head<(int)q.size()){
            int v=q[head++];
            for(int u=1;u<=N;u++) if(ADJ[v][u] && dist[s][u]<0){
                dist[s][u]=dist[s][v]+1; par[s][u]=v; q.push_back(u);
            }
        }
    }
    static int mult[41][41];
    memset(mult, 0, sizeof(mult));
    int deg[41]; memset(deg,0,sizeof(deg));
    for(auto& e: EDGES){ mult[e.first][e.second]=mult[e.second][e.first]=1; deg[e.first]++; deg[e.second]++; }
    vector<int> odd;
    for(int v=1;v<=N;v++) if(deg[v]&1) odd.push_back(v);
    // randomized greedy matchings + 2-swap; objective = total dup after
    // dropping the worst pair (its endpoints become the open-trail ends)
    int P=(int)odd.size();
    vector<pii> bestPairs; long long bestCost=1LL<<60;
    for(int rs=0; rs<10 && P>0; rs++){
        vector<pii> pairs;
        vector<char> used(P,0);
        vector<int> ord(P); for(int i=0;i<P;i++) ord[i]=i;
        shuffle(ord.begin(), ord.end(), rng);
        for(int oi=0; oi<P; oi++){
            int i=ord[oi];
            if(used[i]) continue;
            int bj=-1, bd=1<<29, ties=0;
            for(int j=0;j<P;j++) if(j!=i && !used[j]){
                int d=dist[odd[i]][odd[j]];
                if(d<0) continue;
                if(d<bd){ bd=d; bj=j; ties=1; }
                else if(d==bd){ ties++; if((int)(rng()%ties)==0) bj=j; }
            }
            if(bj<0) return {};
            used[i]=used[bj]=1;
            pairs.push_back({odd[i],odd[bj]});
        }
        bool improved=true; int guard2=0;
        while(improved && ++guard2<60){
            improved=false;
            for(int i=0;i<(int)pairs.size();i++) for(int j=i+1;j<(int)pairs.size();j++){
                int a=pairs[i].first,b=pairs[i].second,c=pairs[j].first,d=pairs[j].second;
                int cur=dist[a][b]+dist[c][d];
                int alt1=(dist[a][c]<0||dist[b][d]<0)? (1<<29) : dist[a][c]+dist[b][d];
                int alt2=(dist[a][d]<0||dist[b][c]<0)? (1<<29) : dist[a][d]+dist[b][c];
                if(alt1<cur && alt1<=alt2){ pairs[i]={a,c}; pairs[j]={b,d}; improved=true; }
                else if(alt2<cur){ pairs[i]={a,d}; pairs[j]={b,c}; improved=true; }
            }
        }
        long long tot=0; int mx=0;
        for(auto& pr: pairs){ int d=dist[pr.first][pr.second]; tot+=d; mx=max(mx,d); }
        long long cost = tot - mx;
        if(cost < bestCost){ bestCost=cost; bestPairs=pairs; }
    }
    vector<pii> pairs = bestPairs;
    int su=-1, sv=-1; // trail endpoints
    if(!pairs.empty()){
        int worst=0;
        for(int i=1;i<(int)pairs.size();i++)
            if(dist[pairs[i].first][pairs[i].second] > dist[pairs[worst].first][pairs[worst].second]) worst=i;
        su=pairs[worst].first; sv=pairs[worst].second;
        pairs.erase(pairs.begin()+worst);
    }
    for(auto& pr: pairs){
        int u=pr.first, v=pr.second;
        if(dist[u][v]<0) return {}; // unreachable pair: give up safely
        int x=v, guard3=0;
        while(x!=u){ if(++guard3>N+2 || x<1) return {}; int p=par[u][x]; mult[p][x]++; mult[x][p]++; x=p; }
    }
    int start = (su>0)? su : EDGES[0].first;
    // Hierholzer (iterative)
    vector<int> st, walk;
    st.push_back(start);
    long long guard = 0;
    while(!st.empty()){
        if(++guard > 200000) return {};
        int v=st.back();
        int u=-1;
        for(int w=1;w<=N;w++) if(mult[v][w]>0){ u=w; break; }
        if(u<0){ walk.push_back(v); st.pop_back(); }
        else { mult[v][u]--; mult[u][v]--; st.push_back(u); }
    }
    // if any edge unconsumed (disconnected), fail
    for(int a=1;a<=N;a++) for(int b=1;b<=N;b++) if(mult[a][b]>0) return {};
    return walk;
}

// ---------------- min-conflicts local search for a fixed K ----------------
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
    long long iters=0;
    auto now=[&](){ return chrono::steady_clock::now(); };
    while(true){
        iters++;
        if((iters & 255)==0 && now() >= dl) break;
        if((iters & 511)==0 && forbidden==0 && missing==0 && missingColors==0){
            result.assign(K, vector<int>(K));
            for(int r=0;r<K;r++) for(int c=0;c<K;c++) result[r][c]=g[r*K+c];
            return true;
        }
        int p=-1, forcedX=-1;
        if((int)(rng()%1000) < 20){
            p = rng() % (K*K);
        } else if(forbidden>0){
            for(int tries=0; tries<12 && !badVec.empty(); tries++){ int idx=rng()%badVec.size(); int cand=badVec[idx]; if(inBad[cand]){ p=cand; break; } else { badVec[idx]=badVec.back(); badVec.pop_back(); } }
            if(p<0) p = rng()%(K*K);
        } else if(missing>0){
            int eidx=-1;
            for(int tries=0; tries<20; tries++){ int e=rng()%EDGES.size(); if(cnt[EDGES[e].first][EDGES[e].second]==0){ eidx=e; break; } }
            if(eidx>=0){
                int a=EDGES[eidx].first,b=EDGES[eidx].second;
                int ep = (rng()&1)?a:b; int other=(ep==a)?b:a;
                if(!cellsOf[ep].empty()){
                    int baseCell = cellsOf[ep][rng()%cellsOf[ep].size()];
                    if(g[baseCell]==ep){
                        int r=baseCell/K, c=baseCell%K;
                        int d=rng()%4; int nr=r+DR[d],nc=c+DC[d];
                        if(nr>=0&&nr<K&&nc>=0&&nc<K){
                            p=nr*K+nc;
                            if((int)(rng()%100) < 35) forcedX=other; // force the contact; repairs follow
                        }
                    }
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
        int xlo=1, xhi=N;
        if(forcedX>0 && forcedX!=o){ xlo=xhi=forcedX; }
        for(int x=xlo;x<=xhi;x++){
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
    vector<int> totalC(N+2,0);
    vector<vector<int>> rowC(K, vector<int>(N+2,0)), colC(K, vector<int>(N+2,0));
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){ int v=G[r][c]; totalC[v]++; rowC[r][v]++; colC[c][v]++; }
    auto rowLoss=[&](int r)->long long{
        long long loss=0;
        if(r>0) for(int c=0;c<K;c++) if(G[r-1][c]!=G[r][c]) loss++;
        if(r+1<K) for(int c=0;c<K;c++) if(G[r][c]!=G[r+1][c]) loss++;
        for(int v=1;v<=N;v++) if(rowC[r][v]==totalC[v] && totalC[v]>0) loss+=100000;
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
    vector<vector<int>> out(K-1, vector<int>(K-1));
    int rr=0;
    for(int r=0;r<K;r++){ if(r==bestR) continue; int cc=0; for(int c=0;c<K;c++){ if(c==bestC) continue; out[rr][cc]=G[r][c]; cc++; } rr++; }
    return out;
}

// nearest-neighbor rescale of a grid to K2 x K2
static vector<vector<int>> rescaleTo(const vector<vector<int>>& G, int K2){
    int K=(int)G.size();
    vector<vector<int>> out(K2, vector<int>(K2));
    for(int r=0;r<K2;r++) for(int c=0;c<K2;c++)
        out[r][c]=G[(int)((ll)r*K/K2)][(int)((ll)c*K/K2)];
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
        auto tStart = chrono::steady_clock::now();
        HARD_DL = tStart + chrono::milliseconds(900);
        if(N==1){ printf("1\n1\n1\n"); continue; }
        // ---- gather walks ----
        vector<vector<int>> walks;
        bool largeGraph = (N >= 30);
        double firstHpTime = largeGraph ? 0.04 : 0.10;
        double extraHpTime = largeGraph ? 0.015 : 0.02;
        int hpWindowMs = largeGraph ? 80 : 160;
        int hpAttempts = largeGraph ? 2 : 3;
        if(findHP(firstHpTime)) walks.push_back(bestHP);
        auto hpEnd = tStart + chrono::milliseconds(hpWindowMs);
        for(int attempt=0; attempt<hpAttempts && chrono::steady_clock::now()<hpEnd; attempt++){
            if(findHP(extraHpTime)) walks.push_back(bestHP);
        }
        vector<int> dfsw = buildWalk();
        {   // short covering walks from a few starts, keep the 3 shortest
            vector<vector<int>> gws;
            vector<int> starts;
            int mindeg=1<<29, mdv=1;
            for(int v=1;v<=N;v++){ int d=0; for(int u=1;u<=N;u++) d+=ADJ[v][u]; if(d<mindeg){mindeg=d;mdv=v;} }
            starts.push_back(mdv);
            for(int i=0;i<5;i++) starts.push_back(1+(int)(rng()%N));
            for(int s: starts){ auto w=greedyCoverWalk(s); if(!w.empty()) gws.push_back(w); }
            sort(gws.begin(), gws.end(), [](const vector<int>&a, const vector<int>&b){ return a.size()<b.size(); });
            for(int i=0;i<(int)gws.size() && i<3;i++) walks.push_back(gws[i]);
        }
        // ---- constructive candidates: keep min-K valid ----
        int bestK = 1000000; vector<vector<int>> bestGrid;
        auto consider=[&](vector<vector<int>>& g){
            if(g.empty()) return;
            int K=(int)g.size();
            if(K<bestK && verifyGrid(g)){ bestK=K; bestGrid=g; }
        };
        auto considerWalk=[&](const vector<int>& w){
            if(w.empty()) return;
            { vector<vector<int>> g; if(constructDiagV2(w,g)<=240) consider(g); }
            { vector<int> rw(w.rbegin(), w.rend()); vector<vector<int>> g; if(constructDiagV2(rw,g)<=240) consider(g); }
            { vector<vector<int>> g; if(constructDiag(w,g)<=240) consider(g); }
        };
        for(auto& w: walks) considerWalk(w);
        considerWalk(dfsw);
        { vector<int> pw = postmanWalk(); considerWalk(pw); }
        { vector<vector<int>> g; constructFrom(walks.empty()? dfsw : walks[0], g); consider(g); }
        if(bestGrid.empty()){
            // graph likely infeasible (disconnected); emit well-formed stripes anyway
            vector<vector<int>> g; constructFrom(dfsw, g); bestGrid=g; bestK=(int)g.size();
        }
        // ---- local search: progressively shrink K until deadline ----
        int lb = 2; while(lb*lb < N) lb++;
        { int lb2=2; while(2*lb2*(lb2-1) < M) lb2++; lb=max(lb,lb2); }
        int targetK = bestK-1, variant = 0;
        while(targetK >= max(2,lb)){
            auto now = chrono::steady_clock::now();
            if(now >= HARD_DL) break;
            long long rem = chrono::duration_cast<chrono::milliseconds>(HARD_DL - now).count();
            if(rem < 20) break;
            long long slice = min(min(220LL, 70LL + 50LL*variant), max(30LL, rem/2));
            vector<vector<int>> seed;
            int which = variant % 3;
            if(which==0){
                int off = bestK - targetK;
                int dr = (int)(rng()%(off+1)), dc = (int)(rng()%(off+1));
                seed.assign(targetK, vector<int>(targetK));
                for(int r=0;r<targetK;r++) for(int c=0;c<targetK;c++) seed[r][c]=bestGrid[r+dr][c+dc];
            } else if(which==1){
                seed = shrinkByOne(bestGrid);
                if((int)seed.size()!=targetK) seed = rescaleTo(bestGrid, targetK);
            } else {
                seed = rescaleTo(bestGrid, targetK);
            }
            if(variant>=3){
                // perturb to escape repeated failure
                int KK=targetK;
                for(int i=0;i<KK*KK/24+1;i++) seed[rng()%KK][rng()%KK]=1+rng()%N;
            }
            auto sl = now + chrono::milliseconds(slice);
            if(sl > HARD_DL) sl = HARD_DL;
            vector<vector<int>> res;
            if(localSearch(seed, targetK, sl, res) && verifyGrid(res)){
                bestGrid=res; bestK=targetK; targetK--; variant=0;
            } else variant++;
        }
        if(!verifyGrid(bestGrid) && !bestGrid.empty()){ /* infeasible input: emit best-effort grid */ }
        // ---- output ----
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
