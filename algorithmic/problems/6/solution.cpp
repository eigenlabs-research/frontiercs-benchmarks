// World Map (IOI 2024) — own solution.
// Core idea: lay a Hamiltonian path (each country once) along anti-diagonals so
// the sequence length is ~N and K ~ N/2; realize non-path edges as isolated
// single "blob" cells on a tripled anchor diagonal; then min-conflicts local
// search shrinks K further. Dense graphs use a random greedy fill instead.
// A DFS Euler-tour stripe construction is the guaranteed-valid fallback.
#include <cstdio>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>
#include <cstring>
#include <cmath>
using namespace std;
typedef pair<int,int> pii;

static int N, M;
static bool ADJ[45][45];
static unsigned long long adjmask[45]; // bit u set in adjmask[v] iff v~u or v==u (N<=40 fits one word)
static vector<pii> EDGES;
static mt19937 rng(20240915u);
static chrono::steady_clock::time_point HARD_DL;
static inline bool past(const chrono::steady_clock::time_point& dl){ return chrono::steady_clock::now() >= dl; }

// ---- full checker-equivalent verification -------------------------------
static bool verifyGrid(const vector<vector<int>>& g){
    int K=(int)g.size();
    if(K<1||K>240) return false;
    for(auto& row:g){ if((int)row.size()!=K) return false; for(int v:row) if(v<1||v>N) return false; }
    vector<char> present(N+1,0);
    for(auto& row:g) for(int v:row) present[v]=1;
    for(int c=1;c<=N;c++) if(!present[c]) return false;
    static int cnt[45][45];
    for(int a=0;a<45;a++) for(int b=0;b<45;b++) cnt[a][b]=0;
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
        int v=g[r][c];
        if(c+1<K){ int u=g[r][c+1]; if(v!=u){ if(!ADJ[v][u]) return false; cnt[min(v,u)][max(v,u)]++; } }
        if(r+1<K){ int u=g[r+1][c]; if(v!=u){ if(!ADJ[v][u]) return false; cnt[min(v,u)][max(v,u)]++; } }
    }
    for(auto& e:EDGES) if(cnt[e.first][e.second]==0) return false;
    return true;
}

// ---- Hamiltonian path (Warnsdorff-ordered backtracking, time-limited) ---
static vector<int> hpBest, hpCur; static vector<char> hpIn; static bool hpOk;
static chrono::steady_clock::time_point hpDL;
static bool hpDFS(int v,int d){
    if(d==N){ hpOk=true; hpBest=hpCur; return true; }
    if((d&7)==0 && past(hpDL)) return false;
    int nb[45],c=0;
    for(int u=1;u<=N;u++) if(ADJ[v][u]&&!hpIn[u]) nb[c++]=u;
    sort(nb,nb+c,[&](int a,int b){
        int ca=0,cb=0;
        for(int w=1;w<=N;w++){ if(ADJ[a][w]&&!hpIn[w])ca++; if(ADJ[b][w]&&!hpIn[w])cb++; }
        return ca<cb;
    });
    for(int i=0;i<c;i++){
        int u=nb[i]; hpIn[u]=1; hpCur.push_back(u);
        if(hpDFS(u,d+1)) return true;
        hpCur.pop_back(); hpIn[u]=0;
        if(hpOk||past(hpDL)) return hpOk;
    }
    return false;
}
static bool findHP(double sec){
    hpOk=false; hpBest.clear();
    hpDL=chrono::steady_clock::now()+chrono::milliseconds((int)(sec*1000));
    if(hpDL>HARD_DL) hpDL=HARD_DL;
    hpIn.assign(N+1,0);
    vector<int> starts(N); for(int i=0;i<N;i++) starts[i]=i+1;
    shuffle(starts.begin(),starts.end(),rng);
    for(int s:starts){
        if(hpOk||past(hpDL)) break;
        hpCur.clear(); fill(hpIn.begin(),hpIn.end(),0);
        hpIn[s]=1; hpCur.push_back(s);
        hpDFS(s,1);
    }
    return hpOk;
}

// ---- greedy covering walk (visits every vertex, repeats allowed) --------
static vector<int> coverWalk(int start){
    vector<char> vis(N+1,0);
    vector<int> walk={start}; vis[start]=1; int cur=start,cnt=1;
    while(cnt<N){
        vector<int> par(N+1,0),dist(N+1,-1),q={cur}; dist[cur]=0; int head=0,fd=-1;
        vector<int> found;
        while(head<(int)q.size()){
            int v=q[head++];
            if(!vis[v]){ if(fd<0) fd=dist[v]; if(dist[v]>fd) break; found.push_back(v); continue; }
            if(fd>=0 && dist[v]>=fd) continue;
            for(int u=1;u<=N;u++) if(ADJ[v][u]&&dist[u]<0){ dist[u]=dist[v]+1; par[u]=v; q.push_back(u); }
        }
        if(found.empty()) return {};
        int bt=-1,bs=1<<29,ties=0;
        for(int v:found){ int s=0; for(int u=1;u<=N;u++) if(ADJ[v][u]&&!vis[u]) s++;
            if(s<bs){bs=s;bt=v;ties=1;} else if(s==bs){ties++; if((int)(rng()%ties)==0)bt=v;} }
        vector<int> path; for(int x=bt;x!=cur;x=par[x]) path.push_back(x);
        reverse(path.begin(),path.end());
        for(int v:path){ walk.push_back(v); if(!vis[v]){vis[v]=1;cnt++;} }
        cur=bt;
    }
    return walk;
}

// ---- Chinese-postman covering walk (short walk over ALL edges) ----------
// Every edge lands on the walk, so the diagonal layout has no chords -> small K
// on sparse graphs. Duplicates a min-weight matching of odd-degree vertices,
// then extracts an open Euler trail (Hierholzer).
static vector<int> postmanWalk(){
    if(M==0) return {};
    static int dist[41][41], par[41][41];
    for(int s=1;s<=N;s++){ for(int v=1;v<=N;v++){dist[s][v]=-1;par[s][v]=0;} vector<int> q={s}; dist[s][s]=0; int h=0;
        while(h<(int)q.size()){ int v=q[h++]; for(int u=1;u<=N;u++) if(ADJ[v][u]&&dist[s][u]<0){dist[s][u]=dist[s][v]+1;par[s][u]=v;q.push_back(u);} } }
    static int mult[41][41]; memset(mult,0,sizeof(mult)); int deg[41]; memset(deg,0,sizeof(deg));
    for(auto&e:EDGES){ mult[e.first][e.second]=mult[e.second][e.first]=1; deg[e.first]++; deg[e.second]++; }
    vector<int> odd; for(int v=1;v<=N;v++) if(deg[v]&1) odd.push_back(v);
    int P=(int)odd.size(); vector<pii> bestPairs; long long bestCost=1LL<<60;
    for(int rs=0;rs<10&&P>0;rs++){
        vector<pii> pairs; vector<char> used(P,0); vector<int> ord(P); for(int i=0;i<P;i++)ord[i]=i; shuffle(ord.begin(),ord.end(),rng);
        for(int oi=0;oi<P;oi++){ int i=ord[oi]; if(used[i])continue; int bj=-1,bd=1<<29,ties=0;
            for(int j=0;j<P;j++) if(j!=i&&!used[j]){ int d=dist[odd[i]][odd[j]]; if(d<0)continue; if(d<bd){bd=d;bj=j;ties=1;} else if(d==bd){ties++; if((int)(rng()%ties)==0)bj=j;} }
            if(bj<0) return {}; used[i]=used[bj]=1; pairs.push_back({odd[i],odd[bj]}); }
        bool imp=true; int g2=0;
        while(imp&&++g2<60){ imp=false; for(int i=0;i<(int)pairs.size();i++)for(int j=i+1;j<(int)pairs.size();j++){
            int a=pairs[i].first,b=pairs[i].second,c=pairs[j].first,d=pairs[j].second; int cur=dist[a][b]+dist[c][d];
            int a1=(dist[a][c]<0||dist[b][d]<0)?(1<<29):dist[a][c]+dist[b][d];
            int a2=(dist[a][d]<0||dist[b][c]<0)?(1<<29):dist[a][d]+dist[b][c];
            if(a1<cur&&a1<=a2){pairs[i]={a,c};pairs[j]={b,d};imp=true;} else if(a2<cur){pairs[i]={a,d};pairs[j]={b,c};imp=true;} } }
        long long tot=0; int mx=0; for(auto&pr:pairs){int d=dist[pr.first][pr.second];tot+=d;mx=max(mx,d);}
        long long cost=tot-mx; if(cost<bestCost){bestCost=cost;bestPairs=pairs;}
    }
    vector<pii> pairs=bestPairs; int su=-1;
    if(!pairs.empty()){ int w=0; for(int i=1;i<(int)pairs.size();i++) if(dist[pairs[i].first][pairs[i].second]>dist[pairs[w].first][pairs[w].second])w=i;
        su=pairs[w].first; pairs.erase(pairs.begin()+w); }
    for(auto&pr:pairs){ int u=pr.first,v=pr.second; if(dist[u][v]<0)return {}; int x=v,g3=0;
        while(x!=u){ if(++g3>N+2||x<1)return {}; int p=par[u][x]; mult[p][x]++; mult[x][p]++; x=p; } }
    int start=(su>0)?su:EDGES[0].first; vector<int> st={start},walk; long long g=0;
    while(!st.empty()){ if(++g>200000)return {}; int v=st.back(),u=-1; for(int w=1;w<=N;w++) if(mult[v][w]>0){u=w;break;}
        if(u<0){walk.push_back(v);st.pop_back();} else {mult[v][u]--;mult[u][v]--;st.push_back(u);} }
    for(int a=1;a<=N;a++)for(int b=1;b<=N;b++) if(mult[a][b]>0) return {};
    return walk;
}

// ---- DFS Euler-tour stripe (guaranteed-valid fallback) ------------------
static vector<int> dfsTour(){
    vector<int> par(N+1,0); vector<char> vis(N+1,0); vector<int> q={1}; vis[1]=1; int head=0;
    while(head<(int)q.size()){ int v=q[head++]; for(int u=1;u<=N;u++) if(ADJ[v][u]&&!vis[u]){vis[u]=1;par[u]=v;q.push_back(u);} }
    vector<vector<int>> ch(N+1); for(int v=1;v<=N;v++) if(par[v]) ch[par[v]].push_back(v);
    vector<int> w={1}; vector<pii> st={{1,0}};
    while(!st.empty()){ auto&b=st.back(); int v=b.first,&ci=b.second;
        if(ci<(int)ch[v].size()){ int c=ch[v][ci++]; w.push_back(c); st.push_back({c,0}); }
        else { st.pop_back(); if(!st.empty()) w.push_back(st.back().first); } }
    return w;
}

// randomized greedy vertex cover of a chord set
static vector<char> vertexCover(const vector<pii>& ch){
    int C=(int)ch.size(); vector<char> best(N+1,0); int bestSz=1<<29;
    vector<vector<int>> nbr(N+1);
    for(auto&c:ch){ nbr[c.first].push_back(c.second); nbr[c.second].push_back(c.first); }
    for(int pass=0;pass<30;pass++){
        vector<char> inS(N+1,0); vector<char> cov(C,0); int rem=C;
        while(rem>0){
            vector<int> d(N+1,0);
            for(int i=0;i<C;i++) if(!cov[i]){ d[ch[i].first]++; d[ch[i].second]++; }
            int bd=0; for(int v=1;v<=N;v++) bd=max(bd,d[v]); if(bd==0) break;
            vector<int> cand; for(int v=1;v<=N;v++) if(d[v]==bd) cand.push_back(v);
            int bv=cand[rng()%cand.size()]; inS[bv]=1;
            for(int i=0;i<C;i++) if(!cov[i]&&(ch[i].first==bv||ch[i].second==bv)){cov[i]=1;rem--;}
        }
        bool ch2=true; while(ch2){ ch2=false; for(int v=1;v<=N;v++) if(inS[v]){ bool all=true; for(int u:nbr[v]) if(!inS[u]){all=false;break;} if(all){inS[v]=0;ch2=true;} } }
        int sz=0; for(int v=1;v<=N;v++) sz+=inS[v];
        if(sz<bestSz){bestSz=sz;best=inS;}
    }
    return best;
}

// ---- anti-diagonal layout with tripled-anchor blobs ---------------------
// walk: vertex sequence, consecutive entries must be adjacent (or equal).
static int constructDiag(const vector<int>& walk, vector<vector<int>>& grid){
    int W=(int)walk.size(); if(W==0) return 1<<29;
    for(int i=0;i+1<W;i++) if(walk[i]!=walk[i+1] && !ADJ[walk[i]][walk[i+1]]) return 1<<29;
    set<pii> wE; for(int i=0;i+1<W;i++){ int a=walk[i],b=walk[i+1]; if(a!=b) wE.insert({min(a,b),max(a,b)}); }
    vector<pii> chords; for(auto&e:EDGES) if(!wE.count(e)) chords.push_back(e);
    vector<char> inS(N+1,0); if(!chords.empty()) inS=vertexCover(chords);
    vector<int> blobCnt(N+1,0); vector<vector<int>> blobs(N+1);
    for(auto&c:chords){ int a=c.first,b=c.second,chosen;
        if(inS[a]&&inS[b]) chosen=(blobCnt[a]<=blobCnt[b])?a:b; else if(inS[a]) chosen=a; else chosen=b;
        int other=(chosen==a)?b:a; blobs[chosen].push_back(other); blobCnt[chosen]++; }
    auto build=[&](bool rev, vector<int>& seq, vector<int>& mid)->int{
        vector<int> w=walk; if(rev) reverse(w.begin(),w.end());
        seq.clear(); mid.assign(N+1,-1); vector<char> done(N+1,0);
        for(int v:w){
            if(blobCnt[v]>0 && !done[v]){ done[v]=1; seq.push_back(v); seq.push_back(v); seq.push_back(v); mid[v]=(int)seq.size()-2; }
            else if(!(!seq.empty() && seq.back()==v)) seq.push_back(v);   // collapse equal repeats
            else seq.push_back(v);
        }
        int pad=0; for(int v=1;v<=N;v++) if(blobCnt[v]>0&&mid[v]>=0) pad=max(pad,blobCnt[v]-1-mid[v]);
        return pad;
    };
    vector<int> seqA,midA,seqB,midB;
    int pA=build(false,seqA,midA), pB=build(true,seqB,midB);
    vector<int>& seq=(pB<pA)?seqB:seqA; vector<int>& mid=(pB<pA)?midB:midA; int pad=min(pA,pB);
    if(pad>0){ vector<int> ns((size_t)pad,seq[0]); ns.insert(ns.end(),seq.begin(),seq.end()); seq=ns;
        for(int v=1;v<=N;v++) if(mid[v]>=0) mid[v]+=pad; }
    int L=(int)seq.size(); int K=(L+2)/2; if(K<1) K=1;
    auto capOK=[&](int k)->bool{ for(int v=1;v<=N;v++) if(blobCnt[v]>0){ int t=mid[v]; int len=min(min(t+1,k),2*k-1-t); if(len<blobCnt[v]) return false; } return true; };
    while(K<=240 && !capOK(K)) K++;
    if(K>240) return 1<<29;
    grid.assign(K,vector<int>(K));
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){ int t=r+c; grid[r][c]=seq[min(t,L-1)]; }
    for(int v=1;v<=N;v++) if(blobCnt[v]>0){ int t=mid[v]; int i0=max(0,t-K+1),i1=min(t,K-1),idx=0;
        for(int i=i0;i<=i1&&idx<(int)blobs[v].size();i++) grid[i][t-i]=blobs[v][idx++];
        if(idx<(int)blobs[v].size()) return 1<<29; }
    return K;
}

// ---- generalized anti-diagonal layout -----------------------------------
// Base band: cell (r,c)=seq[r+c]. A chord (u,w) is realized by overriding one
// cell on some diagonal t with the blob color, when the surrounding bands
// (seq[t-1], seq[t+1], seq[t]) are all adjacent-or-equal to the blob and one
// neighbor band equals the chord partner. No vertex tripling, so L stays ~N.
static int constructDiagGen(const vector<int>& walk, vector<vector<int>>& grid){
    int W=(int)walk.size(); if(W==0) return 1<<29;
    for(int i=0;i+1<W;i++) if(walk[i]!=walk[i+1] && !ADJ[walk[i]][walk[i+1]]) return 1<<29;
    // collapse consecutive duplicates in the walk to form the band sequence
    vector<int> seq; for(int v:walk){ if(seq.empty()||seq.back()!=v) seq.push_back(v); }
    set<pii> wE; for(int i=0;i+1<(int)seq.size();i++){ int a=seq[i],b=seq[i+1]; wE.insert({min(a,b),max(a,b)}); }
    vector<pii> chords; for(auto&e:EDGES) if(!wE.count(e)) chords.push_back(e);
    int C=(int)chords.size();

    vector<char> tripled(N+1,0); int pad=0;
    for(int round=0; round<12; round++){
        // build band sequence with front pad and any tripled vertices
        vector<int> s((size_t)pad, seq.empty()?1:seq[0]); vector<char> done(N+1,0);
        for(int v:seq){ if(tripled[v]&&!done[v]){ done[v]=1; s.push_back(v); s.push_back(v); s.push_back(v);} else s.push_back(v); }
        int L=(int)s.size(); if(L>470) return 1<<29;
        vector<vector<int>> posOf(N+1); for(int i=0;i<L;i++) posOf[s[i]].push_back(i);
        // candidate diagonals (t, blobColor) per chord
        vector<vector<pii>> cand(C);
        for(int ci=0;ci<C;ci++){ set<pii> cs;
            for(int dir=0;dir<2;dir++){
                int anchor=dir?chords[ci].second:chords[ci].first, blob=dir?chords[ci].first:chords[ci].second;
                for(int idx:posOf[anchor]) for(int dt=-1;dt<=1;dt+=2){
                    int t=idx+dt; if(t<1||t>L-2) continue;
                    int s1=s[t-1],s2=s[t+1];
                    bool ok1=(s1==blob)||ADJ[min(s1,blob)][max(s1,blob)];
                    bool ok2=(s2==blob)||ADJ[min(s2,blob)][max(s2,blob)];
                    if(ok1&&ok2) cs.insert({t,blob});
                }
            }
            cand[ci].assign(cs.begin(),cs.end());
        }
        bool anyEmpty=false; for(int i=0;i<C;i++) if(cand[i].empty()) anyEmpty=true;
        vector<int> order(C); for(int i=0;i<C;i++)order[i]=i;
        sort(order.begin(),order.end(),[&](int a,int b){return cand[a].size()<cand[b].size();});
        if(!anyEmpty){
            int K0=(L+2)/2;
            for(int K=K0;K<=min(240,K0+10);K++){
                auto I0=[&](int t){return max(0,t-K+1);}; auto I1=[&](int t){return min(t,K-1);};
                auto dlen=[&](int t){return I1(t)-I0(t)+1;};
                vector<int> load(L,0);
                // margin left on diag t if one more blob added: keep >=1 pure cell and
                // keep a pure contact across each diagonal boundary
                auto margin=[&](int t)->int{ int ld=load[t]+1; int m=(dlen(t)-1)-ld;
                    for(int dt=-1;dt<=1;dt+=2){ int u=t+dt; if(u<0||u>=L)continue;
                        int e=(t%2==0)?t:u, o=(t%2==0)?u:t; int le=(e==t)?ld:load[e], lo=(o==t)?ld:load[o];
                        int lim=(e<o)?I1(o)-I0(e)-1:I1(o)-I0(e); m=min(m,lim-(le+lo)); }
                    return m; };
                vector<pii> assign(C,{-1,-1}); bool ok=true;
                for(int oi=0;oi<C&&ok;oi++){ int ci=order[oi]; int bt=-1,bb=-1,bs=-1;
                    for(auto&tc:cand[ci]){ int sl=margin(tc.first); if(sl>bs){bs=sl;bt=tc.first;bb=tc.second;} }
                    if(bt<0||bs<0){ok=false;break;} load[bt]++; assign[ci]={bt,bb}; }
                if(!ok) continue;
                grid.assign(K,vector<int>(K));
                for(int r=0;r<K;r++)for(int c=0;c<K;c++){int t=r+c; grid[r][c]=s[min(t,L-1)];}
                vector<int> used(L,0);
                for(int ci=0;ci<C;ci++){ int t=assign[ci].first,blob=assign[ci].second;
                    int i=(t%2==0)?I0(t)+used[t]:I1(t)-used[t]; used[t]++; grid[i][t-i]=blob; }
                return K;
            }
        }
        // repair: triple a vertex cover of anchorless chords, or pad for capacity
        if(anyEmpty){
            vector<pii> bad; for(int i=0;i<C;i++) if(cand[i].empty()) bad.push_back(chords[i]);
            vector<int> d(N+1,0); for(auto&b:bad){d[b.first]++;d[b.second]++;} vector<char> chosen(N+1,0);
            while(true){ int bv=-1,bd=0; for(int v=1;v<=N;v++) if(!tripled[v]&&!chosen[v]&&d[v]>bd){bd=d[v];bv=v;}
                if(bv<0)break; chosen[bv]=1; vector<int> d2(N+1,0); bool any=false;
                for(auto&b:bad){ if(chosen[b.first]||chosen[b.second])continue; d2[b.first]++;d2[b.second]++;any=true; } d=d2; if(!any)break; }
            bool prog=false; for(int v=1;v<=N;v++) if(chosen[v]){tripled[v]=1;prog=true;}
            if(!prog) return 1<<29;
        } else pad+=2;
    }
    return 1<<29;
}

// ---- dense random greedy fill (row-major) -------------------------------
static bool denseFill(int K, int attempts, chrono::steady_clock::time_point dl, vector<vector<int>>& out){
    out.clear(); if(K*K<N || 2*K*(K-1)<M) return false;
    int deg[45]; memset(deg,0,sizeof(deg)); for(auto&e:EDGES){deg[e.first]++;deg[e.second]++;}
    for(int at=0;at<attempts;at++){
        if((at&3)==0 && past(dl)) break;   // check often: each attempt fills a full KxK grid
        vector<int> g(K*K,0),cc(N+1,0); static bool seen[45][45]; memset(seen,0,sizeof(seen));
        int miss=N,ec=0; bool ok=true;
        for(int p=0;p<K*K&&ok;p++){
            int r=p/K,c=p%K, up=r?g[(r-1)*K+c]:0, lf=c?g[p-1]:0, rem=K*K-p;
            int cand[45],wt[45],ncc=0,tot=0;
            for(int x=1;x<=N;x++){
                if(up&&up!=x&&!ADJ[up][x]) continue;
                if(lf&&lf!=x&&!ADJ[lf][x]) continue;
                if(miss>=rem && cc[x]>0) continue;
                int add=0; if(up&&up!=x&&!seen[min(up,x)][max(up,x)])add++; if(lf&&lf!=x&&!seen[min(lf,x)][max(lf,x)])add++;
                int w=1+80*add+deg[x]; if(cc[x]==0) w+=200; cand[ncc]=x; wt[ncc]=w; tot+=w; ncc++;
            }
            if(ncc==0){ok=false;break;}
            int pick=(int)(rng()%tot),x=cand[0];
            for(int i=0;i<ncc;i++){ if(pick<wt[i]){x=cand[i];break;} pick-=wt[i]; }
            g[p]=x; if(cc[x]++==0) miss--;
            if(up&&up!=x){int a=min(up,x),b=max(up,x); if(!seen[a][b]){seen[a][b]=1;ec++;}}
            if(lf&&lf!=x){int a=min(lf,x),b=max(lf,x); if(!seen[a][b]){seen[a][b]=1;ec++;}}
        }
        if(ok&&miss==0&&ec==M){ out.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)out[r][c]=g[r*K+c]; return true; }
    }
    return false;
}

// Like denseFill but returns the BEST fully-filled grid found (all colors
// present, maximal edges realized) even if not perfect — a strong local-search
// seed for K just below the perfect-fill floor.
static bool denseFillBest(int K, int attempts, chrono::steady_clock::time_point dl, vector<vector<int>>& out){
    out.clear(); if(K*K<N || 2*K*(K-1)<M) return false;
    int deg[45]; memset(deg,0,sizeof(deg)); for(auto&e:EDGES){deg[e.first]++;deg[e.second]++;}
    vector<int> bestG; int bestMiss=1<<29,bestEc=-1;
    for(int at=0;at<attempts;at++){
        if((at&3)==0 && past(dl)) break;
        vector<int> g(K*K,0),cc(N+1,0); static bool seen[45][45]; memset(seen,0,sizeof(seen));
        int miss=N,ec=0; bool ok=true;
        for(int p=0;p<K*K&&ok;p++){
            int r=p/K,c=p%K, up=r?g[(r-1)*K+c]:0, lf=c?g[p-1]:0, rem=K*K-p;
            int cand[45],wt[45],ncc=0,tot=0;
            for(int x=1;x<=N;x++){
                if(up&&up!=x&&!ADJ[up][x]) continue;
                if(lf&&lf!=x&&!ADJ[lf][x]) continue;
                if(miss>=rem && cc[x]>0) continue;
                int add=0; if(up&&up!=x&&!seen[min(up,x)][max(up,x)])add++; if(lf&&lf!=x&&!seen[min(lf,x)][max(lf,x)])add++;
                int w=1+80*add+deg[x]; if(cc[x]==0) w+=200; cand[ncc]=x; wt[ncc]=w; tot+=w; ncc++;
            }
            if(ncc==0){ok=false;break;}
            int pick=(int)(rng()%tot),x=cand[0];
            for(int i=0;i<ncc;i++){ if(pick<wt[i]){x=cand[i];break;} pick-=wt[i]; }
            g[p]=x; if(cc[x]++==0) miss--;
            if(up&&up!=x){int a=min(up,x),b=max(up,x); if(!seen[a][b]){seen[a][b]=1;ec++;}}
            if(lf&&lf!=x){int a=min(lf,x),b=max(lf,x); if(!seen[a][b]){seen[a][b]=1;ec++;}}
        }
        if(!ok || miss>0) continue;              // require full fill + all colors present
        int missE=M-ec;
        if(missE<bestMiss){ bestMiss=missE; bestEc=ec; bestG=g; if(missE==0) break; }
    }
    if(bestG.empty()) return false;
    out.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)out[r][c]=bestG[r*K+c];
    return bestMiss==0;
}

// ---- grid analysis + targeted missing-edge repair -----------------------
struct Stats{ int forb=0, missE=0, missC=0; vector<int> missIdx; };
static Stats analyze(const vector<vector<int>>& g){
    Stats s; int K=(int)g.size(); static int cnt[45][45];
    for(int a=0;a<45;a++)for(int b=0;b<45;b++)cnt[a][b]=0;
    vector<int> cc(N+1,0);
    for(int r=0;r<K;r++)for(int c=0;c<K;c++){ int v=g[r][c]; if(v>=1&&v<=N)cc[v]++;
        if(c+1<K){int u=g[r][c+1]; if(v!=u){int a=min(v,u),b=max(v,u); if(a>=1&&b<=N&&ADJ[a][b])cnt[a][b]++; else s.forb++;}}
        if(r+1<K){int u=g[r+1][c]; if(v!=u){int a=min(v,u),b=max(v,u); if(a>=1&&b<=N&&ADJ[a][b])cnt[a][b]++; else s.forb++;}} }
    for(int c=1;c<=N;c++) if(cc[c]==0)s.missC++;
    for(int i=0;i<M;i++){ auto e=EDGES[i]; if(cnt[e.first][e.second]==0){s.missE++;s.missIdx.push_back(i);} }
    return s;
}
static long long repScore(const Stats& s){ return 1000000000LL*s.forb+1000000LL*s.missC+s.missE; }
// only valid when missing colors already 0; greedily fixes missing edges by
// forcing an endpoint pair onto an adjacent cell pair that improves the score.
static bool patchRepair(vector<vector<int>> g, int K, chrono::steady_clock::time_point dl, vector<vector<int>>& res){
    if((int)g.size()!=K||K<=0||(int)g[0].size()!=K) return false;
    Stats cur=analyze(g);
    if(cur.forb==0&&cur.missC==0&&cur.missE==0){ res=g; return true; }
    if(cur.missC>0) return false;
    const int DR[2]={0,1},DC[2]={1,0}; int rounds=0;
    while(cur.missE>0 && !past(dl) && rounds++<25){
        long long base=repScore(cur), bestSc=base; int br1=-1,bc1=-1,br2=-1,bc2=-1,bx=0,by=0; Stats bst=cur; bool stop=false;
        for(int eidx:cur.missIdx){ if(stop||past(dl))break; int a=EDGES[eidx].first,b=EDGES[eidx].second;
            for(int ori=0;ori<2&&!stop;ori++){ int x=ori?a:b,y=ori?b:a;
                for(int r=0;r<K&&!stop;r++){ if(past(dl)){stop=true;break;} for(int c=0;c<K&&!stop;c++)for(int d=0;d<2;d++){
                    int nr=r+DR[d],nc=c+DC[d]; if(nr>=K||nc>=K)continue;
                    int ox=g[r][c],oy=g[nr][nc]; if(ox==x&&oy==y)continue;
                    g[r][c]=x; g[nr][nc]=y; Stats s=analyze(g); long long sc=repScore(s);
                    if(sc<bestSc){bestSc=sc;bst=s;br1=r;bc1=c;br2=nr;bc2=nc;bx=x;by=y;}
                    if(sc<base) stop=true;
                    g[r][c]=ox; g[nr][nc]=oy; if(stop)break;
                } } } }
        if(bestSc>=base||br1<0) break;
        g[br1][bc1]=bx; g[br2][bc2]=by; cur=bst;
    }
    if(cur.forb==0&&cur.missC==0&&cur.missE==0&&verifyGrid(g)){ res=g; return true; }
    return false;
}

// Chinese-Postman SNAKE seed: fold an edge-covering walk into a boustrophedon
// path in a KxK grid. Along-snake adjacencies ARE the walk's steps (real edges),
// so every edge on the walk is realized -> miss=0 by construction. Only vertical
// between-row pairs may be forbidden; SA then only repairs those (never hunts for
// missing edges). Wins on sparse/mid graphs where the walk is short (K~sqrt(len)).
// Returns empty if the walk doesn't fit (len > K*K).
static vector<vector<int>> snakeSeed(const vector<int>& walk, int K){
    if(walk.empty() || (int)walk.size() > K*K) return {};
    vector<int> w = walk; while((int)w.size() < K*K) w.push_back(w.back());
    vector<vector<int>> g(K, vector<int>(K)); int oi=0;
    for(int r=0;r<K;r++) for(int cc=0;cc<K;cc++){ int c=(r&1)?K-1-cc:cc; g[r][c]=w[oi++]; }
    return g;
}

// A fresh random-greedy fill (row-major, prefers legal colors w.r.t. up/left
// neighbors, allows an illegal pick when stuck). Not a valid grid by itself —
// it's a good, diverse SA starting point (the random restart that matters).
static vector<vector<int>> greedySeed(int K){
    vector<int> g(K*K);
    for(int p=0;p<K*K;p++){ int r=p/K,c=p%K,up=r?g[(r-1)*K+c]:0,lf=c?g[p-1]:0;
        unsigned long long ok=~0ULL; if(up) ok&=adjmask[up]; if(lf) ok&=adjmask[lf];
        int cand[45],ncc=0; for(int x=1;x<=N;x++) if(ok&(1ULL<<x)) cand[ncc++]=x;
        g[p]= ncc? cand[rng()%ncc] : (1+rng()%N);
    }
    vector<vector<int>> gr(K,vector<int>(K));
    for(int r=0;r<K;r++)for(int c=0;c<K;c++) gr[r][c]=g[r*K+c];
    return gr;
}

// ---- simulated-annealing repair at fixed K ------------------------------
// Accepts worsening moves with prob exp(-delta/T) so it escapes the local
// minima that trap greedy min-conflicts. Incremental delta + a cooling schedule
// with periodic reheat; ~10M steps in ~0.5s for N<=40. Seed may be any grid.
static bool saRepair(vector<vector<int>> seed, int K, double budgetSec,
                     chrono::steady_clock::time_point dl, vector<vector<int>>& res){
    if((int)seed.size()<K||(int)seed[0].size()<K) return false;
    auto t0=chrono::steady_clock::now();
    vector<int> g(K*K); for(int r=0;r<K;r++)for(int c=0;c<K;c++)g[r*K+c]=seed[r][c];
    static int cnt[45][45]; for(int a=0;a<45;a++)for(int b=0;b<45;b++)cnt[a][b]=0;
    vector<int> cc(N+1,0); int forb=0,miss=0,missC=0;
    const int DR[4]={1,-1,0,0},DC[4]={0,0,1,-1};
    for(int r=0;r<K;r++)for(int c=0;c<K;c++){int v=g[r*K+c];cc[v]++;
        if(c+1<K){int u=g[r*K+c+1];if(v!=u){int a=min(v,u),b=max(v,u);cnt[a][b]++;if(!ADJ[a][b])forb++;}}
        if(r+1<K){int u=g[(r+1)*K+c];if(v!=u){int a=min(v,u),b=max(v,u);cnt[a][b]++;if(!ADJ[a][b])forb++;}}}
    for(auto&e:EDGES) if(cnt[e.first][e.second]==0) miss++;
    for(int v=1;v<=N;v++) if(cc[v]==0) missC++;
    auto solved=[&](){ return forb==0&&miss==0&&missC==0; };
    if(solved()){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
    double Tmp=3.0; long long it=0;
    while(true){
        it++;
        if((it&1023)==0){
            if(past(dl)) break;
            double el=chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-t0).count()/1000.0;
            if(el>budgetSec) break;
            Tmp=3.0*exp(-3.0*el/budgetSec); if(((it>>10)&63)==0) Tmp=max(Tmp,1.5); // periodic reheat
        }
        // ~18% forced-edge move: directly create a currently-missing adjacency
        // (single random recolor reaches it only by luck). Same delta/commit path.
        int p, x;
        if(miss>0 && (int)(rng()%100) < 18){
            int e=-1; for(int t=0;t<24;t++){int ee=rng()%(int)EDGES.size(); if(cnt[EDGES[ee].first][EDGES[ee].second]==0){e=ee;break;}}
            if(e<0){ p=rng()%(K*K); x=1+rng()%N; }
            else{ int u=EDGES[e].first,w=EDGES[e].second; if(rng()&1) swap(u,w);
                int nbp=-1; for(int t=0;t<32;t++){ int q=rng()%(K*K); if(g[q]!=u)continue; int qr=q/K,qc=q%K,d=rng()%4,nr=qr+DR[d],nc=qc+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; nbp=nr*K+nc; break; }
                if(nbp<0){ p=rng()%(K*K); x=1+rng()%N; } else { p=nbp; x=w; } }
        } else { p=rng()%(K*K); x=1+rng()%N; }
        int r=p/K,c=p%K, o=g[p];
        if(x==o) continue;
        long long dF=0; int tch[8],de[8],nt=0;
        for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
            if(o!=u){int a=min(o,u),b=max(o,u); if(ADJ[a][b]){int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break; if(j==nt){tch[nt]=e;de[nt]=0;nt++;} de[j]-=1;} else dF--; }
            if(x!=u){int a=min(x,u),b=max(x,u); if(ADJ[a][b]){int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break; if(j==nt){tch[nt]=e;de[nt]=0;nt++;} de[j]+=1;} else dF++; } }
        long long dM=0; for(int j=0;j<nt;j++){int e=tch[j],a=e/45,b=e%45,oc=cnt[a][b],ncnt=oc+de[j]; if(oc>0&&ncnt==0)dM++; else if(oc==0&&ncnt>0)dM--;}
        long long dMC=0; if(cc[o]==1)dMC++; if(cc[x]==0)dMC--;
        long long delta=1000LL*dF+50LL*dM+50LL*dMC;
        bool acc = delta<=0 || (exp(-double(delta)/(Tmp*50.0)) > (double)(rng()&0xffff)/65536.0);
        if(acc){
            for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
                if(o!=u){int a=min(o,u),b=max(o,u),oc=cnt[a][b];cnt[a][b]=oc-1;if(ADJ[a][b]){if(oc==1)miss++;}else forb--;}
                if(x!=u){int a=min(x,u),b=max(x,u),oc=cnt[a][b];cnt[a][b]=oc+1;if(ADJ[a][b]){if(oc==0)miss--;}else forb++;}}
            cc[o]--; if(cc[o]==0)missC++; cc[x]++; if(cc[x]==1)missC--;
            g[p]=x;
            if(solved()){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
        }
    }
    return false;
}

// ---- min-conflicts local search at fixed K ------------------------------
static bool localSearch(vector<vector<int>> seed, int K, chrono::steady_clock::time_point dl, vector<vector<int>>& res){
    if((int)seed.size()<K||(int)seed[0].size()<K) return false;
    vector<int> g(K*K); for(int r=0;r<K;r++)for(int c=0;c<K;c++)g[r*K+c]=seed[r][c];
    static int cnt[45][45]; for(int a=0;a<45;a++)for(int b=0;b<45;b++)cnt[a][b]=0;
    int forb=0,miss=0,missC=0; vector<int> cc(N+2,0);
    for(int r=0;r<K;r++)for(int c=0;c<K;c++){ int v=g[r*K+c]; cc[v]++;
        if(c+1<K){int u=g[r*K+c+1]; if(v!=u){int a=min(v,u),b=max(v,u); cnt[a][b]++; if(!ADJ[a][b])forb++;}}
        if(r+1<K){int u=g[(r+1)*K+c]; if(v!=u){int a=min(v,u),b=max(v,u); cnt[a][b]++; if(!ADJ[a][b])forb++;}} }
    for(auto&e:EDGES) if(cnt[e.first][e.second]==0) miss++;
    for(int c=1;c<=N;c++) if(cc[c]==0) missC++;
    auto done=[&](){ return forb==0&&miss==0&&missC==0; };
    if(done()){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
    const int DR[4]={1,-1,0,0},DC[4]={0,0,1,-1};
    // cellsOf[c] = positions currently colored c, kept up to date incrementally
    // (swap-remove via posInCell index) so we never pay for a full O(K^2) rebuild.
    vector<vector<int>> cellsOf(N+2);
    vector<int> posInCell(K*K,0);
    for(int p=0;p<K*K;p++){ int c=g[p]; posInCell[p]=(int)cellsOf[c].size(); cellsOf[c].push_back(p); }
    auto moveCell=[&](int p,int from,int to){
        auto& vf=cellsOf[from]; int idx=posInCell[p], last=(int)vf.size()-1;
        if(idx!=last){ vf[idx]=vf[last]; posInCell[vf[idx]]=idx; } vf.pop_back();
        posInCell[p]=(int)cellsOf[to].size(); cellsOf[to].push_back(p);
    };
    long long it=0;
    while(true){
        it++;
        if((it&255)==0 && past(dl)) break;
        if((it&511)==0 && done()){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
        int p=-1,forced=-1;
        if((int)(rng()%1000)<20) p=rng()%(K*K);
        else if(forb>0){
            // random forbidden-neighbor cell
            for(int t=0;t<20;t++){ int q=rng()%(K*K),r=q/K,c=q%K,o=g[q]; bool bad=false;
                for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc]; if(o!=u&&!ADJ[min(o,u)][max(o,u)]){bad=true;break;}}
                if(bad){p=q;break;} }
            if(p<0) p=rng()%(K*K);
        } else if(miss>0){
            int e=-1; for(int t=0;t<20;t++){ int ee=rng()%EDGES.size(); if(cnt[EDGES[ee].first][EDGES[ee].second]==0){e=ee;break;} }
            if(e>=0){ int a=EDGES[e].first,b=EDGES[e].second, ep=(rng()&1)?a:b, other=(ep==a)?b:a;
                if(!cellsOf[ep].empty()){ int bc=cellsOf[ep][rng()%cellsOf[ep].size()]; if(g[bc]==ep){ int r=bc/K,c=bc%K,d=rng()%4,nr=r+DR[d],nc=c+DC[d];
                    if(nr>=0&&nr<K&&nc>=0&&nc<K){ p=nr*K+nc; if((int)(rng()%100)<35) forced=other; } } } }
            if(p<0) p=rng()%(K*K);
        } else p=rng()%(K*K);
        if(p<0) continue;
        int r=p/K,c=p%K,o=g[p], nb[4],nn=0;
        for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; nb[nn++]=g[nr*K+nc];}
        long long bestD=0; int bestX=o,ties=1; bool first=true; int xlo=1,xhi=N;
        if(forced>0&&forced!=o){xlo=xhi=forced;}
        for(int x=xlo;x<=xhi;x++){ if(x==o)continue;
            int df=0; int tch[8],de[8],nt=0;
            for(int i=0;i<nn;i++){ int nc=nb[i];
                if(o!=nc){ int a=min(o,nc),b=max(o,nc); if(ADJ[a][b]){int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break; if(j==nt){tch[nt]=e;de[nt]=0;nt++;} de[j]-=1;} else df--; }
                if(x!=nc){ int a=min(x,nc),b=max(x,nc); if(ADJ[a][b]){int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break; if(j==nt){tch[nt]=e;de[nt]=0;nt++;} de[j]+=1;} else df++; } }
            int dm=0; for(int j=0;j<nt;j++){int e=tch[j],a=e/45,b=e%45,oc=cnt[a][b],ncnt=oc+de[j]; if(oc>0&&ncnt==0)dm++; else if(oc==0&&ncnt>0)dm--;}
            int dmc=0; if(cc[o]==1)dmc++; if(cc[x]==0)dmc--;
            long long delta=1000LL*df+50LL*dm+50LL*dmc;
            if(first){bestD=delta;bestX=x;ties=1;first=false;} else if(delta<bestD){bestD=delta;bestX=x;ties=1;} else if(delta==bestD){ties++; if((int)(rng()%ties)==0)bestX=x;}
        }
        int x=bestX;
        if(!first&&x!=o){
            for(int i=0;i<nn;i++){ int nc=nb[i];
                if(o!=nc){int a=min(o,nc),b=max(o,nc),oc=cnt[a][b]; cnt[a][b]=oc-1; if(ADJ[a][b]){if(oc==1)miss++;} else forb--;}
                if(x!=nc){int a=min(x,nc),b=max(x,nc),oc=cnt[a][b]; cnt[a][b]=oc+1; if(ADJ[a][b]){if(oc==0)miss--;} else forb++;} }
            cc[o]--; if(cc[o]==0)missC++; cc[x]++; if(cc[x]==1)missC--;
            moveCell(p,o,x); g[p]=x;
        }
    }
    if(done()){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
    return false;
}

// crop / rescale helpers for reseeding a smaller K
static vector<vector<int>> rescale(const vector<vector<int>>& G,int K2){
    int K=(int)G.size(); vector<vector<int>> o(K2,vector<int>(K2));
    for(int r=0;r<K2;r++)for(int c=0;c<K2;c++) o[r][c]=G[(int)((long long)r*K/K2)][(int)((long long)c*K/K2)];
    return o;
}
// remove the least-valuable row and column (fewest lost adjacencies, and never
// one that holds the only cells of some color) to produce a (K-1)x(K-1) seed.
static vector<vector<int>> shrinkByOne(const vector<vector<int>>& G){
    int K=(int)G.size(); if(K<=1) return G;
    vector<int> tot(N+2,0); vector<vector<int>> rC(K,vector<int>(N+2,0)),cC(K,vector<int>(N+2,0));
    for(int r=0;r<K;r++)for(int c=0;c<K;c++){int v=G[r][c];tot[v]++;rC[r][v]++;cC[c][v]++;}
    auto rowLoss=[&](int r)->long long{ long long l=0;
        if(r>0)for(int c=0;c<K;c++) if(G[r-1][c]!=G[r][c])l++;
        if(r+1<K)for(int c=0;c<K;c++) if(G[r][c]!=G[r+1][c])l++;
        for(int v=1;v<=N;v++) if(rC[r][v]==tot[v]&&tot[v]>0)l+=100000; return l; };
    auto colLoss=[&](int c)->long long{ long long l=0;
        if(c>0)for(int r=0;r<K;r++) if(G[r][c-1]!=G[r][c])l++;
        if(c+1<K)for(int r=0;r<K;r++) if(G[r][c]!=G[r][c+1])l++;
        for(int v=1;v<=N;v++) if(cC[c][v]==tot[v]&&tot[v]>0)l+=100000; return l; };
    int bR=0; long long lr=rowLoss(0); for(int r=1;r<K;r++){long long x=rowLoss(r); if(x<lr){lr=x;bR=r;}}
    int bC=0; long long lc=colLoss(0); for(int c=1;c<K;c++){long long x=colLoss(c); if(x<lc){lc=x;bC=c;}}
    vector<vector<int>> o(K-1,vector<int>(K-1)); int rr=0;
    for(int r=0;r<K;r++){ if(r==bR)continue; int cc=0; for(int c=0;c<K;c++){ if(c==bC)continue; o[rr][cc]=G[r][c]; cc++; } rr++; }
    return o;
}

int main(){
    int T; if(scanf("%d",&T)!=1) return 0;
    while(T--){
        scanf("%d %d",&N,&M);
        memset(ADJ,0,sizeof(ADJ)); memset(adjmask,0,sizeof(adjmask)); EDGES.clear();
        for(int i=0;i<M;i++){int a,b;scanf("%d %d",&a,&b);ADJ[a][b]=ADJ[b][a]=1;EDGES.push_back({min(a,b),max(a,b)});}
        for(int v=1;v<=N;v++){ adjmask[v]=(1ULL<<v); for(int u=1;u<=N;u++) if(ADJ[v][u]) adjmask[v]|=(1ULL<<u); }
        auto t0=chrono::steady_clock::now(); HARD_DL=t0+chrono::milliseconds(820);
        if(N==1){ printf("1\n1\n1\n"); continue; }
        int lb=2; while(lb*lb<N)lb++; { int l2=2; while(2*l2*(l2-1)<M)l2++; lb=max(lb,l2); }

        int bestK=1<<29; vector<vector<int>> best;
        auto consider=[&](vector<vector<int>>& g){ if(g.empty())return; int K=(int)g.size(); if(K<bestK&&verifyGrid(g)){bestK=K;best=g;} };

        // dense path: random fill for edge-dense graphs
        if(8LL*M>=1LL*N*(N-1)){
            auto dl=t0+chrono::milliseconds(N>=35?820:650); if(dl>HARD_DL)dl=HARD_DL;
            int focus=max(lb,(int)(sqrt((double)M)*1.12)), hi=min(N,max(lb+16,focus+6));
            for(int k=lb;k<=hi&&!past(dl);k++){
                int at=(k<focus)?40:(N>=35?800:2200);
                vector<vector<int>> g; if(denseFill(k,at,dl,g)&&verifyGrid(g)){ consider(g); break; }
            }
        }

        // gather covering walks (Hamiltonian preferred)
        vector<vector<int>> walks;
        bool big=(N>=30);
        if(findHP(big?0.05:0.12)) walks.push_back(hpBest);
        auto hpEnd=t0+chrono::milliseconds(big?90:170);
        for(int a=0;a<(big?2:3)&&!past(hpEnd);a++) if(findHP(big?0.02:0.03)) walks.push_back(hpBest);
        {   vector<vector<int>> gws; vector<int> starts;
            int md=1<<29,mv=1; for(int v=1;v<=N;v++){int d=0;for(int u=1;u<=N;u++)d+=ADJ[v][u]; if(d<md){md=d;mv=v;}}
            starts.push_back(mv); for(int i=0;i<5;i++) starts.push_back(1+(int)(rng()%N));
            for(int s:starts){ auto w=coverWalk(s); if(!w.empty()) gws.push_back(w); }
            sort(gws.begin(),gws.end(),[](const vector<int>&a,const vector<int>&b){return a.size()<b.size();});
            for(int i=0;i<(int)gws.size()&&i<3;i++) walks.push_back(gws[i]);
        }
        walks.push_back(dfsTour());
        vector<int> postmanW=postmanWalk();
        { if(!postmanW.empty()) walks.push_back(postmanW); }
        for(auto& w:walks){
            vector<int> rw(w.rbegin(),w.rend());
            { vector<vector<int>> g; if(constructDiagGen(w,g)<=240) consider(g); }
            { vector<vector<int>> g; if(constructDiagGen(rw,g)<=240) consider(g); }
            { vector<vector<int>> g; if(constructDiag(w,g)<=240) consider(g); }
            { vector<vector<int>> g; if(constructDiag(rw,g)<=240) consider(g); }
        }
        // CHINESE-POSTMAN SNAKE candidate. The snake seed realizes every edge on the
        // walk (miss=0), so SA only fixes forbidden verticals — reaching K~sqrt(walkLen),
        // which beats the anti-diagonal on sparse/mid graphs. Try increasing K from the
        // smallest that fits the walk; SA-repair each; feed valid grids to consider()
        // (min-over-methods, so it only helps). Time-boxed to protect the shrink phase.
        {
            vector<int> sw = !postmanW.empty()? postmanW : dfsTour();
            int len=(int)sw.size();
            if(len>0){
                int k0=1; while(k0*k0<len) k0++;           // smallest K fitting the walk
                // Only worth it when the snake can plausibly beat bestK (sparse/mid:
                // k0 < bestK-1). Give each of a FEW K a real slice (~240ms) — spreading
                // thin fails to fix the forbidden verticals in time.
                // Descend, not ascend: the walk fits comfortably a couple above k0 and
                // SA fixes verticals FAST there, then we push down. Start at k0+3 with a
                // short slice; each success re-anchors lower. This lands near the true
                // achievable K (≈k0+2) instead of wasting budget on infeasible k0.
                if(k0 < bestK-1){
                    auto snDl=chrono::steady_clock::now()+chrono::milliseconds(560);
                    for(int K=min(k0+3,bestK-1); K>=max(k0,2) && !past(snDl); K--){
                        vector<vector<int>> seed=snakeSeed(sw,K); if(seed.empty()) continue;
                        long long rem=chrono::duration_cast<chrono::milliseconds>(snDl-chrono::steady_clock::now()).count();
                        if(rem<80) break;
                        long long slice=min(rem,220LL);
                        auto sl=chrono::steady_clock::now()+chrono::milliseconds(slice);
                        vector<vector<int>> res;
                        if(saRepair(seed,K,slice/1000.0,sl,res)&&verifyGrid(res)){ consider(res); }
                        else break; // this K failed fast; lower K only harder — stop
                    }
                }
            }
        }

        // guaranteed fallback: stripe from DFS tour never fails to be well-formed
        if(best.empty()){ vector<vector<int>> g; constructDiag(dfsTour(),g); if(!g.empty()) {best=g;bestK=(int)g.size();} }
        if(best.empty()){ // last resort: single-color-per-row stripes on N-length tour padded
            vector<int> w=dfsTour(); int K=min(240,(int)w.size()); best.assign(K,vector<int>(K,1));
            for(int r=0;r<K;r++) for(int c=0;c<K;c++) best[r][c]=w[min(r, (int)w.size()-1)]; bestK=K;
        }

        // Gated aggressive-SA-jump (density>=0.30): one big SA push a few units below
        // bestK from a fresh greedy seed, time-boxed so the incremental loop keeps its share.
        {
            double density=(N>1)?2.0*M/((double)N*(N-1)):0.0;
            long long remJ=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
            if(density>=0.30 && remJ>300){
                auto jd=chrono::steady_clock::now()+chrono::milliseconds((long long)(remJ*0.45));
                for(int jump=3;jump>=1 && !past(jd);jump--){
                    int tgt=bestK-jump; if(tgt<max(2,lb))continue;
                    long long rem=chrono::duration_cast<chrono::milliseconds>(jd-chrono::steady_clock::now()).count();
                    if(rem<120)break; long long js=min(rem,500LL);
                    auto jsl=chrono::steady_clock::now()+chrono::milliseconds(js);
                    vector<vector<int>> jr;
                    if(saRepair(greedySeed(tgt),tgt,js/1000.0,jsl,jr)&&verifyGrid(jr)){best=jr;bestK=tgt;break;}
                }
            }
        }

        // shrink with local search
        int target=bestK-1,variant=0;
        while(target>=max(2,lb)){
            if(past(HARD_DL)) break;
            long long rem=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
            if(rem<20) break;
            // Give each SA attempt a solid slice, but keep several attempts (random
            // restarts) available within the remaining budget.
            long long slice=min(min(300LL,150LL+40LL*variant),max(40LL,rem/2));
            vector<vector<int>> seed;
            int which=variant%3;
            if(which==0){ seed=best; while((int)seed.size()>target) seed=shrinkByOne(seed); if((int)seed.size()!=target) seed=rescale(best,target); }
            else if(which==1){ int off=bestK-target,dr=(int)(rng()%(off+1)),dc=(int)(rng()%(off+1));
                seed.assign(target,vector<int>(target)); for(int r=0;r<target;r++)for(int c=0;c<target;c++)seed[r][c]=best[r+dr][c+dc]; }
            else seed=rescale(best,target);
            if(variant>=3){ int KK=target; for(int i=0;i<KK*KK/24+1;i++) seed[rng()%KK][rng()%KK]=1+rng()%N; }
            // single shared budget for the whole iteration; every sub-call is
            // clamped to `sl` so the iteration never overruns its slice.
            auto sl=chrono::steady_clock::now()+chrono::milliseconds(slice); if(sl>HARD_DL)sl=HARD_DL;
            double budgetSec=slice/1000.0;
            vector<vector<int>> res;
            bool ok=false;
            // Primary engine: simulated annealing from a FRESH random-greedy fill.
            // A diverse random restart is what lets SA reach K values that both
            // dense-fill and shrink-then-repair plateau above. (Reusing the shrunk
            // grid as a seed traps SA in the same basin, so we don't.)
            {
                // Alternate seeds: shrink-of-best (great for sparse graphs where a
                // random fill is mostly forbidden) and a fresh greedy fill (great for
                // dense graphs where SA can wander to the floor). Both feed the same SA.
                vector<vector<int>> saSeed = (variant%2==0) ? seed : greedySeed(target);
                double sb=budgetSec - (chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-(sl-chrono::milliseconds(slice))).count())/1000.0;
                if(sb>0.01 && saRepair(saSeed,target,sb,sl,res)&&verifyGrid(res)) ok=true;
            }
            // Fallbacks: direct dense fill, then the older greedy repairs.
            if(!ok && !past(sl)){ vector<vector<int>> dg; auto dsl=chrono::steady_clock::now()+chrono::milliseconds(max<long long>(15,slice/4)); if(dsl>sl)dsl=sl;
                if(denseFill(target,60000,dsl,dg)&&verifyGrid(dg)){ res=dg; ok=true; } }
            if(!ok && !past(sl) && patchRepair(seed,target,sl,res)&&verifyGrid(res)) ok=true;
            else if(!ok && !past(sl)&&localSearch(seed,target,sl,res)&&verifyGrid(res)) ok=true;
            if(ok){ best=res;bestK=target;target--;variant=0; }
            else variant++;
        }

        int K=(int)best.size();
        string out; out.reserve((size_t)K*K*4+K*8);
        out+=to_string(K); out+='\n';
        for(int j=0;j<K;j++){ out+=to_string(K); if(j+1<K)out+=' '; } out+='\n';
        for(int i=0;i<K;i++){ for(int j=0;j<K;j++){ out+=to_string(best[i][j]); if(j+1<K)out+=' '; } out+='\n'; }
        fwrite(out.data(),1,out.size(),stdout);
    }
    return 0;
}
