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
        if(past(HARD_DL)) return 1<<29;
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
        if((it&255)==0){
            if(past(dl)) break;
            double el=chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-t0).count()/1000.0;
            if(el>budgetSec) break;
            Tmp=3.0*exp(-3.0*el/budgetSec); if(((it>>8)&63)==0) Tmp=max(Tmp,1.5); // periodic reheat
        }
        int p=rng()%(K*K), r=p/K,c=p%K, o=g[p];
        int x=1+rng()%N; if(x==o) continue;
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


// ---- constrained SA: only recolor to neighbor-legal colors ----------------
// Keeps the contact graph valid (forb=0) and optimizes missE/missC. Critical
// for sparse graphs where unconstrained recoloring almost always creates
// forbidden edges and never leaves the illegal basin.
static bool constrainedSA(vector<vector<int>> seed, int K, double budgetSec,
                          chrono::steady_clock::time_point dl, vector<vector<int>>& res){
    if((int)seed.size()<K||(int)seed[0].size()<K) return false;
    if(past(dl) || past(HARD_DL)) return false;
    auto t0=chrono::steady_clock::now();
    vector<int> g(K*K); for(int r=0;r<K;r++)for(int c=0;c<K;c++)g[r*K+c]=seed[r][c];
    static int cnt[45][45]; for(int a=0;a<45;a++)for(int b=0;b<45;b++)cnt[a][b]=0;
    vector<int> cc(N+1,0); int miss=0,missC=0;
    const int DR[4]={1,-1,0,0},DC[4]={0,0,1,-1};
    // require seed already contact-legal
    for(int r=0;r<K;r++)for(int c=0;c<K;c++){int v=g[r*K+c]; if(v<1||v>N) return false; cc[v]++;
        if(c+1<K){int u=g[r*K+c+1]; if(v!=u){ if(!ADJ[v][u]) return false; cnt[min(v,u)][max(v,u)]++;}}
        if(r+1<K){int u=g[(r+1)*K+c]; if(v!=u){ if(!ADJ[v][u]) return false; cnt[min(v,u)][max(v,u)]++;}}}
    for(auto&e:EDGES) if(cnt[e.first][e.second]==0) miss++;
    for(int v=1;v<=N;v++) if(cc[v]==0) missC++;
    if(miss==0&&missC==0){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
    vector<vector<int>> cellsOf(N+2); vector<int> posIn(K*K);
    for(int p=0;p<K*K;p++){ posIn[p]=(int)cellsOf[g[p]].size(); cellsOf[g[p]].push_back(p); }
    auto moveCell=[&](int p,int from,int to){
        auto&vf=cellsOf[from]; int idx=posIn[p], last=(int)vf.size()-1;
        if(idx!=last){ vf[idx]=vf[last]; posIn[vf[idx]]=idx; } vf.pop_back();
        posIn[p]=(int)cellsOf[to].size(); cellsOf[to].push_back(p);
    };
    double Tmp=8.0; long long it=0;
    long long bestSc=1000000LL*missC+miss; vector<int> bestG=g;
    while(true){
        it++;
        if((it&255)==0){
            if(past(dl)) break;
            double el=chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-t0).count()/1000.0;
            if(el>budgetSec) break;
            Tmp=8.0*exp(-3.5*el/budgetSec); if(((it>>8)&31)==0) Tmp=max(Tmp,2.5);
        }
        // early endgame kick when almost done
        if(miss+missC<=3 && (it&2047)==0){
            long long sc=1000000LL*missC+miss; if(sc<bestSc){bestSc=sc; bestG=g;}
        }
        int mode=(int)(rng()%100), p=-1, forced=-1;
        if(mode<(miss<=4?75:58) && miss>0){
            int e=-1; for(int t=0;t<40;t++){ int ee=(int)(rng()%EDGES.size()); if(cnt[EDGES[ee].first][EDGES[ee].second]==0){e=ee;break;} }
            if(e>=0){ int a=EDGES[e].first,b=EDGES[e].second; if(rng()&1) swap(a,b);
                if(!cellsOf[a].empty()){ int bp=cellsOf[a][rng()%cellsOf[a].size()];
                    int r=bp/K,c=bp%K,d=(int)(rng()%4),nr=r+DR[d],nc=c+DC[d];
                    if(nr>=0&&nr<K&&nc>=0&&nc<K){ p=nr*K+nc; forced=b; } } }
        } else if(mode<70 && missC>0){
            int col=-1; for(int t=0;t<20&&col<0;t++){ int v=1+(int)(rng()%N); if(cc[v]==0) col=v; }
            if(col>0){ for(int t=0;t<40;t++){ int q=(int)(rng()%(K*K)); int r=q/K,c=q%K; bool ok=true;
                for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                    int u=g[nr*K+nc]; if(u!=col&&!ADJ[u][col]){ok=false;break;}}
                if(ok){p=q;forced=col;break;} } }
        }
        if(p<0) p=(int)(rng()%(K*K));
        int r=p/K,c=p%K,o=g[p];
        int cand[45], ncc=0;
        for(int x=1;x<=N;x++){ if(x==o)continue; if(forced>0&&x!=forced)continue;
            bool ok=true; for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                int u=g[nr*K+nc]; if(u!=x&&!ADJ[u][x]){ok=false;break;}} if(ok) cand[ncc++]=x; }
        if(!ncc) continue;
        int x;
        if((int)(rng()%100)<18) x=cand[rng()%ncc];
        else {
            long long bestD=(1LL<<60); x=cand[0]; int ties=0;
            for(int i=0;i<ncc;i++){ int xx=cand[i]; int tch[8],de[8],nt=0; int dM=0,dMC=0;
                for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
                    if(o!=u){int a=min(o,u),b=max(o,u);int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]--;}
                    if(xx!=u){int a=min(xx,u),b=max(xx,u);int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]++;}}
                for(int j=0;j<nt;j++){int e=tch[j],a=e/45,b=e%45,oc=cnt[a][b],nn=oc+de[j]; if(oc>0&&nn==0)dM++; else if(oc==0&&nn>0)dM--;}
                if(cc[o]==1)dMC++; if(cc[xx]==0)dMC--;
                long long delta=1000LL*dMC + dM;
                if(delta<bestD){bestD=delta;x=xx;ties=1;} else if(delta==bestD){ties++; if((int)(rng()%ties)==0)x=xx;}
            }
        }
        int tch[8],de[8],nt=0; int dM=0,dMC=0;
        for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
            if(o!=u){int a=min(o,u),b=max(o,u);int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]--;}
            if(x!=u){int a=min(x,u),b=max(x,u);int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]++;}}
        for(int j=0;j<nt;j++){int e=tch[j],a=e/45,b=e%45,oc=cnt[a][b],nn=oc+de[j]; if(oc>0&&nn==0)dM++; else if(oc==0&&nn>0)dM--;}
        if(cc[o]==1)dMC++; if(cc[x]==0)dMC--;
        long long delta=1000LL*dMC + dM;
        bool acc = delta<=0 || (exp(-double(delta)/max(0.25,Tmp)) > (double)(rng()&0xffff)/65536.0);
        if(!acc) continue;
        for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
            if(o!=u){int a=min(o,u),b=max(o,u);int oc=cnt[a][b]; cnt[a][b]=oc-1; if(oc==1)miss++;}
            if(x!=u){int a=min(x,u),b=max(x,u);int oc=cnt[a][b]; cnt[a][b]=oc+1; if(oc==0)miss--;}}
        cc[o]--; if(cc[o]==0)missC++; cc[x]++; if(cc[x]==1)missC--;
        moveCell(p,o,x); g[p]=x;
        if(miss==0&&missC==0){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
        // track best partial for endgame repair
        long long sc=1000000LL*missC+miss;
        if(sc<bestSc){ bestSc=sc; bestG=g; }
    }
    // Endgame: from best partial, greedily place remaining missing edges (deadline-aware)
    if(bestSc<(1LL<<60) && bestSc>0 && !past(dl)){
        g=bestG;
        for(int a=0;a<45;a++)for(int b=0;b<45;b++)cnt[a][b]=0;
        fill(cc.begin(),cc.end(),0); miss=0; missC=0;
        for(int r=0;r<K;r++)for(int c=0;c<K;c++){int v=g[r*K+c]; cc[v]++;
            if(c+1<K){int u=g[r*K+c+1]; if(v!=u) cnt[min(v,u)][max(v,u)]++;}
            if(r+1<K){int u=g[(r+1)*K+c]; if(v!=u) cnt[min(v,u)][max(v,u)]++;}}
        for(auto&e:EDGES) if(cnt[e.first][e.second]==0) miss++;
        for(int v=1;v<=N;v++) if(cc[v]==0) missC++;
        for(int v=0;v<=N;v++) cellsOf[v].clear();
        for(int p=0;p<K*K;p++){ posIn[p]=(int)cellsOf[g[p]].size(); cellsOf[g[p]].push_back(p); }
        auto applyRecolor=[&](int p,int col){
            int o=g[p]; if(o==col) return;
            int r=p/K,c=p%K;
            for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
                if(o!=u){int a=min(o,u),b=max(o,u);int oc=cnt[a][b];cnt[a][b]=oc-1;if(oc==1)miss++;}
                if(col!=u){int a=min(col,u),b=max(col,u);int oc=cnt[a][b];cnt[a][b]=oc+1;if(oc==0)miss--;}}
            cc[o]--; if(cc[o]==0)missC++; cc[col]++; if(cc[col]==1)missC--;
            moveCell(p,o,col); g[p]=col;
        };
        auto legalAt=[&](int p,int col)->bool{
            int r=p/K,c=p%K;
            for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                int u=g[nr*K+nc]; if(u!=col&&!ADJ[u][col]) return false;}
            return true;
        };
        // Exhaustive neighbor placement when few edges missing
        if(miss>0 && miss<=6 && missC==0){
            for(int guard=0; guard<12 && miss>0 && !past(dl); guard++){
                bool placed=false;
                for(auto&e:EDGES){ if(cnt[e.first][e.second]>0) continue;
                    int A=e.first,B=e.second;
                    for(int pass=0;pass<2&&!placed;pass++){
                        int a=pass?B:A, b=pass?A:B;
                        for(int bi=0;bi<(int)cellsOf[a].size()&&!placed;bi++){
                            int bp=cellsOf[a][bi]; int r=bp/K,c=bp%K;
                            for(int d=0;d<4&&!placed;d++){
                                int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                                int p=nr*K+nc,o=g[p]; if(o==b){placed=true;break;}
                                if(cc[o]==1 && o!=b) continue;
                                if(!legalAt(p,b)) continue;
                                int before=miss; applyRecolor(p,b);
                                if(miss<before) placed=true;
                                else {
                                    // manual revert
                                    for(int d2=0;d2<4;d2++){int rr=nr+DR[d2],cc2=nc+DC[d2]; if(rr<0||rr>=K||cc2<0||cc2>=K)continue; int u=g[rr*K+cc2];
                                        if(b!=u){int aa=min(b,u),bb=max(b,u);int oc=cnt[aa][bb];cnt[aa][bb]=oc-1;if(oc==1)miss++;}
                                        if(o!=u){int aa=min(o,u),bb=max(o,u);int oc=cnt[aa][bb];cnt[aa][bb]=oc+1;if(oc==0)miss--;}}
                                    cc[b]--; if(cc[b]==0)missC++; cc[o]++; if(cc[o]==1)missC--;
                                    moveCell(p,b,o); g[p]=o;
                                }
                            }
                        }
                    }
                    if(placed) break;
                }
                if(!placed) break;
            }
            if(miss==0&&missC==0){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
        }
        for(int round=0;round<60 && (miss>0||missC>0) && !past(dl);round++){
            bool prog=false;
            if(missC>0){
                for(int col=1;col<=N;col++) if(cc[col]==0){
                    for(int p=0;p<K*K;p++){ int o=g[p]; if(cc[o]==1) continue;
                        if(!legalAt(p,col)) continue;
                        applyRecolor(p,col); prog=true; break;
                    }
                    if(prog) break;
                }
            }
            if(miss>0){
                for(auto&e:EDGES) if(cnt[e.first][e.second]==0){
                    int a=e.first,b=e.second;
                    for(int pass=0;pass<2&&!prog;pass++){
                        if(pass) swap(a,b);
                        if(cellsOf[a].empty()) continue;
                        for(int bi=0;bi<(int)cellsOf[a].size()&&!prog;bi++){
                            int bp=cellsOf[a][bi]; int r=bp/K,c=bp%K;
                            for(int d=0;d<4&&!prog;d++){
                                int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                                int p=nr*K+nc,o=g[p]; if(o==b) continue; if(cc[o]==1 && o!=b) continue;
                                if(!legalAt(p,b)) continue;
                                applyRecolor(p,b); prog=true;
                            }
                        }
                    }
                    if(prog) break;
                }
            }
            // pair force: limited scans (deadline-aware)
            if(!prog && miss>0 && miss<=8 && !past(dl)){
                for(auto&e:EDGES){ if(cnt[e.first][e.second]>0) continue; if(past(dl)) break;
                    int a=e.first,b=e.second;
                    // only scan a random subset of cells when K large
                    int lim = (K*K<=400)? K*K : 400;
                    for(int t=0;t<lim&&!prog;t++){
                        int p=(K*K<=400)? t : (int)(rng()%(K*K));
                        int r=p/K,c=p%K;
                        for(int d=0;d<2&&!prog;d++){
                            int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                            int p2=nr*K+nc; int o1=g[p],o2=g[p2];
                            if(cc[o1]==1||cc[o2]==1) continue;
                            for(int sw=0;sw<2&&!prog;sw++){
                                int x=sw?b:a,y=sw?a:b;
                                g[p]=x; g[p2]=y;
                                bool ok=true;
                                for(int d2=0;d2<4;d2++){int rr=r+DR[d2],cc2=c+DC[d2]; if(rr<0||rr>=K||cc2<0||cc2>=K)continue; int u=(rr*K+cc2==p2)?y:g[rr*K+cc2]; if(u!=x&&!ADJ[u][x]){ok=false;break;}}
                                if(ok) for(int d2=0;d2<4;d2++){int rr=nr+DR[d2],cc2=nc+DC[d2]; if(rr<0||rr>=K||cc2<0||cc2>=K)continue; int u=(rr*K+cc2==p)?x:g[rr*K+cc2]; if(u!=y&&!ADJ[u][y]){ok=false;break;}}
                                g[p]=o1; g[p2]=o2;
                                if(!ok) continue;
                                applyRecolor(p,x); applyRecolor(p2,y); prog=true;
                            }
                        }
                    }
                    if(prog) break;
                }
            }
            if(!prog) break;
            if(miss==0&&missC==0){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
        }
        if(miss==0&&missC==0){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
    }
    return false;
}

// Project a seed onto the contact-legal set by recoloring illegal cells.
static bool makeContactLegal(vector<vector<int>>& g){
    int K=(int)g.size(); if(K==0) return false;
    if(past(HARD_DL)) return false;
    const int DR[4]={1,-1,0,0},DC[4]={0,0,1,-1};
    int maxPass=min(K*K*4, max(40, K*8));
    for(int pass=0;pass<maxPass;pass++){
        if((pass&7)==0 && past(HARD_DL)) return false;
        bool any=false;
        for(int r=0;r<K;r++)for(int c=0;c<K;c++){
            int o=g[r][c]; bool bad=false;
            for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                int u=g[nr][nc]; if(o!=u&&!ADJ[o][u]){bad=true;break;}}
            if(!bad) continue; any=true;
            int cand[45],ncc=0;
            for(int x=1;x<=N;x++){ bool ok=true;
                for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                    int u=g[nr][nc]; if(u!=x&&!ADJ[u][x]){ok=false;break;}} if(ok) cand[ncc++]=x; }
            if(ncc) g[r][c]=cand[rng()%ncc]; else g[r][c]=1+(int)(rng()%N);
        }
        if(!any) return true;
    }
    for(int r=0;r<K;r++)for(int c=0;c<K;c++){ int v=g[r][c];
        if(c+1<K){int u=g[r][c+1]; if(v!=u&&!ADJ[v][u]) return false;}
        if(r+1<K){int u=g[r+1][c]; if(v!=u&&!ADJ[v][u]) return false;} }
    return true;
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
    vector<pair<long long,int>> rows,cols;
    for(int r=0;r<K;r++) rows.push_back({rowLoss(r),r});
    for(int c=0;c<K;c++) cols.push_back({colLoss(c),c});
    sort(rows.begin(),rows.end()); sort(cols.begin(),cols.end());
    int bR=rows[0].second, bC=cols[0].second;
    // try a few combos of top row/col candidates; pick min contact-break among non-critical
    long long bestL=rows[0].first+cols[0].first; int topR=min(3,(int)rows.size()), topC=min(3,(int)cols.size());
    for(int i=0;i<topR;i++) for(int j=0;j<topC;j++){
        long long L=rows[i].first+cols[j].first; if(L<bestL){bestL=L; bR=rows[i].second; bC=cols[j].second;}
    }
    // occasional random among top for diversity
    if((rng()%5)==0){ bR=rows[rng()%topR].second; bC=cols[rng()%topC].second; }
    vector<vector<int>> o(K-1,vector<int>(K-1)); int rr=0;
    for(int r=0;r<K;r++){ if(r==bR)continue; int cc=0; for(int c=0;c<K;c++){ if(c==bC)continue; o[rr][cc]=G[r][c]; cc++; } rr++; }
    return o;
}

// Place a covering walk as a snake on a KxK grid, then fill leftovers legally.
static vector<vector<int>> snakeSeed(const vector<int>& walk, int K){
    if(K<=0 || walk.empty()) return {};
    vector<vector<int>> g(K, vector<int>(K,0));
    int p=0, W=(int)walk.size();
    for(int r=0;r<K;r++){
        if(r%2==0){ for(int c=0;c<K;c++){ g[r][c]=walk[min(p,W-1)]; if(p<W)p++; } }
        else { for(int c=K-1;c>=0;c--){ g[r][c]=walk[min(p,W-1)]; if(p<W)p++; } }
    }
    // fill zeros (shouldn't happen) and repair illegal contacts greedily
    const int DR[4]={1,-1,0,0},DC[4]={0,0,1,-1};
    for(int pass=0;pass<3;pass++){
        for(int r=0;r<K;r++)for(int c=0;c<K;c++){
            int o=g[r][c]; bool bad=(o<1||o>N);
            if(!bad){ for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                int u=g[nr][nc]; if(u>=1&&u<=N&&o!=u&&!ADJ[o][u]){bad=true;break;}} }
            if(!bad) continue;
            int cand[45],ncc=0;
            for(int x=1;x<=N;x++){ bool ok=true;
                for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                    int u=g[nr][nc]; if(u>=1&&u<=N&&u!=x&&!ADJ[u][x]){ok=false;break;}}
                if(ok) cand[ncc++]=x; }
            if(ncc) g[r][c]=cand[rng()%ncc]; else g[r][c]=1+(int)(rng()%N);
        }
    }
    return g;
}

int main(){
    int T; if(scanf("%d",&T)!=1) return 0;
    while(T--){
        scanf("%d %d",&N,&M);
        memset(ADJ,0,sizeof(ADJ)); memset(adjmask,0,sizeof(adjmask)); EDGES.clear();
        for(int i=0;i<M;i++){int a,b;scanf("%d %d",&a,&b);ADJ[a][b]=ADJ[b][a]=1;EDGES.push_back({min(a,b),max(a,b)});}
        for(int v=1;v<=N;v++){ adjmask[v]=(1ULL<<v); for(int u=1;u<=N;u++) if(ADJ[v][u]) adjmask[v]|=(1ULL<<u); }
        // deterministic graph-hash seed (portfolio variants can change salt)
        unsigned long long gh=1469598103934665603ull;
        gh^=(unsigned long long)N; gh*=1099511628211ull; gh^=(unsigned long long)M; gh*=1099511628211ull;
        for(auto&e:EDGES){ gh^=(unsigned long long)e.first; gh*=1099511628211ull; gh^=(unsigned long long)e.second; gh*=1099511628211ull; }
        const unsigned long long SEED_SALT=0x31415926ULL; // tune portfolio
        rng.seed((unsigned)(gh^SEED_SALT^(gh>>32)));
        auto t0=chrono::steady_clock::now(); HARD_DL=t0+chrono::milliseconds(985);
        if(N==1){ printf("1\n1\n1\n"); continue; }
        int lb=2; while(lb*lb<N)lb++; { int l2=2; while(2*l2*(l2-1)<M)l2++; lb=max(lb,l2); }

        int bestK=1<<29; vector<vector<int>> best;
        auto consider=[&](vector<vector<int>>& g){ if(g.empty())return; int K=(int)g.size(); if(K<bestK&&verifyGrid(g)){bestK=K;best=g;} };

        // dense path: random fill. Dense graphs: large budget. Sparse: tiny probe only.
        double dens = (N>1)? (2.0*M/(N*(N-1))) : 0.0;
        bool sparse = dens < 0.22;
        {
            long long densMs;
            if(8LL*M>=1LL*N*(N-1)) densMs = (N>=35?700:580);
            else if(dens>=0.45) densMs = (N>=30?420:360);
            else if(sparse) densMs = (N>=32?35:70);
            else densMs = (N>=30?130:(N<=24?340:200));
            auto dl=t0+chrono::milliseconds(densMs); if(dl>HARD_DL)dl=HARD_DL;
            int focus=max(lb,(int)(sqrt((double)M)*1.12)), hi=min(N,max(lb+12,focus+5));
            for(int k=lb;k<=hi&&!past(dl);k++){
                int at=(k<=lb+1)? (N>=35?80:(N<=24?700:250)) : (k<focus)?60:(N>=35?500:1200);
                if(sparse) at=min(at,80);
                else if(dens>=0.28) at=max(at, (k<=lb+1)?(N<=24?1000:400):250);
                vector<vector<int>> g; if(denseFill(k,at,dl,g)&&verifyGrid(g)){ consider(g); break; }
            }
        }

        // gather covering walks (Hamiltonian preferred)
        vector<vector<int>> walks;
        bool big=(N>=30);
        // Cap HP budget — shrink is the high-ROI phase under the 1s wall.
        double hpSec = sparse ? (big?0.02:0.04) : (big?0.03:0.06);
        if(findHP(hpSec)) walks.push_back(hpBest);
        auto hpEnd=t0+chrono::milliseconds(sparse?(big?35:60):(big?50:90));
        for(int a=0;a<(big?1:2)&&!past(hpEnd);a++) if(findHP(big?0.012:0.02)) walks.push_back(hpBest);
        // longest-path style degree-greedy paths (cheap, good for sparse)
        {
            for(int trial=0; trial<(sparse?8:3) && !past(HARD_DL); trial++){
                vector<char> used(N+1,0); vector<int> path;
                int s=1+(int)(rng()%N);
                // prefer low degree start on sparse
                if(trial<3){ int md=1<<29; for(int v=1;v<=N;v++){int d=0;for(int u=1;u<=N;u++)d+=ADJ[v][u]; if(d<md){md=d;s=v;}} }
                path.push_back(s); used[s]=1;
                while((int)path.size()<N){
                    int v=path.back(); int bestu=-1,bd=1<<29,ties=0;
                    for(int u=1;u<=N;u++) if(ADJ[v][u]&&!used[u]){
                        int d=0; for(int w=1;w<=N;w++) if(ADJ[u][w]&&!used[w]) d++;
                        if(d<bd){bd=d;bestu=u;ties=1;} else if(d==bd){ties++; if((int)(rng()%ties)==0)bestu=u;}
                    }
                    if(bestu<0) break;
                    path.push_back(bestu); used[bestu]=1;
                }
                // if incomplete, splice coverWalk remainder
                if((int)path.size()<N){
                    auto w=coverWalk(path.back());
                    if(!w.empty()){
                        for(int i=1;i<(int)w.size();i++) if(!used[w[i]]){ path.push_back(w[i]); used[w[i]]=1; }
                    }
                }
                if((int)path.size()>=max(3,N*2/3)) walks.push_back(path);
            }
        }
        {   vector<vector<int>> gws; vector<int> starts;
            int md=1<<29,mv=1; for(int v=1;v<=N;v++){int d=0;for(int u=1;u<=N;u++)d+=ADJ[v][u]; if(d<md){md=d;mv=v;}}
            starts.push_back(mv); for(int i=0;i<(sparse?7:5);i++) starts.push_back(1+(int)(rng()%N));
            for(int s:starts){ if(past(HARD_DL))break; auto w=coverWalk(s); if(!w.empty()) gws.push_back(w); }
            sort(gws.begin(),gws.end(),[](const vector<int>&a,const vector<int>&b){return a.size()<b.size();});
            for(int i=0;i<(int)gws.size()&&i<(sparse?4:3);i++) walks.push_back(gws[i]);
        }
        walks.push_back(dfsTour());
        { for(int t=0;t<(sparse?3:4);t++){ if(past(HARD_DL))break; vector<int> pw=postmanWalk(); if(!pw.empty()) walks.push_back(pw);} }
        // Prefer shorter walks first (better anti-diagonal K)
        sort(walks.begin(),walks.end(),[](const vector<int>&a,const vector<int>&b){return a.size()<b.size();});
        vector<vector<int>> snakeWalks;
        for(int i=0;i<(int)walks.size() && i<6;i++) snakeWalks.push_back(walks[i]);
        int wlim = sparse ? min((int)walks.size(), 14) : min((int)walks.size(), 18);
        for(int wi=0; wi<wlim; wi++){
            if(past(HARD_DL)) break;
            auto& w=walks[wi];
            // leave at least 55% time for shrink on sparse large N
            if(sparse && N>=30){
                long long used=chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-t0).count();
                if(used>(N>=36?280:380)) break;
            }
            vector<int> rw(w.rbegin(),w.rend());
            { vector<vector<int>> g; if(constructDiagGen(w,g)<=240) consider(g); }
            if(past(HARD_DL)) break;
            { vector<vector<int>> g; if(constructDiagGen(rw,g)<=240) consider(g); }
            { vector<vector<int>> g; if(constructDiag(w,g)<=240) consider(g); }
            { vector<vector<int>> g; if(constructDiag(rw,g)<=240) consider(g); }
        }

        // guaranteed fallback: stripe from DFS tour never fails to be well-formed
        if(best.empty()){ vector<vector<int>> g; constructDiag(dfsTour(),g); if(!g.empty()) {best=g;bestK=(int)g.size();} }
        if(best.empty()){ // last resort: single-color-per-row stripes on N-length tour padded
            vector<int> w=dfsTour(); int K=min(240,(int)w.size()); best.assign(K,vector<int>(K,1));
            for(int r=0;r<K;r++) for(int c=0;c<K;c++) best[r][c]=w[min(r, (int)w.size()-1)]; bestK=K;
        }

        // shrink: geometric jumps then sequential -1; polish remaining time.
        int failsAt=0;
        auto tryAim=[&](int aim, long long slice)->bool{
            if(aim<lb||aim>=bestK||slice<8) return false;
            if(past(HARD_DL)) return false;
            auto sl=chrono::steady_clock::now()+chrono::milliseconds(slice); if(sl>HARD_DL)sl=HARD_DL;
            if(past(sl)) return false;
            auto mkSeed=[&](int which)->vector<vector<int>>{
                vector<vector<int>> seed;
                if(which==0){ seed=best; while((int)seed.size()>aim) seed=shrinkByOne(seed); if((int)seed.size()!=aim) seed=rescale(best,aim); }
                else if(which==1){ int off=bestK-aim; if(off<0)off=0; int dr=(int)(rng()%(off+1)),dc=(int)(rng()%(off+1));
                    seed.assign(aim,vector<int>(aim));
                    for(int r=0;r<aim;r++)for(int c=0;c<aim;c++) seed[r][c]=best[min(bestK-1,r+dr)][min(bestK-1,c+dc)]; }
                else if(which==2) seed=rescale(best,aim);
                else if(which==3){ seed=greedySeed(aim); }
                else if(which==4){ // multi shrinkByOne with random jitter after each step
                    seed=best; int guard=0;
                    while((int)seed.size()>aim && guard++<80){ seed=shrinkByOne(seed);
                        if((int)seed.size()>aim && (rng()&3)==0){ int KK=(int)seed.size(); seed[rng()%KK][rng()%KK]=1+(int)(rng()%N);} }
                    if((int)seed.size()!=aim) seed=rescale(best,aim);
                } else { // snake walk seed at target K
                    if(!snakeWalks.empty()) seed=snakeSeed(snakeWalks[rng()%snakeWalks.size()], aim);
                    else seed=greedySeed(aim);
                }
                if(which>=1 && which<=3){ int KK=aim; for(int i=0;i<KK*KK/12+4;i++) seed[rng()%KK][rng()%KK]=1+(int)(rng()%N); }
                return seed;
            };
            vector<vector<int>> res;
            int nSeed = slice>120 ? 5 : 4;
            if(N>=36) nSeed = slice>100 ? 6 : 5;
            for(int which=0; which<nSeed && !past(sl); which++){
                auto seed=mkSeed(which);
                if(!makeContactLegal(seed)) continue;
                long long rem2=chrono::duration_cast<chrono::milliseconds>(sl-chrono::steady_clock::now()).count();
                if(rem2<6) break;
                double sb=max(0.015, rem2/1000.0*0.96);
                if(constrainedSA(seed,aim,sb,sl,res)&&verifyGrid(res)){ best=res; bestK=aim; return true; }
            }
            // denseFill / denseFillBest seeds for mid+ density or large N last-mile
            if(!past(sl) && (dens>=0.18 || (N>=34 && aim<=bestK-1))){
                vector<vector<int>> dg; auto dsl=chrono::steady_clock::now()+chrono::milliseconds(max<long long>(8,slice/4)); if(dsl>sl)dsl=sl;
                if(denseFill(aim,40000,dsl,dg)&&verifyGrid(dg)){ best=dg; bestK=aim; return true; }
                if(!past(sl)){
                    vector<vector<int>> partial;
                    auto dsl2=chrono::steady_clock::now()+chrono::milliseconds(max<long long>(8,slice/5)); if(dsl2>sl)dsl2=sl;
                    // use best partial fill as SA seed even if missE>0
                    denseFillBest(aim, 8000, dsl2, partial);
                    if(!partial.empty()){
                        if(!makeContactLegal(partial)) {}
                        else {
                            long long rem2=chrono::duration_cast<chrono::milliseconds>(sl-chrono::steady_clock::now()).count();
                            if(rem2>10){ double sb=rem2/1000.0*0.95;
                                if(constrainedSA(partial,aim,sb,sl,res)&&verifyGrid(res)){ best=res; bestK=aim; return true; }
                            }
                        }
                    }
                }
            }
            if(!past(sl)){
                auto seed=mkSeed(0);
                long long rem2=chrono::duration_cast<chrono::milliseconds>(sl-chrono::steady_clock::now()).count();
                if(rem2>12){ vector<vector<int>> res2;
                    if(saRepair(seed,aim,rem2/1000.0,sl,res2)&&verifyGrid(res2)){ best=res2; bestK=aim; return true; }
                }
            }
            return false;
        };

        // Warm-start: snake seeds + ambitious mid-jump when construct is far from LB
        if(bestK-lb>=6 && !past(HARD_DL) && !snakeWalks.empty()){
            long long rem=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
            // try a ladder of aims from mid down
            vector<int> warmAims;
            int mid=max(lb, (bestK+lb)/2);
            int near=max(lb, (bestK*2+lb)/3);
            warmAims.push_back(near); warmAims.push_back(mid);
            if(N>=34) warmAims.push_back(max(lb, (int)(N*0.42+0.99))); // ~target R~0.42
            for(int aim:warmAims){
                if(past(HARD_DL)||aim>=bestK||aim<lb) continue;
                rem=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
                if(rem<50) break;
                // direct snake + legal + SA
                auto sl=chrono::steady_clock::now()+chrono::milliseconds(min(180LL, rem/3)); if(sl>HARD_DL)sl=HARD_DL;
                for(int si=0; si<(int)snakeWalks.size() && si<3 && !past(sl); si++){
                    auto seed=snakeSeed(snakeWalks[si], aim);
                    if(!makeContactLegal(seed)) continue;
                    vector<vector<int>> res;
                    long long rem2=chrono::duration_cast<chrono::milliseconds>(sl-chrono::steady_clock::now()).count();
                    if(rem2<8) break;
                    if(constrainedSA(seed,aim,rem2/1000.0*0.95,sl,res)&&verifyGrid(res)){ best=res; bestK=aim; break; }
                }
            }
            rem=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
            if(rem>80) tryAim(max(lb,(bestK*2+lb)/3), min(200LL, rem/3));
        }
        // Phase A: while gap large, prefer geometric mid-jump then -1
        while(bestK>max(2,lb) && !past(HARD_DL)){
            long long rem=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
            if(rem<12) break;
            int gap=bestK-lb;
            bool improved=false;
            // build aim list
            vector<int> aims;
            if(sparse && gap>=6){
                // once in the hard zone (K near 0.5N), only sequential -1/-2
                if(bestK*2 <= N+6){
                    aims.push_back(bestK-1);
                    if(gap>=2) aims.push_back(max(lb,bestK-2));
                } else {
                    aims.push_back(max(lb, bestK - max(2, gap/5)));
                    aims.push_back(max(lb, (bestK*3+lb)/4));
                    aims.push_back(bestK-1);
                    if(gap>=10) aims.push_back(max(lb,(bestK+lb)/2));
                }
            } else {
                aims.push_back(bestK-1);
                if(gap>=4) aims.push_back(max(lb,bestK-2));
                if(failsAt>=1 && gap>=6) aims.push_back(max(lb,(bestK*2+lb)/3));
                if(failsAt>=2 && gap>=8) aims.push_back(max(lb,(bestK+lb)/2));
            }
            // unique preserve order
            { vector<int> u; for(int a:aims){ if(a<lb||a>=bestK) continue; bool seen=false; for(int x:u) if(x==a)seen=true; if(!seen)u.push_back(a);} aims.swap(u); }
            for(int ai=0; ai<(int)aims.size() && !improved && !past(HARD_DL); ai++){
                rem=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
                if(rem<12) break;
                int aim=aims[ai];
                long long slice;
                if(aim==bestK-1) slice=min(N>=36?360LL:300LL, max(45LL, rem*3/4));
                else if(aim<=bestK-3) slice=min(240LL, max(45LL, rem/2));
                else slice=min(220LL, max(40LL, rem/2));
                if(tryAim(aim, slice)){ failsAt=0; improved=true; }
            }
            if(!improved){
                failsAt++;
                if(failsAt>=5) break;
            }
        }
        // Phase B polish: many short restarts on bestK-1 (last-mile K shaving)
        int polishFails=0;
        while(bestK>max(2,lb) && !past(HARD_DL)){
            long long rem=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
            if(rem<12) break;
            long long slice = (rem>200) ? min(170LL, rem/2) : max(25LL, rem-5);
            if(tryAim(bestK-1, slice)){ polishFails=0; continue; }
            polishFails++;
            rem=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
            // keep hammering while time remains (case1 previously left ~500ms unused)
            if(rem>40){
                // alternate bestK-1 snake direct and bestK-2
                if((polishFails%2)==0 && !snakeWalks.empty()){
                    auto sl=chrono::steady_clock::now()+chrono::milliseconds(min(rem-5, 120LL)); if(sl>HARD_DL)sl=HARD_DL;
                    int aim=bestK-1;
                    auto seed=snakeSeed(snakeWalks[rng()%snakeWalks.size()], aim);
                    if(makeContactLegal(seed)){
                        vector<vector<int>> res;
                        long long rem2=chrono::duration_cast<chrono::milliseconds>(sl-chrono::steady_clock::now()).count();
                        if(rem2>8 && constrainedSA(seed,aim,rem2/1000.0*0.95,sl,res)&&verifyGrid(res)){
                            best=res; bestK=aim; polishFails=0; continue;
                        }
                    }
                }
                if(bestK-2>=lb && tryAim(bestK-2, min(rem-5, max(40LL, rem/2)))){ polishFails=0; continue; }
                if(rem>40 && dens>=0.15){
                    vector<vector<int>> dg; auto dsl=chrono::steady_clock::now()+chrono::milliseconds(min(rem/3,100LL));
                    if(dsl>HARD_DL)dsl=HARD_DL;
                    if(denseFill(bestK-1,30000,dsl,dg)&&verifyGrid(dg)){ best=dg; bestK=bestK-1; polishFails=0; continue; }
                }
            }
            if(polishFails>=8 && rem<80) break; // keep going while time remains
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
