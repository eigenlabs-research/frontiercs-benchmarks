// World Map (IOI 2024) — Mixed strategy solver
// 1. Edge-covering walk → anti-diagonal grid (always valid, K≈L/2)
// 2. Random greedy fill + two-phase SA repair at target K
// 3. Shrink by crop+repair
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <cstring>
using namespace std;
typedef pair<int,int> pii;

static int N, M;
static bool ADJ[45][45];
static unsigned long long adjmask[45];
static vector<pii> EDGES;
static mt19937 rng(20240915u);
static chrono::steady_clock::time_point HARD_DL;
static inline bool past(chrono::steady_clock::time_point dl){
    return chrono::steady_clock::now() >= dl;
}

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

// Edge-covering DFS walk. Every edge appears as a consecutive pair.
static vector<int> edgeCoverWalk(){
    if(M==0) return {1};
    struct Frame { int v, ci; };
    static char used[45][45];
    memset(used,0,sizeof(used));
    vector<int> walk={1};
    vector<Frame> st={{1, 1}};
    while(!st.empty()){
        auto& b=st.back(); int v=b.v, &ci=b.ci;
        int u=-1;
        for(int w=ci; w<=N; w++){
            if(ADJ[v][w] && !used[min(v,w)][max(v,w)]){ u=w; ci=w+1; break; }
            ci=w+1;
        }
        if(u>0){
            used[min(v,u)][max(v,u)]=1;
            walk.push_back(u);
            st.push_back({u,1});
        } else {
            st.pop_back();
            if(!st.empty()) walk.push_back(st.back().v);
        }
    }
    for(auto& e:EDGES) if(!used[e.first][e.second]) return {};
    return walk;
}

// Anti-diagonal grid: grid[r][c]=walk[r+c]. Always valid for edge-covering walk.
static vector<vector<int>> antiDiagonalGrid(const vector<int>& walk){
    int L=(int)walk.size();
    if(L==0) return {};
    int K=(L+2)/2;
    if(K<1||K>240) return {};
    vector<vector<int>> g(K,vector<int>(K));
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
        int t=r+c;
        g[r][c]=walk[t < L ? t : L-1];
    }
    return g;
}

// Uniform-row grid: every row = walk. K=L. Valid for L<=240.
static vector<vector<int>> uniformRowGrid(const vector<int>& walk){
    int L=(int)walk.size();
    if(L<1||L>240) return {};
    vector<vector<int>> g(L,vector<int>(L));
    for(int r=0;r<L;r++) for(int c=0;c<L;c++) g[r][c]=walk[c];
    return g;
}

// Greedy random fill, legal w.r.t. up/left. Biases toward unseen colors/missing edges.
static vector<vector<int>> greedySeed(int K){
    int deg[45]={0};
    for(auto& e:EDGES){ deg[e.first]++; deg[e.second]++; }
    int total=K*K;
    vector<int> g(total), cc(N+1,0);
    int missC=N;
    static int ecnt[45][45];
    for(int a=0;a<45;a++) for(int b=0;b<45;b++) ecnt[a][b]=0;
    for(int p=0; p<total; p++){
        int r=p/K, c=p%K, up=r?g[(r-1)*K+c]:0, lf=c?g[p-1]:0, rem=total-p;
        int cand[45], wt[45], ncc=0, totW=0;
        for(int x=1;x<=N;x++){
            if(up && up!=x && !ADJ[up][x]) continue;
            if(lf && lf!=x && !ADJ[lf][x]) continue;
            if(missC>=rem && cc[x]>0) continue;
            int add=0;
            if(up&&up!=x&&!ecnt[min(up,x)][max(up,x)]) add++;
            if(lf&&lf!=x&&!ecnt[min(lf,x)][max(lf,x)]) add++;
            int w=1+80*add+deg[x];
            if(cc[x]==0) w+=300;
            cand[ncc]=x; wt[ncc]=w; totW+=w; ncc++;
        }
        if(!ncc){ g[p]=1+(p%N); continue; }
        int pick=(int)(rng()%totW), x=cand[0];
        for(int i=0;i<ncc;i++){ if(pick<wt[i]){ x=cand[i]; break; } pick-=wt[i]; }
        g[p]=x;
        if(cc[x]++==0) missC--;
        if(up&&up!=x){ int a=min(up,x),b=max(up,x); if(!ecnt[a][b]) ecnt[a][b]=1; }
        if(lf&&lf!=x){ int a=min(lf,x),b=max(lf,x); if(!ecnt[a][b]) ecnt[a][b]=1; }
    }
    vector<vector<int>> gg(K,vector<int>(K));
    for(int r=0;r<K;r++) for(int c=0;c<K;c++) gg[r][c]=g[r*K+c];
    return gg;
}

// Fix forbidden adjacencies by recoloring illegal cells.
static bool makeContactLegal(vector<vector<int>>& g){
    int K=(int)g.size();
    const int DR[4]={1,-1,0,0},DC[4]={0,0,1,-1};
    for(int pass=0; pass<K*K*2; pass++){
        bool any=false;
        for(int r=0;r<K;r++) for(int c=0;c<K;c++){
            int o=g[r][c]; bool bad=false;
            for(int d=0;d<4;d++){
                int nr=r+DR[d],nc=c+DC[d];
                if(nr<0||nr>=K||nc<0||nc>=K) continue;
                if(o!=g[nr][nc] && !ADJ[o][g[nr][nc]]){ bad=true; break; }
            }
            if(!bad) continue;
            any=true;
            int cand[45],ncc=0;
            for(int x=1;x<=N;x++){
                bool ok=true;
                for(int d=0;d<4;d++){
                    int nr=r+DR[d],nc=c+DC[d];
                    if(nr<0||nr>=K||nc<0||nc>=K) continue;
                    int u=g[nr][nc];
                    if(u!=x && !ADJ[u][x]){ ok=false; break; }
                }
                if(ok) cand[ncc++]=x;
            }
            if(ncc) g[r][c]=cand[rng()%ncc];
            else g[r][c]=1+(int)(rng()%N);
        }
        if(!any) return true;
    }
    // Final check
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
        int v=g[r][c];
        if(c+1<K){ int u=g[r][c+1]; if(v!=u&&!ADJ[v][u]) return false; }
        if(r+1<K){ int u=g[r+1][c]; if(v!=u&&!ADJ[v][u]) return false; }
    }
    return true;
}

// Unconstrained SA: can have forbidden edges (penalty weighted 1000).
// More flexible, escapes local minima.
static bool saRepair(vector<vector<int>> seed, int K, chrono::steady_clock::time_point dl,
                     vector<vector<int>>& res){
    if((int)seed.size()!=K||(int)seed[0].size()!=K) return false;
    auto t0=chrono::steady_clock::now();
    vector<int> g(K*K);
    for(int r=0;r<K;r++) for(int c=0;c<K;c++) g[r*K+c]=seed[r][c];

    static int cnt[45][45];
    for(int a=0;a<45;a++) for(int b=0;b<45;b++) cnt[a][b]=0;
    vector<int> cc(N+1,0);
    int forb=0, miss=0, missC=0;
    const int DR[4]={1,-1,0,0},DC[4]={0,0,1,-1};

    for(int p=0;p<K*K;p++){ int v=g[p],r=p/K,c=p%K; cc[v]++;
        if(c+1<K){ int u=g[p+1]; if(v!=u){ int a=min(v,u),b=max(v,u); cnt[a][b]++; if(!ADJ[a][b]) forb++; } }
        if(r+1<K){ int u=g[p+K]; if(v!=u){ int a=min(v,u),b=max(v,u); cnt[a][b]++; if(!ADJ[a][b]) forb++; } } }
    for(auto& e:EDGES) if(cnt[e.first][e.second]==0) miss++;
    for(int v=1;v<=N;v++) if(cc[v]==0) missC++;

    auto solved=[&](){ return forb==0&&miss==0&&missC==0; };
    if(solved()){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }

    double budgetMs=max(1.0, (double)chrono::duration_cast<chrono::milliseconds>(dl-t0).count());
    long long bestSc=1000000000LL*forb+1000000LL*missC+miss;
    vector<int> bestG=g;
    double Tmp=5.0;
    long long it=0;

    while(true){
        it++;
        if((it&127)==0){
            if(past(dl)) break;
            double el=chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-t0).count()/budgetMs;
            Tmp=5.0*exp(-4.0*el);
            if(((it>>9)&15)==0) Tmp=max(Tmp, 1.2);
        }
        int p=(int)(rng()%(K*K));
        int r=p/K, c=p%K, o=g[p];
        int x=0;
        // 70% random color, 30% smart (biased toward high-degree colors in neighborhood)
        if((int)(rng()%100)<30){
            // Prefer colors that are legal with some neighbor and reduce violations
            int nb[4],nn=0;
            for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; nb[nn++]=g[nr*K+nc];}
            // Find best candidate by delta
            int cand[45],ncc=0; long long bestD=(1LL<<60); int bestX=o;
            for(int xx=1;xx<=N;xx++){ if(xx==o) continue; cand[ncc++]=xx; }
            for(int i=0;i<ncc;i++){ int xx=cand[i];
                int df=0; int tch[8],de[8],nt=0;
                for(int j=0;j<nn;j++){ int nc=nb[j];
                    if(o!=nc){ int a=min(o,nc),b=max(o,nc);
                        if(ADJ[a][b]){int e=a*45+b,k;for(k=0;k<nt;k++)if(tch[k]==e)break;if(k==nt){tch[nt]=e;de[nt]=0;nt++;}de[k]--;}else df--;}
                    if(xx!=nc){ int a=min(xx,nc),b=max(xx,nc);
                        if(ADJ[a][b]){int e=a*45+b,k;for(k=0;k<nt;k++)if(tch[k]==e)break;if(k==nt){tch[nt]=e;de[nt]=0;nt++;}de[k]++;}else df++;}}
                int dm=0; for(int j=0;j<nt;j++){int e=tch[j],a=e/45,b=e%45,oc=cnt[a][b],ncnt=oc+de[j]; if(oc>0&&ncnt==0)dm++; else if(oc==0&&ncnt>0)dm--;}
                int dmc=0; if(cc[o]==1)dmc++; if(cc[xx]==0)dmc--;
                long long delta=1000LL*df+50LL*dm+50LL*dmc;
                if(delta<bestD){ bestD=delta; bestX=xx; }
            }
            x=bestX;
        } else {
            x=1+(int)(rng()%N);
        }
        if(x==o) continue;

        int nb[4],nn=0;
        for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; nb[nn++]=g[nr*K+nc];}
        int df=0; int tch[8],de[8],nt=0;
        for(int i=0;i<nn;i++){ int nc=nb[i];
            if(o!=nc){ int a=min(o,nc),b=max(o,nc);
                if(ADJ[a][b]){int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]--;}else df--;}
            if(x!=nc){ int a=min(x,nc),b=max(x,nc);
                if(ADJ[a][b]){int e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]++;}else df++;}}
        int dm=0; for(int j=0;j<nt;j++){int e=tch[j],a=e/45,b=e%45,oc=cnt[a][b],ncnt=oc+de[j]; if(oc>0&&ncnt==0)dm++; else if(oc==0&&ncnt>0)dm--;}
        int dmc=0; if(cc[o]==1)dmc++; if(cc[x]==0)dmc--;
        long long delta=1000LL*df+50LL*dm+50LL*dmc;
        bool acc=delta<=0||(exp(-(double)delta/max(0.2,Tmp*50.0))>(double)(rng()&0xffff)/65536.0);
        if(!acc) continue;

        for(int i=0;i<nn;i++){ int nc=nb[i];
            if(o!=nc){int a=min(o,nc),b=max(o,nc),oc=cnt[a][b];cnt[a][b]=oc-1;if(ADJ[a][b]){if(oc==1)miss++;}else forb--;}
            else {forb--;}
            if(x!=nc){int a=min(x,nc),b=max(x,nc),oc=cnt[a][b];cnt[a][b]=oc+1;if(ADJ[a][b]){if(oc==0)miss--;}else forb++;}
            else {forb++;}
        }
        cc[o]--; if(cc[o]==0) missC++;
        cc[x]++; if(cc[x]==1) missC--;
        g[p]=x;

        if(solved()){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }

        long long sc=1000000000LL*forb+1000000LL*missC+miss;
        if(sc<bestSc){ bestSc=sc; bestG=g; }
    }

    // If SA couldn't solve, try makeContactLegal + constrained SA on best partial
    if(bestSc>0){
        g=bestG;
        // Rebuild state
        for(int a=0;a<45;a++) for(int b=0;b<45;b++) cnt[a][b]=0;
        fill(cc.begin(),cc.end(),0); forb=0; miss=0; missC=0;
        for(int p=0;p<K*K;p++){ int v=g[p],r=p/K,c=p%K; cc[v]++;
            if(c+1<K){ int u=g[p+1]; if(v!=u){ int a=min(v,u),b=max(v,u); cnt[a][b]++; if(!ADJ[a][b]){ forb++; return false; } } }
            if(r+1<K){ int u=g[p+K]; if(v!=u){ int a=min(v,u),b=max(v,u); cnt[a][b]++; if(!ADJ[a][b]){ forb++; return false; } } } }
        if(forb>0) return false;
        for(auto& e:EDGES) if(cnt[e.first][e.second]==0) miss++;
        for(int v=1;v<=N;v++) if(cc[v]==0) missC++;
        if(miss==0&&missC==0){
            res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true;
        }
    }
    return false;
}

// Constrained SA: only contact-legal moves. Optimizes missE/missC.
static bool constrainedSA(vector<vector<int>> seed, int K, chrono::steady_clock::time_point dl,
                          vector<vector<int>>& res){
    if((int)seed.size()!=K||(int)seed[0].size()!=K) return false;
    auto t0=chrono::steady_clock::now();
    vector<int> g(K*K);
    for(int r=0;r<K;r++) for(int c=0;c<K;c++) g[r*K+c]=seed[r][c];

    static int cnt[45][45];
    for(int a=0;a<45;a++) for(int b=0;b<45;b++) cnt[a][b]=0;
    vector<int> cc(N+1,0);
    int miss=0, missC=0;
    const int DR[4]={1,-1,0,0},DC[4]={0,0,1,-1};

    for(int r=0;r<K;r++) for(int c=0;c<K;c++){ int v=g[r*K+c];
        if(v<1||v>N) return false; cc[v]++;
        if(c+1<K){ int u=g[r*K+c+1]; if(v!=u){ if(!ADJ[v][u]) return false; cnt[min(v,u)][max(v,u)]++; }}
        if(r+1<K){ int u=g[(r+1)*K+c]; if(v!=u){ if(!ADJ[v][u]) return false; cnt[min(v,u)][max(v,u)]++; }}}
    for(auto& e:EDGES) if(cnt[e.first][e.second]==0) miss++;
    for(int v=1;v<=N;v++) if(cc[v]==0) missC++;

    auto solved=[&](){ return miss==0&&missC==0; };
    if(solved()){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }

    vector<vector<int>> cellsOf(N+2);
    vector<int> posIn(K*K);
    for(int p=0;p<K*K;p++){ posIn[p]=(int)cellsOf[g[p]].size(); cellsOf[g[p]].push_back(p); }
    auto moveCell=[&](int p,int from,int to){
        auto& vf=cellsOf[from];
        int idx=posIn[p], last=(int)vf.size()-1;
        if(idx!=last){ vf[idx]=vf[last]; posIn[vf[idx]]=idx; } vf.pop_back();
        posIn[p]=(int)cellsOf[to].size(); cellsOf[to].push_back(p);
    };

    double budgetMs=max(5.0, (double)chrono::duration_cast<chrono::milliseconds>(dl-t0).count());
    long long bestSc=1000000LL*missC+miss;
    vector<int> bestG=g;
    double Tmp=6.0;
    long long it=0;

    while(true){
        it++;
        if((it&127)==0){
            if(past(dl)) break;
            double el=chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-t0).count()/budgetMs;
            Tmp=6.0*exp(-4.0*el);
            if(((it>>9)&7)==0) Tmp=max(Tmp, 1.5);
        }
        int mode=(int)(rng()%100), p=-1, forced=-1;
        if(mode<50 && miss>0){
            int e=-1;
            for(int t=0;t<25;t++){ int ee=(int)(rng()%EDGES.size()); if(cnt[EDGES[ee].first][EDGES[ee].second]==0){e=ee;break;} }
            if(e>=0){
                int a=EDGES[e].first, b=EDGES[e].second;
                if(rng()&1) swap(a,b);
                if(!cellsOf[a].empty()){
                    int bp=cellsOf[a][rng()%cellsOf[a].size()];
                    int r=bp/K, c=bp%K, d=rng()%4, nr=r+DR[d], nc=c+DC[d];
                    if(nr>=0&&nr<K&&nc>=0&&nc<K){ p=nr*K+nc; forced=b; }
                }
            }
        } else if(mode<70 && missC>0){
            int col=-1; for(int t=0;t<12&&col<0;t++){ int v=1+(rng()%N); if(cc[v]==0) col=v; }
            if(col>0){
                for(int t=0;t<25;t++){
                    int q=(int)(rng()%(K*K)), r=q/K, c=q%K; bool ok=true;
                    for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                        int u=g[nr*K+nc]; if(u!=col&&!ADJ[u][col]){ok=false;break;}}
                    if(ok){ p=q; forced=col; break; }
                }
            }
        }
        if(p<0) p=(int)(rng()%(K*K));
        int r=p/K, c=p%K, o=g[p];

        int cand[45], ncc=0;
        for(int x=1;x<=N;x++){ if(x==o)continue; if(forced>0&&x!=forced)continue;
            bool ok=true; for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                int u=g[nr*K+nc]; if(u!=x&&!ADJ[u][x]){ok=false;break;}} if(ok)cand[ncc++]=x;}
        if(!ncc) continue;

        int x;
        if((int)(rng()%100)<12) x=cand[rng()%ncc];
        else {
            long long bestD=(1LL<<60); x=cand[0]; int ties=0;
            for(int i=0;i<ncc;i++){ int xx=cand[i];
                int tch[8],de[8],nt=0, dM=0, dMC=0;
                for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
                    if(o!=u){int a=min(o,u),b=max(o,u),e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]--;}
                    if(xx!=u){int a=min(xx,u),b=max(xx,u),e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]++;}}
                for(int j=0;j<nt;j++){int e=tch[j],a=e/45,b=e%45,oc=cnt[a][b],nn=oc+de[j]; if(oc>0&&nn==0)dM++;else if(oc==0&&nn>0)dM--;}
                if(cc[o]==1)dMC++; if(cc[xx]==0)dMC--;
                long long delta=1000LL*dMC+dM;
                if(delta<bestD){bestD=delta;x=xx;ties=1;}else if(delta==bestD){ties++;if((int)(rng()%ties)==0)x=xx;}
            }
        }

        int tch[8],de[8],nt=0, dM=0, dMC=0;
        for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
            if(o!=u){int a=min(o,u),b=max(o,u),e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]--;}
            if(x!=u){int a=min(x,u),b=max(x,u),e=a*45+b,j;for(j=0;j<nt;j++)if(tch[j]==e)break;if(j==nt){tch[nt]=e;de[nt]=0;nt++;}de[j]++;}}
        for(int j=0;j<nt;j++){int e=tch[j],a=e/45,b=e%45,oc=cnt[a][b],nn=oc+de[j]; if(oc>0&&nn==0)dM++;else if(oc==0&&nn>0)dM--;}
        if(cc[o]==1)dMC++; if(cc[x]==0)dMC--;
        long long delta=1000LL*dMC+dM;
        bool acc=delta<=0||(exp(-(double)delta/max(0.15,Tmp))>(double)(rng()&0xffff)/65536.0);
        if(!acc) continue;

        for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
            if(o!=u){int a=min(o,u),b=max(o,u);int oc=cnt[a][b];cnt[a][b]=oc-1;if(oc==1)miss++;}
            if(x!=u){int a=min(x,u),b=max(x,u);int oc=cnt[a][b];cnt[a][b]=oc+1;if(oc==0)miss--;}}
        cc[o]--; if(cc[o]==0) missC++;
        cc[x]++; if(cc[x]==1) missC--;
        moveCell(p,o,x); g[p]=x;
        if(solved()){ res.assign(K,vector<int>(K)); for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c]; return true; }
        long long sc=1000000LL*missC+miss;
        if(sc<bestSc){ bestSc=sc; bestG=g; }
    }
    // Endgame greedy repair from best partial
    if(bestSc<(1LL<<60)&&bestSc>0){
        g=bestG;
        for(int a=0;a<45;a++)for(int b=0;b<45;b++)cnt[a][b]=0;
        fill(cc.begin(),cc.end(),0);
        for(int r=0;r<K;r++) for(int c=0;c<K;c++){ int v=g[r*K+c]; cc[v]++;
            if(c+1<K){ int u=g[r*K+c+1]; if(v!=u)cnt[min(v,u)][max(v,u)]++;}
            if(r+1<K){ int u=g[(r+1)*K+c]; if(v!=u)cnt[min(v,u)][max(v,u)]++;}}
        miss=0; for(auto& e:EDGES) if(cnt[e.first][e.second]==0) miss++;
        missC=0; for(int v=1;v<=N;v++) if(cc[v]==0) missC++;
        for(int v=0;v<=N;v++)cellsOf[v].clear();
        for(int p=0;p<K*K;p++){posIn[p]=(int)cellsOf[g[p]].size();cellsOf[g[p]].push_back(p);}
        for(int round=0;round<60&&!past(dl);round++){
            bool prog=false;
            if(missC>0){
                for(int col=1;col<=N&&!past(dl);col++) if(cc[col]==0){
                    for(int p=0;p<K*K&&!past(dl);p++){
                        int r=p/K,c=p%K,o=g[p]; if(cc[o]==1)continue;
                        bool ok=true;
                        for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                            int u=g[nr*K+nc]; if(u!=col&&!ADJ[u][col]){ok=false;break;}}
                        if(!ok) continue;
                        for(int d=0;d<4;d++){int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue; int u=g[nr*K+nc];
                            if(o!=u){int a=min(o,u),b=max(o,u);int oc=cnt[a][b];cnt[a][b]=oc-1;if(oc==1)miss++;}
                            if(col!=u){int a=min(col,u),b=max(col,u);int oc=cnt[a][b];cnt[a][b]=oc+1;if(oc==0)miss--;}}
                        cc[o]--; if(cc[o]==0)missC++; cc[col]++; if(cc[col]==1)missC--;
                        moveCell(p,o,col); g[p]=col; prog=true; break;
                    } if(prog) break;
                }
            }
            if(miss>0){
                for(auto& e:EDGES) if(cnt[e.first][e.second]==0){
                    int a=e.first,b=e.second;
                    for(int pass=0;pass<2&&!prog;pass++){
                        if(pass)swap(a,b);
                        if(cellsOf[a].empty())continue;
                        for(int bi=0;bi<(int)cellsOf[a].size()&&!prog;bi++){
                            int bp=cellsOf[a][bi], r=bp/K, c=bp%K;
                            for(int d=0;d<4&&!prog;d++){
                                int nr=r+DR[d],nc=c+DC[d]; if(nr<0||nr>=K||nc<0||nc>=K)continue;
                                int p=nr*K+nc, o=g[p]; if(o==b)continue; if(cc[o]==1&&o!=b)continue;
                                bool ok=true;
                                for(int d2=0;d2<4;d2++){int rr=nr+DR[d2],cc2=nc+DC[d2]; if(rr<0||rr>=K||cc2<0||cc2>=K)continue;
                                    int u=g[rr*K+cc2]; if(u!=b&&!ADJ[u][b]){ok=false;break;}}
                                if(!ok)continue;
                                for(int d2=0;d2<4;d2++){int rr=nr+DR[d2],cc2=nc+DC[d2]; if(rr<0||rr>=K||cc2<0||cc2>=K)continue;int u=g[rr*K+cc2];
                                    if(o!=u){int aa=min(o,u),bb=max(o,u);int oc=cnt[aa][bb];cnt[aa][bb]=oc-1;if(oc==1)miss++;}
                                    if(b!=u){int aa=min(b,u),bb=max(b,u);int oc=cnt[aa][bb];cnt[aa][bb]=oc+1;if(oc==0)miss--;}}
                                cc[o]--; if(cc[o]==0)missC++; cc[b]++; if(cc[b]==1)missC--;
                                moveCell(p,o,b); g[p]=b; prog=true;
                            }
                        }
                    } if(prog)break;
                }
            }
            if(!prog)break;
            if(miss==0&&missC==0){res.assign(K,vector<int>(K));for(int r=0;r<K;r++)for(int c=0;c<K;c++)res[r][c]=g[r*K+c];return true;}
        }
    }
    return false;
}

// Shrink by removing least-valuable row+column
static vector<vector<int>> shrinkByOne(const vector<vector<int>>& G){
    int K=(int)G.size(); if(K<=1) return G;
    vector<int> tot(N+1,0);
    vector<vector<int>> rc(K,vector<int>(N+1)), cc2(K,vector<int>(N+1));
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){ int v=G[r][c]; tot[v]++; rc[r][v]++; cc2[c][v]++; }
    vector<int> rLoss(K), cLoss(K);
    for(int r=0;r<K;r++) for(int c=0;c<K;c++){
        if(c+1<K && G[r][c]!=G[r][c+1]){ rLoss[r]++; cLoss[c]++; cLoss[c+1]++; }
        if(r+1<K && G[r][c]!=G[r+1][c]){ rLoss[r]++; rLoss[r+1]++; cLoss[c]++; }
    }
    auto rowPen=[&](int r)->long long{ long long l=rLoss[r]; for(int v=1;v<=N;v++) if(rc[r][v]==tot[v]&&tot[v]>0) l+=1000000; return l; };
    auto colPen=[&](int c)->long long{ long long l=cLoss[c]; for(int v=1;v<=N;v++) if(cc2[c][v]==tot[v]&&tot[v]>0) l+=1000000; return l; };
    int bR=0; long long rp=rowPen(0); for(int r=1;r<K;r++){ long long x=rowPen(r); if(x<rp){rp=x;bR=r;} }
    int bC=0; long long cp=colPen(0); for(int c=1;c<K;c++){ long long x=colPen(c); if(x<cp){cp=x;bC=c;} }
    vector<vector<int>> o(K-1,vector<int>(K-1));
    int rr=0;
    for(int r=0;r<K;r++){ if(r==bR) continue; int cc=0; for(int c=0;c<K;c++){ if(c==bC)continue; o[rr][cc]=G[r][c]; cc++; } rr++; }
    return o;
}

static vector<vector<int>> rescale(const vector<vector<int>>& G, int K2){
    int K=(int)G.size();
    vector<vector<int>> o(K2,vector<int>(K2));
    for(int r=0;r<K2;r++) for(int c=0;c<K2;c++)
        o[r][c]=G[(int)((long long)r*K/K2)][(int)((long long)c*K/K2)];
    return o;
}

// Try to find a valid grid at target K via greedy+SA. Returns empty if fails.
static bool tryDirectK(int K, int attempts, chrono::steady_clock::time_point dl,
                       vector<vector<int>>& res){
    for(int at=0; at<attempts && !past(dl); at++){
        auto seed=greedySeed(K);
        // Phase 1: unconstrained SA to reach a good partial
        auto dl1=chrono::steady_clock::now()+chrono::milliseconds(15);
        if(dl1>dl) dl1=dl;
        vector<vector<int>> tmp;
        saRepair(seed, K, dl1, tmp);
        if(tmp.empty()){
            // Try makeContactLegal on seed
            if(!makeContactLegal(seed)) continue;
            tmp=seed;
        }
        // Phase 2: constrained SA to finish
        auto rem=chrono::duration_cast<chrono::milliseconds>(dl-chrono::steady_clock::now()).count();
        if(rem<5) break;
        auto dl2=chrono::steady_clock::now()+chrono::milliseconds(rem);
        if(dl2>dl) dl2=dl;
        if(constrainedSA(tmp, K, dl2, res)){
            if(verifyGrid(res)) return true;
        }
    }
    return false;
}

int main(){
    int T; scanf("%d",&T);
    while(T--){
        scanf("%d %d",&N,&M);
        memset(ADJ,0,sizeof(ADJ));
        memset(adjmask,0,sizeof(adjmask));
        EDGES.clear();
        for(int i=0;i<M;i++){
            int a,b; scanf("%d %d",&a,&b);
            ADJ[a][b]=ADJ[b][a]=1;
            EDGES.push_back({min(a,b),max(a,b)});
        }
        for(int v=1;v<=N;v++){
            adjmask[v]=(1ULL<<v);
            for(int u=1;u<=N;u++) if(ADJ[v][u]) adjmask[v]|=(1ULL<<u);
        }

        auto t0=chrono::steady_clock::now();
        HARD_DL=t0+chrono::milliseconds(940);
        int bestK=999;
        vector<vector<int>> best;

        if(N==1){ printf("1\n1\n1\n"); continue; }
        if(N==2){ printf("2\n2 2\n1 2\n2 1\n"); continue; }

        vector<int> walk=edgeCoverWalk();

        // Compute lower bound for K
        int lb=1; while(lb*lb<N) lb++;
        { int l2=2; while(2*l2*(l2-1)<M) l2++; lb=max(lb,l2); }

        // === Phase 1: Direct construction at target K (spend ~800ms) ===
        auto phase1DL=t0+chrono::milliseconds(800);
        if(phase1DL>HARD_DL) phase1DL=HARD_DL;

        // For each K from lb upward, try to construct a valid grid
        int ub=min(240, N*4);
        for(int K=lb; K<=ub && !past(phase1DL); K++){
            int attempts = (K<=lb+1) ? 120 : (K<=lb+5) ? 60 : (K<=lb+10) ? 30 : 15;
            vector<vector<int>> g;
            if(tryDirectK(K, attempts, phase1DL, g)){
                best=g; bestK=K; break;
            }
        }

        // === Phase 2: FASTER direct at larger K if nothing found yet ===
        if(best.empty()){
            int K=min(240, (int)ceil(sqrt((double)M*2.0)+5));
            if(K<lb) K=lb;
            for(; K<=min(240, lb+20) && !past(HARD_DL); K++){
                vector<vector<int>> g;
                if(tryDirectK(K, 30, HARD_DL, g)){
                    best=g; bestK=K; break;
                }
            }
        }

        // === Phase 3: Anti-diagonal fallback ===
        {
            vector<vector<int>> g=antiDiagonalGrid(walk);
            if(!g.empty() && verifyGrid(g) && (int)g.size()<bestK){
                best=g; bestK=(int)g.size();
            }
        }
        if(best.empty()){
            int L=(int)walk.size();
            if(L<=240){
                vector<vector<int>> g=uniformRowGrid(walk);
                if(!g.empty() && verifyGrid(g)){ best=g; bestK=L; }
            }
        }

        // === Phase 4: Shrink from best (spend remaining time) ===
        int failsAt=0;
        while(bestK>lb && !past(HARD_DL)){
            long long rem=chrono::duration_cast<chrono::milliseconds>(HARD_DL-chrono::steady_clock::now()).count();
            if(rem<20) break;

            bool improved=false;
            // Try bestK-1
            long long slice=min(250LL, max(25LL, rem*1/2));
            auto sl=chrono::steady_clock::now()+chrono::milliseconds(slice);
            if(sl>HARD_DL) sl=HARD_DL;
            int aim=bestK-1;

            for(int which=0; which<5 && !past(sl) && !improved; which++){
                vector<vector<int>> seed;
                if(which==0) seed=shrinkByOne(best);
                else if(which==1){ seed=rescale(best,aim);
                    for(int i=0;i<aim*aim/8;i++) seed[rng()%aim][rng()%aim]=1+(rng()%N); }
                else if(which==2){
                    seed.assign(aim,vector<int>(aim));
                    for(int r=0;r<aim;r++) for(int c=0;c<aim;c++)
                        seed[r][c]=best[min(bestK-1,r+1)][min(bestK-1,c+1)];
                }
                else seed=greedySeed(aim);

                if(!makeContactLegal(seed)) continue;
                long long rem2=chrono::duration_cast<chrono::milliseconds>(sl-chrono::steady_clock::now()).count();
                if(rem2<6) break;
                auto dl2=chrono::steady_clock::now()+chrono::milliseconds(rem2);
                if(dl2>sl) dl2=sl;

                vector<vector<int>> res;
                if(constrainedSA(seed, aim, dl2, res)){
                    if(verifyGrid(res)){ best=res; bestK=aim; improved=true; failsAt=0; }
                }
            }

            if(!improved){ failsAt++; }
            else continue;

            // Multi-step jumps on repeated failure
            if(failsAt>=2 && bestK-lb>=3){
                int aim2=max(lb, bestK-2);
                if(failsAt>=3) aim2=max(lb, (bestK+lb)/2);
                if(aim2==bestK-1) aim2=max(lb,bestK-2);
                slice=min(180LL, max(15LL, rem*1/3));
                sl=chrono::steady_clock::now()+chrono::milliseconds(slice);
                if(sl>HARD_DL) sl=HARD_DL;

                for(int which=0; which<4 && !past(sl) && !improved; which++){
                    vector<vector<int>> seed;
                    if(which==0) seed=rescale(best, aim2);
                    else if(which==1){
                        seed.assign(aim2,vector<int>(aim2));
                        for(int r=0;r<aim2;r++) for(int c=0;c<aim2;c++)
                            seed[r][c]=best[(long long)r*bestK/aim2][(long long)c*bestK/aim2];
                    }
                    else seed=greedySeed(aim2);
                    if(!makeContactLegal(seed)) continue;
                    long long rem2=chrono::duration_cast<chrono::milliseconds>(sl-chrono::steady_clock::now()).count();
                    if(rem2<6) break;
                    auto dl2=chrono::steady_clock::now()+chrono::milliseconds(rem2);
                    if(dl2>sl) dl2=sl;
                    vector<vector<int>> res;
                    if(constrainedSA(seed, aim2, dl2, res)){
                        if(verifyGrid(res)){ best=res; bestK=aim2; improved=true; failsAt=0; }
                    }
                }
            }
            if(!improved){
                if(failsAt>=5) break;
            }
        }

        // === Ultimate fallbacks ===
        if(!verifyGrid(best)){
            if(!walk.empty() && (int)walk.size()<=240){
                int L=(int)walk.size();
                best.assign(L,vector<int>(L));
                for(int r=0;r<L;r++) for(int c=0;c<L;c++) best[r][c]=walk[c];
                bestK=L;
            }
        }
        if(!verifyGrid(best)){
            bestK=min(240,max(1,(int)ceil(sqrt(N))));
            best.assign(bestK,vector<int>(bestK));
            for(int r=0;r<bestK;r++) for(int c=0;c<bestK;c++) best[r][c]=1+(r+c)%N;
        }

        // Output
        int K=(int)best.size();
        printf("%d\n",K);
        for(int j=0;j<K;j++) printf("%d%c",K,j==K-1?'\n':' ');
        for(int i=0;i<K;i++){
            for(int j=0;j<K;j++) printf("%d%c",best[i][j],j==K-1?'\n':' ');
        }
    }
    return 0;
}
