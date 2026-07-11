// IOI 2025 "World Map" — improved candidate v2 (min-conflicts + compression)
#include <bits/stdc++.h>
using namespace std;
static double timeLimitSec=0.88;
static chrono::steady_clock::time_point startClock;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-startClock).count(); }
static const int MAXK = 240;
static int gN;
static vector<vector<char>> gAdj;
static vector<pair<int,int>> gEdges;
static set<pair<int,int>> gEdgeSet;
static bool validate(const vector<vector<int>>& C){
    int K = (int)C.size();
    if(K<1 || K>MAXK) return false;
    for(auto& row: C){ if((int)row.size()!=K) return false; for(int v: row) if(v<1||v>gN) return false; }
    static const int dy[4]={1,0,-1,0}, dx[4]={0,1,0,-1};
    set<pair<int,int>> real;
    for(int i=0;i<K;i++) for(int j=0;j<K;j++){
        int a=C[i][j]-1;
        for(int l=0;l<4;l++){
            int ni=i+dy[l], nj=j+dx[l];
            if(ni<0||nj<0||ni>=K||nj>=K) continue;
            int b=C[ni][nj]-1;
            if(a==b) continue;
            if(!gEdgeSet.count({min(a,b),max(a,b)})) return false;
            real.insert({min(a,b),max(a,b)});
        }
    }
    if(real.size()!=gEdgeSet.size()) return false;
    return true;
}
static vector<vector<int>> buildReference2N(){
    int N=gN;
    if(N==1) return {{1}};
    vector<vector<int>> G(N);
    for(auto&e:gEdges){ G[e.first].push_back(e.second); G[e.second].push_back(e.first); }
    vector<char> vis(N,0); vector<int> depth(N,0), tour, RA, RB;
    function<void(int)> dfs=[&](int x){
        vis[x]=1; tour.push_back(x);
        for(int i:G[x]){
            if(!vis[i]){ depth[i]=depth[x]+1; dfs(i); tour.push_back(x); }
            else if(depth[x]<depth[i]){ RA.push_back(x); RB.push_back(i); }
        }
    };
    dfs(0);
    int L=(int)tour.size();
    vector<int> rnk(N,-1), holder(N,-1);
    for(int i=0;i<L;i++){ int d=min(i,(L-1)-i); if(rnk[tour[i]]<d){ rnk[tour[i]]=d; holder[tour[i]]=i; } }
    vector<vector<int>> H(N);
    for(size_t i=0;i<RA.size();i++){ if(rnk[RA[i]]<rnk[RB[i]]) swap(RA[i],RB[i]); H[RA[i]].push_back(RB[i]); }
    int S=2*N;
    vector<vector<int>> ans(S, vector<int>(S,0));
    int cur=0;
    for(int i=0;i<L;i++){
        if(i==holder[tour[i]]){
            int pos=0;
            for(int j=0;j<S;j++){
                int ya=cur-j;      if(0<=ya&&ya<S) ans[j][ya]=tour[i];
                int yb=cur+1-j;    if(0<=yb&&yb<S){ if(pos<(int)H[tour[i]].size()) ans[j][yb]=H[tour[i]][pos++]; else ans[j][yb]=tour[i]; }
                int yc=cur+2-j;    if(0<=yc&&yc<S) ans[j][yc]=tour[i];
            }
            cur+=3;
        } else {
            for(int j=0;j<S;j++){ int ya=cur-j; if(0<=ya&&ya<S) ans[j][ya]=tour[i]; }
            cur+=1;
        }
    }
    for(int i=0;i<S;i++) for(int j=0;j<S;j++){ if(i<2*N&&j<2*N) continue; if(j>0) ans[i][j]=ans[i][j-1]; }
    for(int i=0;i<S;i++) for(int j=0;j<S;j++) if(ans[i][j]==0) ans[i][j]= (j>0?ans[i][j-1]:(i>0?ans[i-1][j]:0));
    for(int i=0;i<S;i++) for(int j=0;j<S;j++) ans[i][j]++;
    return ans;
}
static vector<vector<int>> buildAntiDiag(const vector<int>& walk, const vector<vector<int>>& hostMap){
    int N=gN, L=(int)walk.size();
    if(L==0) return {};
    for(int K=1;K<=MAXK;K++){
        int cap=K;
        vector<pair<int,vector<int>>> S;
        S.reserve(2*K);
        vector<char> hosted(N,0);
        for(int k=0;k<L;k++){
            int v=walk[k];
            S.push_back({v, {}});
            if(!hosted[v] && !hostMap[v].empty()){
                hosted[v]=1;
                const vector<int>& lst=hostMap[v];
                for(size_t p=0;p<lst.size();p+=cap){
                    vector<int> chunk(lst.begin()+p, lst.begin()+min(lst.size(),p+cap));
                    S.push_back({v, chunk});
                    S.push_back({v, {}});
                }
            }
        }
        int len=(int)S.size();
        if(2*K-1 < len) continue;
        int offset=(2*K-1-len)/2;
        bool ok=true;
        for(int t=0;t<len && ok;t++){
            if(S[t].second.empty()) continue;
            int d=t+offset;
            int Ld=min(d+1, min(K, 2*K-1-d));
            if((int)S[t].second.size()>Ld) ok=false;
        }
        if(!ok) continue;
        vector<int> diagVertex(2*K-1);
        vector<vector<int>> diagChunk(2*K-1);
        int frontV=S.front().first, backV=S.back().first;
        for(int d=0;d<2*K-1;d++){
            int t=d-offset;
            if(t<0){ diagVertex[d]=frontV; }
            else if(t>=len){ diagVertex[d]=backV; }
            else { diagVertex[d]=S[t].first; diagChunk[d]=S[t].second; }
        }
        vector<vector<int>> C(K, vector<int>(K,0));
        for(int d=0;d<2*K-1;d++){
            int i0=max(0,d-K+1), i1=min(d,K-1);
            int v=diagVertex[d];
            const vector<int>& chunk=diagChunk[d];
            int ci=0;
            for(int i=i0;i<=i1;i++){
                int j=d-i;
                if(ci<(int)chunk.size()) C[i][j]=chunk[ci++]+1; else C[i][j]=v+1;
            }
        }
        if(validate(C)) return C;
        return {};
    }
    return {};
}
static vector<vector<int>> gdist, gnext;
static void computeSP(){
    int N=gN;
    gdist.assign(N, vector<int>(N, 1e9));
    gnext.assign(N, vector<int>(N, -1));
    for(int i=0;i<N;i++) gdist[i][i]=0;
    for(auto&e:gEdges){ gdist[e.first][e.second]=gdist[e.second][e.first]=1;
        gnext[e.first][e.second]=e.second; gnext[e.second][e.first]=e.first; }
    for(int k=0;k<N;k++)for(int i=0;i<N;i++)if(gdist[i][k]<1e9)for(int j=0;j<N;j++)
        if(gdist[i][k]+gdist[k][j]<gdist[i][j]){ gdist[i][j]=gdist[i][k]+gdist[k][j]; gnext[i][j]=gnext[i][k]; }
}
static vector<int> expandOrder(const vector<int>& order){
    vector<int> walk;
    for(int idx=0; idx<(int)order.size(); idx++){
        int v=order[idx];
        if(walk.empty()){ walk.push_back(v); continue; }
        int cur=walk.back();
        if(gdist[cur][v]>=(int)1e9) return {};
        while(cur!=v){ cur=gnext[cur][v]; walk.push_back(cur); }
    }
    return walk;
}
static bool buildHostMap(const vector<int>& walk, vector<vector<int>>& hostMap){
    int N=gN;
    vector<char> inWalk(N,0); for(int v:walk) inWalk[v]=1;
    set<pair<int,int>> covered;
    for(int i=0;i+1<(int)walk.size();i++){ int a=walk[i],b=walk[i+1]; covered.insert({min(a,b),max(a,b)}); }
    vector<pair<int,int>> unc;
    for(auto&e:gEdges) if(!covered.count(e)){
        if(!inWalk[e.first] && !inWalk[e.second]) return false;
        unc.push_back(e);
    }
    hostMap.assign(N, {});
    vector<char> removed(unc.size(),0);
    int remaining=(int)unc.size();
    while(remaining>0){
        vector<int> deg(N,0);
        for(size_t i=0;i<unc.size();i++) if(!removed[i]){
            if(inWalk[unc[i].first])  deg[unc[i].first]++;
            if(inWalk[unc[i].second]) deg[unc[i].second]++;
        }
        int v=max_element(deg.begin(),deg.end())-deg.begin();
        for(size_t i=0;i<unc.size();i++) if(!removed[i] && ((inWalk[unc[i].first]&&unc[i].first==v)||(inWalk[unc[i].second]&&unc[i].second==v))){
            int y=unc[i].first^unc[i].second^v; hostMap[v].push_back(y); removed[i]=1; remaining--; }
    }
    return true;
}
static vector<vector<int>> greedyFill(int K, mt19937& rng, const set<pair<int,int>>& want){
    int N=gN;
    vector<vector<int>> C(K, vector<int>(K,-1));
    auto code=[&](int a,int b){ if(a>b) swap(a,b); return a*N+b; };
    unordered_set<int> uncov;
    for(auto&e:want) uncov.insert(code(e.first,e.second));
    vector<char> present(N,0);
    for(int i=0;i<K;i++) for(int j=0;j<K;j++){
        int a = (i>0)? C[i-1][j] : -1;
        int b = (j>0)? C[i][j-1] : -1;
        static vector<int> cand; cand.clear();
        for(int c=0;c<N;c++){
            bool okA = (a<0) || c==a || gAdj[a][c];
            if(!okA) continue;
            bool okB = (b<0) || c==b || gAdj[b][c];
            if(!okB) continue;
            cand.push_back(c);
        }
        if(cand.empty()) return {};
        int bestScore=INT_MIN; static vector<int> best; best.clear();
        for(int c:cand){
            int sc=0;
            if(a>=0 && c!=a && uncov.count(code(a,c))) sc+=4;
            if(b>=0 && c!=b && uncov.count(code(b,c))) sc+=4;
            if(!present[c]) sc+=3;
            if(sc>bestScore){ bestScore=sc; best.clear(); best.push_back(c);}
            else if(sc==bestScore) best.push_back(c);
        }
        int c=best[rng()%best.size()];
        C[i][j]=c; present[c]=1;
        if(a>=0&&c!=a) uncov.erase(code(a,c));
        if(b>=0&&c!=b) uncov.erase(code(b,c));
    }
    for(auto&row:C) for(auto&v:row) v++;
    if(validate(C)) return C;
    return {};
}
static vector<vector<int>> annealFill(int K, mt19937& rng, double deadline, vector<vector<int>> C){
    int N=gN;
    if((int)C.size()!=K){ C.assign(K, vector<int>(K,0)); for(auto&r:C) for(auto&v:r) v=rng()%N; }
    static const int dy[4]={1,0,-1,0}, dx[4]={0,1,0,-1};
    vector<vector<int>> er(N, vector<int>(N,0));
    vector<int> colc(N,0);
    long forb=0; int missE=gEdgeSet.size(), missC=0;
    auto addPair=[&](int a,int b,int s){
        if(a==b) return;
        if(gAdj[a][b]){
            int before=er[a][b];
            er[a][b]+=s; er[b][a]+=s;
            if(before==0 && er[a][b]>0) missE--;
            if(before>0 && er[a][b]==0) missE++;
        } else { forb+=s; }
    };
    auto cost=[&](){ return forb*1000L + (long)missE*25 + (long)missC*25; };
    for(int c=0;c<N;c++) colc[c]=0;
    for(int i=0;i<K;i++)for(int j=0;j<K;j++) colc[C[i][j]]++;
    for(int c=0;c<N;c++) if(colc[c]==0) missC++;
    for(int i=0;i<K;i++)for(int j=0;j<K;j++){
        if(i+1<K) addPair(C[i][j],C[i+1][j],+1);
        if(j+1<K) addPair(C[i][j],C[i][j+1],+1);
    }
    vector<pair<int,int>> bad;
    auto refreshBad=[&](){ bad.clear();
        for(int i=0;i<K;i++)for(int j=0;j<K;j++){ int a=C[i][j];
            for(int l=0;l<4;l++){ int ni=i+dy[l],nj=j+dx[l]; if(ni<0||nj<0||ni>=K||nj>=K)continue;
                int b=C[ni][nj]; if(a!=b&&!gAdj[a][b]){ bad.push_back({i,j}); break; } } } };
    double T=1.5; long iters=0;
    long bestCost=cost();
    while(cost()>0){
        if((iters&255)==0){ if(elapsed()>=deadline) break; refreshBad(); T*=0.999; if(T<0.05) T=0.05; }
        iters++;
        int i,j;
        int mode=rng()%100;
        if(!bad.empty() && mode<80){ auto&p=bad[rng()%bad.size()]; i=p.first; j=p.second; }
        else { i=rng()%K; j=rng()%K; }
        int old=C[i][j];
        int nw;
        int rr=rng()%100;
        if(rr<75){
            int bestc=old, bestv=INT_MAX;
            for(int c=0;c<N;c++){
                int lf=0;
                for(int l=0;l<4;l++){ int ni=i+dy[l],nj=j+dx[l]; if(ni<0||nj<0||ni>=K||nj>=K)continue;
                    int b=C[ni][nj]; if(c!=b&&!gAdj[c][b]) lf++; }
                int score=lf*8;
                if(colc[c]==0) score-=1;
                if(score<bestv || (score==bestv && (rng()&1))){ bestv=score; bestc=c; }
            }
            nw=bestc;
        } else if(rr<88 && !gEdges.empty()){
            auto&e=gEdges[rng()%gEdges.size()]; nw = (rng()&1)? e.first : e.second;
        } else nw=rng()%N;
        if(nw==old) continue;
        long f0=forb; int me0=missE, mc0=missC;
        for(int l=0;l<4;l++){ int ni=i+dy[l],nj=j+dx[l]; if(ni<0||nj<0||ni>=K||nj>=K)continue; addPair(old,C[ni][nj],-1); }
        if(--colc[old]==0) missC++;
        if(colc[nw]==0){ missC--; } colc[nw]++;
        for(int l=0;l<4;l++){ int ni=i+dy[l],nj=j+dx[l]; if(ni<0||nj<0||ni>=K||nj>=K)continue; addPair(nw,C[ni][nj],+1); }
        long newCost=cost();
        long oldCost=f0*1000L+(long)me0*25+(long)mc0*25;
        long d=newCost-oldCost;
        bool accept = d<=0 || (exp(-d/(T*(double)K))*4294967296.0 > rng());
        if(accept){ C[i][j]=nw; if(newCost<bestCost) bestCost=newCost; }
        else {
            for(int l=0;l<4;l++){ int ni=i+dy[l],nj=j+dx[l]; if(ni<0||nj<0||ni>=K||nj>=K)continue; addPair(nw,C[ni][nj],-1); }
            if(--colc[nw]==0) missC++;
            if(colc[old]==0){ missC--; } colc[old]++;
            for(int l=0;l<4;l++){ int ni=i+dy[l],nj=j+dx[l]; if(ni<0||nj<0||ni>=K||nj>=K)continue; addPair(old,C[ni][nj],+1); }
        }
    }
    if(cost()!=0) return {};
    vector<vector<int>> out(K, vector<int>(K));
    for(int i=0;i<K;i++)for(int j=0;j<K;j++) out[i][j]=C[i][j]+1;
    if(validate(out)) return out;
    return {};
}
static vector<vector<int>> tryCompress(const vector<vector<int>>& best1, mt19937& rng, double deadline){
    int K=(int)best1.size(); if(K<=2) return {};
    vector<vector<int>> g(K, vector<int>(K));
    for(int i=0;i<K;i++)for(int j=0;j<K;j++) g[i][j]=best1[i][j]-1;
    auto seamRow=[&](int r)->int{ if(r==0||r==K-1) return 0; int s=0;
        for(int j=0;j<K;j++){ int a=g[r-1][j],b=g[r+1][j]; if(a!=b&&!gAdj[a][b]) s++; } return s; };
    auto seamCol=[&](int c)->int{ if(c==0||c==K-1) return 0; int s=0;
        for(int i=0;i<K;i++){ int a=g[i][c-1],b=g[i][c+1]; if(a!=b&&!gAdj[a][b]) s++; } return s; };
    vector<pair<int,int>> rows,cols;
    for(int r=0;r<K;r++) rows.push_back({seamRow(r),r});
    for(int c=0;c<K;c++) cols.push_back({seamCol(c),c});
    sort(rows.begin(),rows.end()); sort(cols.begin(),cols.end());
    int tries=0;
    while(elapsed()<deadline && tries<40){
        int r,c;
        if(tries<6){ r=rows[tries%rows.size()].second; c=cols[tries%cols.size()].second; }
        else { r=rows[rng()%min((size_t)6,rows.size())].second; c=cols[rng()%min((size_t)6,cols.size())].second; }
        tries++;
        vector<vector<int>> seed; seed.reserve(K-1);
        for(int i=0;i<K;i++){ if(i==r) continue; vector<int> row; row.reserve(K-1);
            for(int j=0;j<K;j++){ if(j==c) continue; row.push_back(g[i][j]); } seed.push_back(row); }
        double dl=min(deadline, elapsed()+0.06);
        auto res=annealFill(K-1, rng, dl, seed);
        if(!res.empty()) return res;
    }
    return {};
}
vector<vector<int>> create_map(int N, int M, vector<int> A, vector<int> B){
    startClock=chrono::steady_clock::now();
    gN=N;
    gAdj.assign(N, vector<char>(N,0));
    gEdges.clear(); gEdgeSet.clear();
    for(int i=0;i<M;i++){
        int a=A[i]-1, b=B[i]-1;
        if(a==b) continue;
        if(!gAdj[a][b]){ gAdj[a][b]=gAdj[b][a]=1; gEdges.push_back({min(a,b),max(a,b)}); gEdgeSet.insert({min(a,b),max(a,b)}); }
    }
    if(N==1) return {{1}};
    vector<vector<int>> best=buildReference2N();
    if(!validate(best)){ best.assign(2*N, vector<int>(2*N,1)); }
    int bestK=(int)best.size();
    computeSP();
    mt19937 rng(12345);
    auto tryOrderOn=[&](vector<int> order){
        vector<int> walk=expandOrder(order);
        if((int)walk.size()==0) return;
        vector<vector<int>> hostMap;
        if(!buildHostMap(walk, hostMap)) return;
        vector<vector<int>> C=buildAntiDiag(walk, hostMap);
        if(!C.empty() && (int)C.size()<bestK){ best=C; bestK=(int)C.size(); }
    };
    auto nearestNeighbourOn=[&](int start, const vector<int>& dom)->vector<int>{
        vector<int> order; vector<char> inDom(N,0); for(int v:dom) inDom[v]=1;
        vector<char> used(N,0); int cur=start; order.push_back(cur); used[cur]=1;
        for(size_t step=1;step<dom.size();step++){
            int bn=-1,bd=INT_MAX;
            for(int v:dom) if(!used[v] && gdist[cur][v]<bd){ bd=gdist[cur][v]; bn=v; }
            if(bn<0) return {};
            used[bn]=1; order.push_back(bn); cur=bn;
        }
        return order;
    };
    auto twoOpt=[&](vector<int> o)->vector<int>{
        bool imp=true; int guard=0;
        while(imp && guard++<20 && elapsed()<timeLimitSec){
            imp=false;
            for(int i=0;i+1<(int)o.size();i++) for(int k=i+1;k<(int)o.size();k++){
                int c=o[k], d=(k+1<(int)o.size())?o[k+1]:-1;
                long before=gdist[o[i]][o[i+1]] + (d>=0?gdist[c][d]:0);
                long after =gdist[o[i]][c] + (d>=0?gdist[o[i+1]][d]:0);
                if(after<before){ reverse(o.begin()+i+1, o.begin()+k+1); imp=true; }
            }
        }
        return o;
    };
    vector<int> deg(N,0); for(auto&e:gEdges){ deg[e.first]++; deg[e.second]++; }
    vector<vector<int>> domains;
    { vector<int> all(N); iota(all.begin(),all.end(),0); domains.push_back(all); }
    for(int t=2;t<=4;t++){
        vector<int> dom; for(int v=0;v<N;v++) if(deg[v]>=t) dom.push_back(v);
        if(dom.size()>=1 && dom.size()<(size_t)N) domains.push_back(dom);
    }
    double antiBudget = timeLimitSec*0.28;
    for(size_t di=0; di<domains.size() && elapsed()<antiBudget; di++){
        vector<int>& dom=domains[di];
        if(dom.empty()) continue;
        vector<int> starts=dom; shuffle(starts.begin(), starts.end(), rng);
        for(int si=0; si<(int)starts.size() && elapsed()<antiBudget; si++){
            vector<int> o=nearestNeighbourOn(starts[si], dom);
            if(o.empty()) continue;
            tryOrderOn(o);
            o=twoOpt(o);
            tryOrderOn(o);
        }
        for(int r=0;r<15 && elapsed()<antiBudget;r++){
            vector<int> o=dom; shuffle(o.begin(),o.end(),rng);
            o=twoOpt(o); tryOrderOn(o);
        }
    }
    int lb=max(1,(int)ceil(sqrt((double)N)));
    int cnt=(int)gEdges.size();
    while((long)2*lb*(lb-1) < cnt) lb++;
    double greedyBudget = timeLimitSec*0.30;
    int lo=lb-1, hi=bestK;
    vector<vector<int>> hiGrid=best;
    while(hi-lo>1 && elapsed()<greedyBudget){
        int mid=(lo+hi)/2;
        double slice = elapsed() + max(0.02, (greedyBudget-elapsed())/(log2(max(2,hi-lo))+1));
        bool ok=false; vector<vector<int>> got;
        while(elapsed()<slice && elapsed()<greedyBudget){
            vector<vector<int>> C=greedyFill(mid, rng, gEdgeSet);
            if(!C.empty()){ ok=true; got=C; break; }
        }
        if(ok){ hi=mid; hiGrid=got; } else lo=mid;
    }
    if((int)hiGrid.size()<bestK){ best=hiGrid; bestK=(int)hiGrid.size(); }
    int floorFails=0;
    while(elapsed()<timeLimitSec-0.02 && bestK>lb){
        int target=bestK-1;
        bool improved=false;
        for(int t=0;t<45 && !improved && elapsed()<timeLimitSec;t++){
            vector<vector<int>> C=greedyFill(target, rng, gEdgeSet);
            if(!C.empty()){ best=C; bestK=target; improved=true; }
        }
        if(improved){ floorFails=0; continue; }
        {
            double dl=min(timeLimitSec, elapsed()+0.32);
            vector<vector<int>> C=annealFill(target, rng, dl, {});
            if(!C.empty()){ best=C; bestK=target; improved=true; }
        }
        if(improved){ floorFails=0; continue; }
        {
            vector<vector<int>> C=tryCompress(best, rng, min(timeLimitSec, elapsed()+0.20));
            if(!C.empty() && (int)C.size()==target){ best=C; bestK=target; improved=true; }
        }
        if(improved){ floorFails=0; continue; }
        if(++floorFails>=6) break;
    }
    return best;
}
#ifndef NO_MAIN
int main(){
    int T;
    if(scanf("%d",&T)!=1) return 0;
    for(int t=0;t<T;t++){
        int N,M; if(scanf("%d %d",&N,&M)!=2) return 0;
        vector<int> A(M),B(M);
        for(int i=0;i<M;i++) scanf("%d %d",&A[i],&B[i]);
        vector<vector<int>> C=create_map(N,M,A,B);
        int P=(int)C.size();
        printf("%d\n",P);
        for(int i=0;i<P;i++) printf("%d%c", P, " \n"[i+1==P]);
        printf("\n");
        for(int i=0;i<P;i++) for(int j=0;j<P;j++) printf("%d%c", C[i][j], " \n"[j+1==P]);
        if(t<T-1) printf("\n");
    }
    return 0;
}
#endif
