#include <bits/stdc++.h>
using namespace std;

struct Inst {
    int J, M, N;
    vector<vector<int>> mach;
    vector<vector<long long>> p;
    vector<vector<int>> pos;
    vector<vector<long long>> rem, tail, head;
    vector<long long> mload, jload;
};

struct OpCand { int j,k,m; long long p,s,c; };

long long eval_makespan(const Inst& I, const vector<vector<int>>& seq) {
    int J=I.J, M=I.M, N=I.N;
    vector<vector<int>> adj(N);
    vector<int> indeg(N,0);
    vector<long long> w(N), dist(N,0);
    auto id=[&](int j,int k){return j*M+k;};
    for(int j=0;j<J;j++) for(int k=0;k<M;k++) {
        int u=id(j,k); w[u]=I.p[j][k];
        if(k+1<M){ adj[u].push_back(id(j,k+1)); indeg[id(j,k+1)]++; }
    }
    for(int m=0;m<M;m++) for(int a=0;a+1<J;a++) {
        int j1=seq[m][a], j2=seq[m][a+1];
        int u=id(j1,I.pos[j1][m]), v=id(j2,I.pos[j2][m]);
        adj[u].push_back(v); indeg[v]++;
    }
    deque<int> q;
    for(int u=0;u<N;u++) if(!indeg[u]) { dist[u]=w[u]; q.push_back(u); }
    int cnt=0; long long C=0;
    while(!q.empty()){
        int u=q.front(); q.pop_front(); cnt++; C=max(C,dist[u]);
        for(int v: adj[u]){ if(dist[v]<dist[u]+w[v]) dist[v]=dist[u]+w[v]; if(--indeg[v]==0) q.push_back(v); }
    }
    if(cnt!=N) return LLONG_MAX/4;
    return C;
}

struct Rule {
    // Smaller score is preferred.  If gtMode is true, the Giffler-Thompson conflict machine is used.
    array<double,10> w;
    bool gtMode;
};

vector<vector<int>> build_schedule(const Inst& I, const Rule& R, mt19937_64 &rng) {
    int J=I.J, M=I.M;
    vector<int> next(J,0);
    vector<long long> jr(J,0), mr(M,0);
    vector<vector<int>> seq(M); for(int m=0;m<M;m++) seq[m].reserve(J);
    vector<double> noise(J);
    for(int j=0;j<J;j++) noise[j]=uniform_real_distribution<double>(0.0,1.0)(rng);

    for(int step=0; step<I.N; step++) {
        vector<OpCand> avail; avail.reserve(J);
        long long bestC = LLONG_MAX; int critM = -1;
        for(int j=0;j<J;j++) if(next[j] < M) {
            int k=next[j], m=I.mach[j][k]; long long pp=I.p[j][k];
            long long s=max(jr[j], mr[m]), c=s+pp;
            avail.push_back({j,k,m,pp,s,c});
            if(c < bestC || (c==bestC && pp > (critM<0? -1:0))) { bestC=c; critM=m; }
        }
        const vector<OpCand>* candPtr = &avail;
        vector<OpCand> conflict;
        if(R.gtMode) {
            for(auto &o: avail) if(o.m==critM && o.s < bestC) conflict.push_back(o);
            if(!conflict.empty()) candPtr=&conflict;
        }
        const auto &cand=*candPtr;
        int bi=0; double bs=1e300;
        for(int i=0;i<(int)cand.size();i++) {
            const auto &o=cand[i];
            double sc=0;
            sc += R.w[0] * (double)o.p;
            sc += R.w[1] * (double)I.rem[o.j][o.k];
            sc += R.w[2] * (double)I.tail[o.j][o.k];
            sc += R.w[3] * (double)o.s;
            sc += R.w[4] * (double)o.c;
            sc += R.w[5] * (double)I.mload[o.m];
            sc += R.w[6] * (double)(M-o.k);
            sc += R.w[7] * (double)I.head[o.j][o.k];
            sc += R.w[8] * (double)I.jload[o.j];
            sc += R.w[9] * noise[o.j];
            // stable deterministic tie breakers
            sc += 1e-9 * o.j + 1e-12 * o.k;
            if(sc < bs) { bs=sc; bi=i; }
        }
        auto o=cand[bi];
        long long st=max(jr[o.j], mr[o.m]), en=st+o.p;
        jr[o.j]=en; mr[o.m]=en; seq[o.m].push_back(o.j); next[o.j]++;
    }
    return seq;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Inst I; if(!(cin>>I.J>>I.M)) return 0; I.N=I.J*I.M;
    I.mach.assign(I.J, vector<int>(I.M)); I.p.assign(I.J, vector<long long>(I.M));
    I.pos.assign(I.J, vector<int>(I.M)); I.mload.assign(I.M,0); I.jload.assign(I.J,0);
    for(int j=0;j<I.J;j++) for(int k=0;k<I.M;k++){
        cin>>I.mach[j][k]>>I.p[j][k]; I.pos[j][I.mach[j][k]]=k;
        I.mload[I.mach[j][k]] += I.p[j][k]; I.jload[j] += I.p[j][k];
    }
    I.rem.assign(I.J, vector<long long>(I.M+1,0));
    I.tail.assign(I.J, vector<long long>(I.M,0)); I.head.assign(I.J, vector<long long>(I.M,0));
    for(int j=0;j<I.J;j++){
        long long h=0; for(int k=0;k<I.M;k++){ I.head[j][k]=h; h+=I.p[j][k]; }
        for(int k=I.M-1;k>=0;k--) I.rem[j][k]=I.rem[j][k+1]+I.p[j][k];
        for(int k=0;k<I.M;k++) I.tail[j][k]=I.rem[j][k+1];
    }

    vector<Rule> rules;
    auto add=[&](initializer_list<double> ws, bool gt){ Rule r; r.w.fill(0); int i=0; for(double x:ws) r.w[i++]=x; r.gtMode=gt; rules.push_back(r); };
    // Classic active-schedule dispatch priorities and bottleneck-friendly head/tail variants.
    add({ 1,0,0,0,0,0,0,0,0,0}, true);    // shortest processing time
    add({-1,0,0,0,0,0,0,0,0,0}, true);    // longest processing time
    add({0,-1,0,0,0,0,0,0,0,0}, true);    // most work remaining
    add({0,0,-1,0,0,0,0,0,0,0}, true);    // largest tail
    add({0,0,0,0,1,0,0,0,0,0}, true);     // earliest completion
    add({1,0,-1,0,0,0,0,0,0,0}, true);
    add({0,-1,-1,0,0,0,0,0,0,0}, true);
    add({1,0,-2,0,0,0,0,0,0,0}, true);
    add({0,0,-1,0,0,0,0,1,0,0}, true);    // tail with small head preference
    add({1,-0.3,-1,0,0,0,0,0,0,0}, true);
    // Serial generation alternatives; sometimes better on random instances.
    add({0,-1,0,0,1,0,0,0,0,0}, false);
    add({1,0,-1,0,1,0,0,0,0,0}, false);
    add({0,0,-1,0,1,0,0,0,0,0}, false);

    mt19937_64 rng(123456789ULL + I.J*1009ULL + I.M*9176ULL);
    vector<vector<int>> bestSeq; long long best=LLONG_MAX/4;
    for(const auto &r: rules){
        auto seq=build_schedule(I,r,rng); long long v=eval_makespan(I,seq);
        if(v<best){ best=v; bestSeq=move(seq); }
    }

    auto start=chrono::steady_clock::now();
    long long ops = max(1, I.N);
    double limit = (ops<=700 ? 0.92 : 0.88);
    int iter=0;
    while(true){
        if((iter & 31)==0){
            double t=chrono::duration<double>(chrono::steady_clock::now()-start).count();
            if(t>limit) break;
        }
        Rule r; r.gtMode = (rng()%100) < 82; r.w.fill(0);
        // Random weighted priority rules around useful JSSP features.
        double scaleP = 1.0 / max(1.0, (double)(accumulate(I.jload.begin(), I.jload.end(), 0LL) / max(1,I.N)));
        (void)scaleP;
        uniform_real_distribution<double> U(-1.0,1.0);
        r.w[0] = U(rng) * 1.5;        // processing time
        r.w[1] = U(rng) * 1.2 - 0.35; // remaining work, biased to prefer large
        r.w[2] = U(rng) * 1.6 - 0.45; // tail, useful on bottlenecks
        r.w[3] = U(rng) * 0.35;
        r.w[4] = U(rng) * 0.55 + 0.15; // mild earliest-completion bias
        r.w[5] = U(rng) * 0.12;
        r.w[6] = U(rng) * 0.4 - 0.15;
        r.w[7] = U(rng) * 0.35;
        r.w[8] = U(rng) * 0.25;
        r.w[9] = U(rng) * 1e-3;
        auto seq=build_schedule(I,r,rng); long long v=eval_makespan(I,seq);
        if(v<best){ best=v; bestSeq=move(seq); }
        iter++;
    }

    if(bestSeq.empty()) { // impossible fallback, but keep output valid
        bestSeq.assign(I.M, vector<int>());
        for(int m=0;m<I.M;m++) for(int j=0;j<I.J;j++) bestSeq[m].push_back(j);
    }
    for(int m=0;m<I.M;m++){
        for(int i=0;i<I.J;i++){ if(i) cout << ' '; cout << bestSeq[m][i]; }
        cout << '\n';
    }
    return 0;
}
