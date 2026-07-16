#include <bits/stdc++.h>
using namespace std;

struct Inst {
    int J, M;
    vector<vector<int>> mach;
    vector<vector<long long>> p;
    vector<vector<long long>> rem; // from k to end
    vector<long long> mload;
};

struct Candidate {
    vector<vector<int>> order;
    long long cmax;
};

static uint64_t splitmix64(uint64_t x){
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static inline double priority_value(const Inst& in, int j, int k, const vector<long long>& jr,
                                    const vector<long long>& mr, const array<double,10>& w,
                                    double noise) {
    int m=in.mach[j][k];
    long long p=in.p[j][k];
    long long est=max(jr[j], mr[m]);
    long long waitM = max(0LL, jr[j]-mr[m]);
    long long waitJ = max(0LL, mr[m]-jr[j]);
    long long tail = in.rem[j][k] - p;
    double val = 0.0;
    val += w[0] * (double)est;
    val += w[1] * (double)(est + p);
    val += w[2] * (double)p;
    val += w[3] * (double)in.rem[j][k];
    val += w[4] * (double)tail;
    val += w[5] * (double)jr[j];
    val += w[6] * (double)mr[m];
    val += w[7] * (double)waitM;
    val += w[8] * (double)waitJ;
    val += w[9] * (double)in.mload[m];
    return val + noise;
}

Candidate build_schedule(const Inst& in, const array<double,10>& w, bool rnd, mt19937_64 &rng, bool giffler) {
    int J=in.J, M=in.M, N=J*M;
    vector<int> nextOp(J,0);
    vector<long long> jr(J,0), mr(M,0);
    vector<vector<int>> ord(M); for(int m=0;m<M;m++) ord[m].reserve(J);
    long long cmax=0;
    uniform_real_distribution<double> ud(-1.0,1.0);
    for(int step=0; step<N; ++step){
        int conflictMachine=-1;
        long long bound=LLONG_MAX/4;
        if(giffler){
            for(int j=0;j<J;j++) if(nextOp[j]<M){
                int k=nextOp[j], m=in.mach[j][k];
                long long ft=max(jr[j], mr[m])+in.p[j][k];
                if(ft<bound){ bound=ft; conflictMachine=m; }
            }
        }
        int best=-1;
        double bestVal=1e300;
        for(int j=0;j<J;j++) if(nextOp[j]<M){
            int k=nextOp[j], m=in.mach[j][k];
            if(giffler && !(m==conflictMachine && max(jr[j],mr[m]) < bound)) continue;
            double noise = rnd ? ud(rng) * (0.0005 * (double)(1 + in.rem[j][0])) : 0.0;
            double val = priority_value(in,j,k,jr,mr,w,noise);
            if(val < bestVal - 1e-9 || (fabs(val-bestVal)<=1e-9 && (best<0 || j<best))){
                bestVal=val; best=j;
            }
        }
        int j=best, k=nextOp[j], m=in.mach[j][k];
        long long st=max(jr[j], mr[m]);
        long long ft=st+in.p[j][k];
        jr[j]=mr[m]=ft;
        cmax=max(cmax,ft);
        ord[m].push_back(j);
        nextOp[j]++;
    }
    return {ord,cmax};
}

// Recompute the longest path for the printed machine orders.  This is also a
// safety net against accidental cyclic output (serial SGS orders are acyclic).
long long eval_order(const Inst& in, const vector<vector<int>>& ord){
    int J=in.J,M=in.M,N=J*M;
    vector<vector<int>> id(J, vector<int>(M));
    vector<int> nodeJob(N), nodeK(N), nodeMach(N);
    for(int j=0,c=0;j<J;j++) for(int k=0;k<M;k++,c++){
        id[j][k]=c; nodeJob[c]=j; nodeK[c]=k; nodeMach[c]=in.mach[j][k];
    }
    vector<vector<int>> adj(N);
    vector<int> indeg(N,0);
    for(int j=0;j<J;j++) for(int k=0;k+1<M;k++){
        int a=id[j][k], b=id[j][k+1]; adj[a].push_back(b); indeg[b]++;
    }
    for(int m=0;m<M;m++) for(int t=0;t+1<J;t++){
        int a=ord[m][t], b=ord[m][t+1];
        int ka=-1,kb=-1;
        for(int k=0;k<M;k++){ if(in.mach[a][k]==m) ka=k; if(in.mach[b][k]==m) kb=k; }
        int u=id[a][ka], v=id[b][kb]; adj[u].push_back(v); indeg[v]++;
    }
    queue<int> q; vector<long long> dist(N,0);
    for(int i=0;i<N;i++) if(!indeg[i]) q.push(i);
    int seen=0; long long ans=0;
    while(!q.empty()){
        int u=q.front(); q.pop(); seen++;
        int j=nodeJob[u], k=nodeK[u];
        long long fin=dist[u]+in.p[j][k]; ans=max(ans,fin);
        for(int v:adj[u]){ if(dist[v]<fin) dist[v]=fin; if(--indeg[v]==0) q.push(v); }
    }
    if(seen<N) return LLONG_MAX/4;
    return ans;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Inst in;
    if(!(cin>>in.J>>in.M)) return 0;
    int J=in.J,M=in.M;
    in.mach.assign(J, vector<int>(M));
    in.p.assign(J, vector<long long>(M));
    in.mload.assign(M,0);
    uint64_t seed=1469598103934665603ULL;
    for(int j=0;j<J;j++) for(int k=0;k<M;k++){
        cin>>in.mach[j][k]>>in.p[j][k];
        in.mload[in.mach[j][k]] += in.p[j][k];
        seed ^= (uint64_t)(in.mach[j][k]+1009*(j+1)+9176*(k+1)) + splitmix64((uint64_t)in.p[j][k]);
        seed *= 1099511628211ULL;
    }
    in.rem.assign(J, vector<long long>(M+1,0));
    for(int j=0;j<J;j++) for(int k=M-1;k>=0;k--) in.rem[j][k]=in.rem[j][k+1]+in.p[j][k];
    mt19937_64 rng(splitmix64(seed));

    vector<array<double,10>> rules;
    auto add=[&](initializer_list<double> l){ array<double,10> a{}; int i=0; for(double x:l) if(i<10) a[i++]=x; rules.push_back(a); };
    // Classical dispatching rules expressed as linear priorities.
    add({1,0,0,0,0,0,0,0,0,0});                 // earliest start
    add({0,1,0,0,0,0,0,0,0,0});                 // earliest finish
    add({0,0,1,0,0,0,0,0,0,0});                 // SPT
    add({0,0,-1,0,0,0,0,0,0,0});                // LPT
    add({0,0,0,-1,0,0,0,0,0,0});                // longest remaining work
    add({0,0,0,1,0,0,0,0,0,0});                 // shortest remaining work
    add({1,0,0,-0.7,0,0,0,0,0,0});
    add({1,0,0,0,-1,0,0,0,0,0});
    add({0,1,0,-0.8,0,0,0,0,0,0});
    add({1,0,-0.3,-0.7,0,0,0,0,0,0});
    add({1,0,0.4,-1.0,0,0,0,0,0,0});
    add({0.7,0,0,-1.0,0,0,0,0,0,-0.02});
    add({0.7,0,0,-1.0,0,0,0,0,0,0.02});
    add({1,0,0,0,0,0,0,1,0,0});                 // avoid machine idle caused by future release
    add({1,0,0,-0.5,0,0,0,1,0,0});
    add({1,0,0,-0.5,0,0,0,0,1,0});
    add({0,1,0,-0.4,0,0,0,0.5,0,0});
    add({0.5,0.5,0.2,-0.8,0,0,0,0.3,0,0});

    Candidate best; best.cmax = LLONG_MAX/4;
    for(auto &r: rules){
        for(bool gt: {false,true}){
            Candidate c=build_schedule(in,r,false,rng,gt);
            c.cmax=eval_order(in,c.order);
            if(c.cmax<best.cmax) best=std::move(c);
        }
    }

    auto start=chrono::steady_clock::now();
    int iter=0;
    while(true){
        if((iter & 15)==0){
            double elapsed=chrono::duration<double>(chrono::steady_clock::now()-start).count();
            if(elapsed>0.82) break;
        }
        array<double,10> w{};
        // Biased random family around effective JSSP priorities.
        w[0]=uniform_real_distribution<double>(0.0,1.5)(rng);
        w[1]=uniform_real_distribution<double>(0.0,1.0)(rng);
        w[2]=uniform_real_distribution<double>(-1.2,1.2)(rng);
        w[3]=uniform_real_distribution<double>(-1.8,0.6)(rng);
        w[4]=uniform_real_distribution<double>(-1.0,0.5)(rng);
        w[5]=uniform_real_distribution<double>(-0.2,0.5)(rng);
        w[6]=uniform_real_distribution<double>(-0.2,0.5)(rng);
        w[7]=uniform_real_distribution<double>(-0.2,1.2)(rng);
        w[8]=uniform_real_distribution<double>(-0.2,0.8)(rng);
        w[9]=uniform_real_distribution<double>(-0.05,0.05)(rng);
        Candidate c=build_schedule(in,w,true,rng,(iter&1));
        c.cmax=eval_order(in,c.order);
        if(c.cmax<best.cmax) best=std::move(c);
        iter++;
    }

    // Fallback should never trigger, but ensure exact output shape.
    if(best.order.empty()){
        best.order.assign(M, vector<int>());
        for(int m=0;m<M;m++) for(int j=0;j<J;j++) best.order[m].push_back(j);
    }
    for(int m=0;m<M;m++){
        for(int i=0;i<J;i++){
            if(i) cout << ' ';
            cout << best.order[m][i];
        }
        cout << '\n';
    }
    return 0;
}
