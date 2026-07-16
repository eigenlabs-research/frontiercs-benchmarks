#include <bits/stdc++.h>
using namespace std;

struct Op { int m; long long p; };
int J, M, N;
vector<vector<Op>> ops;
vector<vector<int>> posOf;
vector<long long> totalWork, machLoad;

struct Rule {
    double est, ect, p, rem, tail, ready, mload, k, rnd;
    int mode; // 0 = Giffler-Thompson conflict set, 1 = all currently available ops
};

static uint64_t rng_state;
static inline uint64_t splitmix64() {
    uint64_t z = (rng_state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline double urand() { return (splitmix64() >> 11) * (1.0 / 9007199254740992.0); }

long long makespan_of(const vector<vector<int>>& seq) {
    vector<array<int,2>> to(N);
    vector<unsigned char> deg(N,0), cnt(N,0);
    vector<long long> w(N), dist(N,0);
    auto id = [&](int j, int k){ return j*M + k; };
    for (int j=0;j<J;j++) for (int k=0;k<M;k++) {
        int u=id(j,k); w[u]=ops[j][k].p; cnt[u]=0;
        if (k+1<M) { to[u][cnt[u]++] = id(j,k+1); deg[id(j,k+1)]++; }
    }
    for (int m=0;m<M;m++) for (int i=0;i+1<J;i++) {
        int j1=seq[m][i], j2=seq[m][i+1];
        int u=id(j1,posOf[j1][m]), v=id(j2,posOf[j2][m]);
        to[u][cnt[u]++] = v; deg[v]++;
    }
    deque<int> q;
    for (int u=0;u<N;u++) if (!deg[u]) { dist[u]=w[u]; q.push_back(u); }
    int seen=0; long long ans=0;
    while(!q.empty()) {
        int u=q.front(); q.pop_front(); seen++; ans=max(ans, dist[u]);
        for (int e=0;e<cnt[u];e++) {
            int v=to[u][e];
            long long nd=dist[u]+w[v]; if (nd>dist[v]) dist[v]=nd;
            if (--deg[v]==0) q.push_back(v);
        }
    }
    if (seen != N) return LLONG_MAX/4;
    return ans;
}

vector<vector<int>> construct(const Rule& r) {
    vector<int> k(J,0);
    vector<long long> jobReady(J,0), machReady(M,0), rem=totalWork;
    vector<vector<int>> seq(M); for (int m=0;m<M;m++) seq[m].reserve(J);
    for (int step=0; step<N; step++) {
        long long tstar = LLONG_MAX; int mst = -1;
        vector<long long> ests(J,0), ects(J,0);
        for (int j=0;j<J;j++) if (k[j] < M) {
            int m=ops[j][k[j]].m; long long p=ops[j][k[j]].p;
            long long est=max(jobReady[j], machReady[m]);
            long long ect=est+p; ests[j]=est; ects[j]=ect;
            if (ect < tstar) { tstar=ect; mst=m; }
        }
        int best=-1; double bv=1e300;
        for (int j=0;j<J;j++) if (k[j] < M) {
            int m=ops[j][k[j]].m;
            if (r.mode==0 && !(m==mst && ests[j] < tstar)) continue;
            long long p=ops[j][k[j]].p;
            long long tail = rem[j] - p;
            double val = r.est*ests[j] + r.ect*ects[j] + r.p*p + r.rem*rem[j] +
                         r.tail*tail + r.ready*jobReady[j] + r.mload*machLoad[m] +
                         r.k*k[j] + r.rnd*urand()*max(1LL,totalWork[j]);
            if (val < bv - 1e-9 || (fabs(val-bv)<=1e-9 && j < best)) { bv=val; best=j; }
        }
        if (best < 0) { // should not happen; fall back to earliest completion op
            for (int j=0;j<J;j++) if (k[j]<M && (best<0 || ects[j]<ects[best])) best=j;
        }
        int j=best, kk=k[j], m=ops[j][kk].m; long long p=ops[j][kk].p;
        long long st=max(jobReady[j], machReady[m]), ft=st+p;
        seq[m].push_back(j);
        jobReady[j]=ft; machReady[m]=ft; rem[j]-=p; k[j]++;
    }
    return seq;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if (!(cin>>J>>M)) return 0; N=J*M;
    ops.assign(J, vector<Op>(M)); posOf.assign(J, vector<int>(M));
    totalWork.assign(J,0); machLoad.assign(M,0);
    rng_state = 123456789ULL ^ (uint64_t)J*1000003ULL ^ (uint64_t)M*9176ULL;
    for (int j=0;j<J;j++) for (int k=0;k<M;k++) {
        int m; long long p; cin>>m>>p; ops[j][k]={m,p}; posOf[j][m]=k;
        totalWork[j]+=p; machLoad[m]+=p; rng_state ^= (uint64_t)(m+1)*0x9e3779b97f4a7c15ULL + (uint64_t)p + (rng_state<<6) + (rng_state>>2);
    }

    vector<Rule> rules;
    auto add=[&](double a,double b,double c,double d,double e,double f,double g,double h,int mode){rules.push_back({a,b,c,d,e,f,g,h,0.0,mode});};
    // Active dispatch rules: choose inside the Giffler-Thompson conflict set.
    for (int mode=0; mode<=1; mode++) {
        add(0,1,0,0,0,0,0,0,mode);          // earliest completion
        add(1,0,0,0,0,0,0,0,mode);          // earliest start
        add(0,0,1,0,0,0,0,0,mode);          // SPT
        add(0,0,-1,0,0,0,0,0,mode);         // LPT
        add(0,0,0,-1,0,0,0,0,mode);         // most work remaining
        add(0,0,0,1,0,0,0,0,mode);          // least work remaining
        add(0,0,-1,-0.4,0,0,0,0,mode);
        add(0,0,1,-0.6,0,0,0,0,mode);
        add(0.2,0.3,-1,-0.5,0,0,0,0,mode);
        add(0.4,0,1,-1,0,0,0,0,mode);
        add(0,0,-1,0,0,0,-0.15,0,mode);
        add(0,0,1,0,0,0,-0.15,0,mode);
        add(0,0,-1,-1,0,0,-0.25,0,mode);
        add(0.5,0,0,-1,0,0,0,0,mode);
        add(0,0,0,-1,-0.5,0,0,0,mode);
        add(0,0,0,-1,0,0,0,-1000,mode);    // many ops remaining
    }

    vector<vector<int>> bestSeq; long long best = LLONG_MAX/4;
    auto consider = [&](const Rule& r){
        auto s = construct(r);
        long long v = makespan_of(s);
        if (v < best) { best=v; bestSeq.swap(s); }
    };
    for (const auto& r: rules) consider(r);

    auto start = chrono::steady_clock::now();
    int iter=0;
    while (chrono::duration<double>(chrono::steady_clock::now()-start).count() < 0.78) {
        Rule r;
        r.mode = (splitmix64() & 3) ? 0 : 1;
        auto coef=[&](double scale){ return (urand()*2.0-1.0)*scale; };
        r.est=coef(1.2); r.ect=coef(1.2); r.p=coef(3.0); r.rem=coef(2.0);
        r.tail=coef(1.0); r.ready=coef(0.8); r.mload=coef(0.5); r.k=coef(20000.0);
        r.rnd = (iter%5==0) ? 0.08 : 0.0;
        consider(r); iter++;
    }

    if (bestSeq.empty()) { bestSeq.assign(M, vector<int>()); for(int m=0;m<M;m++) for(int j=0;j<J;j++) bestSeq[m].push_back(j); }
    for (int m=0;m<M;m++) {
        for (int i=0;i<J;i++) { if (i) cout << ' '; cout << bestSeq[m][i]; }
        cout << '\n';
    }
    return 0;
}
