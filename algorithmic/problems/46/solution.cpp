#include <cstdio>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <chrono>
using namespace std;

// Job Shop Scheduling (JSSP) solver.
//
// Feasible machine processing orders are produced by a serial schedule
// generation scheme (SGS) that dispatches ready operations by priority. Five
// deterministic dispatching rules are sampled (longest/shortest remaining
// route time, shortest/longest current processing time, and a first-machine
// SPT rule that helps near-flow-shop instances); their priorities are given a
// small random jitter to break ties, and purely random priorities are also
// sampled for diversification. The jitter is scaled so that restarts remain
// meaningfully different even when processing times are very large. Every
// schedule built by the serial SGS is acyclic, so the output is always valid.
//
// Orders are held in flat arrays with positions maintained incrementally, and
// the longest-path evaluation computes each node's indegree from a closed-form
// expression, which keeps the makespan kernel fast enough to run many
// independent local optimizations. A multi-start wrapper pairs every
// construction with a first-improvement adjacent-transposition descent to a
// local optimum; swaps that would create a cyclic disjunctive graph are
// detected and reverted, preserving validity. After multi-start, an iterated
// local search perturbs the global best and re-optimizes, injecting fresh
// constructions on stagnation. The global best is monotone, so the result
// never degrades. A wall-clock budget keeps every case within the time limit.
// Processing times use 64-bit integers to tolerate large values.

static int J, M, N;
static vector<int> mmat;          // mmat[j*M+k] : machine of operation (j,k)
static vector<int> inv;           // inv[j*M+m]  : step k where job j visits machine m
static vector<int> deg0;          // deg0[j*M+k]: static job-chain indegree (k>0 ? 1 : 0)
static vector<int> deg;           // working indegree, size N
static vector<int> qbuf;          // topological queue, size <= N
static vector<int> gOrders;       // global best orders, size M*J
static vector<int> gPos;          // global best positions, size M*J
static vector<long long> pmat;   // pmat[j*M+k] : processing time
static vector<long long> remain; // remain[j*M+k]: remaining route time from step k
static vector<long long> ef;     // earliest finish, size N
static long long gBest;

using Deadline = chrono::steady_clock::time_point;
static mt19937 rng;

static inline bool past(Deadline d){ return chrono::steady_clock::now() >= d; }

// Swap the jobs at positions i and i+1 on machine m, keeping pos consistent.
static inline void swap_adj(vector<int>& o, vector<int>& pv, int m, int i){
    int base = m * J;
    int x = o[base + i], y = o[base + i + 1];
    o[base + i] = y; o[base + i + 1] = x;
    pv[m * J + x] += 1;   // x moved one step later
    pv[m * J + y] -= 1;   // y moved one step earlier
}

// Rebuild positions from orders (called after every construction).
static void sync_pos(const vector<int>& o, vector<int>& pv){
    for(int m=0;m<M;++m){
        int base=m*J;
        for(int i=0;i<J;++i) pv[m*J + o[base+i]] = i;
    }
}

// Build priority values for a dispatching rule. Larger values dispatch first.
// mode 0: longest remaining processing time (LRPT)
// mode 1: shortest processing time on current op (SPT)
// mode 2: longest processing time on current op (LPT)
// mode 3: shortest remaining processing time (MRPT)
// mode 4: shortest processing time on the first machine (flow-shop aware)
// mode 5: purely random priorities
static vector<long long> build_prio(int mode){
    vector<long long> prio((size_t)N);
    for(int j=0;j<J;++j){
        int jb=j*M;
        long long f0 = pmat[jb];
        for(int k=0;k<M;++k){
            int id=jb+k;
            long long v;
            if(mode==0)      v = remain[id];
            else if(mode==1) v = -pmat[id];
            else if(mode==2) v = pmat[id];
            else if(mode==3) v = -remain[id];
            else if(mode==4) v = -f0;
            else             v = (long long)rng();
            // Scale by 64 and add small jitter so ties break diversely even
            // when processing times are very large.
            prio[id] = v * 64 + (long long)(rng() & 63);
        }
    }
    return prio;
}

// Serial SGS: dispatch ready operations by priority, append jobs to machine
// orders. Returns an upper-bound makespan from the generated schedule.
static long long serial_sgs(const vector<long long>& prio, vector<int>& o){
    vector<long long> mfree(M,0), jfin(J,0);
    vector<int> cnt(M,0);
    using S = pair<long long,int>;
    priority_queue<S> pq;                 // max-heap by priority
    for(int j=0;j<J;++j) pq.push({prio[j*M], j*M});
    long long ms=0;
    while(!pq.empty()){
        auto t=pq.top(); pq.pop();
        int id=t.second, j=id/M, k=id%M, m=mmat[id];
        long long s = mfree[m] < jfin[j] ? jfin[j] : mfree[m];
        long long f = s + pmat[id];
        if(f>ms) ms=f;
        mfree[m]=f; jfin[j]=f;
        o[m*J + cnt[m]++] = j;
        if(k+1<M) pq.push({prio[j*M+k+1], j*M+k+1});
    }
    return ms;
}

// Earliest-feasible makespan = longest path in the disjunctive graph (job-chain
// arcs + oriented machine arcs), computed by Kahn topological relaxation.
// Node indegree equals (k>0) + (position on its machine >0). Returns -1 cyclic.
static long long evaluate(const vector<int>& o, const vector<int>& pv){
    for(int id=0;id<N;++id){ deg[id]=deg0[id]; ef[id]=0; }
    for(int m=0;m<M;++m)
        for(int i=1;i<J;++i){ int j=o[m*J+i]; ++deg[j*M+inv[j*M+m]]; }
    qbuf.clear();
    for(int j=0;j<J;++j) for(int k=0;k<M;++k){
        int id=j*M+k;
        if(!deg[id]){ ef[id]=pmat[id]; qbuf.push_back(id); }
    }
    long long best=0; int done=0, qh=0;
    while(qh<(int)qbuf.size()){
        int id=qbuf[qh++]; ++done;
        long long f=ef[id];
        if(f>best) best=f;
        int j=id/M, k=id%M, m=mmat[id];
        if(k+1<M){
            int t=j*M+k+1;
            if(f>ef[t]) ef[t]=f;
            if(--deg[t]==0){ ef[t]+=pmat[t]; qbuf.push_back(t); }
        }
        int pp=pv[m*J+j];
        if(pp+1<J){
            int nb=o[m*J+pp+1];
            int t=nb*M+inv[nb*M+m];
            if(f>ef[t]) ef[t]=f;
            if(--deg[t]==0){ ef[t]+=pmat[t]; qbuf.push_back(t); }
        }
    }
    return done<N ? -1 : best;
}

// First-improvement adjacent-transposition descent to a local optimum.
// Updates the global incumbent whenever a better feasible schedule is found.
static void local_search(vector<int>& o, vector<int>& pv, Deadline d){
    long long cur=evaluate(o,pv);
    if(cur<0) return;
    if(cur<gBest){ gBest=cur; gOrders=o; gPos=pv; }
    bool improved=true;
    while(improved && !past(d)){
        improved=false;
        int chk=0;
        for(int m=0;m<M;++m){
            for(int i=0;i+1<J;++i){
                swap_adj(o,pv,m,i);
                long long ns=evaluate(o,pv);
                if(ns>=0 && ns<cur){
                    cur=ns; improved=true;
                    if(ns<gBest){ gBest=ns; gOrders=o; gPos=pv; }
                }else{
                    swap_adj(o,pv,m,i); // revert
                }
                if((++chk & 31)==0 && past(d)) return;
            }
        }
    }
}

// Construct one schedule with the given rule and optimize it in place.
static void construct_and_optimize(int mode, Deadline d){
    vector<int> o((size_t)M*J, 0), pv((size_t)M*J, -1);
    serial_sgs(build_prio(mode), o);
    sync_pos(o,pv);
    local_search(o,pv,d);
}

int main(){
    if(scanf("%d%d",&J,&M)!=2) return 0;
    if(J<=0||M<=0) return 0;
    N=J*M;
    size_t n=(size_t)N, mj=(size_t)M*J;
    mmat.assign(n,0); inv.assign(n,0); deg0.assign(n,0); deg.assign(n,0);
    pmat.assign(n,0); remain.assign(n,0); ef.assign(n,0);
    qbuf.reserve(n);
    gOrders.assign(mj,0); gPos.assign(mj,-1);

    for(int j=0;j<J;++j){
        int jb=j*M;
        for(int k=0;k<M;++k){
            int m; long long p;
            if(scanf("%d%lld",&m,&p)!=2) return 1;
            int id=jb+k;
            mmat[id]=m; pmat[id]=p; inv[jb+m]=k;
        }
    }
    for(int j=0;j<J;++j){
        int jb=j*M; long long acc=0;
        for(int k=M-1;k>=0;--k){ acc+=pmat[jb+k]; remain[jb+k]=acc; }
    }
    for(int j=0;j<J;++j){
        int jb=j*M;
        for(int k=1;k<M;++k) deg0[jb+k]=1;   // job-chain predecessor
    }

    Deadline start=chrono::steady_clock::now();
    Deadline dl_ms = start + chrono::milliseconds(600);
    Deadline dl_full = start + chrono::milliseconds(850);
    rng.seed(20240710);
    gBest = (long long)4e18;

    // Deterministic dispatching rules, each optimized to a local optimum.
    for(int mode=0; mode<5; ++mode) construct_and_optimize(mode, dl_ms);

    // Randomized multi-start: diverse constructions each followed by local
    // search. Bounded iteration count guards degenerate tiny instances.
    for(int r=0; r<20000 && !past(dl_ms); ++r) construct_and_optimize(5, dl_ms);

    // Iterated local search on the global best: perturb with random adjacent
    // transpositions (reverting any that create a cycle) and re-optimize.
    if(J>=2){
        int lastImp=0;
        for(int it=0; it<100000000 && !past(dl_full); ++it){
            vector<int> co=gOrders, cp=gPos;
            long long prev=gBest;
            int strength=J/2;
            if(strength<4) strength=4;
            if(strength>J-1) strength=J-1;
            for(int t=0;t<strength;++t){
                int m=(int)(rng()%M);
                int i=(int)(rng()%(J-1));
                swap_adj(co,cp,m,i);
                if(evaluate(co,cp)<0) swap_adj(co,cp,m,i); // revert if cyclic
            }
            local_search(co,cp,dl_full);
            if(gBest<prev) lastImp=it;
            else if(it-lastImp>60){
                // Inject a fresh randomized construction to escape stagnation.
                construct_and_optimize(5, dl_full);
                lastImp=it;
            }
        }
    }

    // Fallback: a feasible LRPT schedule if no search produced one.
    if(gBest>= (long long)4e18){
        vector<int> o(mj,0), pv(mj,-1);
        serial_sgs(build_prio(0), o);
        sync_pos(o,pv);
        gOrders=o; gPos=pv; gBest=evaluate(o,pv);
    }

    for(int m=0;m<M;++m){
        const int* o=&gOrders[m*J];
        for(int i=0;i<J;++i){
            if(i+1<J) printf("%d ", o[i]);
            else printf("%d\n", o[i]);
        }
    }
    return 0;
}
