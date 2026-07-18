#include <bits/stdc++.h>
using namespace std;

struct Op { int m; long long p; };
struct Cand { long long ms; vector<vector<int>> order; };

int J, M;
vector<vector<Op>> ops;
vector<vector<long long>> suffixJob;
vector<long long> machLoad;

long long evalOrder(const vector<vector<int>>& order){
    int N = J*M;
    vector<vector<int>> id(J, vector<int>(M));
    vector<int> jid(N), pos(N), mach(N);
    vector<long long> dur(N);
    for(int j=0;j<J;j++) for(int k=0;k<M;k++){
        int x = j*M+k; id[j][k]=x; jid[x]=j; pos[x]=k; mach[x]=ops[j][k].m; dur[x]=ops[j][k].p;
    }
    vector<vector<int>> adj(N);
    vector<int> indeg(N,0);
    adj.assign(N, {});
    for(int j=0;j<J;j++) for(int k=0;k+1<M;k++) { adj[id[j][k]].push_back(id[j][k+1]); indeg[id[j][k+1]]++; }
    for(int m=0;m<M;m++) for(int a=0;a+1<J;a++){
        int j1=order[m][a], j2=order[m][a+1];
        int k1=-1,k2=-1;
        // use precomputed route search small enough
        for(int k=0;k<M;k++){ if(ops[j1][k].m==m) k1=k; if(ops[j2][k].m==m) k2=k; }
        int u=id[j1][k1], v=id[j2][k2]; adj[u].push_back(v); indeg[v]++;
    }
    deque<int> q; vector<long long> dist(N,0);
    for(int i=0;i<N;i++) if(!indeg[i]) q.push_back(i);
    int seen=0; long long ans=0;
    while(!q.empty()){
        int u=q.front(); q.pop_front(); seen++; ans=max(ans, dist[u]+dur[u]);
        long long nd=dist[u]+dur[u];
        for(int v: adj[u]){ if(dist[v]<nd) dist[v]=nd; if(--indeg[v]==0) q.push_back(v); }
    }
    if(seen<N) return LLONG_MAX/4;
    return ans;
}

struct RNG { unsigned long long x=88172645463325252ull; unsigned next(){ x^=x<<7; x^=x>>9; return (unsigned)x; } double uni(){ return (next()&0xffffff)/double(0x1000000); } };

vector<vector<int>> buildSchedule(const array<double,7>& w, RNG &rng, bool randomize){
    vector<int> pos(J,0);
    vector<long long> jr(J,0), mr(M,0), machDone(M,0);
    vector<vector<int>> order(M); for(int m=0;m<M;m++) order[m].reserve(J);
    int total=J*M;
    for(int step=0; step<total; ++step){
        int best=-1; double bestKey=1e100;
        for(int j=0;j<J;j++) if(pos[j]<M){
            Op o=ops[j][pos[j]];
            long long st=max(jr[j], mr[o.m]);
            long long ct=st+o.p;
            long long rem=suffixJob[j][pos[j]];
            long long slackMachine = machLoad[o.m] - machDone[o.m];
            double key = 0.0;
            key += w[0]*st;
            key += w[1]*ct;
            key += w[2]*o.p;
            key += w[3]*rem;
            key += w[4]*slackMachine;
            key += w[5]*jr[j];
            key += w[6]*mr[o.m];
            if(randomize) key *= (0.92 + 0.16*rng.uni());
            key += 1e-7*j;
            if(key < bestKey){ bestKey=key; best=j; }
        }
        int j=best; Op o=ops[j][pos[j]];
        long long st=max(jr[j], mr[o.m]);
        long long ct=st+o.p;
        jr[j]=ct; mr[o.m]=ct; machDone[o.m]+=o.p; order[o.m].push_back(j); pos[j]++;
    }
    return order;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>J>>M)) return 0;
    ops.assign(J, vector<Op>(M)); machLoad.assign(M,0);
    for(int j=0;j<J;j++) for(int k=0;k<M;k++){ cin>>ops[j][k].m>>ops[j][k].p; machLoad[ops[j][k].m]+=ops[j][k].p; }
    suffixJob.assign(J, vector<long long>(M+1,0));
    for(int j=0;j<J;j++) for(int k=M-1;k>=0;k--) suffixJob[j][k]=suffixJob[j][k+1]+ops[j][k].p;

    vector<array<double,7>> rules = {
        {0,1,0,0,0,0,0},          // earliest completion
        {1,0,0,0,0,0,0},          // earliest start
        {0,1,1,0,0,0,0},          // SPT among early completions
        {0,1,-1,0,0,0,0},         // LPT among early completions
        {0,1,0,-0.35,0,0,0},      // most work remaining
        {0,1,0,-0.8,-0.15,0,0},   // bottleneck/work remaining
        {0.5,0.5,0,-0.6,-0.25,0,0},
        {0,0,0,-1,0,0,0},         // pure MWKR
        {0,0,-1,-0.3,-0.2,0,0},
        {0,1,0,0,-0.6,0,0},
        {0.2,0.8,0.2,-0.5,-0.1,0,0},
        {0.2,0.8,-0.2,-0.5,-0.3,0,0}
    };

    RNG rng;
    Cand best; best.ms=LLONG_MAX/4;
    auto start = chrono::steady_clock::now();
    auto consider = [&](const vector<vector<int>>& ord){
        long long ms=evalOrder(ord);
        if(ms < best.ms){ best.ms=ms; best.order=ord; }
    };
    for(auto &r: rules){ auto ord=buildSchedule(r, rng, false); consider(ord); }

    int iter=0;
    while(true){
        if((iter & 15)==0){
            double elapsed=chrono::duration<double>(chrono::steady_clock::now()-start).count();
            if(elapsed>0.855) break;
        }
        array<double,7> w;
        // sample around useful dispatch families: early finish plus negative remaining/bottleneck terms
        w[0] = rng.uni()*0.8;
        w[1] = 0.3 + rng.uni()*1.4;
        w[2] = -0.8 + rng.uni()*1.6;
        w[3] = -1.2 * rng.uni();
        w[4] = -0.8 * rng.uni();
        w[5] = -0.2 + rng.uni()*0.4;
        w[6] = -0.2 + rng.uni()*0.4;
        auto ord=buildSchedule(w, rng, true);
        consider(ord);
        iter++;
    }
    if(best.order.empty()){
        array<double,7> r={0,1,0,0,0,0,0}; best.order=buildSchedule(r,rng,false);
        best.ms=evalOrder(best.order);
    }

    // Targeted repair/improvement over the incumbent: spend the final slice on
    // first-improvement adjacent swaps in the selected machine permutations.
    // This keeps the representation general and lets evalOrder reject cyclic moves.
    vector<int> mids(M); iota(mids.begin(), mids.end(), 0);
    int noImprove=0, lsIter=0;
    while(J>1){
        if((lsIter & 7)==0){
            double elapsed=chrono::duration<double>(chrono::steady_clock::now()-start).count();
            if(elapsed>0.94) break;
        }
        int m = mids[rng.next()%M];
        int a = rng.next()%(J-1);
        swap(best.order[m][a], best.order[m][a+1]);
        long long ms=evalOrder(best.order);
        if(ms < best.ms){
            best.ms=ms;
            noImprove=0;
        }else{
            swap(best.order[m][a], best.order[m][a+1]);
            noImprove++;
        }
        // Occasionally scan the currently heaviest machine by remaining load; random
        // attempts alone can miss simple adjacent inversions on bottlenecks.
        if(noImprove>64){
            int bm=max_element(machLoad.begin(), machLoad.end())-machLoad.begin();
            bool got=false;
            for(int a2=0; a2+1<J; ++a2){
                swap(best.order[bm][a2], best.order[bm][a2+1]);
                long long ms2=evalOrder(best.order);
                if(ms2 < best.ms){ best.ms=ms2; got=true; break; }
                swap(best.order[bm][a2], best.order[bm][a2+1]);
                if(chrono::duration<double>(chrono::steady_clock::now()-start).count()>0.94) break;
            }
            noImprove = got ? 0 : 0;
        }
        lsIter++;
    }

    for(int m=0;m<M;m++){
        // Defensive repair should not be needed, but ensure each line is a permutation.
        if((int)best.order[m].size()!=J){
            vector<int> used(J,0), fixed;
            for(int x: best.order[m]) if(0<=x && x<J && !used[x]) used[x]=1, fixed.push_back(x);
            for(int j=0;j<J;j++) if(!used[j]) fixed.push_back(j);
            best.order[m]=fixed;
        }
        for(int i=0;i<J;i++){ if(i) cout << ' '; cout << best.order[m][i]; }
        cout << '\n';
    }
    return 0;
}
