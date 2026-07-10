#include <cstdio>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <chrono>
#include <utility>
using namespace std;
typedef long long LL;

int main(){
    int J,M;
    if(scanf("%d%d",&J,&M)!=2) return 0;
    if(J<=0||M<=0) return 0;
    int N=J*M;

    // Operation o = (job j, rank k) encoded as o = j*M + k.
    vector<int> macOf(N), rankK(N), jobOf(N), opIdx(N);
    vector<LL> proc(N);
    vector<vector<int>> rankInJob(J, vector<int>(M)); // rankInJob[j][m] = rank of job j's op on machine m
    vector<LL> jobTotal(J,0);                        // total processing time per job
    for(int j=0;j<J;j++) for(int k=0;k<M;k++){
        int m; long long p;
        if(scanf("%d%lld",&m,&p)!=2) return 1;
        int o=j*M+k;
        macOf[o]=m; rankK[o]=k; jobOf[o]=j; proc[o]=p;
        rankInJob[j][m]=k;
        opIdx[j*M+m]=o;
        jobTotal[j]+=p;
    }

    // Machine orders (permutations of jobs) and their inverses.
    vector<vector<int>> ord(M), pos(M);
    for(int m=0;m<M;m++){
        ord[m].resize(J); pos[m].resize(J);
        for(int j=0;j<J;j++){ ord[m][j]=j; pos[m][j]=j; }
    }

    const LL SENT=4000000000000000000LL;

    // Longest-path makespan via Kahn's algorithm over the disjunctive graph.
    // Returns false if the current machine ordering induces a cycle (i.e. the
    // implied schedule is infeasible). This is the ground-truth evaluation:
    // every accept/reject decision and the final output are validated by it, so
    // feasibility and validity of every permutation are always guaranteed.
    vector<int> inDeg(N), q(N);
    vector<LL> fin(N);
    auto evalFeas=[&](LL &mk)->bool{
        for(int o=0;o<N;o++){
            int m=macOf[o];
            inDeg[o]=(rankK[o]>0?1:0)+(pos[m][jobOf[o]]>0?1:0);
        }
        int qs=0;
        for(int o=0;o<N;o++) if(inDeg[o]==0) q[qs++]=o;
        int processed=0; mk=0;
        for(int f=0;f<qs;f++){
            int o=q[f]; processed++;
            int k=rankK[o], m=macOf[o], j=jobOf[o], p=pos[m][j];
            LL mx=0;
            if(k>0){ int pr=o-1; if(fin[pr]>mx) mx=fin[pr]; }                 // job-chain predecessor
            if(p>0){ int nj=ord[m][p-1]; int pr=opIdx[nj*M+m]; if(fin[pr]>mx) mx=fin[pr]; } // machine predecessor
            LL nv=mx+proc[o]; fin[o]=nv; if(nv>mk) mk=nv;
            if(k+1<M){ int s=o+1; if(--inDeg[s]==0) q[qs++]=s; }               // job successor
            if(p+1<J){ int nj=ord[m][p+1]; int s=opIdx[nj*M+m]; if(--inDeg[s]==0) q[qs++]=s; } // machine successor
        }
        return processed==N;
    };

    LL bestMk=SENT;
    vector<vector<int>> bestOrd;
    auto consider=[&](){ LL mk; if(evalFeas(mk) && mk<bestMk){ bestMk=mk; bestOrd=ord; } };

    // Deterministic xorshift PRNG seeded by instance dimensions (consistent per
    // input, distinct sequences across sizes).
    uint64_t rng=88172645463325252ULL;
    rng^=(uint64_t)J*1099511628211ULL;
    rng^=(uint64_t)M*4294967311ULL;
    rng^=(uint64_t)N*1099511628211ULL;
    auto rnd=[&]()->unsigned{
        rng^=rng<<13; rng^=rng>>7; rng^=rng<<17;
        return (unsigned)(rng>>32);
    };

    auto t0=chrono::steady_clock::now();
    const double LIMIT=0.82;
    auto timeUp=[&]{ return chrono::duration<double>(chrono::steady_clock::now()-t0).count()>LIMIT; };

    // Self-inverse adjacent transposition: swap jobs at positions i and i+1 on machine m.
    auto swapAdj=[&](int m,int i){
        int a=ord[m][i], b=ord[m][i+1];
        swap(ord[m][i],ord[m][i+1]);
        pos[m][a]=i+1; pos[m][b]=i;
    };

    // Serial schedule generation scheme. Each step dispatches the next operation
    // of exactly one job (advancing nextRank), so job precedence is satisfied by
    // construction. Machine order = dispatch order, which yields an active,
    // hence acyclic, schedule. `mode` selects the priority rule:
    //   0=SPT, 1=LPT, 2=LRT(remaining job total), 3=EFT(earliest finish time).
    // When `randomized`, occasionally the second-best eligible job is chosen
    // (restricted candidate list of size 2) to diversify the construction.
    auto build=[&](int mode, bool randomized){
        vector<int> nextRank(J,0);
        vector<LL> jobFin(J,0), machFree(M,0), doneProc(J,0);
        for(int m=0;m<M;m++){ ord[m].clear(); pos[m].assign(J,-1); }
        for(int step=0;step<N;step++){
            int bj=-1, bj2=-1; LL bkey=0, bkey2=0;
            for(int j=0;j<J;j++){
                if(nextRank[j]>=M) continue;
                int k=nextRank[j], m=macOf[j*M+k];
                LL start = jobFin[j] > machFree[m] ? jobFin[j] : machFree[m];
                LL key;
                if(mode==0) key = proc[j*M+k];                       // SPT
                else if(mode==1) key = -proc[j*M+k];                 // LPT
                else if(mode==2) key = -(jobTotal[j]-doneProc[j]);   // LRT
                else key = start + proc[j*M+k];                      // EFT
                if(bj<0 || key<bkey){ bkey2=bkey; bj2=bj; bkey=key; bj=j; }
                else if(bj2<0 || key<bkey2){ bkey2=key; bj2=j; }
            }
            int j = (randomized && bj2>=0 && (rnd()&1u)) ? bj2 : bj;
            int k=nextRank[j], o=j*M+k, m=macOf[o];
            LL start = jobFin[j] > machFree[m] ? jobFin[j] : machFree[m];
            LL finV = start + proc[o];
            jobFin[j]=finV; machFree[m]=finV; doneProc[j]+=proc[o];
            pos[m][j]=(int)ord[m].size(); ord[m].push_back(j);
            nextRank[j]++;
        }
    };

    // Per-machine processing-time ordering. This directly targets bottleneck
    // machines (SPT shortens their cycle; LPT can balance heavy jobs). Independent
    // per-machine sorts may in general induce cycles, so each is validated by
    // evalFeas via consider() before being kept.
    auto seedSort=[&](bool lpt){
        for(int m=0;m<M;m++){
            for(int j=0;j<J;j++) ord[m][j]=j;
            sort(ord[m].begin(), ord[m].end(), [&](int x,int y){
                LL px=proc[(LL)x*M+rankInJob[x][m]];
                LL py=proc[(LL)y*M+rankInJob[y][m]];
                return lpt ? (px>py) : (px<py);
            });
            for(int j=0;j<J;j++) pos[m][ord[m][j]]=j;
        }
    };

    // Seed 0: identity machine order. Always acyclic because every machine arc
    // then points from a lower job index to a higher one, so job index is
    // monotonically non-decreasing along any directed walk.
    consider();

    // Diverse feasible constructions; keep the best by makespan.
    if(J>1){
        build(0,false); consider();          // SPT
        build(1,false); consider();          // LPT
        build(2,false); consider();          // LRT
        build(3,false); consider();          // EFT
        seedSort(false); consider();         // per-machine SPT
        seedSort(true);  consider();         // per-machine LPT
        build(0,true);  consider();          // randomized SPT
        build(3,true);  consider();          // randomized EFT
        if(!timeUp()){ build(0,true); consider(); }
        if(!timeUp()){ build(3,true); consider(); }
        if(!timeUp()){ build(2,true); consider(); }
        if(!timeUp()){ build(0,true); consider(); }
        if(!timeUp()){ build(1,true); consider(); }
    }
    ord=bestOrd;
    LL cur=bestMk;

    // Adjacent transposition neighborhood.
    vector<pair<int,int>> moves;
    if(J>1) for(int m=0;m<M;m++) for(int i=0;i+1<J;i++) moves.push_back({m,i});
    int Nmoves=(int)moves.size();

    // Iterated local search. Each outer iteration restarts from the incumbent
    // best schedule (so perturbation never drifts away from the best found),
    // applies a few feasibility-preserving random adjacent swaps, then performs
    // a first-improvement descent augmented with sampled non-adjacent job swaps
    // to escape local optima of the adjacent-only neighborhood.
    if(Nmoves>0 && M>1){
        while(!timeUp()){
            ord=bestOrd; cur=bestMk;

            // Small perturbation from the incumbent: a few random adjacent
            // transpositions, reverting any that would create a cycle.
            int P = 4 + N/250; if(P>40) P=40;
            for(int z=0;z<P;z++){
                if(timeUp()) break;
                unsigned idx=(unsigned)(rnd()%(unsigned)Nmoves);
                if(idx>(unsigned)(Nmoves-1)) idx=(unsigned)(Nmoves-1);
                int m=moves[idx].first, i=moves[idx].second;
                swapAdj(m,i);
                LL mk; if(!evalFeas(mk)){ swapAdj(m,i); continue; }
                cur=mk;
            }

            // Variable-neighborhood descent until a local optimum under both the
            // adjacent and the (sampled) non-adjacent swap neighborhoods.
            bool redo=true;
            while(redo && !timeUp()){
                redo=false;
                // First-improvement descent over adjacent transpositions.
                bool improved=true;
                while(improved && !timeUp()){
                    improved=false;
                    for(int x=Nmoves-1;x>0;--x){
                        unsigned y=rnd()%(unsigned)(x+1);
                        if(y>(unsigned)x) y=(unsigned)x;
                        swap(moves[x],moves[y]);
                    }
                    for(int xi=0;xi<Nmoves;xi++){
                        if(timeUp()) break;
                        int m=moves[xi].first, i=moves[xi].second;
                        swapAdj(m,i);
                        LL mk; bool ok=evalFeas(mk);
                        if(ok && mk<cur){ cur=mk; improved=true; break; }
                        else swapAdj(m,i); // revert (self-inverse)
                    }
                }
                if(timeUp()) break;

                // Escape moves: sampled non-adjacent job swaps on a machine
                // (larger steps that the adjacent neighborhood cannot make).
                // Each is validated by evalFeas; only improving, feasible moves
                // are accepted, preserving permutation validity throughout.
                if(J>=3){
                    int tries = 4 + N/200; if(tries>40) tries=40;
                    for(int t=0;t<tries;t++){
                        if(timeUp()) break;
                        int m=(int)(rnd()%(unsigned)M);
                        int i=(int)(rnd()%(unsigned)J);
                        int k=(int)(rnd()%(unsigned)J);
                        if(i==k) continue;
                        if(i>k) swap(i,k);
                        int ai=ord[m][i], ak=ord[m][k];
                        swap(ord[m][i],ord[m][k]);
                        pos[m][ai]=k; pos[m][ak]=i;
                        LL mk; bool ok=evalFeas(mk);
                        if(ok && mk<cur){ cur=mk; redo=true; }
                        else { swap(ord[m][i],ord[m][k]); pos[m][ai]=i; pos[m][ak]=k; }
                    }
                }
            }
            if(cur<bestMk){ bestMk=cur; bestOrd=ord; }
        }
    }
    if(cur<bestMk){ bestMk=cur; bestOrd=ord; }
    ord=bestOrd;

    // Emit exactly M lines, each a permutation of {0..J-1}.
    for(int m=0;m<M;m++){
        for(int j=0;j<J;j++){
            printf("%d", ord[m][j]);
            if(j+1<J) putchar(' '); else putchar('\n');
        }
    }
    return 0;
}
