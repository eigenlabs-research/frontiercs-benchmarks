#include <cstdio>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <climits>
#include <unistd.h>
using namespace std;

static int J, M, N;
static vector<vector<int>> m_of;       // [j][k] machine of op k of job j
static vector<vector<long long>> p_of; // [j][k] processing time
static vector<vector<int>> pos;        // [j][m] = k (position of machine m in job j route)
static vector<int> posF;               // posF[j*M+m] = k (flattened pos)
static vector<long long> pnode;        // pnode[u] = p_of[j][k] for u=j*M+k

// Machine-successor op of each node (rebuilt each eval) and work buffers.
static vector<int> msucc, mpred, indeg, qbuf;
static vector<int> jsuc;               // jsuc[u] = u+1 if job successor exists else -1 (fixed)
static vector<char> jpre;              // jpre[u] = 1 if job predecessor exists (fixed)
static vector<long long> dist_, tail_;
static vector<char> crit;

// Compute makespan of seq. Returns -1 if cycle. Optionally fill critical-node flags
// and tails. A node u is critical if dist_[u] + tail_[u] - pnode[u] == C.
// Successors are addressed directly: job successor of u is u+1 (if not last op of
// its job); machine successor comes from msucc[], rebuilt from seq in O(N).
static long long evalSeq(const vector<vector<int>>& seq, bool fillCrit = false){
    const int n = N;
    int* __restrict ind = indeg.data();
    int* __restrict msu = msucc.data();
    const int* __restrict jsu = jsuc.data();
    const char* __restrict jpr = jpre.data();
    const long long* __restrict pn = pnode.data();
    long long* __restrict ds = dist_.data();
    for(int u=0;u<n;++u) ind[u] = jpr[u];
    for(int m=0;m<M;++m){
        const int* s = seq[m].data();
        int prev = s[0]*M + posF[s[0]*M + m];
        mpred[prev] = -1;
        for(int i=1;i<J;++i){
            int v = s[i]*M + posF[s[i]*M + m];
            msu[prev] = v;
            mpred[v] = prev;
            ind[v]++;
            prev = v;
        }
        msu[prev] = -1;
    }
    fill(dist_.begin(), dist_.begin()+n, 0);
    qbuf.clear();
    int qh=0;
    for(int u=0;u<n;++u) if(ind[u]==0){ ds[u]=pn[u]; qbuf.push_back(u); }
    while(qh < (int)qbuf.size()){
        int u = qbuf[qh++];
        long long du = ds[u];
        int v = jsu[u];
        if(v >= 0){
            long long nd = du + pn[v];
            if(nd > ds[v]) ds[v] = nd;
            if(--ind[v]==0) qbuf.push_back(v);
        }
        v = msu[u];
        if(v >= 0){
            long long nd = du + pn[v];
            if(nd > ds[v]) ds[v] = nd;
            if(--ind[v]==0) qbuf.push_back(v);
        }
    }
    if(qh != n) return -1; // cycle
    long long C = 0;
    for(int u=0;u<n;++u) if(ds[u] > C) C = ds[u];
    if(fillCrit){
        long long* __restrict tl = tail_.data();
        const int* __restrict qb = qbuf.data();
        for(int idx=n-1; idx>=0; --idx){
            int u = qb[idx];
            long long mx = 0;
            int v = jsu[u];
            if(v >= 0 && tl[v] > mx) mx = tl[v];
            v = msu[u];
            if(v >= 0 && tl[v] > mx) mx = tl[v];
            tl[u] = pn[u] + mx;
        }
        fill(crit.begin(), crit.begin()+n, 0);
        for(int u=0;u<n;++u) if(ds[u] + tl[u] - pn[u] == C) crit[u] = 1;
    }
    return C;
}

// Giffler-Thompson dispatch with a priority rule (active schedule generation).
// mode: 0=MWR (most work remaining), 1=LPT (current op), 2=SPT, 3=random.
static vector<vector<int>> seedGT(int mode, mt19937& rng){
    vector<int> jp(J, 0);
    vector<long long> jr(J, 0), mf(M, 0), wrem(J, 0);
    for(int j=0;j<J;++j) for(int k=0;k<M;++k) wrem[j] += p_of[j][k];
    vector<vector<int>> seq(M);
    int remaining = N;
    while(remaining > 0){
        long long bf = LLONG_MAX;
        for(int j=0;j<J;++j){
            if(jp[j] >= M) continue;
            int k = jp[j], m = m_of[j][k];
            long long s = max(jr[j], mf[m]);
            long long f = s + p_of[j][k];
            if(f < bf) bf = f;
        }
        int cm = -1;
        for(int j=0;j<J;++j){
            if(jp[j] >= M) continue;
            int k = jp[j], m = m_of[j][k];
            long long s = max(jr[j], mf[m]);
            if(s + p_of[j][k] == bf){ cm = m; break; }
        }
        int cj = -1; long long cp = 0;
        for(int j=0;j<J;++j){
            if(jp[j] >= M) continue;
            int k = jp[j], m = m_of[j][k];
            if(m != cm) continue;
            long long s = max(jr[j], mf[m]);
            if(s < bf){
                long long pr;
                if(mode==0) pr = wrem[j];
                else if(mode==1) pr = p_of[j][k];
                else if(mode==2) pr = -p_of[j][k];
                else pr = (long long)rng();
                if(cj==-1 || pr > cp){ cp = pr; cj = j; }
            }
        }
        if(cj == -1){
            for(int j=0;j<J;++j) if(jp[j]<M && m_of[j][jp[j]]==cm){ cj=j; break; }
        }
        int k = jp[cj], m = m_of[cj][k];
        long long s = max(jr[cj], mf[m]);
        long long f = s + p_of[cj][k];
        seq[m].push_back(cj);
        jr[cj] = f; mf[m] = f; wrem[cj] -= p_of[cj][k];
        jp[cj]++; remaining--;
    }
    return seq;
}

// ---- N7 / Balas-Vazacopoulos critical-block insertion moves ----
// A move takes the op at position i of machine m and inserts it at the front
// (position b) or back (position e) of its critical block [b..e].
struct Mv { int m, b, e, i; bool front; };

static vector<Mv> gmoves;
static vector<pair<long long,int>> gcand;
static vector<int> gord;           // scratch: new job order of the reordered segment
static vector<long long> gestC;    // scratch: estimated completions
static vector<int> tabuTB;         // TB[uop*J + j] : "job j ordered before op uop" tabu until iter

static inline int opOf(int job, int m){ return job*M + posF[job*M + m]; }

// Generate all N7 moves over critical blocks of the current solution (uses crit[]).
static void genMoves(const vector<vector<int>>& cur){
    gmoves.clear();
    for(int m=0;m<M;++m){
        const auto& s = cur[m];
        int i = 0;
        while(i < J){
            if(!crit[opOf(s[i], m)]){ i++; continue; }
            int b = i;
            while(i+1 < J && crit[opOf(s[i+1], m)]) i++;
            int e = i; i++;
            if(e == b) continue;
            for(int t=b+1; t<=e; ++t) gmoves.push_back({m,b,e,t,true});
            for(int t=b; t<e; ++t)
                if(!(t==b && e==b+1))          // avoid duplicate of the L==2 swap
                    gmoves.push_back({m,b,e,t,false});
        }
    }
}

// Taillard/BV-style estimate of the makespan after the move: recompute heads
// forward and tails backward along the reordered machine segment, holding
// everything outside fixed. Exact for the local path; a good estimate overall.
static long long estMove(const vector<vector<int>>& cur, const Mv& mv){
    const auto& s = cur[mv.m];
    int lo, hi;
    if(mv.front){
        lo = mv.b; hi = mv.i;
        gord[0] = s[mv.i];
        for(int t=lo; t<hi; ++t) gord[t-lo+1] = s[t];
    } else {
        lo = mv.i; hi = mv.e;
        for(int t=lo+1; t<=hi; ++t) gord[t-lo-1] = s[t];
        gord[hi-lo] = s[mv.i];
    }
    int L = hi - lo + 1;
    long long prevC = 0;
    if(lo > 0) prevC = dist_[opOf(s[lo-1], mv.m)];
    for(int t=0; t<L; ++t){
        int v = gord[t]; int u = opOf(v, mv.m); int k = pos[v][mv.m];
        long long jp = (k>0) ? dist_[u-1] : 0;
        long long st = prevC > jp ? prevC : jp;
        gestC[t] = st + pnode[u];
        prevC = gestC[t];
    }
    long long prevT = 0;
    if(hi+1 < J) prevT = tail_[opOf(s[hi+1], mv.m)];
    long long bestLen = 0;
    for(int t=L-1; t>=0; --t){
        int v = gord[t]; int u = opOf(v, mv.m); int k = pos[v][mv.m];
        long long js = (k<M-1) ? tail_[u+1] : 0;
        long long tl = pnode[u] + (prevT > js ? prevT : js);
        long long len = gestC[t] - pnode[u] + tl;
        if(len > bestLen) bestLen = len;
        prevT = tl;
    }
    return bestLen;
}

static inline void applyMove(vector<vector<int>>& cur, const Mv& mv){
    auto& s = cur[mv.m];
    if(mv.front) rotate(s.begin()+mv.b, s.begin()+mv.i, s.begin()+mv.i+1);
    else         rotate(s.begin()+mv.i, s.begin()+mv.i+1, s.begin()+mv.e+1);
}
static inline void undoMove(vector<vector<int>>& cur, const Mv& mv){
    auto& s = cur[mv.m];
    if(mv.front) rotate(s.begin()+mv.b, s.begin()+mv.b+1, s.begin()+mv.i+1);
    else         rotate(s.begin()+mv.i, s.begin()+mv.e,   s.begin()+mv.e+1);
}

// ---- Incremental longest-path update after a move ----
// Propagates exact head (dist_) and tail (tail_) changes from the reordered
// machine segment through the disjunctive graph via change-propagation
// worklists. Returns the new makespan, or -2 if the pop cap was hit (possible
// cycle or huge cascade) -- caller must then fall back to full evalSeq.
static vector<int> gwl;
static vector<char> ginq;
static long long incAfterMove(const vector<vector<int>>& cur, const Mv& mv){
    const int m = mv.m;
    const auto& s = cur[m];
    const int lo = mv.front ? mv.b : mv.i;
    const int hi = mv.front ? mv.i : mv.e;
    // Relink msucc/mpred across the modified span.
    {
        int from = lo>0 ? lo-1 : 0;
        int to   = hi<J-1 ? hi+1 : J-1;
        for(int i=from;i<=to;++i){
            int u = opOf(s[i], m);
            msucc[u] = (i<J-1) ? opOf(s[i+1], m) : -1;
            mpred[u] = (i>0)   ? opOf(s[i-1], m) : -1;
        }
    }
    const int cap = 16*N;
    long long* __restrict ds = dist_.data();
    long long* __restrict tl = tail_.data();
    const long long* __restrict pn = pnode.data();
    // Heads (forward). Seed the segment plus the op at hi+1 (its machine pred changed).
    gwl.clear();
    auto pushH = [&](int v){ if(v>=0 && !ginq[v]){ ginq[v]=1; gwl.push_back(v); } };
    int hiH = hi < J-1 ? hi+1 : hi;
    for(int i=lo;i<=hiH;++i) pushH(opOf(s[i], m));
    int wh = 0, pops = 0;
    while(wh < (int)gwl.size()){
        int v = gwl[wh++]; ginq[v] = 0;
        if(++pops > cap){
            for(int t=wh;t<(int)gwl.size();++t) ginq[gwl[t]] = 0;
            return -2;
        }
        long long b = 0;
        if(jpre[v] && ds[v-1] > b) b = ds[v-1];
        int mp = mpred[v];
        if(mp >= 0 && ds[mp] > b) b = ds[mp];
        long long nd = pn[v] + b;
        if(nd != ds[v]){ ds[v] = nd; pushH(jsuc[v]); pushH(msucc[v]); }
    }
    // Tails (backward). Seed the segment plus the op at lo-1 (its machine succ changed).
    gwl.clear();
    int loT = lo > 0 ? lo-1 : lo;
    for(int i=hi;i>=loT;--i){ int v = opOf(s[i], m); if(!ginq[v]){ ginq[v]=1; gwl.push_back(v); } }
    int headPops = pops;
    wh = 0; pops = 0;
    while(wh < (int)gwl.size()){
        int v = gwl[wh++]; ginq[v] = 0;
        if(++pops > cap){
            for(int t=wh;t<(int)gwl.size();++t) ginq[gwl[t]] = 0;
            return -2;
        }
        long long b = 0;
        int js = jsuc[v];
        if(js >= 0 && tl[js] > b) b = tl[js];
        int ms = msucc[v];
        if(ms >= 0 && tl[ms] > b) b = tl[ms];
        long long nt = pn[v] + b;
        if(nt != tl[v]){
            tl[v] = nt;
            if(jpre[v] && !ginq[v-1]){ ginq[v-1]=1; gwl.push_back(v-1); }
            int mp = mpred[v];
            if(mp >= 0 && !ginq[mp]){ ginq[mp]=1; gwl.push_back(mp); }
        }
    }
    long long C = 0;
    for(int u=0;u<N;++u) if(ds[u] > C) C = ds[u];
    for(int u=0;u<N;++u) crit[u] = (ds[u] + tl[u] - pn[u] == C);
#ifdef DIAG
    extern long long g_pops, g_calls;
    g_pops += pops + headPops; g_calls++;
#endif
    return C;
}
#ifdef DIAG
long long g_pops = 0, g_calls = 0;
#endif

// Tabu test: the move is tabu if any ordering it creates is currently forbidden.
static bool isTabu(const vector<vector<int>>& cur, const Mv& mv, int iter){
    const auto& s = cur[mv.m];
    int uj = s[mv.i];
    if(mv.front){
        // creates (u before x) for x in [b..i-1] -> forbidden if TB[xop][uj]
        for(int t=mv.b; t<mv.i; ++t){
            int xop = opOf(s[t], mv.m);
            if(tabuTB[(size_t)xop*J + uj] > iter) return true;
        }
    } else {
        int uop = opOf(uj, mv.m);
        // creates (x before u) for x in [i+1..e] -> forbidden if TB[uop][x]
        for(int t=mv.i+1; t<=mv.e; ++t)
            if(tabuTB[(size_t)uop*J + s[t]] > iter) return true;
    }
    return false;
}

// Collect the tabu entries for the orderings reversed by a move.
// NOTE: call BEFORE applying (uses pre-move positions); commit on success.
static vector<size_t> gpend;
static void collectTabu(const vector<vector<int>>& cur, const Mv& mv){
    gpend.clear();
    const auto& s = cur[mv.m];
    int uj = s[mv.i];
    int uop = opOf(uj, mv.m);
    if(mv.front){
        // reversed: (x before u) -> forbid TB[uop][x]
        for(int t=mv.b; t<mv.i; ++t) gpend.push_back((size_t)uop*J + s[t]);
    } else {
        // reversed: (u before x) -> forbid TB[xop][uj]
        for(int t=mv.i+1; t<=mv.e; ++t)
            gpend.push_back((size_t)opOf(s[t], mv.m)*J + uj);
    }
}

// Random adjacent swaps on a few machines (always tries to stay feasible).
static void perturbRandom(vector<vector<int>>& s, mt19937& rng, int nSwaps){
    for(int t=0;t<nSwaps;++t){
        int m = (int)(rng() % M);
        if(J < 2) break;
        int i = (int)(rng() % (J-1));
        swap(s[m][i], s[m][i+1]);
    }
}

// Path-relink-ish: for each machine, copy order from A or B with coin flip,
// then repair via a fresh GT-like rebuild is too heavy; instead just pick
// whole-machine orders from A/B (always feasible if A,B are).
static vector<vector<int>> mixElites(const vector<vector<int>>& A,
                                     const vector<vector<int>>& B,
                                     mt19937& rng){
    vector<vector<int>> s = A;
    for(int m=0;m<M;++m) if(rng() & 1) s[m] = B[m];
    return s;
}

int main(){
    auto T0 = chrono::steady_clock::now();
    // Leave ~60-80ms for final I/O + margin under the 1s TL.
    const auto budget = chrono::milliseconds(920);

    if(scanf("%d %d", &J, &M) != 2) return 0;
    N = J*M;
    m_of.assign(J, vector<int>(M));
    p_of.assign(J, vector<long long>(M));
    pos.assign(J, vector<int>(M));
    for(int j=0;j<J;++j)
        for(int k=0;k<M;++k)
            if(scanf("%d %lld", &m_of[j][k], &p_of[j][k]) != 2) return 0;
    for(int j=0;j<J;++j) for(int k=0;k<M;++k) pos[j][m_of[j][k]] = k;
    posF.assign(N, 0);
    for(int j=0;j<J;++j) for(int m=0;m<M;++m) posF[j*M+m] = pos[j][m];
    pnode.assign(N, 0);
    for(int j=0;j<J;++j) for(int k=0;k<M;++k) pnode[j*M+k] = p_of[j][k];

    msucc.assign(N, -1);
    mpred.assign(N, -1);
    jsuc.assign(N, -1); jpre.assign(N, 0);
    for(int u=0;u<N;++u){
        if(u % M != M-1) jsuc[u] = u+1;
        if(u % M != 0) jpre[u] = 1;
    }
    indeg.resize(N); dist_.resize(N); tail_.resize(N); crit.resize(N);
    gwl.reserve(4*N); ginq.assign(N, 0);
    qbuf.reserve(N);
    gord.resize(J); gestC.resize(J);
    tabuTB.assign((size_t)N*J, 0);

    // Feasibility fallback: job-index order on every machine (always acyclic).
    vector<vector<int>> best(M, vector<int>(J));
    for(int m=0;m<M;++m) for(int j=0;j<J;++j) best[m][j] = j;
    long long bestC = evalSeq(best);

    // Trivial case: a single job has no ordering choices -- output is optimal.
    if(J <= 1){
        for(int m=0;m<M;++m){ printf("0\n"); }
        fflush(stdout);
        _exit(0);
    }

    vector<vector<int>> cur = best;
    long long curC = bestC;

    // Elite pool (up to 6 diverse good solutions) for multi-start / path relinking.
    struct Elite { long long C; vector<vector<int>> s; };
    vector<Elite> elite;
    auto pushElite = [&](const vector<vector<int>>& s, long long C){
        if(C <= 0) return;
        for(auto& e : elite) if(e.C == C){
            // cheap diversity: reject exact same makespan twin
            return;
        }
        elite.push_back({C, s});
        sort(elite.begin(), elite.end(), [](const Elite& a, const Elite& b){ return a.C < b.C; });
        if((int)elite.size() > 6) elite.resize(6);
    };

    auto trySeed = [&](const vector<vector<int>>& s){
        long long c = evalSeq(s);
        if(c > 0){
            pushElite(s, c);
            if(c < curC){
                cur = s; curC = c;
                if(c < bestC){ best = s; bestC = c; }
            }
        }
    };

    // Seed from multiple GT priorities + a few pure-random GT runs.
    mt19937 rng((unsigned)chrono::steady_clock::now().time_since_epoch().count() ^ 0x9e3779b9u);
    trySeed(seedGT(0, rng)); // MWR
    trySeed(seedGT(1, rng)); // LPT
    trySeed(seedGT(2, rng)); // SPT
    for(int r=0; r<4 && chrono::steady_clock::now() - T0 < budget/8; ++r)
        trySeed(seedGT(3, rng)); // random priority GT

    // Trivial lower bound: max(machine load, job length). If reached, we are
    // provably optimal and can stop immediately.
    long long LB = 0;
    {
        vector<long long> mload(M, 0);
        for(int j=0;j<J;++j){
            long long jl = 0;
            for(int k=0;k<M;++k){ jl += p_of[j][k]; mload[m_of[j][k]] += p_of[j][k]; }
            if(jl > LB) LB = jl;
        }
        for(int m=0;m<M;++m) if(mload[m] > LB) LB = mload[m];
    }

    auto T_end = T0 + budget;

    // Diversify: jump to a new start (GT / elite mix / best+kick) and clear tabu.
    int restartCnt = 0;
    auto diversify = [&](bool strong){
        restartCnt++;
        fill(tabuTB.begin(), tabuTB.end(), 0);
        int mode = (int)(rng() % (strong ? 5 : 3));
        if(mode == 0 && !elite.empty()){
            // Restart from a random elite (prefer better ones).
            int ei = (int)(rng() % elite.size());
            if((rng() & 3) == 0) ei = 0; // bias toward global-best elite
            cur = elite[ei].s;
            // Strong: also mix with another elite (path relink).
            if(strong && elite.size() >= 2 && (rng() & 1)){
                int ej = (int)(rng() % elite.size());
                if(ej != ei) cur = mixElites(cur, elite[ej].s, rng);
            }
        } else if(mode == 1){
            cur = seedGT((int)(rng() % 4), rng);
        } else if(mode == 2){
            cur = best;
        } else if(mode == 3){
            cur = seedGT(0, rng); // MWR
        } else {
            cur = best;
            perturbRandom(cur, rng, 2 + (int)(rng() % max(2, J/4)));
        }
        long long c = evalSeq(cur, true);
        if(c < 0){
            cur = best;
            c = evalSeq(cur, true);
        }
        curC = c;
        if(c > 0 && c < bestC){ best = cur; bestC = c; pushElite(best, bestC); }
        // Extra critical-block kicks from the chosen start.
        int kicks = (strong ? 3 : 1) + (int)(rng() % (strong ? 5 : 3));
        for(int r=0; r<kicks; ++r){
            if(chrono::steady_clock::now() >= T_end) break;
            genMoves(cur);
            if(gmoves.empty()){
                perturbRandom(cur, rng, 1 + (int)(rng() % 3));
                long long nc = evalSeq(cur, true);
                if(nc < 0){ cur = best; evalSeq(cur, true); curC = bestC; break; }
                curC = nc;
                continue;
            }
            const Mv& mv = gmoves[rng() % gmoves.size()];
            applyMove(cur, mv);
            long long nc = evalSeq(cur, true);
            if(nc < 0){ undoMove(cur, mv); evalSeq(cur, true); }
            else {
                curC = nc;
                if(nc < bestC){ best = cur; bestC = nc; pushElite(best, bestC); }
            }
        }
    };

    // ---- Tabu search over the N7 neighborhood with multi-start ILS ----
    // One exact/incremental eval per iteration; candidates ranked by O(segment)
    // head/tail re-estimation. Diversify on stuck / empty neighborhood / failure.
    if(bestC > LB){
        long long c = evalSeq(cur, true);
        if(c > 0) curC = c; else { cur = best; curC = evalSeq(cur, true); }
        pushElite(best, bestC);

        int iter = 0, sinceImp = 0;
        // Adaptive stuck limit: smaller instances get more diversification;
        // large N can burn many cheap iters before kicking.
        const int stuckLim = max(8000, min(40000, 250000 / max(1, N)));
        const int TENURE_MIN = max(5, J / 6);
        const int TENURE_SPAN = max(6, J / 2);
        bool timeUp = false;
        int topK = (N <= 400 ? 32 : 20); // rank more candidates on smaller graphs

        while(!timeUp && chrono::steady_clock::now() < T_end){
            iter++;
            if((iter & 16383) == 0){ // periodic exact refresh (drift insurance)
                long long fc = evalSeq(cur, true);
                if(fc >= 0) curC = fc;
            }
            genMoves(cur);
            int nmv = (int)gmoves.size();
            if(nmv == 0){
                // No critical-block moves: diversify rather than halt.
                if(chrono::steady_clock::now() >= T_end) break;
                diversify(true);
                sinceImp = 0;
                continue;
            }
            gcand.clear();
            for(int idx=0; idx<nmv; ++idx)
                gcand.push_back({estMove(cur, gmoves[idx]), idx});
            int K = nmv < topK ? nmv : topK;
            partial_sort(gcand.begin(), gcand.begin()+K, gcand.end());
            bool sorted_all = (K == nmv);

            bool applied = false;
            for(int pass=0; pass<2 && !applied && !timeUp; ++pass){
                for(int t=0; t<nmv; ++t){
                    if((t & 7)==7 && chrono::steady_clock::now() >= T_end){ timeUp = true; break; }
                    if(t >= K && !sorted_all){
                        sort(gcand.begin(), gcand.end());
                        sorted_all = true;
                    }
                    const Mv& mv = gmoves[gcand[t].second];
                    if(pass==0){
                        bool tb = isTabu(cur, mv, iter);
                        bool asp = gcand[t].first < bestC;
                        if(tb && !asp) continue;
                    }
                    collectTabu(cur, mv);
                    applyMove(cur, mv);
                    long long nc = incAfterMove(cur, mv);
                    if(nc == -2) nc = evalSeq(cur, true); // cap hit: exact recompute
                    if(nc < 0){
                        undoMove(cur, mv);
                        evalSeq(cur, true); // restore dist/tail/crit for cur
                        continue;
                    }
                    int tenure = TENURE_MIN + (int)(rng() % TENURE_SPAN);
                    for(size_t id : gpend) tabuTB[id] = iter + tenure;
                    curC = nc;
                    applied = true;
                    break;
                }
            }
            if(!applied){
                // All moves infeasible or timed out: diversify if time remains.
                if(timeUp || chrono::steady_clock::now() >= T_end) break;
                diversify(true);
                sinceImp = 0;
                continue;
            }

            if(curC < bestC){
                // Confirm with an exact evaluation before recording a new best
                // (also refreshes dist/tail/crit exactly, washing out any drift).
                long long exact = evalSeq(cur, true);
                curC = exact;
                if(exact >= 0 && exact < bestC){
                    best = cur; bestC = exact; sinceImp = 0;
                    pushElite(best, bestC);
                    if(bestC <= LB) break; // provably optimal
                }
            }
            else if(++sinceImp > stuckLim){
                // Escalating diversification: mild kick first, then strong multi-start.
                bool strong = (sinceImp > stuckLim * 2) || ((restartCnt & 3) == 3);
                diversify(strong);
                sinceImp = 0;
            }
        }
    }

    // Fast buffered output (single fwrite) then _exit to skip static-vector
    // teardown, which otherwise adds ~100-200ms of wall time on some systems.
    {
        vector<char> buf;
        buf.reserve((size_t)N*8 + M + 16);
        for(int m=0;m<M;++m){
            for(int j=0;j<J;++j){
                int x = best[m][j];
                if(x == 0){ buf.push_back('0'); }
                else { char tmp[12]; int t = 0; while(x > 0){ tmp[t++] = char('0' + x%10); x /= 10; } while(t > 0) buf.push_back(tmp[--t]); }
                buf.push_back(j+1<J ? ' ' : '\n');
            }
        }
        fwrite(buf.data(), 1, buf.size(), stdout);
        fflush(stdout);
    }
    _exit(0);
}
