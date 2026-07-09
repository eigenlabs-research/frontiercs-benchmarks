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

int main(){
    auto T0 = chrono::steady_clock::now();
    const auto budget = chrono::milliseconds(995); // increased budget for more iterations with longer tenure

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

    auto trySeed = [&](const vector<vector<int>>& s){
        long long c = evalSeq(s);
        if(c > 0 && c < curC){
            cur = s; curC = c;
            if(c < bestC){ best = s; bestC = c; }
        }
    };

    mt19937 rng(777u);
    trySeed(seedGT(0, rng)); // MWR
    trySeed(seedGT(1, rng)); // LPT
    if(chrono::steady_clock::now() - T0 < budget)
        trySeed(seedGT(2, rng)); // SPT

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

    // ---- Tabu search over the N7 neighborhood with ILS kicks ----
    // One exact eval per iteration (on the accepted move); all candidate moves
    // ranked by an O(segment) head/tail re-estimation.
    if(bestC > LB){
        long long c = evalSeq(cur, true);
        if(c > 0) curC = c; else { cur = best; curC = evalSeq(cur, true); }

        int iter = 0, sinceImp = 0;
#ifndef STUCK_LIM
#define STUCK_LIM 1000000000
#endif
#ifndef TEN_MIN
#define TEN_MIN 15
#endif
#ifndef TEN_SPAN_DIV
#define TEN_SPAN_DIV 2
#endif
        const int stuckLim = STUCK_LIM;
        const int TENURE_MIN = TEN_MIN;
        const int TENURE_SPAN = max(4, J/TEN_SPAN_DIV);
        bool timeUp = false;

        // Stagnation-reactive oscillating tenure band: exploit with a short
        // tenure while global bests keep arriving; widen toward the long
        // regime as time since the last improvement grows (anti-cycling),
        // snapping back on any new best.
#ifndef DYN_S1
#define DYN_S1 30
#endif
#ifndef DYN_S2
#define DYN_S2 120
#endif
#ifndef DYN_ESC
#define DYN_ESC 250
#endif
        const int spanShort = max(4, J/3);
        const int spanLong  = TENURE_SPAN;
        auto lastImpT = chrono::steady_clock::now();
        auto dynTen = [&](chrono::steady_clock::time_point nowT)->int{
            long long stag = chrono::duration_cast<chrono::milliseconds>(nowT - lastImpT).count();
            int tmin, span;
            if(stag <= DYN_S1){ tmin = 8; span = spanShort; }
            else if(stag >= DYN_S2){
                tmin = TENURE_MIN + (int)min(6LL, (stag - DYN_S2)/DYN_ESC);
                span = spanLong;
            } else {
                int f = (int)((stag - DYN_S1)*100/(DYN_S2 - DYN_S1)); // 0..100
                tmin = 8 + (TENURE_MIN - 8)*f/100;
                span = spanShort + (spanLong - spanShort)*f/100;
            }
            return tmin + (int)(rng() % (unsigned)span);
        };

#ifdef DIAG
        int iterLastImp = 0;
#endif
        while(!timeUp){
            auto nowT = chrono::steady_clock::now();
            if(nowT >= T_end) break;
            iter++;
            if((iter & 16383) == 0){ // periodic exact refresh (drift insurance)
                long long fc = evalSeq(cur, true);
                if(fc >= 0) curC = fc;
            }
            genMoves(cur);
            int nmv = (int)gmoves.size();
#ifdef DIAG
            static long long totMv = 0; totMv += nmv;
            if(iter % 50000 == 0) fprintf(stderr, "avg nmv=%.1f\n", (double)totMv/iter);
#endif
            if(nmv == 0) break;
            gcand.clear();
            for(int idx=0; idx<nmv; ++idx)
                gcand.push_back({estMove(cur, gmoves[idx]), idx});
            int K = nmv < 24 ? nmv : 24;
            partial_sort(gcand.begin(), gcand.begin()+K, gcand.end());
            bool sorted_all = (K == nmv);

            bool applied = false;
#ifndef EVAL_TOP
#define EVAL_TOP 1
#endif
#if EVAL_TOP > 1
            // Exact-evaluate the EVAL_TOP best admissible candidates, keep the best.
            {
                int bestIdx = -1; long long bestNC = -1; int evald = 0;
                for(int t=0; t<nmv && evald<EVAL_TOP; ++t){
                    if((t & 7)==7 && chrono::steady_clock::now() >= T_end){ timeUp = true; break; }
                    if(t >= K && !sorted_all){ sort(gcand.begin(), gcand.end()); sorted_all = true; }
                    const Mv& mv = gmoves[gcand[t].second];
                    bool tb = isTabu(cur, mv, iter);
                    bool asp = gcand[t].first < bestC;
                    if(tb && !asp) continue;
                    applyMove(cur, mv);
                    long long nc = evalSeq(cur, false); // trashes dist_ only
                    undoMove(cur, mv);
                    if(nc < 0) continue;
                    evald++;
                    if(bestIdx < 0 || nc < bestNC){ bestNC = nc; bestIdx = gcand[t].second; }
                }
                if(bestIdx >= 0){
                    const Mv& mv = gmoves[bestIdx];
                    collectTabu(cur, mv);
                    applyMove(cur, mv);
                    long long nc = evalSeq(cur, true);
                    if(nc >= 0){
                        int tenure = dynTen(nowT);
                        for(size_t id : gpend) tabuTB[id] = iter + tenure;
                        curC = nc; applied = true;
                    } else { undoMove(cur, mv); evalSeq(cur, true); }
                }
            }
#endif
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
#ifdef VERIFY
                    {
                        static vector<long long> vd, vt; static long long mism = 0; static long long checks = 0;
                        vd = dist_; vt = tail_;
                        long long fc = evalSeq(cur, true);
                        checks++;
                        if(nc != -2){
                            if(fc != nc || vd != dist_ || vt != tail_){
                                mism++;
                                fprintf(stderr, "MISMATCH iter=%d inc=%lld full=%lld distOK=%d tailOK=%d\n",
                                        iter, nc, fc, (int)(vd==dist_), (int)(vt==tail_));
                            }
                        }
                        if(checks % 20000 == 0) fprintf(stderr, "verify checks=%lld mism=%lld\n", checks, mism);
                        nc = fc;
                    }
#endif
                    if(nc == -2) nc = evalSeq(cur, true); // cap hit: exact recompute
                    if(nc < 0){
                        undoMove(cur, mv);
                        evalSeq(cur, true); // restore dist/tail/crit for cur
                        continue;
                    }
                    int tenure = dynTen(nowT);
                    for(size_t id : gpend) tabuTB[id] = iter + tenure;
                    curC = nc;
                    applied = true;
                    break;
                }
            }
            if(!applied) break; // timed out or no feasible move at all

            if(curC < bestC){
                // Confirm with an exact evaluation before recording a new best
                // (also refreshes dist/tail/crit exactly, washing out any drift).
                long long exact = evalSeq(cur, true);
                curC = exact;
                if(exact >= 0 && exact < bestC){
                    best = cur; bestC = exact; sinceImp = 0;
                    lastImpT = chrono::steady_clock::now();
#ifdef DIAG
                    iterLastImp = iter;
#endif
                    if(bestC <= LB) break; // provably optimal
                }
            }
            else if(++sinceImp > stuckLim){
                // ILS kick: restart from best with a few random critical-block moves.
                cur = best;
                long long cc = evalSeq(cur, true);
                curC = cc;
#ifndef KICK_BASE
#define KICK_BASE 1
#define KICK_RAND 2
#endif
                int kicks = KICK_BASE + (KICK_RAND ? (int)(rng() % KICK_RAND) : 0);
                for(int r=0; r<kicks; ++r){
                    if(chrono::steady_clock::now() >= T_end){ timeUp = true; break; }
                    genMoves(cur);
                    if(gmoves.empty()) break;
                    const Mv& mv = gmoves[rng() % gmoves.size()];
                    applyMove(cur, mv);
                    long long nc = evalSeq(cur, true);
                    if(nc < 0){ undoMove(cur, mv); evalSeq(cur, true); }
                    else curC = nc;
                }
                fill(tabuTB.begin(), tabuTB.end(), 0);
                sinceImp = 0;
            }
        }
#ifdef DIAG
        extern long long g_pops, g_calls;
        fprintf(stderr, "iters=%d lastImp=%d bestC=%lld avgPops=%.1f (N=%d)\n", iter, iterLastImp, bestC, g_calls? (double)g_pops/g_calls : 0.0, N);
#endif
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
