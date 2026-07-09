#include <cstdio>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <climits>
#include <cmath>
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

    // ---- Tabu-primary ILS with simulated-annealing meta-acceptance (v5_sa) --
    // The workhorse is the incumbent's proven N7 critical-block TABU descent
    // (strong exploitation). We wrap it in an iterated-local-search shell whose
    // acceptance criterion is SIMULATED ANNEALING: after each tabu burst lands
    // on a local optimum, we accept it as the new walk anchor if it improves,
    // or with Metropolis probability exp(-delta/T) if it is worse -- letting the
    // search DRIFT across the landscape and climb out of the basins a pure
    // best-only tabu (the incumbent) stays trapped in. A geometric cooling
    // schedule tightens acceptance over time; a reheat + stronger kick fires on
    // deep stagnation. Perturbation strength ("kick") scales with how long we
    // have gone without a global best (adaptive ILS).
    if(bestC > LB){
        long long c = evalSeq(cur, true);
        if(c > 0) curC = c; else { cur = best; curC = evalSeq(cur, true); }

        // Fast exp for the Metropolis test (argument <= 0).
        auto fastExp = [](double x)->double{
            if(x < -20.0) return 0.0;
            double y = x * 1.4426950408889634; // x / ln2
            double fl = floor(y);
            double fr = y - fl;
            double p = 1.0 + fr*(0.6931472 + fr*(0.2402265 + fr*0.0555041));
            int e = (int)fl;
            union { double d; unsigned long long u; } bits;
            bits.d = p;
            long long expo = (long long)((bits.u >> 52) & 0x7ff) + e;
            if(expo <= 0) return 0.0;
            if(expo >= 2047) return 1e300;
            bits.u = (bits.u & 0x800fffffffffffffULL) | ((unsigned long long)expo << 52);
            return bits.d;
        };

        bool timeUp = false;

        // ---- Stagnation-reactive tabu tenure (mirrors the incumbent) --------
        // Short tenure while improving (exploit); widen as iterations since the
        // last burst-local improvement grow (anti-cycling), snap back on gain.
        const int spanShort = max(4, J/3);
        const int spanLong  = max(4, J/2);
        auto dynTen = [&](int stag)->int{
            int tmin, span;
            if(stag <= 20){ tmin = 8; span = spanShort; }
            else if(stag >= 120){ tmin = 15 + (stag-120)/250; if(tmin>21) tmin=21; span = spanLong; }
            else {
                int f = (stag-20)*100/100;
                tmin = 8 + (15-8)*f/100;
                span = spanShort + (spanLong-spanShort)*f/100;
            }
            return tmin + (int)(rng() % (unsigned)span);
        };

        // ---- A single tabu descent burst -----------------------------------
        // Runs from `cur` for up to `budgetIters` accepted moves (or until it
        // stalls / time is up). Updates global best/bestC. Leaves cur at the
        // final walked point and curC exact. Returns the makespan reached.
        auto tabuBurst = [&](int budgetIters, chrono::steady_clock::time_point dl)->long long{
            fill(tabuTB.begin(), tabuTB.end(), 0);
            long long burstBest = curC; int stag = 0;
            for(int bi=0; bi<budgetIters; ++bi){
                if((bi & 7)==7){
                    auto tn = chrono::steady_clock::now();
                    if(tn >= T_end){ timeUp = true; break; }
                    if(tn >= dl) break;
                }
                genMoves(cur);
                int nmv = (int)gmoves.size();
                if(nmv == 0) break;
                gcand.clear();
                for(int idx=0; idx<nmv; ++idx)
                    gcand.push_back({estMove(cur, gmoves[idx]), idx});
                int K = nmv < 24 ? nmv : 24;
                partial_sort(gcand.begin(), gcand.begin()+K, gcand.end());
                bool sorted_all = (K == nmv);
                bool applied = false;
                for(int pass=0; pass<2 && !applied && !timeUp; ++pass){
                    for(int t=0; t<nmv; ++t){
                        if((t & 7)==7 && chrono::steady_clock::now() >= T_end){ timeUp = true; break; }
                        if(t >= K && !sorted_all){ sort(gcand.begin(), gcand.end()); sorted_all = true; }
                        const Mv& mv = gmoves[gcand[t].second];
                        if(pass==0){
                            bool tb = isTabu(cur, mv, bi);
                            bool asp = gcand[t].first < bestC;
                            if(tb && !asp) continue;
                        }
                        collectTabu(cur, mv);
                        applyMove(cur, mv);
                        long long nc = incAfterMove(cur, mv);
                        if(nc == -2) nc = evalSeq(cur, true);
                        if(nc < 0){ undoMove(cur, mv); evalSeq(cur, true); continue; }
                        int tenure = dynTen(stag);
                        for(size_t id : gpend) tabuTB[id] = bi + tenure;
                        curC = nc; applied = true; break;
                    }
                }
                if(!applied) break;
                if(curC < burstBest){ burstBest = curC; stag = 0; } else stag++;
                if(curC < bestC){
                    long long exact = evalSeq(cur, true);
                    curC = exact;
                    if(exact >= 0 && exact < bestC){
                        best = cur; bestC = exact;
                        if(exact < burstBest) burstBest = exact;
                        if(bestC <= LB){ timeUp = true; break; }
                    }
                }
            }
            return curC;
        };

        // ---- Perturbation: random critical-block moves from `from` ----------
        auto kickFrom = [&](vector<vector<int>>& from, int kicks){
            cur = from;
            long long cc = evalSeq(cur, true);
            curC = cc;
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
        };

        // ---- SA meta-schedule ---------------------------------------------
        double avgP = 0.0;
        for(int u=0;u<N;++u) avgP += (double)pnode[u];
        avgP = avgP / (double)N;
        if(avgP < 1.0) avgP = 1.0;
#ifndef T0_FRAC
#define T0_FRAC 0.35
#endif
#ifndef TMIN_FRAC
#define TMIN_FRAC 0.02
#endif
#ifndef ALPHA
#define ALPHA 0.90
#endif
#ifndef BURST_ITERS
#define BURST_ITERS 1000
#endif
#ifndef INIT_FRAC
#define INIT_FRAC 70
#endif
        // Reserve the tail of the time budget for ILS/SA rounds; the initial
        // deep tabu descent gets the first INIT_FRAC% of the budget.
        auto T_ils = T0 + chrono::milliseconds((long long)995 * INIT_FRAC / 100);
        double Thi = avgP * (double)T0_FRAC;
        double Tmin = avgP * (double)TMIN_FRAC;
        double T = Thi;
        const double alpha = (double)ALPHA;

        // Anchor of the SA walk (the accepted point we perturb from).
        vector<vector<int>> anchor = best;
        long long anchorC = bestC;
        int noGlobalImp = 0;   // ILS rounds since last global best

        // First: a deep tabu descent (time-bounded) to reach a strong local
        // optimum -- this is the exploitation phase that matches the incumbent.
        tabuBurst(2000000000, T_ils);
        if(curC < anchorC){ anchor = cur; anchorC = curC; }

        // ---- Main ILS/SA loop ---------------------------------------------
        while(!timeUp && bestC > LB){
            if(chrono::steady_clock::now() >= T_end) break;

            // Perturb the anchor; harder kicks the longer we stagnate.
            int kicks = 2 + noGlobalImp/2;
            if(kicks > 6) kicks = 6;
            kickFrom(anchor, kicks);
            if(timeUp) break;

            long long before = bestC;
            long long reached = tabuBurst(BURST_ITERS, T_end);
            if(timeUp) break;

            // SA acceptance at the ILS level: compare the local optimum reached
            // to the current anchor. Accept if better, else with Metropolis
            // probability -- this is the annealing that lets the walk drift.
            long long delta = reached - anchorC;
            bool accept;
            if(delta <= 0) accept = true;
            else {
                double pr = fastExp(-(double)delta / (T > 1e-9 ? T : 1e-9));
                accept = ((double)(rng() & 0xffffff) / 16777216.0) < pr;
            }
            if(accept){ anchor = cur; anchorC = reached; }

            if(bestC < before) noGlobalImp = 0; else noGlobalImp++;

            // Cool; reheat + re-anchor on the incumbent best after deep stalls.
            T *= alpha;
            if(T < Tmin){
                T = Thi;
                anchor = best; anchorC = bestC; // re-center exploration on best
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
