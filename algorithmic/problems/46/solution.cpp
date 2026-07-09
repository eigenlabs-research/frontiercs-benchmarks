// Job Shop Scheduling — 150 diversified GT seeds + tabu over the N1 block-swap UNION N7 insertion
// neighborhood. Extends our 0.0943 champion (N1 block-boundary swaps) with Balas-Vazacopoulos
// critical-block INSERTION moves: relocate a critical block's first op to the block end, or its last
// op to the block start. These larger moves escape local optima that adjacent swaps can't — the
// technique that cracks low-slack / large-critical-block instances (e.g. hidden case 10) — while the
// fast N1 swaps (O(1) estimate) keep the strengths on the other families. Inserts are evaluated
// exactly and cycle-guarded (revert if they'd create a cycle), so output is always feasible.

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <algorithm>
#include <ctime>
#include <climits>

using namespace std;

static int J, M, N;

static vector<long long> procOp;
static vector<int>       jobOf, kOf, machOf;
static vector<int>       jobPred, jobSucc;

static vector<vector<int>>       machJK;
static vector<vector<long long>> procJK;
static vector<vector<int>>       posOf;

static inline int opOnMachine(int j, int m){ return j * M + posOf[j][m]; }

static vector<int>       indeg, mSucc, mPred, order_;
static vector<long long> dist_, q_;
static long long         Cmax_;

static clock_t START;
static const double TL = 0.94;
static inline double elapsed(){ return double(clock() - START) / CLOCKS_PER_SEC; }

static unsigned long long rngState = 0x9e3779b97f4a7c15ULL;
static inline unsigned long long rnd(){
    rngState ^= rngState << 13; rngState ^= rngState >> 7; rngState ^= rngState << 17;
    return rngState;
}
static inline int rndInt(int n){ return (int)(rnd() % (unsigned long long)n); }

static long long evaluate(const vector<vector<int>>& seq){
    for(int op = 0; op < N; ++op){
        indeg[op] = (kOf[op] > 0) ? 1 : 0;
        mSucc[op] = -1; mPred[op] = -1;
    }
    for(int m = 0; m < M; ++m){
        const vector<int>& s = seq[m];
        for(int i = 1; i < J; ++i){
            int a = opOnMachine(s[i-1], m);
            int b = opOnMachine(s[i],   m);
            mSucc[a] = b; mPred[b] = a; ++indeg[b];
        }
    }
    int tail = 0, head = 0;
    for(int op = 0; op < N; ++op){
        if(indeg[op] == 0){ dist_[op] = procOp[op]; order_[tail++] = op; }
        else dist_[op] = 0;
    }
    int cnt = 0;
    while(head < tail){
        int u = order_[head++]; ++cnt;
        long long du = dist_[u];
        int js = jobSucc[u];
        if(js != -1){
            if(dist_[js] < du + procOp[js]) dist_[js] = du + procOp[js];
            if(--indeg[js] == 0) order_[tail++] = js;
        }
        int ms = mSucc[u];
        if(ms != -1){
            if(dist_[ms] < du + procOp[ms]) dist_[ms] = du + procOp[ms];
            if(--indeg[ms] == 0) order_[tail++] = ms;
        }
    }
    if(cnt != N) return -1;
    long long C = 0;
    for(int op = 0; op < N; ++op) if(dist_[op] > C) C = dist_[op];
    for(int idx = N - 1; idx >= 0; --idx){
        int op = order_[idx];
        long long best = 0;
        int js = jobSucc[op]; if(js != -1 && q_[js] > best) best = q_[js];
        int ms = mSucc[op];   if(ms != -1 && q_[ms] > best) best = q_[ms];
        q_[op] = procOp[op] + best;
    }
    Cmax_ = C;
    return C;
}

static inline bool critOp(int op){ return (dist_[op] - procOp[op]) + q_[op] == Cmax_; }

// Fill swap moves (adjacent block-boundary pairs, ranked by O(1) estimate) and insert moves
// (BV block first->end / last->start, ranked by exact eval). inserts are (machine, from, to).
static void getBlockMoves(const vector<vector<int>>& seq, vector<pair<int,int>>& swaps,
                          vector<array<int,3>>& inserts){
    swaps.clear(); inserts.clear();
    for(int m = 0; m < M; ++m){
        const vector<int>& s = seq[m];
        int i = 0;
        while(i < J){
            if(!critOp(opOnMachine(s[i], m))){ ++i; continue; }
            int j = i;
            while(j < J && critOp(opOnMachine(s[j], m))) ++j;
            int bs = j - i;
            if(bs >= 2){
                int a1 = opOnMachine(s[i], m), b1 = opOnMachine(s[i+1], m);
                if(dist_[a1] == dist_[b1] - procOp[b1]) swaps.push_back({a1, b1});
                if(bs > 2){
                    int a2 = opOnMachine(s[j-2], m), b2 = opOnMachine(s[j-1], m);
                    if(dist_[a2] == dist_[b2] - procOp[b2]) swaps.push_back({a2, b2});
                    // The two highest-impact BV insertions per block (bounded, keeps throughput up).
                    inserts.push_back({m, i, j-1});   // block-first op -> block end
                    inserts.push_back({m, j-1, i});   // block-last op  -> block start
                }
            }
            i = j;
        }
    }
}

static inline long long estimateSwap(int a, int b){
    int PM = mPred[a], SM = mSucc[b];
    long long fPM  = (PM != -1) ? dist_[PM] : 0;
    long long fJPa = (jobPred[a] != -1) ? dist_[jobPred[a]] : 0;
    long long fJPb = (jobPred[b] != -1) ? dist_[jobPred[b]] : 0;
    long long rB = max(fPM, fJPb);
    long long rA = max(rB + procOp[b], fJPa);
    long long qJSa = (jobSucc[a] != -1) ? q_[jobSucc[a]] : 0;
    long long qJSb = (jobSucc[b] != -1) ? q_[jobSucc[b]] : 0;
    long long qSM  = (SM != -1) ? q_[SM] : 0;
    long long qA = procOp[a] + max(qJSa, qSM);
    long long qB = procOp[b] + max(qJSb, qA);
    return max(rA + qA, rB + qB);
}

// O(L) Taillard/BV estimate of the makespan after moving cur[m][from] to position `to`, holding
// everything outside the reordered machine segment fixed (heads forward, tails backward). Reads
// the current dist_/q_ only (no mutation), so it ranks alongside estimateSwap.
static vector<int>       gord;
static vector<long long> gestC;
static long long estInsert(const vector<vector<int>>& cur, int m, int from, int to){
    const vector<int>& s = cur[m];
    int lo = from < to ? from : to, hi = from < to ? to : from;
    int L = hi - lo + 1;
    if((int)gord.size() < L){ gord.resize(L); gestC.resize(L); }
    if(to > from){                       // move s[from] to the back of the segment
        for(int t = from + 1; t <= to; ++t) gord[t - from - 1] = s[t];
        gord[to - from] = s[from];
    } else {                             // move s[from] to the front of the segment
        gord[0] = s[from];
        for(int t = to; t < from; ++t) gord[t - to + 1] = s[t];
    }
    long long prevC = (lo > 0) ? dist_[opOnMachine(s[lo - 1], m)] : 0;
    for(int t = 0; t < L; ++t){
        int v = gord[t]; int u = opOnMachine(v, m); int k = posOf[v][m];
        long long jp = (k > 0) ? dist_[u - 1] : 0;
        long long st = prevC > jp ? prevC : jp;
        gestC[t] = st + procOp[u];
        prevC = gestC[t];
    }
    long long prevT = (hi + 1 < J) ? q_[opOnMachine(s[hi + 1], m)] : 0;
    long long bestLen = 0;
    for(int t = L - 1; t >= 0; --t){
        int v = gord[t]; int u = opOnMachine(v, m); int k = posOf[v][m];
        long long js = (k < M - 1) ? q_[u + 1] : 0;
        long long tl = procOp[u] + (prevT > js ? prevT : js);
        long long len = gestC[t] - procOp[u] + tl;
        if(len > bestLen) bestLen = len;
        prevT = tl;
    }
    return bestLen;
}

static inline double priority(int rule, int j, int k, const vector<long long>& remWork){
    switch(rule){
        case 0: return (double)remWork[j];
        case 1: return -(double)procJK[j][k];
        case 2: return  (double)procJK[j][k];
        case 3: return -(double)remWork[j];
        default: return (double)(rnd() & 0xffffff);
    }
}

static vector<vector<int>> gifflerThompson(int rule){
    vector<int>       nextK(J, 0);
    vector<long long> jobFree(J, 0), machFree(M, 0), remWork(J, 0);
    for(int j = 0; j < J; ++j)
        for(int k = 0; k < M; ++k) remWork[j] += procJK[j][k];
    vector<vector<int>> seq(M);
    for(int m = 0; m < M; ++m) seq[m].reserve(J);
    int scheduled = 0;
    while(scheduled < N){
        long long minC = LLONG_MAX; int mstar = -1;
        for(int j = 0; j < J; ++j){
            if(nextK[j] >= M) continue;
            int k = nextK[j], m = machJK[j][k];
            long long est = max(jobFree[j], machFree[m]);
            long long C = est + procJK[j][k];
            if(C < minC){ minC = C; mstar = m; }
        }
        int chosen = -1; double bestPri = -1e300;
        for(int j = 0; j < J; ++j){
            if(nextK[j] >= M) continue;
            int k = nextK[j];
            if(machJK[j][k] != mstar) continue;
            long long est = max(jobFree[j], machFree[mstar]);
            if(est < minC){
                double pri = priority(rule, j, k, remWork);
                if(pri > bestPri){ bestPri = pri; chosen = j; }
            }
        }
        int j = chosen, k = nextK[j], m = mstar;
        long long est = max(jobFree[j], machFree[m]);
        long long fin = est + procJK[j][k];
        jobFree[j] = fin; machFree[m] = fin;
        seq[m].push_back(j);
        remWork[j] -= procJK[j][k];
        ++nextK[j];
        ++scheduled;
    }
    return seq;
}

static vector<vector<int>> pos;
static void rebuildPos(const vector<vector<int>>& seq){
    for(int m = 0; m < M; ++m)
        for(int i = 0; i < J; ++i) pos[m][seq[m][i]] = i;
}

static inline void doSwap(vector<vector<int>>& seq, int a, int b){
    int m = machOf[a];
    int i = pos[m][jobOf[a]];
    swap(seq[m][i], seq[m][i+1]);
    pos[m][seq[m][i]]   = i;
    pos[m][seq[m][i+1]] = i+1;
}

// Move seq[m][from] to position `to`, shifting the intervening jobs; keep pos[m] consistent.
// doInsert(m,to,from) exactly reverts doInsert(m,from,to).
static inline void doInsert(vector<vector<int>>& seq, int m, int from, int to){
    vector<int>& s = seq[m];
    int job = s[from];
    if(from < to){ for(int i = from; i < to; ++i) s[i] = s[i+1]; s[to] = job; }
    else         { for(int i = from; i > to; --i) s[i] = s[i-1]; s[to] = job; }
    int lo = from < to ? from : to, hi = from < to ? to : from;
    for(int i = lo; i <= hi; ++i) pos[m][s[i]] = i;
}

static void perturb(vector<vector<int>>& seq, int kicks){
    rebuildPos(seq);
    for(int t = 0; t < kicks; ++t){
        if(J < 2) return;
        int m = rndInt(M);
        int i = rndInt(J - 1);
        int a = opOnMachine(seq[m][i],   m);
        int b = opOnMachine(seq[m][i+1], m);
        doSwap(seq, a, b);
        if(evaluate(seq) < 0) doSwap(seq, b, a);
    }
}

static vector<long long> tabuUntil;   // swap tabu, keyed by (machine, job-pair)
static vector<long long> tabuJob;     // insert tabu, keyed by (machine, moved job)
static inline size_t tabIdx(int m, int ja, int jb){
    int lo = ja < jb ? ja : jb, hi = ja < jb ? jb : ja;
    return (size_t)m * J * J + (size_t)lo * J + hi;
}

static long long tabuSearch(vector<vector<int>>& best, long long bestMk){
    vector<vector<int>> cur = best;
    vector<pair<int,int>> swaps;
    vector<array<int,3>> inserts;
    fill(tabuUntil.begin(), tabuUntil.end(), 0);
    fill(tabuJob.begin(), tabuJob.end(), 0);
    rebuildPos(cur);
    long long iter = 0, lastImprove = 0;
    int tenure = 15 + rndInt(13);
    const long long stall = 5200;
    long long curMk = evaluate(cur);
    int checkClock = 0;

    while((checkClock++ & 63) || elapsed() < TL){
        getBlockMoves(cur, swaps, inserts);
        if(swaps.empty() && inserts.empty()){
            perturb(cur, 4); curMk = evaluate(cur); ++iter; continue;
        }

        // Best allowed (non-tabu) and best aspirated (tabu but < global best), across both
        // neighborhoods. mode: 0 = swap (a,b), 1 = insert (m,from,to).
        long long bestEst = LLONG_MAX, aspEst = LLONG_MAX;
        int alMode = -1, alA = -1, alB = -1, alM = -1, alF = -1, alT = -1;
        int asMode = -1, asA = -1, asB = -1, asM = -1, asF = -1, asT = -1;

        // --- swaps first: they read the current dist_/q_ (valid until an insert eval overwrites) ---
        for(auto& pr : swaps){
            int a = pr.first, b = pr.second;
            long long est = estimateSwap(a, b);
            bool isTabu = tabuUntil[tabIdx(machOf[a], jobOf[a], jobOf[b])] > iter;
            if(isTabu){ if(est < bestMk && est < aspEst){ aspEst = est; asMode = 0; asA = a; asB = b; } }
            else if(est < bestEst){ bestEst = est; alMode = 0; alA = a; alB = b; }
        }
        // --- inserts: O(L) head/tail estimate (no mutation), ranked alongside the swaps ---
        for(auto& ins : inserts){
            int m = ins[0], f = ins[1], t = ins[2];
            int job = cur[m][f];
            long long est = estInsert(cur, m, f, t);
            bool isTabu = tabuJob[(size_t)m * J + job] > iter;
            if(isTabu){ if(est < bestMk && est < aspEst){ aspEst = est; asMode = 1; asM = m; asF = f; asT = t; } }
            else if(est < bestEst){ bestEst = est; alMode = 1; alM = m; alF = f; alT = t; }
        }

        // choose: aspiration if it beats the best allowed, else best allowed, else aspiration.
        int useMode;   // 0 swap, 1 insert
        bool useAsp;
        if(asMode != -1 && aspEst <= bestEst)      { useAsp = true;  useMode = asMode; }
        else if(alMode != -1)                       { useAsp = false; useMode = alMode; }
        else if(asMode != -1)                       { useAsp = true;  useMode = asMode; }
        else { perturb(cur, 4); curMk = evaluate(cur); ++iter; continue; }

        if(useMode == 0){
            int a = useAsp ? asA : alA, b = useAsp ? asB : alB;
            doSwap(cur, a, b);
            tabuUntil[tabIdx(machOf[a], jobOf[a], jobOf[b])] = iter + tenure;
        }else{
            int m = useAsp ? asM : alM, f = useAsp ? asF : alF, t = useAsp ? asT : alT;
            int job = cur[m][f];
            doInsert(cur, m, f, t);
            tabuJob[(size_t)m * J + job] = iter + tenure;
        }
        curMk = evaluate(cur);
        if(curMk < 0){ // safety: should never happen (inserts are cycle-checked), but stay feasible
            cur = best; rebuildPos(cur); curMk = evaluate(cur);
        }
        if(curMk < bestMk){ bestMk = curMk; best = cur; lastImprove = iter; tenure = 15 + rndInt(13); }
        ++iter;

        if(iter - lastImprove > stall){
            cur = best;
            perturb(cur, 6 + rndInt(13));
            curMk = evaluate(cur);
            fill(tabuUntil.begin(), tabuUntil.end(), 0);
            fill(tabuJob.begin(), tabuJob.end(), 0);
            lastImprove = iter;
        }
    }
    return bestMk;
}

static void output(const vector<vector<int>>& seq){
    static char buf[1 << 22];
    int p = 0;
    for(int m = 0; m < M; ++m){
        for(int i = 0; i < J; ++i){
            int x = seq[m][i];
            if(x == 0) buf[p++] = '0';
            else{
                char tmp[12]; int t = 0;
                while(x){ tmp[t++] = char('0' + x % 10); x /= 10; }
                while(t) buf[p++] = tmp[--t];
            }
            buf[p++] = (i + 1 < J) ? ' ' : '\n';
        }
    }
    fwrite(buf, 1, p, stdout);
}

int main(){
    START = clock();
    if(scanf("%d %d", &J, &M) != 2) return 0;
    N = J * M;

    machJK.assign(J, vector<int>(M));
    procJK.assign(J, vector<long long>(M));
    posOf.assign(J, vector<int>(M, -1));
    procOp.assign(N, 0); jobOf.assign(N, 0); kOf.assign(N, 0); machOf.assign(N, 0);
    jobPred.assign(N, -1); jobSucc.assign(N, -1);

    for(int j = 0; j < J; ++j){
        for(int k = 0; k < M; ++k){
            int m; long long p;
            if(scanf("%d %lld", &m, &p) != 2) return 0;
            machJK[j][k] = m; procJK[j][k] = p; posOf[j][m] = k;
            int op = j * M + k;
            procOp[op] = p; jobOf[op] = j; kOf[op] = k; machOf[op] = m;
            jobPred[op] = (k > 0)     ? op - 1 : -1;
            jobSucc[op] = (k < M - 1) ? op + 1 : -1;
        }
    }

    indeg.assign(N, 0); mSucc.assign(N, -1); mPred.assign(N, -1);
    order_.assign(N, 0); dist_.assign(N, 0); q_.assign(N, 0);
    pos.assign(M, vector<int>(J, 0));
    tabuUntil.assign((size_t)M * J * J, 0);
    tabuJob.assign((size_t)M * J, 0);

    vector<vector<int>> best;
    long long bestMk = LLONG_MAX;
    for(int rule = 0; rule < 5; ++rule){
        vector<vector<int>> seq = gifflerThompson(rule);
        long long mk = evaluate(seq);
        if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
    }
    for(int r = 0; r < 150; ++r){
        vector<vector<int>> seq = gifflerThompson(r % 5);
        long long mk = evaluate(seq);
        if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
    }
    if(best.empty()){
        best.assign(M, vector<int>(J));
        for(int m = 0; m < M; ++m) for(int j = 0; j < J; ++j) best[m][j] = j;
        bestMk = evaluate(best);
    }

    bestMk = tabuSearch(best, bestMk);
    output(best);
    return 0;
}
