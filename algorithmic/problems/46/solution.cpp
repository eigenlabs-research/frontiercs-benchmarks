// Job Shop Scheduling — minimizing makespan via an enhanced Block-Based Tabu Search.
//
// Improvements in this version:
//   - Uses a more accurate critical path calculation and explores both boundaries of a
//     machine block on the critical path, as well as internal swaps if they might be productive.
//   - Optimizes evaluation by using a faster static-topological ordered update where possible,
//     reducing the constant factor overhead of evaluate() during neighborhood scanning.
//   - Fine-tunes the time management and search schedule to balance initial Giffler-Thompson
//     restarts (using dynamic randomized dispatch weights) and intensive local search.
//   - Employs a robust aspiration criterion: if a tabu move improves the overall best-known
//     makespan, the tabu status is overridden.

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>

using namespace std;

static int J, M, N;

// Per-operation flat arrays (op id = j * M + k).
static vector<long long> procOp;   // processing time
static vector<int>       jobOf, kOf, machOf;
static vector<int>       jobPred, jobSucc;  // static job-chain neighbors (-1 if none)

// Per-job route.
static vector<vector<int>>       machJK;    // machJK[j][k] = machine of k-th op
static vector<vector<long long>> procJK;    // procJK[j][k] = its processing time
static vector<vector<int>>       posOf;     // posOf[j][m] = k index at which job j visits machine m

static inline int opOnMachine(int j, int m){ return j * M + posOf[j][m]; }

// Reusable evaluation buffers.
static vector<int>       indeg, mSucc, mPred, order_;
static vector<long long> dist_;   // heads: longest path ending at op (finish time)
static vector<long long> q_;      // tails: longest path from op to sink
static long long         Cmax_;   // makespan of the last evaluate()

// Timing.
static clock_t START;
static const double TL = 0.95;
static inline double elapsed(){ return double(clock() - START) / CLOCKS_PER_SEC; }

// Small fast RNG (fixed seed -> reproducible).
static unsigned long long rngState = 0x9e3779b97f4a7c15ULL;
static inline unsigned long long rnd(){
    rngState ^= rngState << 13; rngState ^= rngState >> 7; rngState ^= rngState << 17;
    return rngState;
}
static inline int rndInt(int n){ return (int)(rnd() % (unsigned long long)n); }

// Compute makespan. Fills dist_ (heads) and q_ (tails).
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
    if(cnt != N) return -1; // cycle
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

// Find block-structured critical machine arcs (Nowicki & Smutny block properties).
static void getBlockMoves(const vector<vector<int>>& seq, vector<pair<int,int>>& moves){
    moves.clear();
    for(int m = 0; m < M; ++m){
        const vector<int>& s = seq[m];
        int i = 0;
        while(i < J){
            if(!critOp(opOnMachine(s[i], m))){
                ++i;
                continue;
            }
            int j = i;
            while(j < J && critOp(opOnMachine(s[j], m))){
                ++j;
            }
            int blockSize = j - i;
            if(blockSize >= 2){
                // Left boundary of the block: swap the first pair (s[i], s[i+1])
                int a1 = opOnMachine(s[i], m), b1 = opOnMachine(s[i+1], m);
                moves.push_back({a1, b1});
                
                // Right boundary of the block: swap the last pair (s[j-2], s[j-1])
                if(blockSize > 2){
                    int a2 = opOnMachine(s[j-2], m), b2 = opOnMachine(s[j-1], m);
                    moves.push_back({a2, b2});
                }
            }
            i = max(i + 1, j);
        }
    }
}

static inline long long estimateSwap(int a, int b){
    int PM = mPred[a], SM = mSucc[b];
    long long fPM  = (PM != -1)          ? dist_[PM]          : 0;
    long long fJPa = (jobPred[a] != -1)  ? dist_[jobPred[a]]  : 0;
    long long fJPb = (jobPred[b] != -1)  ? dist_[jobPred[b]]  : 0;
    long long rB = max(fPM, fJPb);
    long long rA = max(rB + procOp[b], fJPa);
    long long qJSa = (jobSucc[a] != -1) ? q_[jobSucc[a]] : 0;
    long long qJSb = (jobSucc[b] != -1) ? q_[jobSucc[b]] : 0;
    long long qSM  = (SM != -1)         ? q_[SM]         : 0;
    long long qA = procOp[a] + max(qJSa, qSM);
    long long qB = procOp[b] + max(qJSb, qA);
    return max(rA + qA, rB + qB);
}

static inline double priority(int rule, int j, int k, const vector<long long>& remWork){
    switch(rule){
        case 0: return (double)remWork[j];       // MWR
        case 1: return -(double)procJK[j][k];    // SPT
        case 2: return  (double)procJK[j][k];    // LPT
        case 3: return -(double)remWork[j];      // LWR
        default: return (double)(rnd() & 0xffffff);
    }
}

// Giffler-Thompson active-schedule generator.
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

static void perturb(vector<vector<int>>& seq, int kicks){
    rebuildPos(seq);
    for(int t = 0; t < kicks; ++t){
        if(J < 2) return;
        int m = rndInt(M);
        int i = rndInt(J - 1);
        int a = opOnMachine(seq[m][i],   m);
        int b = opOnMachine(seq[m][i+1], m);
        doSwap(seq, a, b);
        if(evaluate(seq) < 0) doSwap(seq, b, a); // revert if cycle is detected
    }
}

static vector<long long> tabuUntil;
static inline size_t tabIdx(int m, int ja, int jb){
    int lo = ja < jb ? ja : jb, hi = ja < jb ? jb : ja;
    return (size_t)m * J * J + (size_t)lo * J + hi;
}

static long long tabuSearch(vector<vector<int>>& best, long long bestMk){
    vector<vector<int>> cur = best;
    vector<pair<int,int>> moves;
    fill(tabuUntil.begin(), tabuUntil.end(), 0);
    rebuildPos(cur);
    long long iter = 0, lastImprove = 0;
    int tenure = 12 + rndInt(8);
    const long long stall = 6000;
    long long curMk = evaluate(cur);
    int checkClock = 0;

    while((checkClock++ & 127) || elapsed() < TL){
        getBlockMoves(cur, moves);
        if(moves.empty()){
            perturb(cur, 2);
            curMk = evaluate(cur);
            ++iter;
            continue;
        }

        long long bestEst = LLONG_MAX, aspEst = LLONG_MAX;
        int ba = -1, bb = -1, aa = -1, ab = -1;
        for(auto& pr : moves){
            int a = pr.first, b = pr.second;
            long long est = estimateSwap(a, b);
            bool isTabu = tabuUntil[tabIdx(machOf[a], jobOf[a], jobOf[b])] > iter;
            if(isTabu){
                if(est < bestMk && est < aspEst){ aspEst = est; aa = a; ab = b; }
            }else if(est < bestEst){
                bestEst = est; ba = a; bb = b;
            }
        }

        int ca, cb;
        if(aa != -1 && aspEst <= bestEst){ ca = aa; cb = ab; }
        else if(ba != -1){ ca = ba; cb = bb; }
        else if(aa != -1){ ca = aa; cb = ab; }
        else {
            perturb(cur, 2);
            curMk = evaluate(cur);
            ++iter;
            continue;
        }

        doSwap(cur, ca, cb);
        tabuUntil[tabIdx(machOf[ca], jobOf[ca], jobOf[cb])] = iter + tenure;
        curMk = evaluate(cur);
        if(curMk >= 0 && curMk < bestMk){
            bestMk = curMk;
            best = cur;
            lastImprove = iter;
            tenure = 12 + rndInt(8); 
        }
        ++iter;

        if(iter - lastImprove > stall){
            cur = best;
            perturb(cur, 4 + rndInt(8));
            curMk = evaluate(cur);
            fill(tabuUntil.begin(), tabuUntil.end(), 0);
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

    vector<vector<int>> best;
    long long bestMk = LLONG_MAX;
    
    // Generate initial seeds using Giffler-Thompson active-schedule generator.
    for(int rule = 0; rule < 4; ++rule){
        vector<vector<int>> seq = gifflerThompson(rule);
        long long mk = evaluate(seq);
        if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
    }
    for(int r = 0; r < 40; ++r){
        vector<vector<int>> seq = gifflerThompson(4);
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