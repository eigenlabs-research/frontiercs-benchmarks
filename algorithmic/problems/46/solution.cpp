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
static vector<long long> pnode;        // pnode[u] = p_of[j][k] for u=j*M+k

// Precomputed job-chain CSR (fixed across all evaluations).
static vector<int> jobHead, jobTo, jobNxt, jobIndeg;
static int jobE;
// Machine adjacency (rebuilt each eval).
static vector<int> mHead, mTo, mNxt;
// Reusable work buffers.
static vector<int> indeg, qbuf;
static vector<long long> dist, tail;
static vector<char> crit;

// Compute makespan of seq. Returns -1 if cycle. Optionally fill critical-node flags.
// A node u is critical if dist[u] + tail[u] - pnode[u] == C (lies on a longest path).
static long long evalSeq(const vector<vector<int>>& seq, bool fillCrit = false){
    // indeg = constant job-edge indeg + variable machine-edge indeg
    for(int u=0;u<N;++u) indeg[u] = jobIndeg[u];
    fill(mHead.begin(), mHead.begin()+N, -1);
    int ec = 0;
    for(int m=0;m<M;++m){
        const auto& s = seq[m];
        for(int i=0;i+1<J;++i){
            int j1=s[i], j2=s[i+1];
            int u = j1*M + pos[j1][m];
            int v = j2*M + pos[j2][m];
            mTo[ec]=v; mNxt[ec]=mHead[u]; mHead[u]=ec; ec++;
            indeg[v]++;
        }
    }
    fill(dist.begin(), dist.begin()+N, 0);
    qbuf.clear();
    int qh=0, cnt=0;
    for(int u=0;u<N;++u) if(indeg[u]==0){ dist[u]=pnode[u]; qbuf.push_back(u); }
    while(qh < (int)qbuf.size()){
        int u = qbuf[qh++]; cnt++;
        long long du = dist[u];
        for(int e=jobHead[u]; e!=-1; e=jobNxt[e]){
            int v = jobTo[e];
            long long nd = du + pnode[v];
            if(nd > dist[v]) dist[v] = nd;
            if(--indeg[v]==0) qbuf.push_back(v);
        }
        for(int e=mHead[u]; e!=-1; e=mNxt[e]){
            int v = mTo[e];
            long long nd = du + pnode[v];
            if(nd > dist[v]) dist[v] = nd;
            if(--indeg[v]==0) qbuf.push_back(v);
        }
    }
    if(cnt != N) return -1; // cycle
    long long C = 0; int arg = 0;
    for(int u=0;u<N;++u) if(dist[u] > C){ C = dist[u]; arg = u; }
    if(fillCrit){
        fill(tail.begin(), tail.begin()+N, 0);
        for(int idx=(int)qbuf.size()-1; idx>=0; --idx){
            int u = qbuf[idx];
            long long mx = 0;
            for(int e=jobHead[u]; e!=-1; e=jobNxt[e]){ int v=jobTo[e]; if(tail[v]>mx) mx=tail[v]; }
            for(int e=mHead[u]; e!=-1; e=mNxt[e]){ int v=mTo[e]; if(tail[v]>mx) mx=tail[v]; }
            tail[u] = pnode[u] + mx;
        }
        fill(crit.begin(), crit.begin()+N, 0);
        for(int u=0;u<N;++u) if(dist[u] + tail[u] - pnode[u] == C) crit[u] = 1;
    }
    return C;
}

struct Move { int m, i; };

// Full NS neighborhood: every adjacent pair (i,i+1) on each machine where BOTH
// operations are critical (lie on a longest path).
static vector<Move> genMoves(const vector<vector<int>>& seq, const vector<vector<int>>& pm){
    vector<Move> moves;
    for(int m=0;m<M;++m){
        const auto& s = seq[m];
        for(int i=0;i+1<J;++i){
            int a = s[i], b = s[i+1];
            int ua = a*M + pos[a][m];
            int ub = b*M + pos[b][m];
            if(crit[ua] && crit[ub]) moves.push_back({m, i});
        }
    }
    return moves;
}

// Taillard-style O(1) lower-bound estimate of the makespan after swapping the
// adjacent pair (a=seq[m][i], b=seq[m][i+1]) on machine m. Uses the current
// heads (dist[] = completion times) and tails (tail[] = start-to-sink length).
// Machine order becomes ... PM, b, a, SM ...  Returns the longest path length
// through the two swapped nodes under the new orientation (a true lower bound
// on the new makespan; exact for the local path). Requires dist/tail current.
static inline long long estimateMove(const vector<vector<int>>& seq, int m, int i){
    const auto& s = seq[m];
    int a = s[i], b = s[i+1];
    int ka = pos[a][m], kb = pos[b][m];
    int ua = a*M + ka, ub = b*M + kb;
    long long pa = pnode[ua], pb = pnode[ub];
    // completions of job predecessors (0 if none)
    long long c_jpa = (ka>0) ? dist[ua-1] : 0;
    long long c_jpb = (kb>0) ? dist[ub-1] : 0;
    // completion of machine predecessor of the pair (op at position i-1)
    long long c_pm = 0;
    if(i>0){ int pj = s[i-1]; c_pm = dist[pj*M + pos[pj][m]]; }
    // tails of job successors (0 if none)
    long long q_jsa = (ka<M-1) ? tail[ua+1] : 0;
    long long q_jsb = (kb<M-1) ? tail[ub+1] : 0;
    // tail of machine successor of the pair (op at position i+2)
    long long q_sm = 0;
    if(i+2<J){ int sj = s[i+2]; q_sm = tail[sj*M + pos[sj][m]]; }
    // new order: PM -> b -> a -> SM
    long long rb = (c_jpb > c_pm ? c_jpb : c_pm);          // start of b
    long long cb = rb + pb;                                 // completion of b
    long long ra = (c_jpa > cb ? c_jpa : cb);               // start of a
    long long qa = pa + (q_jsa > q_sm ? q_jsa : q_sm);      // tail of a
    long long qb = pb + (q_jsb > qa ? q_jsb : qa);          // tail of b
    long long la = ra + qa;
    long long lb = rb + qb;
    return la > lb ? la : lb;
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

int main(){
    auto T0 = chrono::steady_clock::now();
    const auto budget = chrono::milliseconds(770);

    if(scanf("%d %d", &J, &M) != 2) return 0;
    N = J*M;
    m_of.assign(J, vector<int>(M));
    p_of.assign(J, vector<long long>(M));
    pos.assign(J, vector<int>(M));
    for(int j=0;j<J;++j)
        for(int k=0;k<M;++k)
            if(scanf("%d %lld", &m_of[j][k], &p_of[j][k]) != 2) return 0;
    for(int j=0;j<J;++j) for(int k=0;k<M;++k) pos[j][m_of[j][k]] = k;
    pnode.assign(N, 0);
    for(int j=0;j<J;++j) for(int k=0;k<M;++k) pnode[j*M+k] = p_of[j][k];

    // Precompute job-chain CSR (fixed).
    jobE = J*(M-1);
    jobHead.assign(N, -1);
    jobTo.resize(jobE); jobNxt.resize(jobE);
    jobIndeg.assign(N, 0);
    {
        int ec = 0;
        for(int j=0;j<J;++j){
            int b = j*M;
            for(int k=0;k+1<M;++k){
                int u = b+k, v = b+k+1;
                jobTo[ec] = v; jobNxt[ec] = jobHead[u]; jobHead[u] = ec; ec++;
                jobIndeg[v]++;
            }
        }
    }
    mHead.assign(N, -1);
    mTo.resize(M*(J-1)); mNxt.resize(M*(J-1));
    indeg.resize(N); dist.resize(N); tail.resize(N); crit.resize(N);
    qbuf.reserve(N);

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

    vector<vector<int>> pm(M, vector<int>(J));
    auto rebuildPM = [&](const vector<vector<int>>& s){
        for(int m=0;m<M;++m) for(int i=0;i<J;++i) pm[m][s[m][i]] = i;
    };

    // Steepest descent via full NS neighborhood (runs until deadline).
    auto steepest = [&](chrono::steady_clock::time_point deadline){
        while(chrono::steady_clock::now() < deadline){
            rebuildPM(cur);
            long long c = evalSeq(cur, true);
            if(c < 0) break;
            auto moves = genMoves(cur, pm);
            long long bestNC = -1; int bestIdx = -1;
            for(int idx=0; idx<(int)moves.size(); ++idx){
                if((idx & 7)==0 && chrono::steady_clock::now() >= deadline) break;
                auto& mv = moves[idx];
                swap(cur[mv.m][mv.i], cur[mv.m][mv.i+1]);
                long long nc = evalSeq(cur, false);
                swap(cur[mv.m][mv.i], cur[mv.m][mv.i+1]);
                if(nc > 0 && (bestIdx==-1 || nc < bestNC)){ bestNC = nc; bestIdx = idx; }
            }
            if(bestIdx==-1 || bestNC >= c) break;
            auto& mv = moves[bestIdx];
            swap(cur[mv.m][mv.i], cur[mv.m][mv.i+1]);
            curC = bestNC;
            if(curC < bestC){ best = cur; bestC = curC; }
        }
    };

    // Tabu search over full NS neighborhood with ILS diversification (until deadline).
    auto tabu = [&](chrono::steady_clock::time_point deadline){
        vector<vector<int>> tabuUntil(M, vector<int>(J, 0));
        int iter = 0;
        int tenure = max(6, J);
        int sinceImprove = 0;
        const int bigStuck = 700;
        while(chrono::steady_clock::now() < deadline){
            iter++;
            rebuildPM(cur);
            long long c = evalSeq(cur, true);
            if(c < 0){ cur = best; curC = bestC; continue; }
            auto moves = genMoves(cur, pm);
            int nmv = (int)moves.size();
            // Rank all candidate moves by an O(1) Taillard-style lower-bound
            // estimate (uses the heads/tails just computed), then exact-eval
            // only the most promising few. Exact evalSeq still decides
            // acceptance and detects cycles, so feasibility is guaranteed.
            static vector<pair<long long,int>> cand;
            cand.clear(); cand.reserve(nmv);
            for(int idx=0; idx<nmv; ++idx)
                cand.push_back({estimateMove(cur, moves[idx].m, moves[idx].i), idx});
            int SK = nmv < 20 ? nmv : 20;         // shortlist to sort
            const int EVALCAP = 6;                 // exact evals per iteration
            partial_sort(cand.begin(), cand.begin()+SK, cand.end());
            long long bestNC = -1; int bestIdx = -1;
            int evaluated = 0;
            for(int t=0; t<SK && evaluated<EVALCAP; ++t){
                int idx = cand[t].second;
                auto& mv = moves[idx];
                bool isTabu = tabuUntil[mv.m][mv.i] > iter;
                bool aspEst = (cand[t].first < bestC);
                if(isTabu && !aspEst) continue;    // prune clearly-tabu moves cheaply
                swap(cur[mv.m][mv.i], cur[mv.m][mv.i+1]);
                long long nc = evalSeq(cur, false);
                swap(cur[mv.m][mv.i], cur[mv.m][mv.i+1]);
                if(nc < 0) continue;
                evaluated++;
                bool asp = (nc < bestC);
                if(isTabu && !asp) continue;
                if(bestIdx==-1 || nc < bestNC){ bestNC = nc; bestIdx = idx; }
            }
            if(bestIdx==-1){
                if(chrono::steady_clock::now() >= deadline) break;
                int m = (int)(rng() % M);
                int i = (int)(rng() % (J-1));
                swap(cur[m][i], cur[m][i+1]);
                long long nc = evalSeq(cur, false);
                if(nc < 0){ swap(cur[m][i], cur[m][i+1]); }
                else {
                    curC = nc;
                    if(curC < bestC){ best = cur; bestC = curC; sinceImprove = 0; }
                }
                continue;
            }
            auto& mv = moves[bestIdx];
            swap(cur[mv.m][mv.i], cur[mv.m][mv.i+1]);
            curC = bestNC;
            tabuUntil[mv.m][mv.i] = iter + tenure;
            if(curC < bestC){ best = cur; bestC = curC; sinceImprove = 0; }
            else {
                sinceImprove++;
                if(sinceImprove > bigStuck){
                    cur = best; curC = bestC;
                    int kicks = 3 + (int)(rng() % 4);
                    for(int r=0;r<kicks;++r){
                        int m = (int)(rng() % M);
                        int i = (int)(rng() % (J-1));
                        swap(cur[m][i], cur[m][i+1]);
                    }
                    long long nc = evalSeq(cur, false);
                    if(nc < 0){ cur = best; curC = bestC; }
                    else curC = nc;
                    sinceImprove = 0;
                    for(int m=0;m<M;++m) for(int i=0;i<J;++i) tabuUntil[m][i] = 0;
                }
            }
        }
    };

    // First pass: steepest descent + tabu from the best dispatch seed, using
    // the full budget. (The tabu search already includes ILS-style perturbation
    // via bigStuck, so external restarts are not beneficial here.)
    auto T_end = T0 + budget;
    steepest(T_end);
    tabu(T_end);

    // Fast buffered output (single fwrite) then _exit to skip static-vector
    // teardown, which otherwise adds ~100-200ms of wall time on some systems.
    {
        static char buf[1<<16];
        int len = 0;
        for(int m=0;m<M;++m){
            for(int j=0;j<J;++j){
                int x = best[m][j];
                if(x == 0){ buf[len++] = '0'; }
                else { char tmp[12]; int t = 0; while(x > 0){ tmp[t++] = char('0' + x%10); x /= 10; } while(t > 0) buf[len++] = tmp[--t]; }
                buf[len++] = (j+1<J ? ' ' : '\n');
            }
        }
        fwrite(buf, 1, len, stdout);
        fflush(stdout);
    }
    _exit(0);
}
