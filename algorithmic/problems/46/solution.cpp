// Job Shop Scheduling — FrontierCS problem 46 — solver v2
//
// score = clamp(1 - P/B, 0, 1), B a per-case constant => maximizing score is
// EXACTLY minimizing makespan P. So: shortest schedule in 1 second, nothing else.
//
// v1 -> v2, the one change that matters: neighbor evaluation went from a full
// O(N) topological decode to Taillard's O(1) head/tail estimate. v1 burned one
// full decode per neighbor (~5-30 per iteration) and managed only ~5-20k
// iterations. v2 ranks moves by the O(1) estimate and fully recomputes only the
// winner, buying ~20-30x more search in the same second.
//
// Design:
//   heads r[v]  = earliest start of op v      (longest path source -> v)
//   tails q[v]  = longest path from end of v -> sink
//   makespan    = max_v ( r[v] + p[v] + q[v] )
//   critical op = r[v] + p[v] + q[v] == makespan
//   N5 (Nowicki-Smutnicki): swap the first / last adjacent pair of each critical
//   block. Provably cannot create a cycle, so every move stays feasible.
//   i-TSAB-style elite pool: on stagnation, backtrack to an elite and perturb.
//
// Self-timed: stops on the wall clock, never on an iteration count, so it ports
// from our ARM dev box to their (probably x86) judge without a TLE or leaving
// budget unspent.

#include <bits/stdc++.h>
using namespace std;
using Clock = chrono::steady_clock;

static Clock::time_point T0;
static inline double el() { return chrono::duration<double>(Clock::now() - T0).count(); }
static double TIME_LIMIT = 0.85;   // 1.0s judge limit
static int P_TENURE_A = 6, P_TENURE_B = 5, P_STALL_A = 600, P_STALL_B = 8, P_KICK = 3, P_ELITE = 4;

static int J, M, N;
static vector<int> MACH, PT, KOF;          // flat [j*M + k]
static inline int OP(int j, int k) { return j * M + k; }
static inline int pt(int v) { return PT[v]; }

// schedule: seq[m][i] = job at position i on machine m ; pos[m][j] = its index
static vector<vector<int>> seqM, posM;

// heads / tails / order
static vector<int> R, Q, ORD, indeg_, deg_;

static inline int jobPred(int v) { return (v % M) ? v - 1 : -1; }
static inline int jobSucc(int v) { return ((v % M) + 1 < M) ? v + 1 : -1; }
static inline int macPred(const vector<vector<int>>& sq, int v) {
    int j = v / M, m = MACH[v], i = posM[m][j];
    if (i == 0) return -1;
    int j2 = sq[m][i - 1];
    return OP(j2, KOF[j2 * M + m]);
}
static inline int macSucc(const vector<vector<int>>& sq, int v) {
    int j = v / M, m = MACH[v], i = posM[m][j];
    if (i + 1 >= J) return -1;
    int j2 = sq[m][i + 1];
    return OP(j2, KOF[j2 * M + m]);
}

static void rebuildPos(const vector<vector<int>>& sq) {
    for (int m = 0; m < M; m++)
        for (int i = 0; i < J; i++) posM[m][sq[m][i]] = i;
}

// full decode: heads, tails, topological order. returns makespan, -1 if cyclic.
static int decode(const vector<vector<int>>& sq) {
    rebuildPos(sq);
    fill(indeg_.begin(), indeg_.end(), 0);
    for (int v = 0; v < N; v++) {
        if (jobPred(v) >= 0) indeg_[v]++;
        int j = v / M, m = MACH[v];
        if (posM[m][j] > 0) indeg_[v]++;
    }
    fill(R.begin(), R.end(), 0);
    ORD.clear();
    static vector<int> stk; stk.clear();
    deg_ = indeg_;
    for (int v = 0; v < N; v++) if (!deg_[v]) stk.push_back(v);
    while (!stk.empty()) {
        int v = stk.back(); stk.pop_back();
        ORD.push_back(v);
        int f = R[v] + pt(v);
        int a = jobSucc(v), b = macSucc(sq, v);
        if (a >= 0) { R[a] = max(R[a], f); if (!--deg_[a]) stk.push_back(a); }
        if (b >= 0) { R[b] = max(R[b], f); if (!--deg_[b]) stk.push_back(b); }
    }
    if ((int)ORD.size() != N) return -1;           // cycle
    fill(Q.begin(), Q.end(), 0);
    int mk = 0;
    for (int i = N - 1; i >= 0; i--) {
        int v = ORD[i];
        int t = 0;
        int a = jobSucc(v), b = macSucc(sq, v);
        if (a >= 0) t = max(t, Q[a] + pt(a));
        if (b >= 0) t = max(t, Q[b] + pt(b));
        Q[v] = t;
        mk = max(mk, R[v] + pt(v) + Q[v]);
    }
    return mk;
}

struct Move { int m, i, j; };   // j<0: swap i,i+1 (N5). j>=0: move op at i to after j (N7 insertion)

// Taillard O(1) estimate of the makespan after swapping u=sq[m][i], v=sq[m][i+1].
// Lower bound on the true new makespan; exact when the critical path stays local.
static inline int estimate(const vector<vector<int>>& sq, int m, int i, int curMk) {
    int ju = sq[m][i], jv = sq[m][i + 1];
    int u = OP(ju, KOF[ju * M + m]);
    int v = OP(jv, KOF[jv * M + m]);

    int mp = (i > 0) ? OP(sq[m][i-1], KOF[sq[m][i-1] * M + m]) : -1;
    int ms = (i + 2 < J) ? OP(sq[m][i+2], KOF[sq[m][i+2] * M + m]) : -1;

    int jpu = jobPred(u), jsu = jobSucc(u);
    int jpv = jobPred(v), jsv = jobSucc(v);

    // after swap: v then u
    int rv = 0;
    if (jpv >= 0) rv = max(rv, R[jpv] + pt(jpv));
    if (mp  >= 0) rv = max(rv, R[mp]  + pt(mp));
    int ru = 0;
    if (jpu >= 0) ru = max(ru, R[jpu] + pt(jpu));
    ru = max(ru, rv + pt(v));

    int qu = 0;
    if (jsu >= 0) qu = max(qu, Q[jsu] + pt(jsu));
    if (ms  >= 0) qu = max(qu, Q[ms]  + pt(ms));
    int qv = 0;
    if (jsv >= 0) qv = max(qv, Q[jsv] + pt(jsv));
    qv = max(qv, qu + pt(u));

    return max(ru + pt(u) + qu, rv + pt(v) + qv);
}

// Estimate makespan after moving the op at position a to just after position b
// (a < b) on machine m, or to just before position b (b < a). O(block length):
// we recompute heads forward and tails backward only across the shifted window,
// reusing R/Q outside it. Cheap enough to afford a much richer neighborhood.
static int estimateIns(const vector<vector<int>>& sq, int m, int a, int b) {
    int lo = min(a, b), hi = max(a, b);
    static vector<int> ops; ops.clear();
    for (int i = lo; i <= hi; i++) { int j = sq[m][i]; ops.push_back(OP(j, KOF[j*M+m])); }
    // new order of the window after the move
    static vector<int> nw; nw.clear();
    if (a < b) { for (int i = 1; i < (int)ops.size(); i++) nw.push_back(ops[i]); nw.push_back(ops[0]); }
    else       { nw.push_back(ops.back()); for (int i = 0; i + 1 < (int)ops.size(); i++) nw.push_back(ops[i]); }

    int prevOp = (lo > 0) ? OP(sq[m][lo-1], KOF[sq[m][lo-1]*M+m]) : -1;
    int nextOp = (hi + 1 < J) ? OP(sq[m][hi+1], KOF[sq[m][hi+1]*M+m]) : -1;

    static vector<int> rr, qq;
    rr.assign(nw.size(), 0); qq.assign(nw.size(), 0);
    for (int i = 0; i < (int)nw.size(); i++) {
        int v = nw[i], r = 0;
        int jp = jobPred(v);
        if (jp >= 0) r = max(r, R[jp] + pt(jp));
        if (i == 0) { if (prevOp >= 0) r = max(r, R[prevOp] + pt(prevOp)); }
        else        r = max(r, rr[i-1] + pt(nw[i-1]));
        rr[i] = r;
    }
    for (int i = (int)nw.size() - 1; i >= 0; i--) {
        int v = nw[i], q = 0;
        int js = jobSucc(v);
        if (js >= 0) q = max(q, Q[js] + pt(js));
        if (i + 1 == (int)nw.size()) { if (nextOp >= 0) q = max(q, Q[nextOp] + pt(nextOp)); }
        else q = max(q, qq[i+1] + pt(nw[i+1]));
        qq[i] = q;
    }
    int est = 0;
    for (int i = 0; i < (int)nw.size(); i++) est = max(est, rr[i] + pt(nw[i]) + qq[i]);
    return est;
}

// critical blocks -> N5 swaps + N7 insertion moves
static int USE_N7 = 1;
static void genMoves(const vector<vector<int>>& sq, int mk, vector<Move>& out) {
    out.clear();
    for (int m = 0; m < M; m++) {
        int i = 0;
        while (i < J) {
            int j = sq[m][i];
            int v = OP(j, KOF[j * M + m]);
            if (R[v] + pt(v) + Q[v] != mk) { i++; continue; }
            int b = i;
            while (b + 1 < J) {
                int j2 = sq[m][b + 1];
                int v2 = OP(j2, KOF[j2 * M + m]);
                if (R[v2] + pt(v2) + Q[v2] != mk) break;
                b++;
            }
            if (b > i) {
                out.push_back({m, i, -1});                       // N5: first pair
                if (b - 1 != i) out.push_back({m, b - 1, -1});   // N5: last pair
                if (USE_N7 && b - i >= 2) {
                    // N7: move each interior/edge op to the front or back of the block
                    for (int t = i + 1; t <= b; t++) out.push_back({m, t, i});   // -> front
                    for (int t = i; t < b; t++)      out.push_back({m, t, b});   // -> back
                }
            }
            i = b + 1;
        }
    }
}

// apply a move to a sequence
static inline void applyMove(vector<vector<int>>& sq, const Move& mv) {
    if (mv.j < 0) { swap(sq[mv.m][mv.i], sq[mv.m][mv.i + 1]); return; }
    int a = mv.i, b = mv.j;
    int job = sq[mv.m][a];
    if (a < b) { for (int i = a; i < b; i++) sq[mv.m][i] = sq[mv.m][i+1]; sq[mv.m][b] = job; }
    else       { for (int i = a; i > b; i--) sq[mv.m][i] = sq[mv.m][i-1]; sq[mv.m][b] = job; }
}
static inline void undoMove(vector<vector<int>>& sq, const Move& mv) {
    if (mv.j < 0) { swap(sq[mv.m][mv.i], sq[mv.m][mv.i + 1]); return; }
    Move inv{mv.m, mv.j, mv.i};
    applyMove(sq, inv);
}

// Giffler-Thompson active-schedule construction
static vector<vector<int>> gt(int rule, mt19937& rng) {
    vector<int> nk(J, 0), jr(J, 0), mr(M, 0), wl(J, 0);
    for (int j = 0; j < J; j++)
        for (int k = 0; k < M; k++) wl[j] += PT[OP(j, k)];
    vector<vector<int>> sq(M);
    for (int m = 0; m < M; m++) sq[m].reserve(J);

    for (int done = 0; done < N; done++) {
        int bf = INT_MAX, bm = -1;
        for (int j = 0; j < J; j++) {
            if (nk[j] >= M) continue;
            int v = OP(j, nk[j]), m = MACH[v];
            int f = max(jr[j], mr[m]) + pt(v);
            if (f < bf) { bf = f; bm = m; }
        }
        int pick = -1; long long bk = LLONG_MIN;
        for (int j = 0; j < J; j++) {
            if (nk[j] >= M) continue;
            int v = OP(j, nk[j]), m = MACH[v];
            if (m != bm) continue;
            int st = max(jr[j], mr[m]);
            if (st >= bf) continue;
            long long key;
            switch (rule) {
                case 0: key = wl[j]; break;                       // MWR
                case 1: key = -pt(v); break;                      // SPT
                case 2: key = pt(v); break;                       // LPT
                case 3: key = M - nk[j]; break;                   // MOR
                case 4: key = -(long long)st; break;              // ECT
                default: key = (long long)(rng() & 0xffff);       // random
            }
            if (key > bk) { bk = key; pick = j; }
        }
        if (pick < 0) for (int j = 0; j < J; j++)
            if (nk[j] < M && MACH[OP(j, nk[j])] == bm) { pick = j; break; }
        int j = pick, v = OP(j, nk[j]), m = MACH[v];
        int st = max(jr[j], mr[m]), f = st + pt(v);
        sq[m].push_back(j);
        jr[j] = f; mr[m] = f; wl[j] -= pt(v); nk[j]++;
    }
    return sq;
}

int main() {
    T0 = Clock::now();
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if (!(cin >> J >> M)) return 0;
    N = J * M;
    MACH.assign(N, 0); PT.assign(N, 0); KOF.assign(J * M, 0);
    for (int j = 0; j < J; j++)
        for (int k = 0; k < M; k++) {
            int m, p; cin >> m >> p;
            MACH[OP(j, k)] = m; PT[OP(j, k)] = p;
            KOF[j * M + m] = k;
        }
    R.assign(N, 0); Q.assign(N, 0); indeg_.assign(N, 0); deg_.assign(N, 0);
    ORD.reserve(N);
    seqM.assign(M, vector<int>(J)); posM.assign(M, vector<int>(J));

    { auto E=[&](const char*k,int&v){ if(const char*s=getenv(k)) v=atoi(s); };
      E("P_TENURE_A",P_TENURE_A); E("P_TENURE_B",P_TENURE_B);
      E("P_STALL_A",P_STALL_A);   E("P_STALL_B",P_STALL_B);
      E("P_KICK",P_KICK);         E("P_ELITE",P_ELITE);
      if(const char*s=getenv("P_TIME")) TIME_LIMIT=atof(s);
      if(const char*s=getenv("P_N7")) USE_N7=atoi(s); }
    unsigned seed = 0x9e3779b9;
    if (const char* s = getenv("P_SEED")) seed = (unsigned)atoi(s);
    mt19937 rng(seed);

    // ---- construction: all priority rules + randomized GT restarts ----------
    vector<vector<int>> best, cur;
    int bestMk = INT_MAX;
    for (int rule = 0; rule <= 5; rule++) {
        auto s = gt(rule, rng);
        int mk = decode(s);
        if (mk > 0 && mk < bestMk) { bestMk = mk; best = s; }
    }
    while (el() < 0.06) {
        auto s = gt(5, rng);
        int mk = decode(s);
        if (mk > 0 && mk < bestMk) { bestMk = mk; best = s; }
    }

    // ---- tabu search with O(1) move evaluation ------------------------------
    cur = best;
    int curMk = decode(cur);
    vector<vector<int>> tabu(M, vector<int>(J, 0));
    const int TENURE = P_TENURE_A + (J + M) / max(1,P_TENURE_B);

    // elite pool for i-TSAB style backtracking
    vector<pair<int, vector<vector<int>>>> elite;
    elite.emplace_back(bestMk, best);

    vector<Move> mv;
    long long iter = 0, evals = 0;
    int stall = 0;

    while (el() < TIME_LIMIT) {
        iter++;
        genMoves(cur, curMk, mv);
        if (mv.empty()) {                       // no critical block -> perturb
            for (int t = 0; t < 4; t++) {
                int m = rng() % M, i = rng() % (J - 1);
                swap(cur[m][i], cur[m][i + 1]);
            }
            int mk = decode(cur);
            if (mk < 0) { cur = best; curMk = decode(cur); } else curMk = mk;
            continue;
        }

        // rank all candidate moves by their cheap estimate
        static vector<pair<int,int>> ranked; ranked.clear();
        for (int t = 0; t < (int)mv.size(); t++) {
            int e = (mv[t].j < 0) ? estimate(cur, mv[t].m, mv[t].i, curMk)
                                  : estimateIns(cur, mv[t].m, mv[t].i, mv[t].j);
            evals++;
            bool isTabu = tabu[mv[t].m][mv[t].i] > iter;
            bool aspire = e < bestMk;           // aspiration: beats global best
            if (isTabu && !aspire) continue;
            ranked.push_back({e, t});
        }
        if (ranked.empty()) for (int t = 0; t < (int)mv.size(); t++) ranked.push_back({INT_MAX, t});
        sort(ranked.begin(), ranked.end());

        // apply the best estimated move; N7 insertions can be cyclic, so verify
        // and fall through to the next candidate if the decode rejects it.
        int mk = -1; Move mm{0,0,-1};
        for (int c = 0; c < (int)ranked.size() && c < 4; c++) {
            mm = mv[ranked[c].second];
            applyMove(cur, mm);
            mk = decode(cur);
            if (mk > 0) break;
            undoMove(cur, mm);
            mk = -1;
        }
        if (mk < 0) { curMk = decode(cur); stall++; continue; }
        curMk = mk;
        tabu[mm.m][mm.i] = iter + TENURE + (rng() % 5);

        if (curMk < bestMk) {
            bestMk = curMk; best = cur; stall = 0;
            elite.emplace_back(bestMk, best);
            if ((int)elite.size() > P_ELITE) {
                sort(elite.begin(), elite.end(),
                     [](auto& a, auto& b) { return a.first < b.first; });
                elite.resize(P_ELITE);
            }
        } else if (++stall > P_STALL_A + P_STALL_B * (J + M)) {
            // backtrack to a random elite and kick
            auto& e = elite[rng() % elite.size()];
            cur = e.second;
            int kicks = P_KICK + rng() % 6;
            for (int t = 0; t < kicks; t++) {
                int m = rng() % M, i = rng() % (J - 1);
                swap(cur[m][i], cur[m][i + 1]);
            }
            int k = decode(cur);
            if (k < 0) { cur = best; k = decode(cur); }
            curMk = k;
            for (auto& row : tabu) fill(row.begin(), row.end(), 0);
            stall = 0;
        }
    }

    if (getenv("JSSP_STATS"))
        fprintf(stderr, "iters=%lld evals=%lld mk=%d t=%.3f\n", iter, evals, bestMk, el());

    string out; out.reserve(N * 5);
    for (int m = 0; m < M; m++) {
        for (int i = 0; i < J; i++) { if (i) out += ' '; out += to_string(best[m][i]); }
        out += '\n';
    }
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
