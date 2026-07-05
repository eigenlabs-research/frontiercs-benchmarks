// Vertex Cover Challenge — kernelization + NuMVC/FastVC-style local search.
// Score = K*/K per case; 2s CPU, 512MB. N<=10^4, M<=5*10^5.
//
// Pipeline:
//  1. Parse, dedupe edges.
//  2. Reductions to fixpoint: deg-0, deg-1 (pendant), deg-2 triangle, dominance,
//     Nemhauser-Trotter LP (bipartite Hopcroft-Karp + Konig) crown-style reduction.
//  3. Initial cover on kernel: greedy + redundancy pruning.
//  4. NuMVC local search (edge weighting, CC, BMS removal) until ~1.82s.
//  5. Reconstruct + output.
#include <bits/stdc++.h>
using namespace std;

static chrono::steady_clock::time_point T0;
static double TL_MS = 1820.0;
static inline double elapsed_ms() {
    return chrono::duration<double, milli>(chrono::steady_clock::now() - T0).count();
}

// ---------- fast input ----------
static vector<char> inbuf;
static size_t inpos = 0;
static inline int readInt() {
    while (inpos < inbuf.size() && (inbuf[inpos] < '0' || inbuf[inpos] > '9') && inbuf[inpos] != '-') inpos++;
    bool neg = false;
    if (inpos < inbuf.size() && inbuf[inpos] == '-') { neg = true; inpos++; }
    long long v = 0;
    while (inpos < inbuf.size() && inbuf[inpos] >= '0' && inbuf[inpos] <= '9') { v = v * 10 + (inbuf[inpos] - '0'); inpos++; }
    return (int)(neg ? -v : v);
}

struct RNG {
    unsigned long long s;
    RNG(unsigned long long x) : s(x ? x : 88172645463325252ULL) {}
    inline unsigned long long nxt() { s ^= s << 7; s ^= s >> 9; return s; }
    inline int rint(int n) { return (int)(nxt() % (unsigned)n); }
};

static int envInt(const char* name, int def) {
    if (const char* e = getenv(name)) return atoi(e);
    return def;
}
int N, M;
vector<int> eu, ev;                 // deduped edges
vector<int> adjHead, adjNext, adjTo, adjEdge; // linked adjacency over live graph? -> use CSR after dedupe

// final decision per original vertex: -1 unknown, 0 excluded, 1 included
vector<int8_t> decided;

int main() {
    T0 = chrono::steady_clock::now();
    if (const char* e = getenv("VC_TL")) { double v = atof(e); if (v > 50 && v < 10000) TL_MS = v; }
    // read stdin
    {
        size_t cap = 1 << 22; inbuf.resize(cap); size_t len = 0;
        while (true) {
            if (len == cap) { cap <<= 1; inbuf.resize(cap); }
            size_t got = fread(inbuf.data() + len, 1, cap - len, stdin);
            if (got == 0) break;
            len += got;
        }
        inbuf.resize(len);
    }
    N = readInt(); M = readInt();
    if (N <= 0) return 0;
    // dedupe edges
    {
        vector<pair<int,int>> es; es.reserve(M);
        for (int i = 0; i < M; i++) {
            int a = readInt(), b = readInt();
            if (a == b) continue;
            if (a > b) swap(a, b);
            es.push_back({a, b});
        }
        sort(es.begin(), es.end());
        es.erase(unique(es.begin(), es.end()), es.end());
        M = (int)es.size();
        eu.resize(M); ev.resize(M);
        for (int i = 0; i < M; i++) { eu[i] = es[i].first - 1; ev[i] = es[i].second - 1; }
    }
    decided.assign(N, -1);

    // ---------- reduction phase on a mutable graph ----------
    // adjacency sets as sorted vectors + alive flags; small N so ok
    vector<vector<int>> g(N);
    for (int i = 0; i < M; i++) { g[eu[i]].push_back(ev[i]); g[ev[i]].push_back(eu[i]); }
    for (auto& v : g) { sort(v.begin(), v.end()); v.erase(unique(v.begin(), v.end()), v.end()); }
    vector<int> deg(N);
    for (int i = 0; i < N; i++) deg[i] = (int)g[i].size();
    vector<char> alive(N, 1);

    auto removeVertex = [&](int u) { // u decided; detach from neighbors
        alive[u] = 0;
        for (int w : g[u]) if (alive[w]) {
            auto& gw = g[w];
            auto it = lower_bound(gw.begin(), gw.end(), u);
            if (it != gw.end() && *it == u) { gw.erase(it); deg[w]--; }
        }
        g[u].clear(); deg[u] = 0;
    };

    // queue-based degree reductions + dominance (bounded)
    {
        deque<int> q;
        for (int i = 0; i < N; i++) q.push_back(i);
        vector<char> inq(N, 1);
        auto push = [&](int u) { if (alive[u] && !inq[u]) { inq[u] = 1; q.push_back(u); } };
        long long domBudget = 20000000; // cap dominance work
        while (!q.empty()) {
            int u = q.front(); q.pop_front(); inq[u] = 0;
            if (!alive[u]) continue;
            // refresh adjacency (drop dead)
            {
                auto& gu = g[u];
                gu.erase(remove_if(gu.begin(), gu.end(), [&](int w) { return !alive[w]; }), gu.end());
                deg[u] = (int)gu.size();
            }
            if (deg[u] == 0) { decided[u] = 0; alive[u] = 0; continue; }
            if (deg[u] == 1) { // pendant: take the neighbor
                int w = g[u][0];
                decided[u] = 0; removeVertex(u);
                decided[w] = 1;
                vector<int> nb = g[w];
                removeVertex(w);
                for (int x : nb) push(x);
                continue;
            }
            if (deg[u] == 2) {
                int a = g[u][0], b = g[u][1];
                // triangle: take a,b
                if (binary_search(g[a].begin(), g[a].end(), b)) {
                    decided[u] = 0; removeVertex(u);
                    decided[a] = 1; decided[b] = 1;
                    vector<int> nba = g[a]; removeVertex(a);
                    vector<int> nbb = g[b]; removeVertex(b);
                    for (int x : nba) push(x);
                    for (int x : nbb) push(x);
                    continue;
                }
            }
            // dominance: if exists neighbor w with N(u)\{w} ⊆ N[w] ... check u dominated
            if (domBudget > 0 && deg[u] <= 64) {
                for (int w : g[u]) {
                    if (!alive[w] || deg[w] < deg[u]) continue;
                    domBudget -= deg[u];
                    // check N(u) ⊆ N[w]
                    bool sub = true;
                    for (int x : g[u]) {
                        if (x == w) continue;
                        if (!alive[x]) continue;
                        if (!binary_search(g[w].begin(), g[w].end(), x)) { sub = false; break; }
                    }
                    if (sub) { // include w
                        decided[w] = 1;
                        vector<int> nb = g[w];
                        removeVertex(w);
                        for (int x : nb) push(x);
                        push(u);
                        break;
                    }
                    if (domBudget <= 0) break;
                }
            }
        }
    }

    // ---------- Nemhauser-Trotter (LP) on the remaining kernel ----------
    // bipartite double cover; Hopcroft-Karp; Konig; x=1 -> include, x=0 -> exclude
    {
        vector<int> ids; ids.reserve(N);
        for (int i = 0; i < N; i++) if (alive[i]) ids.push_back(i);
        int K = (int)ids.size();
        if (K > 0) {
            vector<int> pos(N, -1);
            for (int i = 0; i < K; i++) pos[ids[i]] = i;
            // bipartite: left i (0..K-1), right j
            vector<vector<int>> bg(K);
            for (int i = 0; i < K; i++) {
                int u = ids[i];
                for (int w : g[u]) if (alive[w]) bg[i].push_back(pos[w]);
            }
            vector<int> mL(K, -1), mR(K, -1), dist(K);
            auto bfs = [&]() {
                deque<int> dq;
                bool found = false;
                for (int i = 0; i < K; i++) { if (mL[i] < 0) { dist[i] = 0; dq.push_back(i); } else dist[i] = -1; }
                while (!dq.empty()) {
                    int i = dq.front(); dq.pop_front();
                    for (int j : bg[i]) {
                        int i2 = mR[j];
                        if (i2 < 0) found = true;
                        else if (dist[i2] < 0) { dist[i2] = dist[i] + 1; dq.push_back(i2); }
                    }
                }
                return found;
            };
            function<bool(int)> dfs = [&](int i) -> bool {
                for (int j : bg[i]) {
                    int i2 = mR[j];
                    if (i2 < 0 || (dist[i2] == dist[i] + 1 && dfs(i2))) {
                        mL[i] = j; mR[j] = i; return true;
                    }
                }
                dist[i] = -1;
                return false;
            };
            while (bfs()) {
                for (int i = 0; i < K; i++) if (mL[i] < 0) dfs(i);
                if (elapsed_ms() > TL_MS * 0.4) break; // safety
            }
            // Konig: Z = unmatched left ∪ reachable via alternating paths
            vector<char> visL(K, 0), visR(K, 0);
            deque<int> dq;
            for (int i = 0; i < K; i++) if (mL[i] < 0) { visL[i] = 1; dq.push_back(i); }
            while (!dq.empty()) {
                int i = dq.front(); dq.pop_front();
                for (int j : bg[i]) {
                    if (visR[j]) continue;
                    visR[j] = 1;
                    int i2 = mR[j];
                    if (i2 >= 0 && !visL[i2]) { visL[i2] = 1; dq.push_back(i2); }
                }
            }
            // bip cover: left not visL, right visR. x_v = (inCover(Lv)+inCover(Rv))/2
            for (int i = 0; i < K; i++) {
                int u = ids[i];
                int c = (visL[i] ? 0 : 1) + (visR[i] ? 1 : 0);
                if (c == 2) { // x=1: include
                    decided[u] = 1;
                    removeVertex(u);
                } else if (c == 0) { // x=0: exclude
                    decided[u] = 0;
                    removeVertex(u); // its edges must be covered by neighbors (all x=1) — they get included
                }
            }
        }
    }

    if (getenv("VC_DEBUG")) fprintf(stderr, "t_reduce_done=%.1f\n", elapsed_ms());
    // ---------- build kernel CSR for local search ----------
    vector<int> kid(N, -1), rid;
    for (int i = 0; i < N; i++) if (alive[i]) { kid[i] = (int)rid.size(); rid.push_back(i); }
    int KN = (int)rid.size();
    vector<pair<int,int>> ke;
    for (int i = 0; i < M; i++) {
        int a = eu[i], b = ev[i];
        if (alive[a] && alive[b]) ke.push_back({kid[a], kid[b]});
    }
    // sanity: any edge with a decided-0 endpoint must have the other endpoint decided-1
    for (int i = 0; i < M; i++) {
        int a = eu[i], b = ev[i];
        if (decided[a] == 0 && decided[b] == 0) { decided[a] = 1; } // failsafe (shouldn't happen)
        if (decided[a] == 0 && decided[b] == -1 && !alive[b]) decided[b] = 1;
        if (decided[b] == 0 && decided[a] == -1 && !alive[a]) decided[a] = 1;
    }
    int KM = (int)ke.size();

    // local search on kernel
    vector<char> best(KN, 1);
    if (KN > 0 && KM > 0) {
        vector<int> head(KN + 1, 0), nxtTo(2 * KM), nxtEdge(2 * KM);
        {
            vector<int> cnt(KN, 0);
            for (auto& e : ke) { cnt[e.first]++; cnt[e.second]++; }
            for (int i = 0; i < KN; i++) head[i + 1] = head[i] + cnt[i];
            vector<int> fill = head;
            for (int i = 0; i < KM; i++) {
                auto [a, b] = ke[i];
                nxtTo[fill[a]] = b; nxtEdge[fill[a]] = i; fill[a]++;
                nxtTo[fill[b]] = a; nxtEdge[fill[b]] = i; fill[b]++;
            }
        }
        vector<int> w(KM, 1);
        vector<char> inC(KN, 0);
        vector<int> cl; cl.reserve(KN);      // cover vertex list
        vector<int> clPos(KN, -1);
        // bucket queue over loss values for exact min-loss removal
        int maxDeg = 0;
        { for (int i = 0; i < KN; i++) maxDeg = max(maxDeg, head[i + 1] - head[i]); }
        vector<vector<int>> bkt(maxDeg + 2);
        vector<int> bktPos(KN, -1);
        vector<int> lossOf(KN, 0);
        int minPtr = 0;
        auto bktInsert = [&](int v2, int L) {
            lossOf[v2] = L; bktPos[v2] = (int)bkt[L].size(); bkt[L].push_back(v2);
            if (L < minPtr) minPtr = L;
        };
        auto bktErase = [&](int v2) {
            int L = lossOf[v2], p2 = bktPos[v2];
            int last = bkt[L].back(); bkt[L][p2] = last; bktPos[last] = p2; bkt[L].pop_back();
            bktPos[v2] = -1;
        };
        auto bktMove = [&](int v2, int newL) {
            if (bktPos[v2] < 0) return;
            bktErase(v2); bktInsert(v2, newL);
        };
        auto bktPopMin = [&](RNG& r2) -> int {
            while (minPtr < (int)bkt.size() && bkt[minPtr].empty()) minPtr++;
            if (minPtr >= (int)bkt.size()) return -1;
            auto& B = bkt[minPtr];
            int idx = (int)(r2.nxt() % B.size());
            int v2 = B[idx];
            return v2;
        };
        vector<long long> score(KN, 0); // for v in C: -(loss); for v not in C: +gain
        vector<char> conf(KN, 1);
        vector<int> uncov; uncov.reserve(KM);
        vector<int> uncovPos(KM, -1);
        long long curW = 0;

        RNG rng(0x1234abcd5678ULL ^ (unsigned long long)KN * 1000003ULL ^ (unsigned long long)KM);

        // initial cover: greedy max-degree
        {
            vector<int> dg(KN, 0);
            for (auto& e : ke) { dg[e.first]++; dg[e.second]++; }
            // process edges; take higher degree endpoint (like baseline but on kernel), then prune
            vector<int> order(KM); iota(order.begin(), order.end(), 0);
            for (int i = 0; i < KM; i++) {
                auto [a, b] = ke[i];
                if (inC[a] || inC[b]) continue;
                inC[dg[a] >= dg[b] ? a : b] = 1;
            }
            // prune redundant (loss 0): vertices whose every edge is co-covered
            vector<int> lossCnt(KN, 0);
            for (int i = 0; i < KM; i++) {
                auto [a, b] = ke[i];
                if (inC[a] && !inC[b]) lossCnt[a]++;
                else if (inC[b] && !inC[a]) lossCnt[b]++;
            }
            vector<int> stk;
            for (int i = 0; i < KN; i++) if (inC[i] && lossCnt[i] == 0) stk.push_back(i);
            vector<char> instk(KN, 0);
            for (int v : stk) instk[v] = 1;
            while (!stk.empty()) {
                int u = stk.back(); stk.pop_back(); instk[u] = 0;
                if (!inC[u] || lossCnt[u] != 0) continue;
                inC[u] = 0;
                for (int it = head[u]; it < head[u + 1]; it++) {
                    int v2 = nxtTo[it];
                    if (inC[v2]) {
                        lossCnt[v2]++;
                    }
                }
                // removing u makes its neighbors' co-coverage change: for each edge (u,v2), v2 now solely covers it
                for (int it = head[u]; it < head[u + 1]; it++) {
                    int v2 = nxtTo[it];
                    if (inC[v2]) lossCnt[v2] += 0; // already handled above conceptually
                }
                // recompute neighbors' loss from scratch is safer:
            }
            // (simple approach: recompute loss and repeat prune until fixpoint)
            bool changed = true;
            while (changed) {
                changed = false;
                fill(lossCnt.begin(), lossCnt.end(), 0);
                for (int i = 0; i < KM; i++) {
                    auto [a, b] = ke[i];
                    if (inC[a] && !inC[b]) lossCnt[a]++;
                    else if (inC[b] && !inC[a]) lossCnt[b]++;
                }
                for (int i = 0; i < KN; i++) if (inC[i] && lossCnt[i] == 0) {
                    // check no uncovered edge appears: all edges of i must have other endpoint in C
                    bool ok = true;
                    for (int it = head[i]; it < head[i + 1]; it++) if (!inC[nxtTo[it]]) { ok = false; break; }
                    if (ok) { inC[i] = 0; changed = true; }
                }
            }
        }

        // init score/uncov from inC
        auto rebuild = [&]() {
            uncov.clear();
            fill(uncovPos.begin(), uncovPos.end(), -1);
            fill(score.begin(), score.end(), 0LL);
            curW = 0;
            for (int i = 0; i < KM; i++) {
                auto [a, b] = ke[i];
                if (!inC[a] && !inC[b]) { uncovPos[i] = (int)uncov.size(); uncov.push_back(i); curW += w[i]; score[a] += w[i]; score[b] += w[i]; }
                else if (inC[a] && !inC[b]) score[a] -= w[i];
                else if (inC[b] && !inC[a]) score[b] -= w[i];
            }
        };
        rebuild();
        cl.clear(); fill(clPos.begin(), clPos.end(), -1);
        for (int i = 0; i < KN; i++) if (inC[i]) { clPos[i] = (int)cl.size(); cl.push_back(i); }
        for (int i = 0; i < KN; i++) if (inC[i]) bktInsert(i, (int)(-score[i]));

        auto coverSize = [&]() { long long c = 0; for (int i = 0; i < KN; i++) c += inC[i]; return c; };
        long long curSize = coverSize();
        long long bestSize = LLONG_MAX;
        auto saveBest = [&]() {
            if (curSize < bestSize) { bestSize = curSize; for (int i = 0; i < KN; i++) best[i] = inC[i]; }
        };
        if (uncov.empty()) saveBest();

        auto addV = [&](int u) {
            inC[u] = 1; curSize++;
            clPos[u] = (int)cl.size(); cl.push_back(u);
            score[u] = -score[u];
            for (int it = head[u]; it < head[u + 1]; it++) {
                int v2 = nxtTo[it], e = nxtEdge[it];
                if (!inC[v2]) {
                    int p = uncovPos[e];
                    if (p >= 0) {
                        int last = uncov.back(); uncov[p] = last; uncovPos[last] = p; uncov.pop_back(); uncovPos[e] = -1;
                    }
                    score[v2] -= w[e];
                    conf[v2] = 1;
                } else {
                    score[v2] += w[e];
                    bktMove(v2, (int)(-score[v2]));
                }
            }
            bktInsert(u, (int)(-score[u]));
        };
        auto remV = [&](int u) {
            inC[u] = 0; curSize--;
            { int p2 = clPos[u]; int last = cl.back(); cl[p2] = last; clPos[last] = p2; cl.pop_back(); clPos[u] = -1; }
            bktErase(u);
            score[u] = -score[u];
            conf[u] = 0;
            for (int it = head[u]; it < head[u + 1]; it++) {
                int v2 = nxtTo[it], e = nxtEdge[it];
                if (!inC[v2]) {
                    uncovPos[e] = (int)uncov.size(); uncov.push_back(e);
                    score[v2] += w[e];
                    conf[v2] = 1;
                } else {
                    score[v2] -= w[e];
                    bktMove(v2, (int)(-score[v2]));
                }
            }
        };

        if (getenv("VC_DEBUG")) fprintf(stderr, "t_init_done=%.1f KN=%d KM=%d curSize=%lld\n", elapsed_ms(), KN, KM, curSize);
        // sorted adjacency (for adjacency tests in 2-improvements)
        vector<int> sadj(2 * KM); vector<int> shead = head;
        {
            for (int i = 0; i < KN; i++) {
                int a = head[i], b = head[i + 1];
                for (int it = a; it < b; it++) sadj[it] = nxtTo[it];
                sort(sadj.begin() + a, sadj.begin() + b);
            }
        }
        auto isAdj = [&](int a, int b) {
            int lo = head[a], hi = head[a + 1];
            return binary_search(sadj.begin() + lo, sadj.begin() + hi, b);
        };
        // ARW-style polish on a feasible cover (complement independent set view).
        // cov: char vector (1=in cover). Returns improved size; modifies cov in place.
        vector<int> tight(KN);
        vector<int> oneTightOf(KN);
        auto oneTightTo = [&](vector<char>& cov, int u) -> int {
            int cnt = 0, last = -1;
            for (int it = head[u]; it < head[u + 1]; it++) {
                int x = nxtTo[it];
                if (!cov[x]) { cnt++; last = x; if (cnt > 1) return -1; }
            }
            return cnt == 1 ? last : -2; // -2: zero non-cover neighbors (loss 0)
        };
        auto polish = [&](vector<char>& cov, double deadline, RNG& r2) -> long long {
            // I = complement of cov. compute tight for non-I (cover) vertices... tight defined for u in cover:
            // tight(u) = # neighbors NOT in cover (i.e., in I).
            long long covSz = 0;
            for (int i = 0; i < KN; i++) covSz += cov[i];
            bool anyImp = true;
            int rounds = 0;
            while (anyImp && elapsed_ms() < deadline && rounds < 200) {
                rounds++;
                anyImp = false;
                // build 1-tight groups: cover vertex u with exactly one I-neighbor v -> group[v]
                for (int i = 0; i < KN; i++) { tight[i] = 0; oneTightOf[i] = -1; }
                for (int u = 0; u < KN; u++) {
                    if (!cov[u]) continue;
                    int cnt = 0, last = -1;
                    for (int it = head[u]; it < head[u + 1]; it++) {
                        int x = nxtTo[it];
                        if (!cov[x]) { cnt++; last = x; if (cnt > 1) break; }
                    }
                    tight[u] = cnt;
                    if (cnt == 1) oneTightOf[u] = last;
                }
                // group by I-vertex
                // (KN small; use buckets)
                static vector<vector<int>> grp; grp.assign(KN, {});
                for (int u = 0; u < KN; u++) if (cov[u] && oneTightOf[u] >= 0) grp[oneTightOf[u]].push_back(u);
                for (int v2 = 0; v2 < KN; v2++) {
                    auto& G2 = grp[v2];
                    if ((int)G2.size() < 2) continue;
                    if (cov[v2]) continue; // must be I vertex
                    // find non-adjacent pair (cap work); validate freshness at apply time
                    int found_u = -1, found_w = -1;
                    int cap = min((int)G2.size(), 24);
                    for (int a2 = 0; a2 < cap && found_u < 0; a2++) {
                        for (int b2 = a2 + 1; b2 < cap; b2++) {
                            int u = G2[a2], w2 = G2[b2];
                            if (!cov[u] || !cov[w2]) continue; // stale membership
                            if (oneTightTo(cov, u) != v2 || oneTightTo(cov, w2) != v2) continue; // FRESH check
                            if (!isAdj(u, w2)) { found_u = u; found_w = w2; break; }
                        }
                    }
                    if (found_u >= 0) {
                        // apply: cover loses u,w gains v2  => K-1
                        cov[found_u] = 0; cov[found_w] = 0; cov[v2] = 1;
                        covSz--;
                        anyImp = true;
                        // local tight updates are complex; lazy: mark round dirty (rebuild next round)
                        // invalidate stale group entries by checking cov flags (done above)
                        // also neighbors' oneTight may change; conservative: continue scanning with cov checks
                    }
                }
                if (elapsed_ms() > deadline) break;
                if (!anyImp) {
                    // plateau (1,1) swaps to shuffle: pick random 1-tight u, swap with its I-neighbor
                    int swaps = 0;
                    for (int attempt = 0; attempt < KN && swaps < max(8, KN / 64); attempt++) {
                        int u = (int)(r2.nxt() % KN);
                        if (!cov[u]) continue;
                        int v2 = oneTightTo(cov, u); // FRESH: exactly one non-cover neighbor
                        if (v2 < 0) continue;
                        cov[u] = 0; cov[v2] = 1;
                        swaps++;
                    }
                    if (swaps > 0) { anyImp = true; } // loop again after shuffle
                    // note: rounds cap prevents infinite plateau cycling
                }
            }
            return covSz;
        };
        // main NuMVC loop
        long long lastImprove = 0;
        const int KICK = envInt("VC_KICK", 50);      // polish after KICK*1000 stale iters
        const int KICKSZ = envInt("VC_KICKSZ", 30);  // vertices removed per kick (per-mille of cover)
        const int BMS = envInt("VC_BMS", 64);
        const int WEIGHTED = envInt("VC_WEIGHTED", 0);
        const int EPS = envInt("VC_EPS", 0); // per-mille random removal
        const long long GAMMA = envInt("VC_GAMMA", 32);
        const int RHOx10 = envInt("VC_RHO", 3);
        long long iters = 0;
        long long feasCnt = 0;
        int checkMask = 1023;
        vector<int> cList; // maintained lazily for BMS sampling: sample random ids
        while (true) {
            if ((iters++ & checkMask) == 0 && elapsed_ms() > TL_MS) break;
            if (uncov.empty()) {
                feasCnt++;
                if (curSize < bestSize) lastImprove = iters;
                saveBest();
                // remove vertex with min loss (exact via buckets)
                if (cl.empty()) break;
                int bestV = bktPopMin(rng);
                if (bestV < 0) break;
                remV(bestV);
                continue;
            }
            // stale? -> ARW polish on best cover, then resume from polished state
            if (KICK > 0 && iters - lastImprove > (long long)KICK * 1000 && bestSize != LLONG_MAX) {
                lastImprove = iters;
                vector<char> cov(best.begin(), best.end());
                double slice = min(TL_MS, elapsed_ms() + 120.0);
                long long newSz = polish(cov, slice, rng);
                if (newSz < bestSize) {
                    bestSize = newSz;
                    for (int i2 = 0; i2 < KN; i2++) best[i2] = cov[i2];
                    // adopt polished cover as current state
                    for (int i2 = 0; i2 < KN; i2++) inC[i2] = cov[i2];
                    rebuild();
                    curSize = coverSize();
                    cl.clear(); fill(clPos.begin(), clPos.end(), -1);
                    for (int i2 = 0; i2 < KN; i2++) if (inC[i2]) { clPos[i2] = (int)cl.size(); cl.push_back(i2); }
                    for (auto& B : bkt) B.clear();
                    fill(bktPos.begin(), bktPos.end(), -1);
                    minPtr = 0;
                    for (int i2 = 0; i2 < KN; i2++) if (inC[i2]) bktInsert(i2, (int)(-score[i2]));
                }
                continue;
            }
            // remove step: exact min loss via buckets (epsilon random diversification)
            if (cl.empty()) break;
            int rv = -1;
            if (EPS > 0 && (int)(rng.nxt() % 1000) < EPS) {
                rv = cl[rng.rint((int)cl.size())];
            } else {
                rv = bktPopMin(rng);
                if (rv < 0) break;
            }
            remV(rv);
            // add step: random uncovered edge, pick better endpoint with conf
            int e = uncov[rng.rint((int)uncov.size())];
            int a = ke[e].first, b = ke[e].second;
            int pick;
            bool ca = conf[a], cb = conf[b];
            if (ca && !cb) pick = a;
            else if (cb && !ca) pick = b;
            else {
                if (score[a] > score[b]) pick = a;
                else if (score[b] > score[a]) pick = b;
                else pick = (rng.nxt() & 1) ? a : b;
            }
            addV(pick);
            if (WEIGHTED) {
                // weight update: bump all uncovered
                for (int ui = 0; ui < (int)uncov.size(); ui++) {
                    int e2 = uncov[ui];
                    w[e2]++;
                    score[ke[e2].first]++;
                    score[ke[e2].second]++;
                }
                curW += (long long)uncov.size();
                // forget weights occasionally
                if ((iters & 0x1FFF) == 0) {
                    long long avg = 0;
                    for (int i2 = 0; i2 < KM; i2++) avg += w[i2];
                    if (avg / max(1, KM) >= GAMMA) {
                        for (int i2 = 0; i2 < KM; i2++) w[i2] = (w[i2] * RHOx10) / 10 + 1;
                        rebuild();
                    }
                }
            }
        }
        if (bestSize == LLONG_MAX) { // never feasible? fall back to all-in-cover kernel
            for (int i = 0; i < KN; i++) best[i] = 1;
        }
        if (getenv("VC_DEBUG")) {
            long long wsum = 0; for (int i2 = 0; i2 < KM; i2++) wsum += w[i2];
            fprintf(stderr, "t_ls_done=%.1f iters=%lld feas=%lld bestSize=%lld meanW=%.1f uncov=%zu\n",
                    elapsed_ms(), iters, feasCnt, bestSize, (double)wsum / max(1, KM), uncov.size());
        }
    } else {
        best.assign(KN, 0);
    }

    // (debug moved into LS scope end)
    // ---------- output ----------
    string out; out.reserve(2 * N + 16);
    for (int i = 0; i < N; i++) {
        int val;
        if (decided[i] == 1) val = 1;
        else if (decided[i] == 0) val = 0;
        else val = best[kid[i]] ? 1 : 0;
        out += val ? "1\n" : "0\n";
    }
    fwrite(out.data(), 1, out.size(), stdout);
    if (getenv("VC_DEBUG")) fprintf(stderr, "t_output_done=%.1f\n", elapsed_ms());
    fflush(stdout);
    _Exit(0);
}
