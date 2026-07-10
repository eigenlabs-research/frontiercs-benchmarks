// Code by: Zain Ul Haq
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

// v19.9: v19.5 + PP_CAPDW default 3 (was 2); local small-band A/B winner
#include <bits/stdc++.h>
using namespace std;



static chrono::steady_clock::time_point T0;
static double TL_MS = 1890.0; // Hard global deadline (env POLYPACK_TL override)
static inline double elapsed_ms() {
    return chrono::duration<double, milli>(chrono::steady_clock::now() - T0).count();
}

struct T { int w, h; vector<pair<int,int>> c; vector<int> lo, hi; int r, f, minx, miny;
           unsigned short rmask[10]; vector<pair<signed char, signed char>> nbr; };
struct P { int id, k; vector<pair<int,int>> b; vector<T> t; int minW = 1e9, minH = 1e9, minA = 1e9; };
struct Pl { int idx, ti, x, y; };
struct R { long long A; int W, H; vector<Pl> pl; bool ok = false; int packW = 0; };
struct RNG {
    unsigned long long s;
    RNG(unsigned long long x) { s = x ? x : 1; }
    inline unsigned long long nxt() { s ^= s << 7; s ^= s >> 9; return s; }
    inline int rint(int n) { return (int)(nxt() % n); }
    inline bool coin() { return nxt() & 1; }
    inline double uni() { return (double)(nxt() >> 11) * (1.0 / 9007199254740992.0); }
};
static inline pair<int,int> rotp(pair<int,int> p, int r) {
    if (r == 0) return p;
    if (r == 1) return make_pair(-p.second, p.first);
    if (r == 2) return make_pair(-p.first, -p.second);
    return make_pair(p.second, -p.first);
}

static vector<char> inbuf;
static size_t inpos = 0;
static inline int readInt() {
    while (inpos < inbuf.size() && (inbuf[inpos] < '0' || inbuf[inpos] > '9') && inbuf[inpos] != '-') inpos++;
    bool neg = false;
    if (inpos < inbuf.size() && inbuf[inpos] == '-') { neg = true; inpos++; }
    int v = 0;
    while (inpos < inbuf.size() && inbuf[inpos] >= '0' && inbuf[inpos] <= '9') { v = v * 10 + (inbuf[inpos] - '0'); inpos++; }
    return neg ? -v : v;
}

int n;
long long S = 0;
vector<P> ps;
static int gFASTFIT = 1; // skyline x0-scan: skip dsum/dr for positions that can't win primary key

static R pack(int W, const vector<int>& o0, RNG& rng, bool randtie, int dynLIM0, bool adaptive,
              double adaptTLms, double deadlineMs, double panicMs = -1.0, int adaptMode = 0) {
    vector<int> h(W, -1);
    long long g = -1;
    vector<Pl> pl; pl.reserve(o0.size());
    vector<int> o = o0;
    int t = 0, nm = (int)o.size();
    int maxBound = max(1, nm / 4);
    int dynLIM = dynLIM0;
    auto tStart = chrono::steady_clock::now();
    auto batchStart = tStart;
    double emaBatch = -1.0;
    int stepCnt = 0;
    int checkCnt = 0;
    bool panic = false;
    while (t < nm) {
        if (!panic && panicMs > 0 && elapsed_ms() > panicMs) { panic = true; dynLIM = 1; adaptive = false; }
        int limCnt = max(1, dynLIM);
        int lim = min(nm, t + limCnt);
        long long bestg = LLONG_MAX; int bti = -1, bx = 0, by = 0, bl = INT_MAX;
        long long bds = LLONG_MAX, bdr = LLONG_MAX; int by0 = INT_MAX, bx0 = INT_MAX;
        int bestpos = t, bestid = -1;
        for (int pos = t; pos < lim; pos++) {
            if (!panic && deadlineMs > 0 && ((++checkCnt & 7) == 0) && elapsed_ms() > deadlineMs) return R{};
            int id = o[pos];
            auto& p = ps[id];
            long long bestg2 = LLONG_MAX; int bti2 = -1, bx2 = 0, by2 = 0, bl2 = INT_MAX;
            long long bds2 = LLONG_MAX, bdr2 = LLONG_MAX; int by02 = INT_MAX, bx02 = INT_MAX;
            for (int ti = 0; ti < (int)p.t.size(); ti++) {
                auto& tsh = p.t[ti];
                if (tsh.w > W) continue;
                int Rpos = W - tsh.w + 1;
                for (int x0 = 0; x0 < Rpos; x0++) {
                    int y0 = 0;
                    for (int j = 0; j < tsh.w; j++) {
                        if (tsh.lo[j] != INT_MAX) { int v = h[x0 + j] - tsh.lo[j] + 1; if (v > y0) y0 = v; }
                    }
                    int nhbuf[32];
                    int l = -1;
                    for (int j = 0; j < tsh.w; j++) {
                        int nh = h[x0 + j];
                        if (tsh.hi[j] != INT_MIN) { int cand = y0 + tsh.hi[j]; if (nh < cand) nh = cand; }
                        nhbuf[j] = nh;
                        if (nh > l) l = nh;
                    }
                    long long gg = g; if (gg < l) gg = l;
                    if (gFASTFIT && bti2 != -1 && gg > bestg2) continue;
                    long long dsum = 0;
                    for (int j = 0; j < tsh.w; j++) { int inc = nhbuf[j] - h[x0 + j]; if (inc > 0) dsum += inc; }
                    long long dr = 0;
                    if (x0 > 0) {
                        long long old = llabs((long long)h[x0] - h[x0 - 1]);
                        long long nw = llabs((long long)nhbuf[0] - h[x0 - 1]);
                        dr += nw - old;
                    }
                    for (int j = 0; j < tsh.w - 1; j++) {
                        long long old = llabs((long long)h[x0 + j + 1] - h[x0 + j]);
                        long long nw = llabs((long long)nhbuf[j + 1] - nhbuf[j]);
                        dr += nw - old;
                    }
                    if (x0 + tsh.w < W) {
                        long long old = llabs((long long)h[x0 + tsh.w] - h[x0 + tsh.w - 1]);
                        long long nw = llabs((long long)h[x0 + tsh.w] - nhbuf[tsh.w - 1]);
                        dr += nw - old;
                    }
                    bool take = false;
                    if (gg < bestg2) take = true;
                    else if (gg == bestg2) {
                        if (dsum < bds2) take = true;
                        else if (dsum == bds2) {
                            if (l < bl2) take = true;
                            else if (l == bl2) {
                                if (dr < bdr2) take = true;
                                else if (dr == bdr2) {
                                    if (y0 < by02) take = true;
                                    else if (y0 == by02) {
                                        if (x0 < bx02) take = true;
                                        else if (x0 == bx02 && randtie && rng.coin()) take = true;
                                    }
                                }
                            }
                        }
                    }
                    if (take) { bestg2 = gg; bti2 = ti; bx2 = x0; by2 = y0; bl2 = l; bds2 = dsum; bdr2 = dr; by02 = y0; bx02 = x0; }
                }
            }
            if (bti2 == -1) continue;
            bool take = false;
            if (bestg2 < bestg) take = true;
            else if (bestg2 == bestg) {
                if (bds2 < bds) take = true;
                else if (bds2 == bds) {
                    if (bl2 < bl) take = true;
                    else if (bl2 == bl) {
                        if (bdr2 < bdr) take = true;
                        else if (bdr2 == bdr) {
                            if (by02 < by0) take = true;
                            else if (by02 == by0) {
                                if (bx02 < bx0) take = true;
                                else if (bx02 == bx0 && randtie && rng.coin()) take = true;
                            }
                        }
                    }
                }
            }
            if (take) { bestg = bestg2; bti = bti2; bx = bx2; by = by2; bl = bl2; bds = bds2; bdr = bdr2; by0 = by02; bx0 = bx02; bestpos = pos; bestid = id; }
        }
        auto adapt = [&]() {
            if (adaptive && stepCnt == 5) {
                auto now = chrono::steady_clock::now();
                double elapsed = chrono::duration<double, milli>(now - tStart).count();
                double batch = chrono::duration<double, milli>(now - batchStart).count();
                double remT = max(0.0, adaptTLms - elapsed);
                int remSteps = max(1, nm - t);
                double budget = remT * 5.0 / remSteps;
                if (adaptMode == 1) {
                    emaBatch = (emaBatch < 0.0) ? batch : (0.6 * emaBatch + 0.4 * batch);
                    int floorL = max(1, dynLIM0 / 4);
                    if (emaBatch < budget * 0.85) dynLIM = min(maxBound, dynLIM + max(1, dynLIM / 8));
                    else if (emaBatch > budget * 1.15) dynLIM = max(floorL, dynLIM - max(1, dynLIM / 6));
                } else {
                    if (batch < budget) dynLIM = min(maxBound, dynLIM + 1); else dynLIM = max(1, dynLIM - 1);
                }
                batchStart = now;
                stepCnt = 0;
            }
        };
        if (bti == -1) {
            t++; stepCnt++;
            adapt();
            if (!panic && deadlineMs > 0 && elapsed_ms() > deadlineMs) return R{};
            continue;
        }
        auto& tsh = ps[bestid].t[bti];
        for (int j = 0; j < tsh.w; j++) {
            int nh = h[bx + j];
            if (tsh.hi[j] != INT_MIN) { int cand = by + tsh.hi[j]; if (nh < cand) nh = cand; }
            h[bx + j] = nh;
        }
        if (g < bl) g = bl;
        pl.push_back({bestid, bti, bx, by});
        if (bestpos != t) swap(o[t], o[bestpos]);
        t++;
        stepCnt++;
        adapt();
        if (!panic && deadlineMs > 0 && elapsed_ms() > deadlineMs) return R{};
    }
    if ((int)pl.size() != nm) return R{};  // validate all pieces placed
    int H = (int)g + 1;
    int maxX = -1;
    for (auto& pp : pl) {
        auto& t2 = ps[pp.idx].t[pp.ti];
        for (auto& q : t2.c) { int x = pp.x + q.first; if (x > maxX) maxX = x; }
    }
    int Wused = 0;
    if (maxX >= 0) {
        vector<char> used(maxX + 1, false);
        for (auto& pp : pl) { auto& t2 = ps[pp.idx].t[pp.ti]; for (auto& q : t2.c) used[pp.x + q.first] = true; }
        for (int x = 0; x <= maxX; x++) if (used[x]) Wused++;
    }
    int Wfinal = max(1, Wused);
    long long A = 1LL * H * Wfinal;
    R res; res.A = A; res.W = Wfinal; res.H = H; res.pl = move(pl); res.ok = true; res.packW = W;
    return res;
}

static vector<uint64_t> g_grid;
static R pack_blf(int W, const vector<int>& order, int policy, double deadlineMs) {
    if (W > 64 || W <= 0) return R{};
    int nm = (int)order.size();
    auto& grid = g_grid;
    size_t estH = (size_t)(S / max(1, W)) + 48;
    grid.assign(max<size_t>(64, estH), 0ULL);
    const uint64_t FULLROW = (W == 64) ? ~0ULL : ((1ULL << W) - 1);
    int y0full = 0;
    int topY = 0;
    vector<Pl> pl; pl.reserve(nm);
    int checkCnt = 0;
    for (int oi = 0; oi < nm; oi++) {
        if (deadlineMs > 0 && ((++checkCnt & 15) == 0) && elapsed_ms() > deadlineMs) return R{};
        int id = order[oi];
        auto& p = ps[id];
        int bestY = INT_MAX, bestX = 0, bestTi = -1; long long bestKey = LLONG_MAX;
        if ((int)grid.size() < topY + 12) grid.resize(topY + 64, 0ULL);
        for (int ti = 0; ti < (int)p.t.size(); ti++) {
            auto& t = p.t[ti];
            if (t.w > W) continue;
            uint64_t limMask = (W - t.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - t.w + 1)) - 1);
            int yFound = -1, xFound = -1;
            for (int y = y0full; y <= topY; y++) {
                uint64_t conflict = 0;
                for (int dy = 0; dy < t.h; dy++) {
                    uint64_t g = grid[y + dy];
                    if (!g) continue;
                    unsigned m = t.rmask[dy];
                    while (m) { int b = __builtin_ctz(m); conflict |= (g >> b); m &= m - 1; }
                }
                uint64_t allowed = ~conflict & limMask;
                if (allowed) { yFound = y; xFound = (int)__builtin_ctzll(allowed); break; }
            }
            if (yFound < 0) continue;
            long long key = (policy == 0) ? ((long long)yFound << 8) | xFound
                                          : ((long long)(yFound + t.h) << 24) | ((long long)yFound << 8) | xFound;
            if (key < bestKey) { bestKey = key; bestY = yFound; bestX = xFound; bestTi = ti; }
        }
        if (bestTi < 0) return R{};
        auto& t = p.t[bestTi];
        if ((int)grid.size() < bestY + t.h + 1) grid.resize(bestY + t.h + 64, 0ULL);
        for (int dy = 0; dy < t.h; dy++) grid[bestY + dy] |= ((uint64_t)t.rmask[dy]) << bestX;
        if (bestY + t.h > topY) topY = bestY + t.h;
        while (y0full < topY && grid[y0full] == FULLROW) y0full++;
        pl.push_back({id, bestTi, bestX, bestY});
    }
    uint64_t colU = 0; int rows = 0;
    for (int y = 0; y < topY; y++) if (grid[y]) { colU |= grid[y]; rows++; }
    int Wu = __builtin_popcountll(colU);
    R res; res.W = Wu ? Wu : 1; res.H = rows ? rows : 1;
    res.A = (long long)res.W * res.H; res.pl = move(pl); res.ok = true; res.packW = W;
    return res;
}

static int envInt(const char* name, int def) {
    if (const char* e = getenv(name)) return atoi(e);
    return def;
}

static vector<uint64_t> c_grid;
static bool crownRepack(R& r, double deadlineMs) {
    if (!r.ok || r.packW <= 0 || r.packW > 64) return false;
    int W = r.packW;
    int H0 = 0;
    for (auto& p : r.pl) {
        auto& t = ps[p.idx].t[p.ti];
        if (p.y + t.h > H0) H0 = p.y + t.h;
    }
    if (H0 <= 1) return false;
    auto& grid = c_grid;
    grid.assign(H0, 0ULL);
    for (auto& p : r.pl) {
        auto& t = ps[p.idx].t[p.ti];
        for (int dy = 0; dy < t.h; dy++) grid[p.y + dy] |= ((uint64_t)t.rmask[dy]) << p.x;
    }
    bool improvedAny = false;
    long long origA = r.A;
    int origW = r.W, origH = r.H;
    vector<Pl> origPl = r.pl;
    int rounds = 0;
    while (rounds++ < 64) {
        if (deadlineMs > 0 && elapsed_ms() > deadlineMs) break;
        int H = (int)grid.size();
        while (H > 0 && grid[H - 1] == 0) H--;
        if (H <= 1) break;
        bool anyDepthOk = false;
        for (int crownDepth = 1; crownDepth <= 3; crownDepth++) {
            if (H - crownDepth < 1) break;
            vector<int> crown;
            for (int i = 0; i < (int)r.pl.size(); i++) {
                auto& p = r.pl[i];
                auto& t = ps[p.idx].t[p.ti];
                if (p.y + t.h > H - crownDepth) crown.push_back(i);
            }
            if (crown.empty()) { anyDepthOk = true; break; } // nothing to relocate, H is fine
            for (int ci : crown) {
                auto& p = r.pl[ci];
                auto& t = ps[p.idx].t[p.ti];
                for (int dy = 0; dy < t.h; dy++) grid[p.y + dy] &= ~(((uint64_t)t.rmask[dy]) << p.x);
            }
            sort(crown.begin(), crown.end(), [&](int a, int b) { return ps[r.pl[a].idx].k > ps[r.pl[b].idx].k; });
            int targetH = H - crownDepth;
            vector<array<int,3>> newPos(crown.size()); // ti, x, y
            bool ok = true;
            for (size_t c = 0; c < crown.size() && ok; c++) {
                auto& p = r.pl[crown[c]];
                auto& piece = ps[p.idx];
                int bTi = -1, bX = 0, bY = INT_MAX;
                for (int ti = 0; ti < (int)piece.t.size(); ti++) {
                    auto& t = piece.t[ti];
                    if (t.w > W || t.h > targetH) continue;
                    uint64_t limMask = (W - t.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - t.w + 1)) - 1);
                    int ymax = targetH - t.h;
                    for (int y = 0; y <= ymax; y++) {
                        if (bY <= y) break;
                        uint64_t conflict = 0;
                        for (int dy = 0; dy < t.h; dy++) {
                            uint64_t g = grid[y + dy];
                            if (!g) continue;
                            unsigned m = t.rmask[dy];
                            while (m) { int b = __builtin_ctz(m); conflict |= (g >> b); m &= m - 1; }
                        }
                        uint64_t allowed = ~conflict & limMask;
                        if (allowed) {
                            int x = (int)__builtin_ctzll(allowed);
                            if (y < bY || (y == bY && x < bX)) { bY = y; bX = x; bTi = ti; }
                            break;
                        }
                    }
                }
                if (bTi < 0) { ok = false; break; }
                auto& t = piece.t[bTi];
                for (int dy = 0; dy < t.h; dy++) grid[bY + dy] |= ((uint64_t)t.rmask[dy]) << bX;
                newPos[c] = {bTi, bX, bY};
            }
            if (ok) {
                uint64_t colU = 0; int rows = 0;
                for (auto& g : grid) if (g) { colU |= g; rows++; }
                int newW = max(1, __builtin_popcountll(colU));
                int newH = max(1, rows);
                long long newA = (long long)newW * newH;
                if (newA < r.A) {
                    for (size_t c = 0; c < crown.size(); c++) {
                        auto& p = r.pl[crown[c]];
                        p.ti = newPos[c][0]; p.x = newPos[c][1]; p.y = newPos[c][2];
                    }
                    r.W = newW; r.H = newH; r.A = newA;
                    improvedAny = true;
                    anyDepthOk = true;
                    if (getenv("PP_DEBUG")) fprintf(stderr, "crown: shaved %d row(s) (crown size %zu) A=%lld->%lld\n", crownDepth, crown.size(), origA, newA);
                    break; // success at this depth, move to next round
                } else {
                    ok = false;
                }
            }
            if (!ok) {
                fill(grid.begin(), grid.end(), 0ULL);
                for (auto& p : r.pl) {
                    auto& t = ps[p.idx].t[p.ti];
                    for (int dy = 0; dy < t.h; dy++) grid[p.y + dy] |= ((uint64_t)t.rmask[dy]) << p.x;
                }
            }
        }
        if (!anyDepthOk) break; // all depths failed, can't improve further
    }
    if (!improvedAny) {
        r.pl = move(origPl);
        r.W = origW; r.H = origH; r.A = origA;
    }
    return improvedAny;
}

static R pack_blf2(int W, const vector<int>& order, int window, double deadlineMs, RNG& rng, bool randtie) {
    if (W > 64 || W <= 0) return R{};
    int nm = (int)order.size();
    auto& grid = g_grid;
    size_t estH = (size_t)(S / max(1, W)) + 48;
    grid.assign(max<size_t>(64, estH), 0ULL);
    const uint64_t FULLROW = (W == 64) ? ~0ULL : ((1ULL << W) - 1);
    int y0full = 0, topY = 0;
    long long g = -1; // current max occupied row index
    vector<int> o = order;
    vector<Pl> pl; pl.reserve(nm);
    int checkCnt = 0;
    for (int t = 0; t < nm; t++) {
        if (deadlineMs > 0 && ((++checkCnt & 3) == 0) && elapsed_ms() > deadlineMs) return R{};
        int lim = min(nm, t + max(1, window));
        long long bG = LLONG_MAX; int bContact = -1, bY = INT_MAX, bX = INT_MAX, bTi = -1, bPos = -1;
        if ((int)grid.size() < topY + 12) grid.resize(topY + 64, 0ULL);
        for (int pos = t; pos < lim; pos++) {
            auto& p = ps[o[pos]];
            for (int ti = 0; ti < (int)p.t.size(); ti++) {
                auto& tt = p.t[ti];
                if (tt.w > W) continue;
                uint64_t limMask = (W - tt.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - tt.w + 1)) - 1);
                int yF = -1, xF = -1;
                for (int y = y0full; y <= topY; y++) {
                    uint64_t conflict = 0;
                    for (int dy = 0; dy < tt.h; dy++) {
                        uint64_t gg = grid[y + dy];
                        if (!gg) continue;
                        unsigned m = tt.rmask[dy];
                        while (m) { int b = __builtin_ctz(m); conflict |= (gg >> b); m &= m - 1; }
                    }
                    uint64_t allowed = ~conflict & limMask;
                    if (allowed) { yF = y; xF = (int)__builtin_ctzll(allowed); break; }
                }
                if (yF < 0) continue;
                long long g2 = max(g, (long long)(yF + tt.h - 1));
                if (g2 > bG) continue;
                int contact = 0;
                for (auto& nb : tt.nbr) {
                    int cx = xF + nb.first, cy = yF + nb.second;
                    if (cx < 0 || cx >= W || cy < 0) { contact++; continue; }
                    if (cy < (int)grid.size() && (grid[cy] >> cx) & 1ULL) contact++;
                }
                bool take = false;
                if (g2 < bG) take = true;
                else { // g2 == bG
                    if (contact > bContact) take = true;
                    else if (contact == bContact) {
                        if (yF < bY) take = true;
                        else if (yF == bY) {
                            if (xF < bX) take = true;
                            else if (xF == bX && randtie && rng.coin()) take = true;
                        }
                    }
                }
                if (take) { bG = g2; bContact = contact; bY = yF; bX = xF; bTi = ti; bPos = pos; }
            }
        }
        if (bPos < 0) return R{};
        auto& p = ps[o[bPos]];
        auto& tt = p.t[bTi];
        if ((int)grid.size() < bY + tt.h + 1) grid.resize(bY + tt.h + 64, 0ULL);
        for (int dy = 0; dy < tt.h; dy++) grid[bY + dy] |= ((uint64_t)tt.rmask[dy]) << bX;
        if (bY + tt.h > topY) topY = bY + tt.h;
        while (y0full < topY && grid[y0full] == FULLROW) y0full++;
        if (bG > g) g = bG;
        pl.push_back({o[bPos], bTi, bX, bY});
        if (bPos != t) swap(o[t], o[bPos]);
    }
    uint64_t colU = 0; int rows = 0;
    for (int y = 0; y < topY; y++) if (grid[y]) { colU |= grid[y]; rows++; }
    int Wu = __builtin_popcountll(colU);
    R res; res.W = Wu ? Wu : 1; res.H = rows ? rows : 1;
    res.A = (long long)res.W * res.H; res.pl = move(pl); res.ok = true; res.packW = W;
    return res;
}

struct B3Cache { int y, x, contact; bool valid; };
static vector<B3Cache> b3cache; // [pieceId * 8 + orientIdx]
static R pack_blf3(int W, const vector<int>& order, int window, double deadlineMs, RNG& rng, bool randtie) {
    if (W > 64 || W <= 0) return R{};
    int nm = (int)order.size();
    auto& grid = g_grid;
    size_t estH = (size_t)(S / max(1, W)) + 48;
    grid.assign(max<size_t>(64, estH), 0ULL);
    const uint64_t FULLROW = (W == 64) ? ~0ULL : ((1ULL << W) - 1);
    int y0full = 0, topY = 0;
    long long g = -1;
    vector<int> o = order;
    vector<Pl> pl; pl.reserve(nm);
    if ((int)b3cache.size() < nm * 8) b3cache.resize(nm * 8);
    for (int i = 0; i < nm * 8; i++) b3cache[i].valid = false;
    int checkCnt = 0;
    for (int t = 0; t < nm; t++) {
        if (deadlineMs > 0 && ((++checkCnt & 15) == 0) && elapsed_ms() > deadlineMs) return R{};
        int lim = min(nm, t + max(1, window));
        long long bG = LLONG_MAX; int bContact = -1, bY = INT_MAX, bX = INT_MAX, bTi = -1, bPos = -1;
        if ((int)grid.size() < topY + 12) grid.resize(topY + 64, 0ULL);
        for (int pos = t; pos < lim; pos++) {
            int pid = o[pos];
            auto& p = ps[pid];
            for (int ti = 0; ti < (int)p.t.size(); ti++) {
                auto& tt = p.t[ti];
                if (tt.w > W) continue;
                B3Cache& cc = b3cache[pid * 8 + ti];
                if (!cc.valid) {
                    uint64_t limMask = (W - tt.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - tt.w + 1)) - 1);
                    int yF = -1, xF = -1;
                    for (int y = y0full; y <= topY; y++) {
                        uint64_t conflict = 0;
                        for (int dy = 0; dy < tt.h; dy++) {
                            uint64_t gg = grid[y + dy];
                            if (!gg) continue;
                            unsigned m = tt.rmask[dy];
                            while (m) { int b = __builtin_ctz(m); conflict |= (gg >> b); m &= m - 1; }
                        }
                        uint64_t allowed = ~conflict & limMask;
                        if (allowed) { yF = y; xF = (int)__builtin_ctzll(allowed); break; }
                    }
                    if (yF < 0) { continue; }
                    int contact = 0;
                    for (auto& nb : tt.nbr) {
                        int cx = xF + nb.first, cy = yF + nb.second;
                        if (cx < 0 || cx >= W || cy < 0) { contact++; continue; }
                        if (cy < (int)grid.size() && (grid[cy] >> cx) & 1ULL) contact++;
                    }
                    cc = {yF, xF, contact, true};
                }
                long long g2 = max(g, (long long)(cc.y + tt.h - 1));
                if (g2 > bG) continue;
                bool take = false;
                if (g2 < bG) take = true;
                else {
                    if (cc.contact > bContact) take = true;
                    else if (cc.contact == bContact) {
                        if (cc.y < bY) take = true;
                        else if (cc.y == bY) {
                            if (cc.x < bX) take = true;
                            else if (cc.x == bX && randtie && rng.coin()) take = true;
                        }
                    }
                }
                if (take) { bG = g2; bContact = cc.contact; bY = cc.y; bX = cc.x; bTi = ti; bPos = pos; }
            }
        }
        if (bPos < 0) return R{};
        int placedId = o[bPos];
        auto& tt = ps[placedId].t[bTi];
        if ((int)grid.size() < bY + tt.h + 1) grid.resize(bY + tt.h + 64, 0ULL);
        for (int dy = 0; dy < tt.h; dy++) grid[bY + dy] |= ((uint64_t)tt.rmask[dy]) << bX;
        if (bY + tt.h > topY) topY = bY + tt.h;
        while (y0full < topY && grid[y0full] == FULLROW) y0full++;
        if (bG > g) g = bG;
        pl.push_back({placedId, bTi, bX, bY});
        if (bPos != t) swap(o[t], o[bPos]);
        int px0 = bX - 1, px1 = bX + tt.w, py0 = bY - 1, py1 = bY + tt.h;
        int lim2 = min(nm, (t + 1) + max(1, window));
        for (int pos = t + 1; pos < lim2; pos++) {
            int pid = o[pos];
            auto& p = ps[pid];
            for (int ti = 0; ti < (int)p.t.size(); ti++) {
                B3Cache& cc = b3cache[pid * 8 + ti];
                if (!cc.valid) continue;
                auto& ct = p.t[ti];
                if (cc.x <= px1 && cc.x + ct.w - 1 >= px0 && cc.y <= py1 && cc.y + ct.h - 1 >= py0)
                    cc.valid = false;
            }
        }
    }
    uint64_t colU = 0; int rows = 0;
    for (int y = 0; y < topY; y++) if (grid[y]) { colU |= grid[y]; rows++; }
    int Wu = __builtin_popcountll(colU);
    R res; res.W = Wu ? Wu : 1; res.H = rows ? rows : 1;
    res.A = (long long)res.W * res.H; res.pl = move(pl); res.ok = true; res.packW = W;
    return res;
}