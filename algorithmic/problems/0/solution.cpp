#include <algorithm>
#include <chrono>
#include <cmath>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
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

static vector<short> g_own;
static vector<int> g_plSlot;
static R pack_capped(int W, int Hcap, const vector<int>& order, int window, double deadlineMs, RNG& rng) {
    if (W > 64 || W <= 0 || Hcap <= 0) return R{};
    int nm = (int)order.size();
    auto& grid = g_grid;
    int gridH = Hcap + 12;
    grid.assign(gridH, 0ULL);
    const uint64_t FULLROW = (W == 64) ? ~0ULL : ((1ULL << W) - 1);
    int y0full = 0, topY = 0;
    vector<int> o = order;
    vector<Pl> pl; pl.reserve(nm);
    if ((int)b3cache.size() < nm * 8) b3cache.resize(nm * 8);
    for (int i = 0; i < nm * 8; i++) b3cache[i].valid = false;
    static int REPAIR = envInt("PP_REPAIR", 25); // ejection-chain depth (0 = off)
    auto& own = g_own;
    auto& plSlot = g_plSlot;
    if (REPAIR > 0) {
        own.assign((size_t)gridH * W, -1);
        if ((int)plSlot.size() < nm) plSlot.resize(nm);
        for (int i = 0; i < nm; i++) plSlot[i] = -1;
    }
    int checkCnt = 0;
    auto place = [&](int pid, int ti, int x, int y) {
        auto& tt = ps[pid].t[ti];
        for (int dy = 0; dy < tt.h; dy++) grid[y + dy] |= ((uint64_t)tt.rmask[dy]) << x;
        if (REPAIR > 0) {
            for (auto& q : tt.c) own[(size_t)(y + q.second) * W + (x + q.first)] = (short)pid;
        }
        if (y + tt.h > topY) topY = min(y + tt.h, Hcap);
        while (y0full < topY && grid[y0full] == FULLROW) y0full++;
        if (plSlot.size() && REPAIR > 0 && plSlot[pid] >= 0) pl[plSlot[pid]] = {pid, ti, x, y};
        else { if (REPAIR > 0) plSlot[pid] = (int)pl.size(); pl.push_back({pid, ti, x, y}); }
    };
    auto eject = [&](int pid) {
        int slot = plSlot[pid];
        Pl p = pl[slot];
        auto& tt = ps[pid].t[p.ti];
        for (int dy = 0; dy < tt.h; dy++) grid[p.y + dy] &= ~(((uint64_t)tt.rmask[dy]) << p.x);
        for (auto& q : tt.c) own[(size_t)(p.y + q.second) * W + (p.x + q.first)] = -1;
        if (p.y < y0full) y0full = p.y;
    };
    auto invalidateWindow = [&](int t) {
        int lim2 = min(nm, t + max(1, window) + 1);
        for (int pos = t; pos < lim2; pos++) {
            int pid = o[pos];
            for (int ti = 0; ti < (int)ps[pid].t.size(); ti++) b3cache[pid * 8 + ti].valid = false;
        }
    };
    auto scanFit = [&](const T& tt) -> pair<int,int> {
        uint64_t limMask = (W - tt.w + 1 >= 64) ? ~0ULL : ((1ULL << (W - tt.w + 1)) - 1);
        int ymax = min(topY, Hcap - tt.h);
        for (int y = y0full; y <= ymax; y++) {
            uint64_t conflict = 0;
            for (int dy = 0; dy < tt.h; dy++) {
                uint64_t gg = grid[y + dy];
                if (!gg) continue;
                unsigned m = tt.rmask[dy];
                while (m) { int b = __builtin_ctz(m); conflict |= (gg >> b); m &= m - 1; }
            }
            uint64_t allowed = ~conflict & limMask;
            if (allowed) return {y, (int)__builtin_ctzll(allowed)};
        }
        return {-1, -1};
    };
    auto tryRepair = [&](int stuckPid, int tCur) -> bool {
        int pending = stuckPid;
        int lastPlaced = -1;
        for (int depth = 0; depth <= REPAIR; depth++) {
            if (deadlineMs > 0 && elapsed_ms() > deadlineMs) return false;
            auto& p = ps[pending];
            if (depth > 0) { // pending may fit normally now
                int bTi = -1, bYy = INT_MAX, bXx = 0;
                for (int ti = 0; ti < (int)p.t.size(); ti++) {
                    auto& tt = p.t[ti];
                    if (tt.w > W || tt.h > Hcap) continue;
                    auto [y, x] = scanFit(tt);
                    if (y >= 0 && (y < bYy || (y == bYy && x < bXx))) { bYy = y; bXx = x; bTi = ti; }
                }
                if (bTi >= 0) { place(pending, bTi, bXx, bYy); return true; }
            }
            if (depth == REPAIR) return false;
            int cTi = -1, cY = -1, cX = -1; int cBlocker = -1; int cCells = INT_MAX;
            for (int ti = 0; ti < (int)p.t.size(); ti++) {
                auto& tt = p.t[ti];
                if (tt.w > W || tt.h > Hcap) continue;
                int ymax = min(topY, Hcap - tt.h);
                for (int y = 0; y <= ymax && (cY < 0 || y <= cY); y++) {
                    for (int x = 0; x + tt.w <= W; x++) {
                        int blocker = -1, bcells = 0; bool bad = false;
                        for (auto& q : tt.c) {
                            size_t idx = (size_t)(y + q.second) * W + (x + q.first);
                            short ow = own[idx];
                            if (ow < 0) continue;
                            if (ow == pending) { bad = true; break; }
                            if (ow == lastPlaced) { bad = true; break; }
                            if (blocker < 0) blocker = ow;
                            else if (blocker != ow) { bad = true; break; }
                            bcells++;
                        }
                        if (bad || blocker < 0) continue;
                        if (cY < 0 || y < cY || (y == cY && bcells < cCells)) {
                            cTi = ti; cY = y; cX = x; cBlocker = blocker; cCells = bcells;
                        }
                    }
                }
            }
            if (cBlocker < 0) return false;
            eject(cBlocker);
            place(pending, cTi, cX, cY);
            lastPlaced = pending;
            pending = cBlocker;
        }
        return false;
    };
    for (int t = 0; t < nm; t++) {
        if (deadlineMs > 0 && ((++checkCnt & 15) == 0) && elapsed_ms() > deadlineMs) return R{};
        int lim = min(nm, t + max(1, window));
        int bContact = -1, bY = INT_MAX, bX = INT_MAX, bTi = -1, bPos = -1;
        bool redo = false;
        for (int pos = t; pos < lim; pos++) {
            int pid = o[pos];
            auto& p = ps[pid];
            bool any = false;
            for (int ti = 0; ti < (int)p.t.size(); ti++) {
                auto& tt = p.t[ti];
                if (tt.w > W || tt.h > Hcap) continue;
                B3Cache& cc = b3cache[pid * 8 + ti];
                if (!cc.valid) {
                    auto [yF, xF] = scanFit(tt);
                    if (yF < 0) { cc.valid = false; continue; }
                    int contact = 0;
                    for (auto& nb : tt.nbr) {
                        int cx = xF + nb.first, cy = yF + nb.second;
                        if (cx < 0 || cx >= W || cy < 0) { contact++; continue; }
                        if (cy >= Hcap) { contact++; continue; } // cap edge counts as wall
                        if (cy < (int)grid.size() && (grid[cy] >> cx) & 1ULL) contact++;
                    }
                    cc = {yF, xF, contact, true};
                }
                any = true;
                static int CAPSEL = envInt("PP_CAPSEL", 1);
                bool take = false;
                if (CAPSEL == 0) {
                    if (cc.y < bY) take = true;
                    else if (cc.y == bY) {
                        if (cc.contact > bContact) take = true;
                        else if (cc.contact == bContact && cc.x < bX) take = true;
                    }
                } else {
                    if (cc.contact > bContact) take = true;
                    else if (cc.contact == bContact) {
                        if (cc.y < bY) take = true;
                        else if (cc.y == bY && cc.x < bX) take = true;
                    }
                }
                if (take) { bContact = cc.contact; bY = cc.y; bX = cc.x; bTi = ti; bPos = pos; }
            }
            if (!any) {
                if (REPAIR > 0 && tryRepair(pid, t)) {
                    if (pos != t) swap(o[t], o[pos]);
                    invalidateWindow(t);
                    redo = true;
                    break;
                }
                return R{}; // permanently infeasible under the cap
            }
        }
        if (redo) continue; // t stays: the placed piece occupies o[t]... advance below
        if (bPos < 0) return R{};
        int placedId = o[bPos];
        place(placedId, bTi, bX, bY);
        if (bPos != t) swap(o[t], o[bPos]);
        int px0 = bX - 1, px1 = bX + ps[placedId].t[bTi].w, py0 = bY - 1, py1 = bY + ps[placedId].t[bTi].h;
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
    for (int y = 0; y < min((int)grid.size(), Hcap); y++) if (grid[y]) { colU |= grid[y]; rows++; }
    int Wu = __builtin_popcountll(colU);
    R res; res.W = Wu ? Wu : 1; res.H = rows ? rows : 1;
    res.A = (long long)res.W * res.H; res.pl = move(pl); res.ok = true; res.packW = W;
    return res;
}

static int BF_HD = 40;
static vector<uint64_t> bf_occ; static vector<int> bf_colH; static int bf_W;
static inline int bf_fit(const T& o, int x) {
    int ysky = 0; for (int j = 0; j < o.w; j++) { int v = bf_colH[x + j] - o.lo[j]; if (v > ysky) ysky = v; }
    int start = ysky - BF_HD; if (start < 0) start = 0;
    for (int y = start; y <= ysky; y++) { bool bad = false;
        for (int dy = 0; dy < o.h && !bad; dy++) { uint64_t rm = ((uint64_t)o.rmask[dy]) << x; uint64_t orow = (y + dy < (int)bf_occ.size()) ? bf_occ[y + dy] : 0ULL; if (rm & orow) bad = true; }
        if (!bad) return y; }
    return ysky;
}
static inline void bf_put(const T& o, int x, int y) {
    if ((int)bf_occ.size() < y + o.h) bf_occ.resize(y + o.h, 0ULL);
    for (int dy = 0; dy < o.h; dy++) bf_occ[y + dy] |= ((uint64_t)o.rmask[dy]) << x;
    for (int j = 0; j < o.w; j++) { int top = y + o.hi[j] + 1; if (top > bf_colH[x + j]) bf_colH[x + j] = top; }
}
static inline uint64_t bf_key(int w, const int* v) { uint64_t k = (uint64_t)w << 52; for (int j = 0; j < w; j++) { int s = v[j] + 16; if (s < 0) s = 0; if (s > 31) s = 31; k |= (uint64_t)s << (5 * j); } return k; }
static int bf_pass(int W, const vector<int>& repr, const unordered_map<uint64_t, vector<pair<int,int>>>& idx,
                   vector<int>& avail, int total, vector<int>& outKind, vector<int>& outOri, vector<int>& outX, vector<int>& outY, int tie) {
    bf_W = W; bf_occ.assign(8, 0ULL); bf_colH.assign(W, 0);
    outKind.clear(); outOri.clear(); outX.clear(); outY.clear();
    int rem = total; int p[12];
    long long bH, bW; int bC, bk, bo, bx, by, curH0;
    auto eval = [&](int ki, int oi, int x) {
        const T& o = ps[repr[ki]].t[oi];
        int y = bf_fit(o, x);
        int rH = curH0; for (int j = 0; j < o.w; j++) { int t = y + o.hi[j] + 1; if (t > rH) rH = t; }
        long long aw = 0; for (int j = 0; j < o.w; j++) { int pb = y + o.lo[j]; if (pb > bf_colH[x + j]) aw += pb - bf_colH[x + j]; }
        int ct = 0; for (auto& nb : o.nbr) { int cx = x + nb.first, cy = y + nb.second; if (cx < 0 || cx >= W || cy < 0) { ct++; continue; } if (cy < (int)bf_occ.size() && ((bf_occ[cy] >> cx) & 1ULL)) ct++; }
        bool cbetter = tie ? (ct >= bC) : (ct > bC);
        if (rH < bH || (rH == bH && (aw < bW || (aw == bW && cbetter)))) { bH = rH; bW = aw; bC = ct; bk = ki; bo = oi; bx = x; by = y; }
    };
    while (rem > 0) {
        int mh = INT_MAX; curH0 = 0; for (int c = 0; c < W; c++) { if (bf_colH[c] < mh) mh = bf_colH[c]; if (bf_colH[c] > curH0) curH0 = bf_colH[c]; }
        int nc = 0; while (nc < W && bf_colH[nc] != mh) nc++;
        bH = LLONG_MAX; bW = LLONG_MAX; bC = -1; bk = -1; bo = -1; bx = -1; by = -1;
        for (int w = 1; w <= 10; w++) {
            int xlo = nc - w + 1; if (xlo < 0) xlo = 0; int xhi = nc; if (xhi > W - w) xhi = W - w;
            for (int x = xlo; x <= xhi; x++) {
                int base = bf_colH[x]; for (int j = 0; j < w; j++) p[j] = bf_colH[x + j] - base;
                auto it = idx.find(bf_key(w, p)); if (it == idx.end()) continue;
                for (auto& pr : it->second) if (avail[pr.first] > 0) eval(pr.first, pr.second, x);
            }
        }
        if (bk < 0) // no indexed candidate: full scan over kinds (rare)
            for (int ki = 0; ki < (int)repr.size(); ki++) { if (avail[ki] <= 0) continue;
                const P& Kp = ps[repr[ki]];
                for (int oi = 0; oi < (int)Kp.t.size(); oi++) { const T& o = Kp.t[oi]; if (o.w > W) continue;
                    int xlo = nc - o.w + 1; if (xlo < 0) xlo = 0; int xhi = nc; if (xhi > W - o.w) xhi = W - o.w;
                    for (int x = xlo; x <= xhi; x++) eval(ki, oi, x);
                }
            }
        if (bk < 0) return -1;
        bf_put(ps[repr[bk]].t[bo], bx, by); avail[bk]--; rem--;
        outKind.push_back(bk); outOri.push_back(bo); outX.push_back(bx); outY.push_back(by);
    }
    int H = 0; for (int c = 0; c < W; c++) if (bf_colH[c] > H) H = bf_colH[c];
    return H;
}
static R bfSolve(int minW, int base, double deadline) {
    unordered_map<string, int> km; vector<int> repr; vector<int> kcnt; vector<int> kindOf(n);
    vector<vector<int>> members;
    string s; s.reserve(64);
    for (int i = 0; i < n; i++) {
        const string* bestk = nullptr; string keep;
        for (auto& t : ps[i].t) { s.clear(); for (auto& q : t.c) { s.push_back((char)('0' + q.first)); s.push_back((char)('0' + q.second)); }
            if (!bestk || s < keep) { keep = s; bestk = &keep; } }
        auto it = km.find(keep);
        int ki; if (it == km.end()) { ki = repr.size(); km[keep] = ki; repr.push_back(i); kcnt.push_back(0); members.push_back({}); }
        else ki = it->second;
        kcnt[ki]++; members[ki].push_back(i); kindOf[i] = ki;
    }
    int K = repr.size();
    unordered_map<uint64_t, vector<pair<int,int>>> idx; idx.reserve(K * 4);
    int rp[12];
    for (int ki = 0; ki < K; ki++) for (int oi = 0; oi < (int)ps[repr[ki]].t.size(); oi++) {
        const T& o = ps[repr[ki]].t[oi]; if (o.w > 63) continue; int b0 = o.lo[0];
        for (int j = 0; j < o.w; j++) rp[j] = o.lo[j] - b0;
        idx[bf_key(o.w, rp)].push_back({ki, oi});
    }
    long long bestA = LLONG_MAX; int bestW = 0, bestH = 0; vector<int> bK, bO, bX, bY, tK, tO, tX, tY;
    double lastPass = 0;
    for (int d = 0; d <= 24; d++) {
        for (int sgn = (d ? -1 : 1); sgn <= 1; sgn += 2) {
            int W = base + sgn * d; if (W < minW || W > 63) continue;
            for (int tie = 1; tie >= 0; tie--) {   // try both tie-break directions, keep best area
                if (elapsed_ms() + lastPass * 1.3 > deadline) goto DONE;
                vector<int> avail = kcnt; double t0 = elapsed_ms();
                int H = bf_pass(W, repr, idx, avail, n, tK, tO, tX, tY, tie);
                lastPass = elapsed_ms() - t0;
                if (H > 0) { long long A = (long long)W * H; if (A < bestA) { bestA = A; bestW = W; bestH = H; swap(bK, tK); swap(bO, tO); swap(bX, tX); swap(bY, tY); } }
            }
        }
    }
    DONE:;
    if (getenv("PP_DEBUG")) fprintf(stderr, "bf: bestW=%d base=%d H=%d fill=%.4f K=%d idx=%zu lastPass=%.0f t=%.0f\n", bestW, base, bestH, (double)S/((double)bestW*bestH), K, idx.size(), lastPass, elapsed_ms());
    R r; if (bestW == 0) { r.ok = false; return r; }
    r.W = bestW; r.H = bestH; r.A = bestA; r.ok = true; r.packW = bestW;
    vector<int> mcur(K, 0);
    r.pl.reserve(n);
    for (size_t z = 0; z < bK.size(); z++) { int ki = bK[z]; int pi = members[ki][mcur[ki]++]; r.pl.push_back({pi, bO[z], bX[z], bY[z]}); }
    return r;
}

int main() {
    T0 = chrono::steady_clock::now();
    if (const char* e = getenv("POLYPACK_TL")) { double v = atof(e); if (v > 50 && v < 10000) TL_MS = v; }
    int ffEnv = envInt("PP_FASTFIT", -1);    // -1 => auto-gate by S below; 0/1 => explicit override
    int BIGBLF = envInt("PP_BIGBLF", 0);     // 1: big cases skip champion pack, BLF-only
    int SMALLBLF = envInt("PP_SMALLBLF", 0); // 1: small cases skip champion sweep
    int JUMP = envInt("PP_JUMP", 15);        // % chance of W jump in restarts
    int BLF2 = envInt("PP_BLF2", 1);         // 1: restarts use blf2 (hole-aware window best-fit)
    int BLF2SWEEP = envInt("PP_BLF2SWEEP", 0); // 1: small-case sweep uses blf2 instead of skyline
    int BLF2ORD = envInt("PP_BLF2ORD", 1);   // 0: big-first order, 1: champion order
    int B3 = envInt("PP_B3", 1);             // 1: use cached blf3 instead of blf2
    int PHASE2 = envInt("PP_PHASE2", 1);     // 1: blf2 pass over best sweep Ws
    double P2FRAC = envInt("PP_P2FRAC", 0) / 100.0;  // 0 => auto by size
    double P2ENDF = envInt("PP_P2END", 0) / 100.0;   // 0 => auto by size
    long long P2MAXS = envInt("PP_P2MAXS", 22000);   // phase2 only when S below this
    long long B2RESTS = envInt("PP_B2RESTS", 50000);  // blf2 restarts when S below this (was 1300)

    {
        size_t cap = 1 << 20; inbuf.resize(cap); size_t len = 0;
        while (true) {
            if (len == cap) { cap <<= 1; inbuf.resize(cap); }
            size_t got = fread(inbuf.data() + len, 1, cap - len, stdin);
            if (got == 0) break;
            len += got;
        }
        inbuf.resize(len);
    }
    n = readInt();
    if (n <= 0) return 0;
    ps.resize(n);
    for (int i = 0; i < n; i++) {
        int k = readInt();
        ps[i].id = i + 1; ps[i].k = k; ps[i].b.resize(k);
        for (int j = 0; j < k; j++) { int x = readInt(), y = readInt(); ps[i].b[j] = {x, y}; }
        S += k;
    }
    gFASTFIT = (ffEnv < 0) ? (S >= 12000 ? 1 : 0) : ffEnv;
    if (P2FRAC <= 0.0) P2FRAC = (S < 6000) ? 0.25 : 0.55;
    if (P2ENDF <= 0.0) P2ENDF = (S < 6000) ? 0.35 : 0.70;
    for (int i = 0; i < n; i++) {
        auto& p = ps[i];
        unordered_set<string> seen; seen.reserve(32);
        for (int rf = 0; rf < 2; rf++) {
            vector<pair<int,int>> src = p.b;
            if (rf) for (auto& q : src) q.first = -q.first;
            for (int r = 0; r < 4; r++) {
                vector<pair<int,int>> v = src;
                for (auto& q : v) q = rotp(q, r);
                int minx = INT_MAX, miny = INT_MAX, maxx = INT_MIN, maxy = INT_MIN;
                for (auto& q : v) { minx = min(minx, q.first); miny = min(miny, q.second); maxx = max(maxx, q.first); maxy = max(maxy, q.second); }
                vector<pair<int,int>> v2 = v;
                for (auto& q : v2) { q.first -= minx; q.second -= miny; }
                sort(v2.begin(), v2.end());
                string key; key.reserve(v2.size() * 8);
                for (auto& q : v2) { key.append(to_string(q.first)); key.push_back(','); key.append(to_string(q.second)); key.push_back(';'); }
                if (seen.insert(key).second) {
                    T t; t.w = maxx - minx + 1; t.h = maxy - miny + 1; t.c = v2; t.r = r; t.f = rf; t.minx = minx; t.miny = miny;
                    t.lo.assign(t.w, INT_MAX); t.hi.assign(t.w, INT_MIN);
                    memset(t.rmask, 0, sizeof(t.rmask));
                    for (auto& q : v2) {
                        int x = q.first, y = q.second;
                        if (t.lo[x] > y) t.lo[x] = y; if (t.hi[x] < y) t.hi[x] = y;
                        t.rmask[y] |= (unsigned short)(1u << x);
                    }
                    {   // neighbor offsets (adjacent to piece, not in piece)
                        set<pair<int,int>> inp(v2.begin(), v2.end()), nbs;
                        const int DX[4] = {1, -1, 0, 0}, DY[4] = {0, 0, 1, -1};
                        for (auto& q : v2)
                            for (int d = 0; d < 4; d++) {
                                pair<int,int> np = {q.first + DX[d], q.second + DY[d]};
                                if (!inp.count(np)) nbs.insert(np);
                            }
                        for (auto& np : nbs) t.nbr.push_back({(signed char)np.first, (signed char)np.second});
                    }
                    p.t.push_back(move(t));
                }
            }
        }
        for (auto& t : p.t) { p.minW = min(p.minW, t.w); p.minH = min(p.minH, t.h); p.minA = min(p.minA, t.w * t.h); }
        if (p.t.empty()) { T t; t.w = 1; t.h = 1; t.c = {{0, 0}}; t.lo = {0}; t.hi = {0}; t.r = 0; t.f = 0; t.minx = 0; t.miny = 0; memset(t.rmask, 0, sizeof(t.rmask)); t.rmask[0] = 1; p.t.push_back(t); p.minW = 1; p.minH = 1; p.minA = 1; }
    }

    vector<int> idx(n); iota(idx.begin(), idx.end(), 0);
    unsigned long long seed = 0x9e3779b97f4a7c15ULL ^ (S << 1) ^ (unsigned long long)(n * 1469598103934665603ULL);
    RNG rng(seed);
    auto ord4 = [&]() {
        vector<int> res = idx;
        stable_sort(res.begin(), res.end(), [&](int a, int b) {
            int da = min(ps[a].minW, ps[a].minH);
            int db = min(ps[b].minW, ps[b].minH);
            if (da != db) return da < db;
            if (ps[a].k != ps[b].k) return ps[a].k < ps[b].k;
            return ps[a].id > ps[b].id;
        });
        return res;
    };

    int minW = 0; for (auto& p : ps) minW = max(minW, p.minW);
    double factor;
    if (S < 1000) factor = 0.4;
    else if (S < 3000) factor = 0.5;
    else if (S < 10000) factor = 0.27;
    else if (S < 30000) factor = 0.08;
    else if (S < 50000) factor = 0.028;  // measured optimum W~33 at S~38k (was 0.01 -> W~19)
    else factor = 0.009;                 // measured optimum W~24-28 at S~58-96k
    int base = max(minW, (int)floor(sqrt((double)S * factor)));
    if (const char* e = getenv("PP_BASEW")) { int v = atoi(e); if (v >= minW && v <= 4000) base = v; }

    vector<int> Ws;
    {
        unordered_set<int> used; used.reserve(512);
        auto addW = [&](int w) { if (w < minW) w = minW; if (used.insert(w).second) Ws.push_back(w); };
        addW(base);
        int span = min(96, max(20, base / 2));
        for (int d = 1; d <= span; d++) { addW(base - d); addW(base + d); }
        addW(minW);
        addW((int)max<long long>(minW, (S + base - 1) / base));
        for (int m = 2; m <= 6; m++) { addW(base * m / 3); addW((int)max<long long>(minW, S / ((base * m / 3) ? (base * m / 3) : 1))); }
        sort(Ws.begin(), Ws.end(), [&](int a, int b) { int da = abs(a - base), db = abs(b - base); if (da != db) return da < db; return a < b; });
    }

    bool big = S > (long long)envInt("PP_BIGTHRESH", 70000);
    vector<int> baseOrder = ord4();

    auto better = [&](const R& a, const R& b) {
        if (!b.ok) return true;
        if (a.A != b.A) return a.A < b.A;
        if (a.H != b.H) return a.H < b.H;
        return a.W < b.W;
    };
    R bestR;
    int bfEnv = envInt("PP_BF", -1);
    long long BF_N = envInt("PP_BFN", 450);
    bool useBF = (bfEnv < 0) ? (n >= BF_N && base <= 63 && minW <= 63) : (bfEnv > 0);
    if (useBF) { bestR = bfSolve(minW, min(base, 63), TL_MS - 120.0); }
    if (!bestR.ok) {
    double tFB0 = elapsed_ms();
    bestR = pack(base, baseOrder, rng, false, 1, false, 0.0, -1.0);
    double tFB = max(0.05, elapsed_ms() - tFB0); // measured window=1 pack cost
    int wprobeEnv = envInt("PP_WPROBE", 0);
    bool doWprobe = (wprobeEnv > 0 && S > 7000);
    if (doWprobe) {
        long long bestProbe = bestR.ok ? bestR.A : LLONG_MAX;
        int bestW = base;
        double probeBudget = TL_MS * 0.12;
        for (int d = 1; d <= 4 && elapsed_ms() < probeBudget; d++) {
            for (int sgn = -1; sgn <= 1 && elapsed_ms() < probeBudget; sgn += 2) {
                int W = base + sgn * d;
                if (W < minW || W > 4000) continue;
                R pr = pack(W, baseOrder, rng, false, 1, false, 0.0, probeBudget + 50.0);
                if (pr.ok) {
                    if (better(pr, bestR)) { bestProbe = pr.A; bestW = W; bestR = move(pr); }
                    else if (pr.A < bestProbe) { bestProbe = pr.A; bestW = W; }
                }
            }
        }
        base = bestW;
        if (base != Ws[0]) {
            vector<int> W2; W2.push_back(base);
            for (int w : Ws) if (w != base) W2.push_back(w);
            Ws = move(W2);
        }
    }
    if (!bestR.ok) {
        R r; r.pl.reserve(n);
        int y = 0; int wmax = 1;
        for (int i = 0; i < n; i++) { auto& t = ps[i].t[0]; r.pl.push_back({i, 0, 0, y}); y += t.h; wmax = max(wmax, t.w); }
        r.W = wmax; r.H = y; r.A = 1LL * wmax * y; r.ok = true;
        bestR = move(r);
    }

    const double SEARCH_END = TL_MS;           // hard abort for any pack
    const double SOFT_END = TL_MS - 15.0;      // don't start new work after this


    vector<int> ordBLF = idx;
    stable_sort(ordBLF.begin(), ordBLF.end(), [&](int a, int b) {
        if (ps[a].k != ps[b].k) return ps[a].k > ps[b].k;
        int ma = max(ps[a].minW, ps[a].minH), mb = max(ps[b].minW, ps[b].minH);
        if (ma != mb) return ma > mb;
        return ps[a].id < ps[b].id;
    });
    const vector<int>& ordB2 = BLF2ORD ? baseOrder : ordBLF;

    int swin = max(1, n / 4);
    if (!big && !((big && BIGBLF) || (!big && SMALLBLF))) {
        static double SW1F = envInt("PP_SW1F", 80) / 100.0;
        double sweepBudget = (PHASE2 && S < P2MAXS) ? min(SOFT_END, TL_MS * P2FRAC) : SOFT_END;
        double budget1 = SW1F * max(50.0, sweepBudget - elapsed_ms());
        double swinPred = budget1 / tFB;
        int cap = max(24, (int)swinPred);
        if (cap < swin) swin = cap;
        if (S > 7000) {
            int expL = min(max(1, n / 4), (int)(350000 / max(1LL, S - 3500)));
            static double RTHRESH = envInt("PP_RTHRESH", 200) / 100.0;
            if (swinPred < RTHRESH * expL) big = true; // env too slow for this size: big path wins
        }
        if (getenv("PP_DEBUG")) fprintf(stderr, "tFB=%.2f swin=%d (n/4=%d) big=%d\n", tFB, swin, max(1, n / 4), (int)big);
    }
    bool skipSweep = (big && BIGBLF) || (!big && SMALLBLF);
    int expLIM = big ? min(max(1, n / 4), (int)(350000 / max(1LL, S - 3500))) : 0;
    double avg = 250.0; int cnt = 0;
    vector<pair<long long,int>> sweepRes; // (area, W)
    bool doPhase2 = (!big && PHASE2 && S < P2MAXS);
    double SWFRAC = envInt("PP_SWFRAC", 100) / 100.0; // sweep cap for the no-phase2 small path
    double sweepEnd = doPhase2 ? min(SOFT_END, TL_MS * P2FRAC)
                               : (!big ? min(SOFT_END, TL_MS * SWFRAC) : SOFT_END);
    for (int wi = 0; wi < (skipSweep ? 0 : (int)Ws.size()); wi++) {
        double used = elapsed_ms();
        if (used + avg * 1.15 > sweepEnd) break;
        int W = Ws[wi];
        double t1 = elapsed_ms();
        R r;
        if (big) {
            double reserve = max(30.0, n * 0.006);
            int BIG2 = envInt("PP_BIG2", 0);
            if (BIG2) {
                double half = (SEARCH_END - reserve - t1) * 0.5;
                double panic1 = t1 + half;
                double adapt1 = max(50.0, half - 15.0);
                r = pack(W, baseOrder, rng, false, expLIM, true, adapt1, SEARCH_END, panic1);
                if (r.ok) { crownRepack(r, SEARCH_END); if (better(r, bestR)) bestR = move(r); }
                double t2 = elapsed_ms();
                double panic2 = SEARCH_END - reserve;
                double adapt2 = max(50.0, panic2 - t2 - 15.0);
                r = pack(W + 1, baseOrder, rng, true, expLIM, true, adapt2, SEARCH_END, panic2);
            } else {
                double panicAt = SEARCH_END - reserve;
                double adaptTL = max(50.0, panicAt - t1 - 15.0);
                r = pack(W, baseOrder, rng, false, expLIM, true, adaptTL, SEARCH_END, panicAt);
            }
        } else if (BLF2SWEEP) {
            r = B3 ? pack_blf3(W, ordB2, max(1, n / 4), SEARCH_END, rng, false)
                   : pack_blf2(W, ordB2, max(1, n / 4), SEARCH_END, rng, false);
        } else {
            double panicAt = min(sweepEnd - 10.0, SEARCH_END - 60.0);
            double adaptTL = max(50.0, panicAt - t1 - 10.0);
            r = pack(W, baseOrder, rng, false, swin, true, adaptTL, SEARCH_END, panicAt, 1);
        }
        double dt = elapsed_ms() - t1;
        cnt++; avg = (avg * (cnt - 1) + dt) / cnt;
        if (!r.ok) break; // deadline hit mid-pack
        crownRepack(r, SEARCH_END);
        if (r.ok) sweepRes.push_back({r.A, W});
        if (better(r, bestR)) bestR = move(r);
    }

    if (doPhase2) {
        sort(sweepRes.begin(), sweepRes.end());
        double p2Stop = min(SOFT_END, TL_MS * P2ENDF);
        double avg2 = 20.0; int cnt2 = 0;
        vector<int> p2W;
        {
            unordered_set<int> seenW; seenW.reserve(64);
            for (auto& sr : sweepRes) if (seenW.insert(sr.second).second) p2W.push_back(sr.second);
            static int P2SPAN = envInt("PP_P2SPAN", 0);
            int c0 = sweepRes.empty() ? base : sweepRes[0].second;
            for (int d = 1; d <= P2SPAN; d++)
                for (int sgn = -1; sgn <= 1; sgn += 2) {
                    int W = c0 + sgn * d;
                    if (W >= minW && W <= 4000 && seenW.insert(W).second) p2W.push_back(W);
                }
        }
        for (size_t i = 0; i < p2W.size(); i++) {
            double used = elapsed_ms();
            if (used + avg2 * 1.2 > p2Stop) break;
            int W = p2W[i];
            double t1 = elapsed_ms();
            static int B3WINDIV = envInt("PP_B3WINDIV", 4);
            int b3win = max(1, B3WINDIV > 0 ? n / B3WINDIV : n);
            R r = B3 ? pack_blf3(W, ordB2, b3win, SEARCH_END, rng, false)
                     : pack_blf2(W, ordB2, b3win, SEARCH_END, rng, false);
            double dt = elapsed_ms() - t1;
            cnt2++; avg2 = (avg2 * (cnt2 - 1) + dt) / cnt2;
            if (!r.ok) break;
            crownRepack(r, SEARCH_END);
            if (better(r, bestR)) bestR = move(r);
        }
    }

    {
        vector<int> ilsOrd; int ilsW = 0; (void)ilsW;
        double avgBLF = 10.0; int cntBLF = 0;
        double avgSky = avg;
        vector<int> obuf;
        while (true) {
            double used = elapsed_ms();
            if (used + 5.0 > SOFT_END) break;
            int bw = bestR.W;
            bool doBLF = true;
            if (!big && (rng.nxt() & 3) == 0) doBLF = false;
            if (doBLF) {
                if (used + avgBLF * 1.3 > SOFT_END) {
                    if (!big && used + avgSky * 1.3 <= SOFT_END) doBLF = false; else break;
                }
            } else if (used + avgSky * 1.3 > SOFT_END) {
                if (used + avgBLF * 1.3 <= SOFT_END) doBLF = true; else break;
            }
            static int CAPSHARE = envInt("PP_CAPSHARE", 100);
            static int CAPBIGN = envInt("PP_CAPBIGN", 6000);
            bool capEligible = !big || (n <= CAPBIGN);
            if (doBLF && capEligible && bestR.packW > 0 && bestR.packW <= 64 && (int)(rng.nxt() % 100) < CAPSHARE) {
                int W;
                {
                    static int dwBase = envInt("PP_CAPDW", 2);
                    static int dwWide = envInt("PP_CAPDWW", 6);
                    static int wideP = envInt("PP_CAPWIDEP", 25);
                    int dw = ((int)(rng.nxt() % 100) < wideP) ? dwWide : dwBase;
                    W = bestR.packW + (dw > 0 ? (rng.rint(2 * dw + 1) - dw) : 0);
                    if (W < minW) W = minW;
                    if (W > 64) W = 64;
                }
                long long capA = bestR.A - 1;
                int Hcap = (int)(capA / W);
                if (Hcap < minW || Hcap <= 1) { continue; } // too squat to be plausible
                bool fromBest = !ilsOrd.empty() && (rng.nxt() & 1);
                obuf = fromBest ? ilsOrd : ((cntBLF & 1) ? ordBLF : ordB2);
                int swaps = 1 + rng.rint(12);
                for (int sswap = 0; sswap < swaps; sswap++) {
                    int a = rng.rint(n), b = rng.rint(n);
                    swap(obuf[a], obuf[b]);
                }
                double t1 = elapsed_ms();
                static int CAPWINDIV = envInt("PP_CAPWINDIV", 1);
                R r = pack_capped(W, Hcap, obuf, max(1, CAPWINDIV > 0 ? n / CAPWINDIV : n), SEARCH_END, rng);
                double dt = elapsed_ms() - t1;
                cntBLF++; avgBLF = (avgBLF * (cntBLF - 1) + dt) / cntBLF;
                if (elapsed_ms() > SEARCH_END) break;
                if (r.ok) {
                    crownRepack(r, SEARCH_END);
                    if (better(r, bestR)) { ilsOrd = obuf; bestR = move(r); }
                }
            } else if (doBLF) {
                int W;
                if ((int)(rng.nxt() % 100) < JUMP) {
                    double mult = 0.6 + 0.8 * rng.uni();
                    W = (int)llround(bw * mult);
                } else {
                    W = bw + (rng.rint(9) - 4);
                }
                if (W < minW) W = minW;
                if (W > 64) W = 64;
                if (W < minW) break;
                bool fromBest = !ilsOrd.empty() && (rng.nxt() & 1);
                obuf = fromBest ? ilsOrd : (BLF2 ? ordB2 : ordBLF);
                if (cntBLF > 0) {
                    int swaps = fromBest ? (2 + rng.rint(10)) : max(1, n / 4);
                    for (int sswap = 0; sswap < swaps; sswap++) {
                        int a = rng.rint(n), b = rng.rint(n);
                        swap(obuf[a], obuf[b]);
                    }
                }
                int policy = (int)(rng.nxt() & 1);
                double t1 = elapsed_ms();
                R r = (BLF2 && S < B2RESTS) ? (B3 ? pack_blf3(W, obuf, max(1, n / 4), SEARCH_END, rng, cntBLF > 0)
                                                  : pack_blf2(W, obuf, max(1, n / 4), SEARCH_END, rng, cntBLF > 0))
                           : pack_blf(W, obuf, policy, SEARCH_END);
                double dt = elapsed_ms() - t1;
                cntBLF++; avgBLF = (avgBLF * (cntBLF - 1) + dt) / cntBLF;
                if (getenv("PP_DEBUG")) fprintf(stderr, "BLF W=%d dt=%.1f ok=%d A=%lld best=%lld\n", W, dt, (int)r.ok, r.ok ? r.A : -1, bestR.A);
                if (!r.ok) break;
                crownRepack(r, SEARCH_END);
                if (better(r, bestR)) { ilsOrd = obuf; ilsW = W; bestR = move(r); }
            } else {
                int W = bw + (rng.rint(7) - 3);
                if (W < minW) W = minW;
                obuf = baseOrder;
                int swaps = max(1, n / 4);
                for (int sswap = 0; sswap < swaps; sswap++) {
                    int a = rng.rint(n), b = rng.rint(n);
                    swap(obuf[a], obuf[b]);
                }
                double t1 = elapsed_ms();
                R r = pack(W, obuf, rng, true, max(1, n / 4), false, 0.0, SEARCH_END);
                double dt = elapsed_ms() - t1;
                cnt++; avgSky = (avgSky * 0.7 + dt * 0.3);
                if (!r.ok) break;
                crownRepack(r, SEARCH_END);
                if (better(r, bestR)) bestR = move(r);
            }
        }
    }
    // late restart: try shuffled orderings at multiple widths near the best, using BLF3 for speed
    if (bestR.ok && bestR.packW > 0 && bestR.packW <= 64 && elapsed_ms() < SOFT_END - 50) {
        int widths[3] = {bestR.packW, bestR.packW - 1, bestR.packW + 1};
        for (int widx = 0; widx < 3; widx++) {
            int lateW = widths[widx];
            if (lateW < minW || lateW > 64) continue;
            int lateHcap = max(1, (int)(bestR.A / lateW) - 1);
            if (lateHcap < minW) continue;
            for (int attempt = 0; attempt < 3 && elapsed_ms() < SOFT_END - 10; attempt++) {
                vector<int> randOrd = idx;
                int ns = max(1, n / 3);
                for (int s = 0; s < ns; s++) { int a = rng.rint(n), b = rng.rint(n); swap(randOrd[a], randOrd[b]); }
                // Try BLF3 (cached, fast) for quick iterations
                R r = pack_blf3(lateW, randOrd, max(1, n / 4), SEARCH_END, rng, attempt > 0);
                if (!r.ok) r = pack_capped(lateW, lateHcap, randOrd, max(1, n / 4), SEARCH_END, rng);
                if (r.ok) { crownRepack(r, SEARCH_END); if (better(r, bestR)) bestR = move(r); }
            }
        }
    }
    } // end champion search block (skipped when best-fit produced the result)
    if (getenv("PP_DEBUG")) fprintf(stderr, "t_search_done=%.1f\n", elapsed_ms());
    if (!useBF) crownRepack(bestR, TL_MS + 10.0);
    if (getenv("PP_DEBUG")) fprintf(stderr, "t_crown_done=%.1f\n", elapsed_ms());

    int maxX = -1, maxY = -1;
    for (auto& p : bestR.pl) {
        auto& t = ps[p.idx].t[p.ti];
        for (auto& q : t.c) { int x = p.x + q.first, y = p.y + q.second; if (x > maxX) maxX = x; if (y > maxY) maxY = y; }
    }
    vector<int> mapx(max(0, maxX + 1), -1), mapy(max(0, maxY + 1), -1);
    int Wc = 0, Hc = 0;
    {
        vector<char> ux(maxX + 1, false), uy(maxY + 1, false);
        for (auto& p : bestR.pl) {
            auto& t = ps[p.idx].t[p.ti];
            for (auto& q : t.c) { ux[p.x + q.first] = true; uy[p.y + q.second] = true; }
        }
        for (int x = 0; x <= maxX; x++) if (ux[x]) mapx[x] = Wc++;
        for (int y = 0; y <= maxY; y++) if (uy[y]) mapy[y] = Hc++;
    }
    if (Wc == 0) Wc = 1;
    if (Hc == 0) Hc = 1;

    vector<array<int,4>> ans(n, {0, 0, 0, 0});
    for (auto& p : bestR.pl) {
        auto& t = ps[p.idx].t[p.ti];
        int bx = mapx[p.x];
        int by = mapy[p.y];
        int Xi = bx - t.minx;
        int Yi = by - t.miny;
        int Ri = (4 - (t.r % 4) + 4) % 4;
        int Fi = t.f;
        ans[p.idx] = {Xi, Yi, Ri, Fi};
    }
    string out;
    out.reserve(16 * (n + 1));
    out += to_string(Wc); out += ' '; out += to_string(Hc); out += '\n';
    for (int i = 0; i < n; i++) {
        out += to_string(ans[i][0]); out += ' ';
        out += to_string(ans[i][1]); out += ' ';
        out += to_string(ans[i][2]); out += ' ';
        out += to_string(ans[i][3]); out += '\n';
    }
    fwrite(out.data(), 1, out.size(), stdout);
    if (getenv("PP_DEBUG")) fprintf(stderr, "t_output_done=%.1f\n", elapsed_ms());
    fflush(stdout);
    _Exit(0); // skip destructor teardown of large heaps
}
