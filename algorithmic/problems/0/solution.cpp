// Endgame prototype for small n: rectangle ladder + cell-driven DFS (exact cover with blanks).
// Reads instance on stdin, writes solution on stdout. Env: POLYPACK_TL (ms), EG_NL0 (initial node limit).
#include <bits/stdc++.h>
using namespace std;

static chrono::steady_clock::time_point T0;
static double TL = 1890.0;
static inline double el() { return chrono::duration<double, milli>(chrono::steady_clock::now() - T0).count(); }

struct Ori { int k, w, h; int dx[10], dy[10]; int r, f, minx, miny, ax, ay; int id; };
struct Grp { int k; vector<Ori> os; vector<int> members; };

int n; long long S = 0;
vector<vector<pair<int,int>>> raw;   // original piece cells
vector<Grp> G;
vector<vector<Ori>> memOri;   // per piece: orientations in canonical (key-sorted) order

static pair<int,int> rot_ccw(pair<int,int> p, int r) {
    for (int i = 0; i < r; i++) p = {-p.second, p.first};
    return p;
}

// build orientations for a cell list; returns list of (key, Ori)
static vector<pair<string, Ori>> build_oris(const vector<pair<int,int>>& cells) {
    vector<pair<string, Ori>> out;
    set<string> seen;
    for (int f = 0; f < 2; f++) for (int r = 0; r < 4; r++) {
        vector<pair<int,int>> v;
        for (auto c : cells) { if (f) c.first = -c.first; v.push_back(rot_ccw(c, r)); }
        int mnx = INT_MAX, mny = INT_MAX, mxx = INT_MIN, mxy = INT_MIN;
        for (auto& c : v) { mnx = min(mnx, c.first); mny = min(mny, c.second); mxx = max(mxx, c.first); mxy = max(mxy, c.second); }
        for (auto& c : v) { c = {c.first - mnx, c.second - mny}; }
        sort(v.begin(), v.end(), [](auto& a, auto& b){ return a.second != b.second ? a.second < b.second : a.first < b.first; }); // row-major (y, then x)
        string key;
        for (auto& c : v) { key += to_string(c.first); key += ','; key += to_string(c.second); key += ';'; }
        if (!seen.insert(key).second) continue;
        Ori o; o.k = v.size(); o.w = mxx - mnx + 1; o.h = mxy - mny + 1;
        o.ax = v[0].first; o.ay = v[0].second;                 // anchor = row-major minimal cell
        for (int i = 0; i < o.k; i++) { o.dx[i] = v[i].first - o.ax; o.dy[i] = v[i].second - o.ay; }
        o.r = r; o.f = f; o.minx = mnx; o.miny = mny;
        out.push_back({key, o});
    }
    return out;
}

// hole census for the floor bound
static int hole_cells(const vector<pair<int,int>>& cells) {
    set<pair<int,int>> cs(cells.begin(), cells.end());
    int mnx = INT_MAX, mny = INT_MAX, mxx = INT_MIN, mxy = INT_MIN;
    for (auto& c : cs) { mnx = min(mnx,c.first); mny = min(mny,c.second); mxx = max(mxx,c.first); mxy = max(mxy,c.second); }
    set<pair<int,int>> seen; deque<pair<int,int>> dq;
    dq.push_back({mnx-1, mny-1}); seen.insert(dq[0]);
    while (!dq.empty()) {
        auto [x, y] = dq.front(); dq.pop_front();
        for (auto [dx, dy] : {pair{1,0},{-1,0},{0,1},{0,-1}}) {
            int nx = x+dx, ny = y+dy;
            if (nx < mnx-1 || nx > mxx+1 || ny < mny-1 || ny > mxy+1) continue;
            if (cs.count({nx,ny}) || seen.count({nx,ny})) continue;
            seen.insert({nx,ny}); dq.push_back({nx,ny});
        }
    }
    int h = 0;
    for (int x = mnx; x <= mxx; x++) for (int y = mny; y <= mxy; y++)
        if (!cs.count({x,y}) && !seen.count({x,y})) h++;
    return h;
}

// ---------------- DFS ----------------
struct Placement { int grp; Ori o; int px, py; };
static int W, H; static long long A;
static vector<uint64_t> grid;
static vector<int> cnt;               // remaining per group
static vector<int> ord;               // group visit order
static vector<Placement> stk, best;
static int slack_;                    // blanks remaining
static long long nodes, nodeLimit;
static double attemptEnd;
static bool aborted;

static bool dfs(int pos) {
    while (pos < A && ((grid[pos / W] >> (pos % W)) & 1ULL)) pos++;
    if (pos >= A) { best = stk; return true; }
    if (++nodes > nodeLimit) { aborted = true; return false; }
    if ((nodes & 2047) == 0 && el() > attemptEnd) { aborted = true; return false; }
    int x = pos % W, y = pos / W;
    for (int gi : ord) {
        if (!cnt[gi]) continue;
        auto& g = G[gi];
        for (int oi = 0; oi < (int)g.os.size(); oi++) {
            auto& o = g.os[oi];
            if (y + o.h - o.ay > H) continue;
            bool ok = true;
            for (int i = 1; i < o.k && ok; i++) {
                int xx = x + o.dx[i], yy = y + o.dy[i];
                if (xx < 0 || xx >= W) { ok = false; break; }
                if ((grid[yy] >> xx) & 1ULL) ok = false;
            }
            if (!ok) continue;
            for (int i = 0; i < o.k; i++) grid[y + o.dy[i]] |= 1ULL << (x + o.dx[i]);
            cnt[gi]--; stk.push_back({gi, o, x, y});
            if (dfs(pos + 1)) return true;
            stk.pop_back(); cnt[gi]++;
            for (int i = 0; i < o.k; i++) grid[y + o.dy[i]] &= ~(1ULL << (x + o.dx[i]));
            if (aborted) return false;
        }
    }
    if (slack_ > 0) {                  // blank this cell, tried last
        slack_--; grid[y] |= 1ULL << x;
        if (dfs(pos + 1)) return true;
        grid[y] &= ~(1ULL << x); slack_++;
    }
    return false;
}


static mt19937_64* g_rng;
static int remPieces;
static long long dfsBudget;

// randomized greedy construction; switches to exact DFS for the last few pieces
static bool construct(int egThresh) {
    int pos = 0;
    struct Cand { int gi, oi; };
    static vector<Cand> cs;
    while (pos < A) {
        if ((grid[pos / W] >> (pos % W)) & 1ULL) { pos++; continue; }
        if (remPieces <= egThresh) {           // exact endgame on the tail
            nodes = 0; nodeLimit = dfsBudget; aborted = false;
            return dfs(pos);
        }
        int x = pos % W, y = pos / W;
        cs.clear();
        int bestK = 0;
        for (int gi : ord) {
            if (!cnt[gi]) continue;
            auto& g = G[gi];
            if (g.k < bestK) break;            // ord is size-desc: nothing better follows
            for (int oi = 0; oi < (int)g.os.size(); oi++) {
                auto& o = g.os[oi];
                if (y + o.h > H) continue;
                bool ok = true;
                for (int i = 1; i < o.k && ok; i++) {
                    int xx = x + o.dx[i], yy = y + o.dy[i];
                    if (xx < 0 || xx >= W || ((grid[yy] >> xx) & 1ULL)) ok = false;
                }
                if (!ok) continue;
                (void)bestK;
                cs.push_back({gi, oi});
            }
        }
        if (cs.empty()) {
            if (slack_ <= 0) return false;
            slack_--; grid[y] |= 1ULL << x; pos++;
            continue;
        }
        // among largest fitting pieces, pick max contact (occupied/wall neighbors), random tie
        int bi = -1; int bc = -1;
        for (auto& c : cs) {
            auto& o = G[c.gi].os[c.oi];
            int ct = 0;
            for (int i = 0; i < o.k; i++) {
                int cx = x + o.dx[i], cy = y + o.dy[i];
                const int DX[4] = {1,-1,0,0}, DY[4] = {0,0,1,-1};
                for (int dd = 0; dd < 4; dd++) {
                    int nx = cx + DX[dd], ny = cy + DY[dd];
                    if (nx < 0 || nx >= W || ny < 0) { ct++; continue; }
                    if (ny >= H) continue;
                    if ((grid[ny] >> nx) & 1ULL) ct++;
                }
            }
            static int WCT = getenv("EG_WCT") ? atoi(getenv("EG_WCT")) : 8;
            static int WK  = getenv("EG_WK")  ? atoi(getenv("EG_WK"))  : 2;
            ct = ct * WCT + G[c.gi].k * WK + (int)((*g_rng)() & 7);
            if (ct > bc) { bc = ct; bi = (int)(&c - &cs[0]); }
        }
        auto& o = G[cs[bi].gi].os[cs[bi].oi];
        for (int i = 0; i < o.k; i++) grid[y + o.dy[i]] |= 1ULL << (x + o.dx[i]);
        cnt[cs[bi].gi]--; remPieces--;
        stk.push_back({cs[bi].gi, o, x, y});
        pos++;
    }
    best = stk;
    return true;
}


static bool attemptLNS_impl(const vector<Placement>& baseSol, int w, int h, long long a,
                            int ruinPct, mt19937_64& rng, int egThresh, long long nl) {
    W = w; H = h; A = a;
    grid.assign(H, 0ULL);
    cnt.assign(G.size(), 0);
    for (size_t g = 0; g < G.size(); g++) cnt[g] = (int)G[g].members.size();
    stk.clear(); remPieces = (int)0; slack_ = (int)(A - S);
    remPieces = 0; for (size_t g = 0; g < G.size(); g++) remPieces += cnt[g];
    for (auto& p : baseSol) {
        auto& o = p.o;
        bool fits = (p.py + o.h <= H);
        for (int i = 0; i < o.k && fits; i++) {
            int cx = p.px + o.dx[i];
            if (cx < 0 || cx >= W) fits = false;
        }
        if (!fits) continue;
        if ((int)(rng() % 100) < ruinPct) continue;
        for (int i = 0; i < o.k; i++) grid[p.py + o.dy[i]] |= 1ULL << (p.px + o.dx[i]);
        cnt[p.grp]--; remPieces--;
        stk.push_back(p);
    }
    g_rng = &rng; dfsBudget = nl; attemptEnd = TL - 60.0;
    return construct(egThresh);
}

int main() {
    T0 = chrono::steady_clock::now();
    if (const char* e = getenv("POLYPACK_TL")) { double v = atof(e); if (v > 50 && v < 20000) TL = v; }
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if (!(cin >> n)) return 0;
    raw.resize(n);
    for (int i = 0; i < n; i++) {
        int k; cin >> k; raw[i].resize(k);
        for (auto& c : raw[i]) cin >> c.first >> c.second;
        S += k;
    }
    // group by canonical form (min orientation key)
    unordered_map<string, int> gid;
    vector<vector<pair<string, Ori>>> oris(n);
    int minW = 1;
    long long holeDef = 0; int m1 = 0, d1 = 0;
    for (int i = 0; i < n; i++) {
        oris[i] = build_oris(raw[i]);
        string canon = oris[i][0].first;
        for (auto& po : oris[i]) canon = min(canon, po.first);
        sort(oris[i].begin(), oris[i].end(), [](auto& a, auto& b){ return a.first < b.first; });
        for (size_t j = 0; j < oris[i].size(); j++) oris[i][j].second.id = (int)j;
        memOri.resize(n);
        for (auto& po : oris[i]) memOri[i].push_back(po.second);
        auto it = gid.find(canon);
        int g;
        if (it == gid.end()) {
            g = G.size(); gid[canon] = g;
            G.push_back({});
            G[g].k = raw[i].size();
            G[g].os = memOri[i];   // canonical order; ids consistent across members
        } else g = it->second;
        G[g].members.push_back(i);
        int mw = INT_MAX; for (auto& po : oris[i]) mw = min(mw, po.second.w);
        minW = max(minW, mw);
        holeDef += hole_cells(raw[i]);
        if ((int)raw[i].size() == 1) m1++;
        if ((int)raw[i].size() == 2) d1++;
    }
    long long floorA = S + max(0LL, holeDef - m1 - 2LL*d1);
    // candidate rectangles sorted by (A, H, W)
    vector<array<long long,3>> cands; // A, H, W
    long long Amax = (long long)ceil(S * 1.14);
    for (int w = minW; w <= 64; w++) {
        for (long long h = max<long long>(1, (floorA + w - 1) / w); w * h <= Amax; h++) {
            // every piece must fit
            cands.push_back({w * h, h, w});
        }
    }
    sort(cands.begin(), cands.end());
    // descending frontier walk + escalating hammer at the wall
    long long bestA = LLONG_MAX; int bestW = 0, bestH = 0;
    vector<Placement> sol;
    unsigned long long seed = 0x9e3779b97f4a7c15ULL ^ (S << 1) ^ n;
    mt19937_64 rng(seed);
    ord.resize(G.size()); iota(ord.begin(), ord.end(), 0);
    double SOFT = TL - 60.0;
    // group candidates by area, descending
    sort(cands.begin(), cands.end(), [](auto& a, auto& b){
        if (a[0] != b[0]) return a[0] > b[0];
        if (a[1] != b[1]) return a[1] < b[1];
        return a[2] < b[2]; });
    vector<vector<array<long long,3>>> levels;
    for (auto& c : cands) {
        if (levels.empty() || levels.back()[0][0] != c[0]) levels.push_back({});
        levels.back().push_back(c);
    }
    long long NL0 = 20000;
    if (const char* e = getenv("EG_NL0")) NL0 = atoll(e);
    auto attempt = [&](int w, int h, long long a, long long nl) -> bool {
        W = w; H = h; A = a;
        grid.assign(H, 0ULL);
        cnt.assign(G.size(), 0);
        for (size_t g = 0; g < G.size(); g++) cnt[g] = (int)G[g].members.size();
        shuffle(ord.begin(), ord.end(), rng);
        stable_sort(ord.begin(), ord.end(), [](int a2, int b2){ return G[a2].k > G[b2].k; });
        for (auto& g : G) shuffle(g.os.begin(), g.os.end(), rng);
        slack_ = (int)(A - S);
        stk.clear(); attemptEnd = SOFT;
        g_rng = &rng; remPieces = n; dfsBudget = nl;
        static int egThresh = getenv("EG_TAIL") ? atoi(getenv("EG_TAIL")) : 10;
        bool okk = construct(egThresh);
        if (getenv("EG_DBG") && ++*(long long*)&nodeLimit != 0 && el() < 1e18) fprintf(stderr, "  attempt A=%lld %dx%d nl=%lld -> ok=%d nodes=%lld aborted=%d t=%.0f\n", A, W, H, nl, (int)okk, nodes, (int)aborted, el());
        if (okk) { bestA = A; bestW = W; bestH = H; sol = best; return true; }
        return false;
    };
    int DBG = getenv("EG_DBG") ? 1 : 0;
    if (DBG) fprintf(stderr, "G=%zu levels=%zu floorA=%lld S=%lld minW=%d Amax=%lld\n", G.size(), levels.size(), floorA, S, minW, Amax);
    long long atCnt = 0;
    size_t li = 0;
    // phase 1: fast descent (one cheap attempt per level)
    while (li < levels.size() && el() < SOFT * 0.5) {
        if (levels[li][0][0] >= bestA) { li++; continue; }
        bool ok = false;
        for (auto& c : levels[li]) if (attempt((int)c[2], (int)c[1], c[0], NL0)) { ok = true; break; }
        li++;
        if (!ok) break;   // hit the wall
    }
    // phase 2: hammer the frontier: window of next few levels below bestA, escalating budgets
    while (el() < SOFT) {
        size_t f = 0;
        while (f < levels.size() && levels[f][0][0] >= bestA) f++;
        if (f >= levels.size()) break;
        bool improved = false;
        long long nl = NL0;
        while (el() < SOFT && !improved) {
            for (size_t wnd = f; wnd < min(f + 3, levels.size()) && !improved; wnd++)
                for (auto& c : levels[wnd]) {
                    if (el() > SOFT) break;
                    static int egT2 = getenv("EG_TAIL") ? atoi(getenv("EG_TAIL")) : 10;
                    bool ok2 = false;
                    if (!sol.empty() && (rng() & 1)) {
                        vector<Placement> baseCopy = sol;
                        ok2 = attemptLNS_impl(baseCopy, (int)c[2], (int)c[1], c[0], 12 + (int)(rng()%15), rng, egT2, nl);
                        if (ok2) { bestA = c[0]; bestW = (int)c[2]; bestH = (int)c[1]; sol = best; }
                    } else ok2 = attempt((int)c[2], (int)c[1], c[0], nl);
                    if (ok2) { improved = true; break; }
                }
            nl = min(nl * 3, (long long)3e9);
        }
        if (!improved) break;
    }
    if (bestW == 0) { // emergency: vertical stack (should never happen)
        printf("10 %lld\n", S); return 0;
    }
    if (getenv("EG_CHECK")) { // self-check in internal frame
        vector<vector<int>> occ2(bestH, vector<int>(bestW, -1));
        long long placed = 0; bool bad = false;
        for (size_t si = 0; si < sol.size(); si++) {
            auto& p = sol[si]; auto& o = p.o;
            for (int i = 0; i < o.k; i++) {
                int cx = p.px + o.dx[i], cy = p.py + o.dy[i];
                if (cx < 0 || cx >= bestW || cy < 0 || cy >= bestH) { fprintf(stderr, "OOB sol[%zu]\n", si); bad = true; continue; }
                if (occ2[cy][cx] >= 0) { fprintf(stderr, "INTERNAL OVERLAP sol[%zu] vs sol[%d] at (%d,%d)\n", si, occ2[cy][cx], cx, cy); bad = true; }
                occ2[cy][cx] = (int)si; placed++;
            }
        }
        fprintf(stderr, "self-check: pieces=%zu cells=%lld S=%lld bad=%d bestA=%lld\n", sol.size(), placed, S, (int)bad, bestA);
    }
    // trim empty rows/cols (safe: pieces span contiguous rows/cols)
    vector<char> ux(bestW, 0), uy(bestH, 0);
    vector<array<int,4>> cellsOf; // decode occupancy from sol
    vector<array<int,4>> ans(n);
    {
        vector<vector<char>> occ(bestH, vector<char>(bestW, 0));
        for (auto& p : sol) {
            auto& o = p.o;
            for (int i = 0; i < o.k; i++) { occ[p.py + o.dy[i]][p.px + o.dx[i]] = 1; }
        }
        for (int y = 0; y < bestH; y++) for (int x = 0; x < bestW; x++) if (occ[y][x]) { ux[x] = 1; uy[y] = 1; }
    }
    vector<int> mx(bestW), my(bestH); int Wc = 0, Hc = 0;
    for (int x = 0; x < bestW; x++) { mx[x] = Wc; if (ux[x]) Wc++; }
    for (int y = 0; y < bestH; y++) { my[y] = Hc; if (uy[y]) Hc++; }
    if (!Wc) Wc = 1; if (!Hc) Hc = 1;
    vector<int> used(G.size(), 0);
    for (auto& p : sol) {
        auto& g = G[p.grp];
        int pi = g.members[used[p.grp]++];
        auto& o = memOri[pi][p.o.id];   // same geometry, member's own (r,f,minx,miny)
        int ox = mx[p.px] - o.ax? 0:0; // placeholder, computed below properly
        (void)ox;
        int px = mx[p.px], py = my[p.py];
        int X = (px - o.ax) - o.minx;
        int Y = (py - o.ay) - o.miny;
        int R = (4 - (o.r % 4)) % 4;
        ans[pi] = {X, Y, R, o.f};
    }
    string out = to_string(Wc) + " " + to_string(Hc) + "\n";
    for (int i = 0; i < n; i++)
        out += to_string(ans[i][0]) + " " + to_string(ans[i][1]) + " " + to_string(ans[i][2]) + " " + to_string(ans[i][3]) + "\n";
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
