// Treasure Packing — exact 2-constraint bounded knapsack (12 types) via branch & bound
// with Lagrangian (LP-equivalent) bounds, multi-greedy + hill-climb primal, hard deadline.
// I/O: JSON in (fixed shape), JSON out (same keys).
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef double dd;

static chrono::steady_clock::time_point T0;
static double TL_MS = 750.0;
static inline double elapsed_ms() {
    return chrono::duration<double, milli>(chrono::steady_clock::now() - T0).count();
}

static const ll MAXM = 20000000LL, MAXL = 25000000LL;

int n = 12;
vector<string> keys;
vector<ll> Q, V, Mw, Lv;

struct RNG {
    unsigned long long s;
    RNG(unsigned long long x) : s(x ? x : 0x9E3779B97F4A7C15ULL) {}
    inline unsigned long long nxt() { s ^= s << 7; s ^= s >> 9; return s; }
    inline int rint(int m) { return (int)(nxt() % (unsigned)m); }
};

// ---------- JSON input ----------
static void parseInput() {
    string all, line;
    {
        std::ios::sync_with_stdio(false);
        std::ostringstream ss; ss << cin.rdbuf();
        all = ss.str();
    }
    size_t i = 0;
    auto skipws = [&]() { while (i < all.size() && (isspace((unsigned char)all[i]) || all[i] == ',' || all[i] == ':')) i++; };
    // expect '{'
    while (i < all.size() && all[i] != '{') i++;
    i++;
    for (int k = 0; k < 12; k++) {
        // key
        while (i < all.size() && all[i] != '"') i++;
        i++;
        string key;
        while (i < all.size() && all[i] != '"') key += all[i++];
        i++;
        // '['
        while (i < all.size() && all[i] != '[') i++;
        i++;
        ll vals[4];
        for (int t = 0; t < 4; t++) {
            skipws();
            bool neg = false;
            if (i < all.size() && all[i] == '-') { neg = true; i++; }
            ll x = 0;
            while (i < all.size() && isdigit((unsigned char)all[i])) x = x * 10 + (all[i++] - '0');
            vals[t] = neg ? -x : x;
        }
        while (i < all.size() && all[i] != ']') i++;
        i++;
        keys.push_back(key);
        Q.push_back(vals[0]); V.push_back(vals[1]); Mw.push_back(vals[2]); Lv.push_back(vals[3]);
    }
}

// ---------- feasibility & value ----------
static inline bool feasible(const array<ll,12>& x) {
    ll m = 0, l = 0;
    for (int i = 0; i < n; i++) { m += x[i] * Mw[i]; l += x[i] * Lv[i]; }
    return m <= MAXM && l <= MAXL;
}
static inline ll valueOf(const array<ll,12>& x) {
    ll v = 0;
    for (int i = 0; i < n; i++) v += x[i] * V[i];
    return v;
}

array<ll,12> bestX;
ll bestVal = -1;
static void tryImprove(const array<ll,12>& x) {
    ll v = valueOf(x);
    if (v > bestVal && feasible(x)) { bestVal = v; bestX = x; }
}

// greedy fill by given score ordering into (x) starting from base usage
static void greedyFill(array<ll,12>& x, const vector<int>& order) {
    ll m = 0, l = 0;
    for (int i = 0; i < n; i++) { m += x[i] * Mw[i]; l += x[i] * Lv[i]; }
    for (int idx : order) {
        ll capm = (MAXM - m) / Mw[idx];
        ll capl = (MAXL - l) / Lv[idx];
        ll add = min(Q[idx] - x[idx], min(capm, capl));
        if (add > 0) { x[idx] += add; m += add * Mw[idx]; l += add * Lv[idx]; }
    }
}

// hill climb: single-type increments after random decrements
static void hillClimb(array<ll,12>& x, RNG& rng, int iters, double deadline) {
    ll m = 0, l = 0, v = 0;
    for (int i = 0; i < n; i++) { m += x[i] * Mw[i]; l += x[i] * Lv[i]; v += x[i] * V[i]; }
    array<ll,12> cur = x;
    ll curV = v; ll curM = m, curL = l;
    array<ll,12> loc = cur; ll locV = curV;
    for (int it = 0; it < iters; it++) {
        if ((it & 255) == 0 && elapsed_ms() > deadline) break;
        // random move: remove d from i, then fill greedily by random density
        int i = rng.rint(n);
        ll d = cur[i] > 0 ? 1 + (ll)rng.rint((int)min<ll>(cur[i], 4)) : 0;
        array<ll,12> nx = cur;
        ll nm = curM, nl = curL, nv = curV;
        if (d > 0) { nx[i] -= d; nm -= d * Mw[i]; nl -= d * Lv[i]; nv -= d * V[i]; }
        // fill: random order among top-value-density types
        for (int t = 0; t < 4; t++) {
            int j = rng.rint(n);
            ll add = min(Q[j] - nx[j], min((MAXM - nm) / Mw[j], (MAXL - nl) / Lv[j]));
            if (add > 0) { nx[j] += add; nm += add * Mw[j]; nl += add * Lv[j]; nv += add * V[j]; }
        }
        if (nv >= curV) { cur = nx; curV = nv; curM = nm; curL = nl; if (curV > locV) { loc = cur; locV = curV; } }
    }
    x = loc;
}

// ---------- Lagrangian bound ----------
// Dualize volume: g(la) = V0 + la*Lrem + max{ sum (V_i - la*Lv_i) y_i : sum Mw_i y_i <= Mrem, 0<=y_i<=r_i } (fractional)
// Returns pair(bound, fractional item index at optimum or -1) and the greedy y (continuous).
struct BoundResult { dd bound; int fracIdx; dd fracVal; array<dd,12> y; dd lastLa = 0; };
static BoundResult evalDualVol(dd la, ll Mrem, ll Lrem, ll V0, const array<ll,12>& r) {
    // sort by density (v - la*l)/m desc
    array<pair<dd,int>,12> ord;
    int cnt = 0;
    for (int i = 0; i < n; i++) {
        dd red = (dd)V[i] - la * (dd)Lv[i];
        if (red > 1e-12 && r[i] > 0) ord[cnt++] = { red / (dd)Mw[i], i };
    }
    for (int a2 = 1; a2 < cnt; a2++) { // insertion sort (n<=12)
        auto key = ord[a2]; int b2 = a2 - 1;
        while (b2 >= 0 && ord[b2].first < key.first) { ord[b2 + 1] = ord[b2]; b2--; }
        ord[b2 + 1] = key;
    }
    dd rem = (dd)Mrem;
    dd val = (dd)V0 + la * (dd)Lrem;
    int fracIdx = -1; dd fracVal = 0;
    array<dd,12> y{};
    for (int t = 0; t < cnt; t++) {
        int i = ord[t].second;
        dd red = (dd)V[i] - la * (dd)Lv[i];
        dd take = min((dd)r[i], rem / (dd)Mw[i]);
        if (take <= 1e-12) continue;
        y[i] = take;
        val += red * take;
        rem -= take * (dd)Mw[i];
        if (take < (dd)r[i] - 1e-9) { fracIdx = i; fracVal = take; break; }
        if (rem < 1e-9) { /* continue to check exact boundary */ }
    }
    return { val, fracIdx, fracVal, y };
}
// symmetric: dualize mass
static BoundResult evalDualMass(dd mu, ll Mrem, ll Lrem, ll V0, const array<ll,12>& r) {
    array<pair<dd,int>,12> ord;
    int cnt = 0;
    for (int i = 0; i < n; i++) {
        dd red = (dd)V[i] - mu * (dd)Mw[i];
        if (red > 1e-12 && r[i] > 0) ord[cnt++] = { red / (dd)Lv[i], i };
    }
    for (int a2 = 1; a2 < cnt; a2++) { // insertion sort (n<=12)
        auto key = ord[a2]; int b2 = a2 - 1;
        while (b2 >= 0 && ord[b2].first < key.first) { ord[b2 + 1] = ord[b2]; b2--; }
        ord[b2 + 1] = key;
    }
    dd rem = (dd)Lrem;
    dd val = (dd)V0 + mu * (dd)Mrem;
    int fracIdx = -1; dd fracVal = 0;
    array<dd,12> y{};
    for (int t = 0; t < cnt; t++) {
        int i = ord[t].second;
        dd red = (dd)V[i] - mu * (dd)Mw[i];
        dd take = min((dd)r[i], rem / (dd)Lv[i]);
        if (take <= 1e-12) continue;
        y[i] = take;
        val += red * take;
        rem -= take * (dd)Lv[i];
        if (take < (dd)r[i] - 1e-9) { fracIdx = i; fracVal = take; break; }
    }
    return { val, fracIdx, fracVal, y };
}

// minimize convex g over la in [lo0, hi0] by ternary search
template <typename F>
static BoundResult minimizeDual(F eval, dd lo0, dd hi0, int iters) {
    dd lo = lo0, hi = hi0;
    BoundResult best = eval(lo);
    BoundResult rHi = eval(hi);
    if (rHi.bound < best.bound) best = rHi;
    for (int it = 0; it < iters; it++) {
        dd a = lo + (hi - lo) / 3.0, b = hi - (hi - lo) / 3.0;
        BoundResult ra = eval(a), rb = eval(b);
        if (ra.bound < rb.bound) hi = b; else lo = a;
        if (ra.bound < best.bound) best = ra;
        if (rb.bound < best.bound) best = rb;
    }
    BoundResult rm = eval((lo + hi) / 2);
    if (rm.bound < best.bound) best = rm;
    best.lastLa = (lo + hi) / 2;
    return best;
}

// ---------- B&B ----------
long long nodesExplored = 0;
bool provedOptimal = true;

struct Node {
    array<ll,12> lo, hi;
    ll ub;
    int depth;
    int bi; ll t; // stored branch decision
    bool operator<(const Node& o) const { return ub < o.ub; } // max-heap by bound
};
static priority_queue<Node> heap;

// evaluate node: compute bound, update incumbent, return (ub, branch info) via out params
static bool evalNode(const array<ll,12>& lo, const array<ll,12>& hi, int depth,
                     ll& ubOut, int& biOut, ll& tOut, dd& fracPartOut) {
    if (elapsed_ms() > TL_MS) { provedOptimal = false; return false; }
    nodesExplored++;
    // committed part
    ll V0 = 0, Mused = 0, Lused = 0;
    for (int i = 0; i < n; i++) { V0 += lo[i] * V[i]; Mused += lo[i] * Mw[i]; Lused += lo[i] * Lv[i]; }
    if (Mused > MAXM || Lused > MAXL) return false;
    ll Mrem = MAXM - Mused, Lrem = MAXL - Lused;
    array<ll,12> r;
    for (int i = 0; i < n; i++) {
        ll cap = min(hi[i] - lo[i], min(Mrem / Mw[i], Lrem / Lv[i]));
        r[i] = max(0LL, cap);
    }
    // quick incumbent from committed + greedy
    {
        array<ll,12> x = lo;
        // fill by combined density
        vector<int> ord(n); iota(ord.begin(), ord.end(), 0);
        sort(ord.begin(), ord.end(), [&](int a, int b) {
            dd da = (dd)V[a] / ((dd)Mw[a] / MAXM + (dd)Lv[a] / MAXL);
            dd db = (dd)V[b] / ((dd)Mw[b] / MAXM + (dd)Lv[b] / MAXL);
            return da > db;
        });
        // respect hi bounds during fill
        ll m = Mused, l = Lused;
        for (int idx : ord) {
            ll add = min(hi[idx] - x[idx], min((MAXM - m) / Mw[idx], (MAXL - l) / Lv[idx]));
            if (add > 0) { x[idx] += add; m += add * Mw[idx]; l += add * Lv[idx]; }
        }
        tryImprove(x);
    }
    // bounds
    dd laMax = 0, muMax = 0;
    for (int i = 0; i < n; i++) {
        laMax = max(laMax, (dd)V[i] / (dd)Lv[i]);
        muMax = max(muMax, (dd)V[i] / (dd)Mw[i]);
    }
    int iters = depth == 0 ? 96 : 72;
    BoundResult b1 = minimizeDual([&](dd la) { return evalDualVol(la, Mrem, Lrem, V0, r); }, 0.0, laMax, iters);
    BoundResult b2 = minimizeDual([&](dd mu) { return evalDualMass(mu, Mrem, Lrem, V0, r); }, 0.0, muMax, iters);
    BoundResult& bb = (b1.bound <= b2.bound) ? b1 : b2;
    ll ub = (ll)floor(min(b1.bound, b2.bound) + 1e-6);
    if (ub <= bestVal) return false;
    // branching variable: fractional item; if none, pick the one with largest remaining range among used
    int bi = bb.fracIdx;
    dd bval = bb.fracVal;
    if (bi < 0) {
        // integral LP: the greedy y is integral; candidate = lo + y
        array<ll,12> x = lo;
        for (int i = 0; i < n; i++) x[i] += (ll)llround(bb.y[i]);
        for (int i = 0; i < n; i++) x[i] = min(x[i], hi[i]);
        tryImprove(x);
        // choose any var with r>0 to branch on (largest reduced value involvement)
        for (int i = 0; i < n; i++) if (r[i] > 0 && bb.y[i] > 0.5) { bi = i; bval = bb.y[i] / 2; break; }
        if (bi < 0) { ubOut = ub; biOut = -1; return false; } // fully fixed / leaf
    }
    ll t = lo[bi] + (ll)floor(bval);
    t = max(lo[bi], min(t, hi[bi] - 1));
    ubOut = ub; biOut = bi; tOut = t; fracPartOut = bval - floor(bval);
    return true;
}

static void bestFirst(array<ll,12> lo0, array<ll,12> hi0) {
    ll ub; int bi; ll t; dd fp;
    if (!evalNode(lo0, hi0, 0, ub, bi, t, fp)) return;
    heap.push({lo0, hi0, ub, 0, bi, t});
    while (!heap.empty()) {
        if (elapsed_ms() > TL_MS) { provedOptimal = false; break; }
        Node nd = heap.top(); heap.pop();
        if (nd.ub <= bestVal) continue; // pruned by newer incumbent
        array<ll,12> hi1 = nd.hi; hi1[nd.bi] = nd.t;
        array<ll,12> lo2 = nd.lo; lo2[nd.bi] = nd.t + 1;
        for (int c = 0; c < 2; c++) {
            const array<ll,12>& clo = (c == 0) ? nd.lo : lo2;
            const array<ll,12>& chi = (c == 0) ? hi1 : nd.hi;
            ll cub; int cbi; ll ct; dd cfp;
            if (!evalNode(clo, chi, nd.depth + 1, cub, cbi, ct, cfp)) continue;
            if (cub > bestVal && cbi >= 0) heap.push({clo, chi, cub, nd.depth + 1, cbi, ct});
        }
        if ((int)heap.size() > 1500000) { provedOptimal = false; break; } // memory guard
    }
}

int main() {
    T0 = chrono::steady_clock::now();
    if (const char* e = getenv("TP_TL")) { double v = atof(e); if (v > 20 && v < 5000) TL_MS = v; }
    parseInput();

    // ---- primal heuristics ----
    bestX.fill(0);
    bestVal = 0;
    {
        vector<vector<int>> orders;
        vector<int> base(n); iota(base.begin(), base.end(), 0);
        auto mk = [&](function<dd(int)> f) {
            vector<int> o = base;
            sort(o.begin(), o.end(), [&](int a, int b) { return f(a) > f(b); });
            orders.push_back(o);
        };
        mk([&](int i) { return (dd)V[i] / (dd)Mw[i]; });
        mk([&](int i) { return (dd)V[i] / (dd)Lv[i]; });
        mk([&](int i) { return (dd)V[i] / ((dd)Mw[i] / MAXM + (dd)Lv[i] / MAXL); });
        mk([&](int i) { return (dd)V[i] / max((dd)Mw[i] / MAXM, (dd)Lv[i] / MAXL); });
        mk([&](int i) { return (dd)V[i]; });
        for (auto& o : orders) {
            array<ll,12> x{}; x.fill(0);
            greedyFill(x, o);
            tryImprove(x);
        }
        RNG rng(0xABCDEF12345ULL);
        array<ll,12> x = bestX;
        hillClimb(x, rng, 20000, TL_MS * 0.25);
        tryImprove(x);
    }

    // ---- exact search ----
    array<ll,12> lo{}, hi{};
    lo.fill(0);
    for (int i = 0; i < n; i++) hi[i] = min(Q[i], min(MAXM / Mw[i], MAXL / Lv[i]));
    bestFirst(lo, hi);

    // ---- output ----
    string out = "{\n";
    for (int i = 0; i < n; i++) {
        out += " \"" + keys[i] + "\": " + to_string(bestX[i]);
        out += (i + 1 < n) ? ",\n" : "\n";
    }
    out += "}\n";
    fwrite(out.data(), 1, out.size(), stdout);
    if (getenv("TP_DEBUG")) fprintf(stderr, "val=%lld nodes=%lld optimal=%d t=%.1f\n", bestVal, nodesExplored, (int)provedOptimal, elapsed_ms());
    fflush(stdout);
    _Exit(0);
}
