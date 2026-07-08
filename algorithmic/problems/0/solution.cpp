// Pack the Polyominoes — bottom-left-fill over bitmask rows.
// Sort pieces by size descending; for each, scan positions bottom-up,
// left-to-right, trying all distinct reflect/rotate orientations, and place
// at the first fit. Repeat for a few candidate widths, keep smallest area.
#include <bits/stdc++.h>
using namespace std;
typedef uint64_t u64;
typedef long long ll;

static chrono::steady_clock::time_point T0;
static inline double elapsedSec() {
    return chrono::duration<double>(chrono::steady_clock::now() - T0).count();
}

static inline pair<int,int> rotcw(int x, int y, int r) {
    switch (r & 3) {
        case 0:  return {x, y};
        case 1:  return {y, -x};
        case 2:  return {-x, -y};
        default: return {-y, x};
    }
}

struct Orient {
    int w, h, tminx, tminy, R, F;
    int b0; // lowest set bit of rows[0]
    uint16_t rows[10];
};

struct Piece {
    int k;
    int bx, by; // base subtracted from original coords
    vector<Orient> os;
};

struct Result {
    ll area = LLONG_MAX;
    int W = 0, H = 0;
    // per piece: ox, oy, orient index
    vector<array<int,3>> placed;
};

int main() {
    T0 = chrono::steady_clock::now();
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    vector<Piece> ps(n);
    ll S = 0;
    for (int i = 0; i < n; i++) {
        int k; cin >> k;
        ps[i].k = k;
        S += k;
        vector<pair<ll,ll>> raw(k);
        ll mnx = LLONG_MAX, mny = LLONG_MAX;
        for (auto &c : raw) { cin >> c.first >> c.second; mnx = min(mnx, c.first); mny = min(mny, c.second); }
        ps[i].bx = (int)mnx; ps[i].by = (int)mny;
        vector<pair<int,int>> cells(k);
        for (int j = 0; j < k; j++) cells[j] = {(int)(raw[j].first - mnx), (int)(raw[j].second - mny)};

        set<vector<pair<int,int>>> seen;
        for (int F = 0; F < 2; F++) {
            for (int R = 0; R < 4; R++) {
                vector<pair<int,int>> t(k);
                int tmnx = INT_MAX, tmny = INT_MAX, tmxx = INT_MIN, tmxy = INT_MIN;
                for (int j = 0; j < k; j++) {
                    int x = F ? -cells[j].first : cells[j].first;
                    t[j] = rotcw(x, cells[j].second, R);
                    tmnx = min(tmnx, t[j].first);  tmxx = max(tmxx, t[j].first);
                    tmny = min(tmny, t[j].second); tmxy = max(tmxy, t[j].second);
                }
                vector<pair<int,int>> nm(k);
                for (int j = 0; j < k; j++) nm[j] = {t[j].first - tmnx, t[j].second - tmny};
                sort(nm.begin(), nm.end());
                if (!seen.insert(nm).second) continue;
                Orient o;
                o.w = tmxx - tmnx + 1; o.h = tmxy - tmny + 1;
                o.tminx = tmnx; o.tminy = tmny; o.R = R; o.F = F;
                memset(o.rows, 0, sizeof(o.rows));
                for (auto &p : nm) o.rows[p.second] |= (uint16_t)(1u << p.first);
                o.b0 = __builtin_ctz((unsigned)o.rows[0]);
                ps[i].os.push_back(o);
            }
        }
        // prefer flat, wide orientations
        sort(ps[i].os.begin(), ps[i].os.end(), [](const Orient &a, const Orient &b) {
            if (a.h != b.h) return a.h < b.h;
            return a.w > b.w;
        });
    }

    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
        if (ps[a].k != ps[b].k) return ps[a].k > ps[b].k;
        return ps[a].os[0].w > ps[b].os[0].w;
    });

    const double HARD = 1.55;   // emergency shelf-fill beyond this
    const double NEXT = 1.10;   // don't start another width beyond this

    auto attempt = [&](int W) -> Result {
        int words = (W + 63) >> 6;
        int rows = 0;
        vector<u64> g;
        vector<int> cnt;
        auto ensureRows = [&](int need) {
            if (need > rows) {
                int nr = max(need, rows + rows / 2 + 64);
                g.resize((size_t)nr * words, 0);
                cnt.resize(nr, 0);
                rows = nr;
            }
        };
        ensureRows((int)(S / W) + 64);
        int topY = 0, fnf = 0;
        Result res;
        res.W = W;
        res.placed.assign(n, {0,0,0});
        bool emergency = false;
        int shX = 0, shY = 0, shH = 0, ctr = 0;

        auto orIn = [&](const Orient &o, int x, int y) {
            int s = x & 63, wi = x >> 6;
            for (int r = 0; r < o.h; r++) {
                u64 m = o.rows[r];
                size_t base = (size_t)(y + r) * words + wi;
                g[base] |= m << s;
                if (s) { u64 hi = m >> (64 - s); if (hi) g[base + 1] |= hi; }
                cnt[y + r] += __builtin_popcount((unsigned)m);
            }
        };

        for (int idx : order) {
            Piece &p = ps[idx];
            if (!emergency && ((++ctr & 31) == 0) && elapsedSec() > HARD) {
                emergency = true; shY = topY; shX = 0; shH = 0;
            }
            if (emergency) {
                const Orient &o = p.os[0];
                if (shX + o.w > W) { shY += shH; shX = 0; shH = 0; }
                ensureRows(shY + o.h);
                orIn(o, shX, shY);
                res.placed[idx] = {shX, shY, 0};
                shX += o.w; shH = max(shH, o.h);
                topY = max(topY, shY + o.h);
                continue;
            }
            int bx = -1, by = -1, bo = -1;
            for (int y = fnf; y <= topY && bo < 0; y++) {
                if (cnt[y] == W) continue;
                ensureRows(y + 11);
                // scan only anchors (free cells): every valid placement puts
                // the leftmost cell of its bottom row on a free cell
                size_t rowBase = (size_t)y * words;
                for (int wi = 0; wi < words && bo < 0; wi++) {
                    u64 freeBits = ~g[rowBase + wi];
                    int lim = W - (wi << 6);
                    if (lim < 64) freeBits &= (lim <= 0) ? 0ULL : ((1ULL << lim) - 1);
                    while (freeBits && bo < 0) {
                        int b = __builtin_ctzll(freeBits);
                        freeBits &= freeBits - 1;
                        int ax = (wi << 6) + b;
                        for (int oi = 0; oi < (int)p.os.size(); oi++) {
                            const Orient &o = p.os[oi];
                            int x = ax - o.b0;
                            if (x < 0 || x + o.w > W) continue;
                            int s = x & 63, xwi = x >> 6;
                            bool ok = true;
                            for (int r = 0; r < o.h; r++) {
                                u64 m = o.rows[r];
                                size_t base = (size_t)(y + r) * words + xwi;
                                if (g[base] & (m << s)) { ok = false; break; }
                                if (s) {
                                    u64 hi = m >> (64 - s);
                                    if (hi && (g[base + 1] & hi)) { ok = false; break; }
                                }
                            }
                            if (ok) { bx = x; by = y; bo = oi; break; }
                        }
                    }
                }
            }
            if (bo < 0) { bx = 0; by = topY; bo = 0; ensureRows(by + 11); }
            const Orient &o = p.os[bo];
            orIn(o, bx, by);
            res.placed[idx] = {bx, by, bo};
            topY = max(topY, by + o.h);
            while (fnf < topY && cnt[fnf] == W) fnf++;
        }
        res.H = max(topY, 1);
        res.area = (ll)W * res.H;
        return res;
    };

    int base = (int)ceil(sqrt((double)S));
    double mults[4] = {0.97, 1.00, 1.06, 1.14};
    vector<int> widths;
    for (double m : mults) {
        int w = max(12, (int)ceil(base * m));
        if (find(widths.begin(), widths.end(), w) == widths.end()) widths.push_back(w);
    }

    Result best;
    for (size_t i = 0; i < widths.size(); i++) {
        if (i > 0 && elapsedSec() > NEXT) break;
        Result r = attempt(widths[i]);
        if (r.area < best.area || (r.area == best.area && r.H < best.H)) best = std::move(r);
    }

    string out;
    out.reserve((size_t)n * 16 + 32);
    out += to_string(best.W); out += ' '; out += to_string(best.H); out += '\n';
    for (int i = 0; i < n; i++) {
        auto [ox, oy, oi] = best.placed[i];
        const Orient &o = ps[i].os[oi];
        int lbx = o.F ? -ps[i].bx : ps[i].bx;
        auto lb = rotcw(lbx, ps[i].by, o.R);
        ll X = (ll)ox - o.tminx - lb.first;
        ll Y = (ll)oy - o.tminy - lb.second;
        out += to_string(X); out += ' ';
        out += to_string(Y); out += ' ';
        out += to_string(o.R); out += ' ';
        out += to_string(o.F); out += '\n';
    }
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
