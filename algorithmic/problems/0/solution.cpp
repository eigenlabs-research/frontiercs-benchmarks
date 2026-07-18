#include <bits/stdc++.h>
using namespace std;

struct Cell { int x, y; };
struct Variant {
    int R, F;
    int minx, miny;
    int w, h;
    vector<Cell> cells;
};
struct Piece {
    vector<pair<int,int>> orig;
    vector<Variant> vars;
    int k;
};
struct Placement { long long X=0, Y=0; int R=0, F=0; };
struct Result { long long area = LLONG_MAX; int W=0, H=0; vector<Placement> pl; };

static pair<int,int> transformCell(int x, int y, int R, int F) {
    if (F) x = -x;
    switch (R & 3) {
        case 0: return {x, y};
        case 1: return {y, -x};
        case 2: return {-x, -y};
        default: return {-y, x};
    }
}

static vector<Variant> makeVariants(const vector<pair<int,int>>& cells) {
    vector<Variant> out;
    set<vector<pair<int,int>>> seen;
    for (int F=0; F<=1; ++F) for (int R=0; R<4; ++R) {
        vector<pair<int,int>> t;
        int mnx = INT_MAX, mny = INT_MAX, mxx = INT_MIN, mxy = INT_MIN;
        for (auto [x,y]: cells) {
            auto [a,b] = transformCell(x,y,R,F);
            t.push_back({a,b});
            mnx = min(mnx,a); mny = min(mny,b); mxx = max(mxx,a); mxy = max(mxy,b);
        }
        vector<pair<int,int>> norm;
        for (auto [a,b]: t) norm.push_back({a-mnx,b-mny});
        sort(norm.begin(), norm.end());
        if (!seen.insert(norm).second) continue;
        Variant v; v.R=R; v.F=F; v.minx=mnx; v.miny=mny; v.w=mxx-mnx+1; v.h=mxy-mny+1;
        for (auto [a,b]: norm) v.cells.push_back({a,b});
        // Put lower cells first; this often breaks ties a little more stably.
        sort(v.cells.begin(), v.cells.end(), [](const Cell& a, const Cell& b){ return a.y==b.y ? a.x<b.x : a.y<b.y; });
        out.push_back(v);
    }
    sort(out.begin(), out.end(), [](const Variant& a, const Variant& b){
        int aa=a.w*a.h, bb=b.w*b.h;
        if (aa!=bb) return aa<bb;
        if (max(a.w,a.h)!=max(b.w,b.h)) return max(a.w,a.h)<max(b.w,b.h);
        if (a.h!=b.h) return a.h<b.h;
        return a.w<b.w;
    });
    return out;
}

static Result packWidth(const vector<Piece>& pieces, const vector<int>& order, int W) {
    const int n = (int)pieces.size();
    vector<int> height(W, 0);
    vector<Placement> ans(n);
    int curH = 0, usedW = 0;

    for (int id: order) {
        const Piece& p = pieces[id];
        bool ok = false;
        const Variant* bestV = nullptr;
        int bestX = 0, bestY = 0, bestNewH = INT_MAX;
        long long bestRise = LLONG_MAX;

        for (const Variant& v: p.vars) {
            if (v.w > W) continue;
            for (int x=0; x+v.w<=W; ++x) {
                int y = 0;
                for (const Cell& c: v.cells) {
                    int need = height[x + c.x] - c.y;
                    if (need > y) y = need;
                }
                int newH = curH;
                long long rise = 0;
                for (const Cell& c: v.cells) {
                    int col = x + c.x;
                    int nh = y + c.y + 1;
                    if (nh > newH) newH = nh;
                    if (nh > height[col]) rise += nh - height[col];
                }
                if (!ok || newH < bestNewH ||
                    (newH == bestNewH && (rise < bestRise ||
                    (rise == bestRise && (y < bestY || (y == bestY && x < bestX)))))) {
                    ok = true; bestV = &v; bestX = x; bestY = y; bestNewH = newH; bestRise = rise;
                }
            }
        }
        if (!ok) { Result bad; return bad; }
        const Variant& v = *bestV;
        for (const Cell& c: v.cells) height[bestX + c.x] = max(height[bestX + c.x], bestY + c.y + 1);
        curH = max(curH, bestNewH);
        usedW = max(usedW, bestX + v.w);
        ans[id] = { (long long)bestX - v.minx, (long long)bestY - v.miny, v.R, v.F };
    }
    Result r;
    r.W = max(1, usedW);
    r.H = max(1, curH);
    r.area = 1LL * r.W * r.H;
    r.pl.swap(ans);
    return r;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    vector<Piece> pieces(n);
    long long total = 0;
    int minFeasibleW = 1;
    for (int i=0;i<n;++i) {
        int k; cin >> k; pieces[i].k = k; total += k;
        pieces[i].orig.resize(k);
        for (int j=0;j<k;++j) cin >> pieces[i].orig[j].first >> pieces[i].orig[j].second;
        pieces[i].vars = makeVariants(pieces[i].orig);
        int mw = INT_MAX;
        for (auto &v: pieces[i].vars) mw = min(mw, v.w);
        minFeasibleW = max(minFeasibleW, mw);
    }

    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        const Piece &A = pieces[a], &B = pieces[b];
        int besta=1000000, bestb=1000000, maxa=0, maxb=0;
        for (auto &v: A.vars) { besta=min(besta, v.w*v.h); maxa=max(maxa, max(v.w,v.h)); }
        for (auto &v: B.vars) { bestb=min(bestb, v.w*v.h); maxb=max(maxb, max(v.w,v.h)); }
        if (A.k != B.k) return A.k > B.k;
        if (besta != bestb) return besta > bestb;
        if (maxa != maxb) return maxa > maxb;
        return a < b;
    });

    double root = sqrt((double)max(1LL,total));
    vector<int> cand;
    auto addW = [&](int w){ if (w >= minFeasibleW) cand.push_back(w); };
    addW(minFeasibleW);
    double factors[] = {0.70, 0.85, 1.00, 1.18, 1.40, 1.70, 2.10, 2.70, 3.50};
    for (double f: factors) addW((int)ceil(root * f));
    // A few rounded widths reduce unlucky off-by-one skyline effects.
    for (int base: vector<int>{16,24,32,48,64,96,128,192,256,384,512,768,1024})
        if (abs(base - (int)root) < (int)root*3/2 + 50) addW(base);
    sort(cand.begin(), cand.end());
    cand.erase(unique(cand.begin(), cand.end()), cand.end());

    Result best;
    auto better = [](const Result& a, const Result& b) {
        return a.area < b.area || (a.area == b.area && (a.H < b.H || (a.H == b.H && a.W < b.W)));
    };
    auto start = chrono::steady_clock::now();
    set<int> tried;
    auto tryWidth = [&](int w) {
        if (w < minFeasibleW || tried.count(w)) return;
        tried.insert(w);
        Result r = packWidth(pieces, order, w);
        if (better(r, best)) best = std::move(r);
    };

    // First run the incumbent coarse width portfolio.
    for (int idx=0; idx<(int)cand.size(); ++idx) {
        if (idx > 0) {
            double elapsed = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            if (elapsed > 1.35) break;
        }
        tryWidth(cand[idx]);
    }

    // Targeted change: once a promising basin is known, spend the remaining time on
    // neighbouring strip widths.  The skyline packer is quite sensitive to width;
    // the previous geometric grid can skip a strictly better W by dozens of cells.
    if (best.area != LLONG_MAX) {
        vector<int> focus;
        for (int c: {best.W, (int)ceil(root), max(minFeasibleW, best.W - 32), best.W + 32}) {
            int rad = (c == (int)ceil(root) ? 28 : 36);
            for (int d=-rad; d<=rad; ++d) focus.push_back(c+d);
        }
        sort(focus.begin(), focus.end(), [&](int a, int b){
            int da = min(abs(a - best.W), abs(a - (int)ceil(root)));
            int db = min(abs(b - best.W), abs(b - (int)ceil(root)));
            if (da != db) return da < db;
            return a < b;
        });
        focus.erase(unique(focus.begin(), focus.end()), focus.end());
        for (int w: focus) {
            double elapsed = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            if (elapsed > 1.82) break;
            tryWidth(w);
        }
    }

    // This should never be needed, but keeps the program total for pathological inputs.
    if (best.area == LLONG_MAX) best = packWidth(pieces, order, max(minFeasibleW, (int)ceil(root)));

    cout << best.W << ' ' << best.H << '\n';
    for (int i=0;i<n;++i) cout << best.pl[i].X << ' ' << best.pl[i].Y << ' ' << best.pl[i].R << ' ' << best.pl[i].F << '\n';
    return 0;
}
