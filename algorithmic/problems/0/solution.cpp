#include <bits/stdc++.h>
using namespace std;

struct Cell { int x, y; };
struct Orient {
    vector<Cell> cells;
    int w = 0, h = 0, r = 0, f = 0;
    int minx = 0, miny = 0; // before normalisation, after reflect+rotate
};
struct Piece {
    vector<Cell> orig;
    vector<Orient> ors;
    int k = 0;
};
struct AnsPlace { long long X=0, Y=0; int r=0, f=0; };
struct PackResult { long long area = LLONG_MAX; int W=0, H=0; vector<AnsPlace> p; bool ok=false; };

static pair<int,int> trans(int x, int y, int r, int f) {
    if (f) x = -x;
    switch (r & 3) {
        case 0: return {x, y};
        case 1: return {y, -x};
        case 2: return {-x, -y};
        default: return {-y, x};
    }
}

static vector<Orient> make_orients(const vector<Cell>& a) {
    vector<Orient> res;
    set<vector<pair<int,int>>> seen;
    for (int f=0; f<2; ++f) for (int r=0; r<4; ++r) {
        vector<pair<int,int>> v;
        int mnx=INT_MAX, mny=INT_MAX, mxx=INT_MIN, mxy=INT_MIN;
        for (auto c: a) {
            auto [x,y] = trans(c.x,c.y,r,f);
            v.push_back({x,y});
            mnx=min(mnx,x); mny=min(mny,y); mxx=max(mxx,x); mxy=max(mxy,y);
        }
        vector<pair<int,int>> norm;
        for (auto [x,y]: v) norm.push_back({x-mnx, y-mny});
        sort(norm.begin(), norm.end());
        norm.erase(unique(norm.begin(), norm.end()), norm.end());
        if ((int)norm.size() != (int)a.size()) continue; // should not happen for valid input
        if (seen.insert(norm).second) {
            Orient o; o.r=r; o.f=f; o.minx=mnx; o.miny=mny; o.w=mxx-mnx+1; o.h=mxy-mny+1;
            for (auto [x,y]: norm) o.cells.push_back({x,y});
            // checking low rows first tends to fill cavities earlier
            sort(o.cells.begin(), o.cells.end(), [](const Cell& A, const Cell& B){
                if (A.y != B.y) return A.y < B.y;
                return A.x < B.x;
            });
            res.push_back(o);
        }
    }
    sort(res.begin(), res.end(), [](const Orient& A, const Orient& B){
        int aa=A.w*A.h, bb=B.w*B.h;
        if (aa!=bb) return aa<bb;
        if (A.h!=B.h) return A.h<B.h;
        return A.w<B.w;
    });
    return res;
}

static PackResult pack_width(const vector<Piece>& pieces, const vector<int>& order, int W) {
    int n = pieces.size();
    PackResult pr; pr.W = W; pr.p.assign(n, {});
    vector<vector<unsigned char>> occ;
    int H = 0;
    long long filled = 0;

    auto ensureH = [&](int nh) {
        while ((int)occ.size() < nh) occ.emplace_back(W, 0);
    };
    auto placeContact = [&](const Orient& o, int x, int y)->int {
        if (x < 0 || x + o.w > W) return -1;
        int contact = 0;
        for (const auto &c: o.cells) {
            int xx = x + c.x, yy = y + c.y;
            if (yy < (int)occ.size() && occ[yy][xx]) return -1;
        }
        for (const auto &c: o.cells) {
            int xx = x + c.x, yy = y + c.y;
            if (xx == 0) ++contact;
            if (xx == W-1) ++contact;
            if (yy == 0) ++contact;
            static const int dx[4] = {1,-1,0,0};
            static const int dy[4] = {0,0,1,-1};
            for (int d=0; d<4; ++d) {
                int nx = xx + dx[d], ny = yy + dy[d];
                if (0 <= nx && nx < W && 0 <= ny && ny < (int)occ.size() && occ[ny][nx]) ++contact;
            }
        }
        return contact;
    };

    for (int id: order) {
        const Piece& pc = pieces[id];
        int bestNewH = INT_MAX, bestY = 0, bestX = 0, bestOi = -1, bestContact = -1;
        int initial = INT_MAX;
        for (int oi=0; oi<(int)pc.ors.size(); ++oi) if (pc.ors[oi].w <= W) initial = min(initial, H + pc.ors[oi].h);
        if (initial == INT_MAX) return pr;
        bestNewH = initial;
        // No placement with y+h >= bestNewH can improve the primary objective for this step.
        for (int oi=0; oi<(int)pc.ors.size(); ++oi) {
            const Orient& o = pc.ors[oi];
            if (o.w > W) continue;
            int maxY = bestNewH - o.h;
            if (maxY < 0) continue;
            for (int y=0; y<=maxY; ++y) {
                int newH = max(H, y + o.h);
                if (newH > bestNewH) continue;
                for (int x=0; x + o.w <= W; ++x) {
                    int contact = placeContact(o,x,y);
                    if (contact >= 0) {
                        if (newH < bestNewH ||
                            (newH == bestNewH && (contact > bestContact ||
                             (contact == bestContact && (y < bestY || (y == bestY && x < bestX)))))) {
                            bestNewH = newH; bestY = y; bestX = x; bestOi = oi; bestContact = contact;
                        }
                    }
                }
            }
        }
        if (bestOi < 0) { // place just above current skyline
            for (int oi=0; oi<(int)pc.ors.size(); ++oi) if (pc.ors[oi].w <= W) { bestOi=oi; bestX=0; bestY=H; bestNewH=H+pc.ors[oi].h; break; }
        }
        const Orient& o = pc.ors[bestOi];
        ensureH(bestY + o.h);
        for (auto c: o.cells) occ[bestY+c.y][bestX+c.x] = 1;
        H = max(H, bestY + o.h);
        filled += pc.k;
        pr.p[id] = {bestX - o.minx, bestY - o.miny, o.r, o.f};
    }
    pr.H = max(1,H); pr.area = 1LL * pr.W * pr.H; pr.ok = true;
    return pr;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n; if (!(cin >> n)) return 0;
    vector<Piece> pieces(n);
    long long total = 0;
    for (int i=0;i<n;++i) {
        int k; cin >> k; pieces[i].k = k; total += k;
        pieces[i].orig.resize(k);
        for (int j=0;j<k;++j) cin >> pieces[i].orig[j].x >> pieces[i].orig[j].y;
        pieces[i].ors = make_orients(pieces[i].orig);
    }
    vector<int> order(n); iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        const Piece &A=pieces[a], &B=pieces[b];
        int besta=A.ors[0].w*A.ors[0].h, bestb=B.ors[0].w*B.ors[0].h;
        int maxda=0, maxdb=0;
        for (auto &o:A.ors) maxda=max(maxda, max(o.w,o.h));
        for (auto &o:B.ors) maxdb=max(maxdb, max(o.w,o.h));
        if (A.k != B.k) return A.k > B.k;
        if (besta != bestb) return besta > bestb;
        if (maxda != maxdb) return maxda > maxdb;
        return a < b;
    });
    int minW = 1;
    for (auto &p: pieces) {
        int mw = INT_MAX;
        for (auto &o: p.ors) mw = min(mw, o.w);
        minW = max(minW, mw);
    }
    double root = sqrt((double)max(1LL,total));
    vector<double> fac = {0.72,0.82,0.92,1.02,1.15,1.30,1.50,1.75,2.10};
    vector<int> widths;
    for (double f: fac) widths.push_back(max(minW, (int)ceil(root*f)));
    for (int extra: {0,1,2,4,8,16}) widths.push_back(minW + extra);
    sort(widths.begin(), widths.end()); widths.erase(unique(widths.begin(), widths.end()), widths.end());

    PackResult best;
    auto start = chrono::steady_clock::now();
    for (int idx=0; idx<(int)widths.size(); ++idx) {
        auto now = chrono::steady_clock::now();
        double elapsed = chrono::duration<double>(now-start).count();
        if (idx > 0 && elapsed > 1.75) break;
        PackResult cur = pack_width(pieces, order, widths[idx]);
        if (cur.ok && (cur.area < best.area || (cur.area == best.area && (cur.H < best.H || (cur.H == best.H && cur.W < best.W))))) best = move(cur);
    }
    if (!best.ok) best = pack_width(pieces, order, max(minW, (int)ceil(root*1.3)));
    cout << best.W << ' ' << best.H << '\n';
    for (int i=0;i<n;++i) cout << best.p[i].X << ' ' << best.p[i].Y << ' ' << best.p[i].r << ' ' << best.p[i].f << '\n';
    return 0;
}
