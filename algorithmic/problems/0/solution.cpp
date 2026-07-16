#include <bits/stdc++.h>
using namespace std;

struct Cell { int x, y; };
struct Orient {
    vector<Cell> cells;       // normalized cells
    int w = 0, h = 0;
    int r = 0, f = 0;
    int minx = 0, miny = 0;   // min coordinates before normalization
};
struct Piece {
    vector<pair<int,int>> raw;
    vector<Orient> ori;
    int k = 0;
};
struct AnsPlace { long long X=0, Y=0; int r=0, f=0; };
struct PackResult { long long area=LLONG_MAX, W=0, H=0; vector<AnsPlace> place; };

static pair<int,int> transform_cell(int x, int y, int r, int f) {
    if (f) x = -x;
    switch (r & 3) {
        case 0: return { x,  y};
        case 1: return { y, -x};
        case 2: return {-x, -y};
        default:return {-y,  x};
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    vector<Piece> pieces(n);
    long long total = 0;
    int globalMinFitW = 1, maxSpan = 1;

    for (int i = 0; i < n; ++i) {
        int k; cin >> k; pieces[i].k = k; total += k;
        pieces[i].raw.resize(k);
        for (int j = 0; j < k; ++j) cin >> pieces[i].raw[j].first >> pieces[i].raw[j].second;

        set<vector<pair<int,int>>> seen;
        int minWidthForPiece = INT_MAX;
        for (int f = 0; f <= 1; ++f) for (int r = 0; r < 4; ++r) {
            vector<pair<int,int>> v;
            v.reserve(k);
            int minx = INT_MAX, miny = INT_MAX, maxx = INT_MIN, maxy = INT_MIN;
            for (auto [x,y] : pieces[i].raw) {
                auto [tx,ty] = transform_cell(x, y, r, f);
                v.push_back({tx,ty});
                minx = min(minx, tx); miny = min(miny, ty);
                maxx = max(maxx, tx); maxy = max(maxy, ty);
            }
            for (auto &p : v) { p.first -= minx; p.second -= miny; }
            sort(v.begin(), v.end());
            if (!seen.insert(v).second) continue;
            Orient o; o.r = r; o.f = f; o.minx = minx; o.miny = miny;
            o.w = maxx - minx + 1; o.h = maxy - miny + 1;
            o.cells.reserve(k);
            for (auto [x,y] : v) o.cells.push_back({x,y});
            pieces[i].ori.push_back(o);
            minWidthForPiece = min(minWidthForPiece, o.w);
            maxSpan = max(maxSpan, max(o.w, o.h));
        }
        globalMinFitW = max(globalMinFitW, minWidthForPiece);
    }

    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
        auto bestBox = [&](int id) {
            int best = INT_MAX, bh = 0;
            for (auto &o: pieces[id].ori) {
                int box = o.w * o.h;
                if (box < best || (box == best && max(o.w,o.h) > bh)) { best = box; bh = max(o.w,o.h); }
            }
            return pair<int,int>{best, bh};
        };
        auto A = bestBox(a), B = bestBox(b);
        int wasteA = A.first - pieces[a].k, wasteB = B.first - pieces[b].k;
        if (pieces[a].k != pieces[b].k) return pieces[a].k > pieces[b].k;
        if (wasteA != wasteB) return wasteA > wasteB;
        return A.second > B.second;
    });

    auto packWidth = [&](int targetW) -> PackResult {
        targetW = max(targetW, globalMinFitW);
        vector<vector<unsigned long long>> rows(32, vector<unsigned long long>((targetW + 63) >> 6, 0));
        vector<int> colH(targetW, 0);
        vector<AnsPlace> place(n);
        int curH = 0, usedW = 0;
        auto ensureRows = [&](int h) {
            if ((int)rows.size() < h) rows.resize(h, vector<unsigned long long>((targetW + 63) >> 6, 0));
        };
        auto occupied = [&](int x, int y)->bool {
            return y >= 0 && y < (int)rows.size() && ((rows[y][x>>6] >> (x&63)) & 1ULL);
        };
        auto setocc = [&](int x, int y) {
            rows[y][x>>6] |= 1ULL << (x&63);
        };

        for (int id : order) {
            long long bestKey = LLONG_MAX;
            int bestX = 0, bestY = 0, bestOi = -1;
            const Piece &pc = pieces[id];
            for (int oi = 0; oi < (int)pc.ori.size(); ++oi) {
                const Orient &o = pc.ori[oi];
                if (o.w > targetW) continue;
                for (int x = 0; x <= targetW - o.w; ++x) {
                    int y = 0;
                    for (const Cell &c : o.cells) y = max(y, colH[x + c.x] - c.y);
                    if (y < 0) y = 0;
                    for (;;) {
                        bool bad = false;
                        for (const Cell &c : o.cells) if (occupied(x+c.x, y+c.y)) { bad = true; break; }
                        if (!bad) break;
                        ++y;
                    }
                    int newH = max(curH, y + o.h);
                    int newW = max(usedW, x + o.w);
                    // Prefer smaller area increase, then low/top-left placements.  The width term discourages needless right spread.
                    long long key = 1000000LL * newH + 1000LL * newW + y * 16LL + x;
                    if (key < bestKey) { bestKey = key; bestX = x; bestY = y; bestOi = oi; }
                }
            }
            if (bestOi < 0) { // should not happen if targetW >= globalMinFitW
                bestOi = 0; bestX = 0; bestY = curH;
            }
            const Orient &o = pc.ori[bestOi];
            ensureRows(bestY + o.h + 1);
            for (const Cell &c : o.cells) {
                int gx = bestX + c.x, gy = bestY + c.y;
                setocc(gx, gy);
                colH[gx] = max(colH[gx], gy + 1);
            }
            curH = max(curH, bestY + o.h);
            usedW = max(usedW, bestX + o.w);
            place[id] = { (long long)bestX - o.minx, (long long)bestY - o.miny, o.r, o.f };
        }
        PackResult res;
        res.W = max(1, usedW); res.H = max(1, curH); res.area = res.W * res.H; res.place = move(place);
        return res;
    };

    double root = sqrt((double)max(1LL,total));
    vector<int> widths;
    auto addW = [&](int w){ w = max(w, globalMinFitW); if (w > 0) widths.push_back(w); };
    addW((int)(root * 0.75)); addW((int)(root * 0.90)); addW((int)(root * 1.05));
    addW((int)(root * 1.25)); addW((int)(root * 1.50)); addW((int)(root * 1.85));
    addW(maxSpan); addW(globalMinFitW);
    sort(widths.begin(), widths.end()); widths.erase(unique(widths.begin(), widths.end()), widths.end());

    // Limit very large cases to the most promising widths.
    if (total > 60000 && widths.size() > 4) {
        vector<int> keep;
        for (int w : widths) if (w >= (int)(root*0.9) && w <= (int)(root*1.55)) keep.push_back(w);
        if (!keep.empty()) widths = keep;
    }

    PackResult best;
    for (int w : widths) {
        PackResult cur = packWidth(w);
        if (cur.area < best.area || (cur.area == best.area && make_pair(cur.H,cur.W) < make_pair(best.H,best.W))) best = move(cur);
    }

    cout << best.W << ' ' << best.H << '\n';
    for (int i = 0; i < n; ++i) {
        cout << best.place[i].X << ' ' << best.place[i].Y << ' ' << best.place[i].r << ' ' << best.place[i].f << '\n';
    }
    return 0;
}
