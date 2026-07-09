#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <climits>
#include <random>
#include <chrono>
#include <cassert>

using namespace std;
using ll = long long;

static chrono::steady_clock::time_point T_START;
static inline double elapsed(){
    return chrono::duration<double>(chrono::steady_clock::now() - T_START).count();
}
static double TIME_LIMIT = 0.96;

struct JsonParser {
    const string& s;
    size_t i = 0;
    explicit JsonParser(const string& src) : s(src) {}
    void skipWS(){ while(i < s.size() && (s[i]==' '||s[i]=='\t'||s[i]=='\n'||s[i]=='\r')) ++i; }
    char peek(){ skipWS(); return i<s.size()? s[i] : 0; }
    char getc(){ skipWS(); return i<s.size()? s[i++] : 0; }
    void expectChar(char c){ skipWS(); if(i<s.size() && s[i]==c) ++i; }
    string readString(){
        skipWS();
        if(i<s.size() && s[i]=='"') ++i;
        string out;
        while(i < s.size() && s[i] != '"'){
            if(s[i]=='\\'){
                ++i;
                if(i<s.size()){ char e = s[i]; ++i;
                    switch(e){
                        case '"': out.push_back('"'); break;
                        case '\\': out.push_back('\\'); break;
                        case '/': out.push_back('/'); break;
                        case 'b': out.push_back('\b'); break;
                        case 'f': out.push_back('\f'); break;
                        case 'n': out.push_back('\n'); break;
                        case 'r': out.push_back('\r'); break;
                        case 't': out.push_back('\t'); break;
                        default: out.push_back(e); break;
                    }
                }
            } else { out.push_back(s[i]); ++i; }
        }
        if(i<s.size()) ++i;
        return out;
    }
    ll readLong(){
        skipWS();
        bool neg = false;
        if(i<s.size() && (s[i]=='+'||s[i]=='-')){ neg = (s[i]=='-'); ++i; }
        ll v = 0;
        while(i<s.size() && s[i]>='0' && s[i]<='9'){ v = v*10 + (s[i]-'0'); ++i; }
        return neg ? -v : v;
    }
    bool readBool(){
        skipWS();
        if(i<s.size() && s[i]=='t'){ i += 4; return true; }
        if(i<s.size() && s[i]=='f'){ i += 5; return false; }
        return false;
    }
};

struct ItemType {
    string name;
    int w, h;
    ll v;
    int limit;
    ll area;
    double density; // v / area
};
struct Bin { int W, H; bool allow_rotate; };

static Bin g_bin;
static vector<ItemType> g_items;

struct Placed {
    int typeId;
    int x, y;
    int rot; // 0 or 1
};

struct Skyline {
    int W, H;
    struct Seg { int x; int height; }; // height applies to [x, next seg x)
    vector<Seg> segs; // sorted by x, last seg.x == W
    void init(int W_, int H_){
        W = W_; H = H_;
        segs.clear();
        segs.push_back({0, 0});
        segs.push_back({W, 0}); // sentinel marking end
    }
    bool findBest(int rw, int rh, int& outX, int& outY) const {
        if(rw <= 0 || rh <= 0 || rw > W || rh > H) return false;
        int bestX = -1, bestY = INT32_MAX;
        int n = (int)segs.size() - 1; // last index is sentinel
        for(int i = 0; i < n; ++i){
            int x = segs[i].x;
            if(x + rw > W) break;
            int mx = 0;
            int j = i;
            int curX = x;
            int covered = 0;
            int k = i;
            while(covered < rw){
                if(k >= n) { mx = INT32_MAX; break; } // ran out
                int segLeft = segs[k].x;
                int segRight = segs[k+1].x;
                int segStart = max(segLeft, x);
                int segEnd = min(segRight, x + rw);
                if(segEnd <= segStart){ ++k; continue; }
                if(segs[k].height > mx) mx = segs[k].height;
                covered += segEnd - segStart;
                ++k;
            }
            if(mx == INT32_MAX) continue;
            int y = mx;
            if(y + rh <= H){
                if(y < bestY || (y == bestY && x < bestX)){
                    bestY = y; bestX = x;
                }
            }
        }
        if(bestX < 0) return false;
        outX = bestX; outY = bestY;
        return true;
    }
    void place(int x, int y, int rw, int rh){
        int x2 = x + rw;
        int newH = y + rh;
        vector<Seg> out;
        out.reserve(segs.size() + 4);
        int i = 0;
        while(i < (int)segs.size() && segs[i].x < x){
            out.push_back(segs[i]);
            ++i;
        }
        if(!out.empty() && out.back().x < x){
        }
        if(!out.empty() && out.back().height == newH){
        } else {
            out.push_back({x, newH});
        }
        while(i < (int)segs.size() && segs[i].x < x2){
            ++i;
        }
        int hAfter = 0;
        {
            int lo = 0, hi = (int)segs.size() - 1;
            while(lo < hi){
                int mid = (lo + hi + 1) / 2;
                if(segs[mid].x <= x2) lo = mid; else hi = mid - 1;
            }
            hAfter = segs[lo].height;
        }
        if(newH != hAfter){
            out.push_back({x2, hAfter});
        }
        while(i < (int)segs.size()){
            if(segs[i].x > x2){
                out.push_back(segs[i]);
            } else if(segs[i].x == x2){
            }
            ++i;
        }
        if(out.empty() || out.back().x != W){
            out.push_back({W, 0});
        }
        vector<Seg> compact;
        compact.reserve(out.size());
        for(auto& s : out){
            if(!compact.empty() && compact.back().height == s.height){
                continue;
            }
            compact.push_back(s);
        }
        if(compact.empty() || compact[0].x != 0){
            compact.insert(compact.begin(), {0, 0});
        }
        if(compact.back().x != W){
            compact.push_back({W, 0});
        }
        segs = std::move(compact);
    }
};

struct MaxRects {
    int W, H;
    vector<int> fx, fy, fw, fh;
    void init(int W_, int H_){
        W = W_; H = H_;
        fx.clear(); fy.clear(); fw.clear(); fh.clear();
        fx.push_back(0); fy.push_back(0); fw.push_back(W_); fh.push_back(H_);
    }
    bool findBest(int rw, int rh, int& outX, int& outY, int& outIdx) const {
        if(rw <= 0 || rh <= 0) return false;
        int bestShort = INT32_MAX, bestLong = INT32_MAX; outIdx = -1;
        int n = (int)fx.size();
        for(int i = 0; i < n; ++i){
            if(fw[i] >= rw && fh[i] >= rh){
                int lw = fw[i] - rw, lh = fh[i] - rh;
                int s = lw < lh ? lw : lh;
                int l = lw < lh ? lh : lw;
                if(s < bestShort || (s == bestShort && l < bestLong)){
                    bestShort = s; bestLong = l; outX = fx[i]; outY = fy[i]; outIdx = i;
                }
            }
        }
        return outIdx >= 0;
    }
    static bool contains(int ax,int ay,int aw,int ah,int bx,int by,int bw,int bh){
        return ax <= bx && ay <= by && ax + aw >= bx + bw && ay + ah >= by + bh;
    }
    void place(int x, int y, int rw, int rh){
        int px1 = x, py1 = y, px2 = x + rw, py2 = y + rh;
        int n = (int)fx.size();
        static vector<int> ax, ay, aw, ah;
        ax.clear(); ay.clear(); aw.clear(); ah.clear();
        for(int i = 0; i < n; ){
            int gx1 = fx[i], gy1 = fy[i], gx2 = fx[i] + fw[i], gy2 = fy[i] + fh[i];
            if(px1 < gx2 && px2 > gx1 && py1 < gy2 && py2 > gy1){
                if(px1 > gx1){ ax.push_back(gx1); ay.push_back(gy1); aw.push_back(px1 - gx1); ah.push_back(gy2 - gy1); }
                if(px2 < gx2){ ax.push_back(px2); ay.push_back(gy1); aw.push_back(gx2 - px2); ah.push_back(gy2 - gy1); }
                if(py1 > gy1){ ax.push_back(gx1); ay.push_back(gy1); aw.push_back(gx2 - gx1); ah.push_back(py1 - gy1); }
                if(py2 < gy2){ ax.push_back(gx1); ay.push_back(py2); aw.push_back(gx2 - gx1); ah.push_back(gy2 - py2); }
                fx[i] = fx[n-1]; fy[i] = fy[n-1]; fw[i] = fw[n-1]; fh[i] = fh[n-1];
                fx.pop_back(); fy.pop_back(); fw.pop_back(); fh.pop_back();
                --n;
            } else ++i;
        }
        for(size_t k = 0; k < ax.size(); ++k){
            if(aw[k] <= 0 || ah[k] <= 0) continue;
            fx.push_back(ax[k]); fy.push_back(ay[k]); fw.push_back(aw[k]); fh.push_back(ah[k]);
        }
        prune();
    }
    void prune(){
        int n = (int)fx.size();
        for(int i = 0; i < n; ++i){
            for(int j = i + 1; j < n; ){
                if(contains(fx[i],fy[i],fw[i],fh[i], fx[j],fy[j],fw[j],fh[j])){
                    fx[j]=fx[n-1]; fy[j]=fy[n-1]; fw[j]=fw[n-1]; fh[j]=fh[n-1];
                    fx.pop_back(); fy.pop_back(); fw.pop_back(); fh.pop_back(); --n;
                } else if(contains(fx[j],fy[j],fw[j],fh[j], fx[i],fy[i],fw[i],fh[i])){
                    fx[i]=fx[n-1]; fy[i]=fy[n-1]; fw[i]=fw[n-1]; fh[i]=fh[n-1];
                    fx.pop_back(); fy.pop_back(); fw.pop_back(); fh.pop_back(); --n;
                    --i; break;
                } else ++j;
            }
        }
    }
};

struct PackResult {
    vector<Placed> placements;
    ll totalValue;
    vector<int> used;
};

static PackResult greedyFill(const vector<int>& order, bool allowRotate, bool tryBothOrient){
    PackResult res;
    res.totalValue = 0;
    res.used.assign(g_items.size(), 0);
    Skyline sk;
    sk.init(g_bin.W, g_bin.H);

    vector<int> rem(g_items.size());
    for(size_t i = 0; i < g_items.size(); ++i) rem[i] = g_items[i].limit;

    bool progress = true;
    while(progress){
        progress = false;
        for(int t : order){
            if(rem[t] <= 0) continue;
            const ItemType& it = g_items[t];
            int x0=-1,y0=-1, x1=-1,y1=-1;
            if(it.w <= sk.W && it.h <= sk.H){
                sk.findBest(it.w, it.h, x0, y0);
            }
            int chosenRot = -1, ox=0, oy=0, rw=0, rh=0;
            if(tryBothOrient && allowRotate){
                if(it.h <= sk.W && it.w <= sk.H){
                    sk.findBest(it.h, it.w, x1, y1);
                }
                if(x0>=0 && x1>=0){
                    if(y0 <= y1){ chosenRot=0; ox=x0; oy=y0; rw=it.w; rh=it.h; }
                    else { chosenRot=1; ox=x1; oy=y1; rw=it.h; rh=it.w; }
                } else if(x0>=0){ chosenRot=0; ox=x0; oy=y0; rw=it.w; rh=it.h; }
                else if(x1>=0){ chosenRot=1; ox=x1; oy=y1; rw=it.h; rh=it.w; }
            } else {
                if(x0>=0){ chosenRot=0; ox=x0; oy=y0; rw=it.w; rh=it.h; }
            }
            if(chosenRot >= 0){
                Placed p; p.typeId = t; p.x = ox; p.y = oy; p.rot = chosenRot;
                res.placements.push_back(p);
                res.totalValue += it.v;
                res.used[t]++;
                rem[t]--;
                sk.place(ox, oy, rw, rh);
                progress = true;
            }
        }
    }
    return res;
}

static void mrGreedyLoop(MaxRects& mr, vector<int>& rem, const vector<int>& order,
                         bool allowRotate, bool tryBothOrient, PackResult& res){
    bool progress = true;
    while(progress){
        progress = false;
        for(int t : order){
            if(rem[t] <= 0) continue;
            const ItemType& it = g_items[t];
            int x0,y0,i0, x1,y1,i1;
            bool f0 = false, f1 = false;
            if(it.w <= mr.W && it.h <= mr.H) f0 = mr.findBest(it.w, it.h, x0, y0, i0);
            if(tryBothOrient && allowRotate && it.h <= mr.W && it.w <= mr.H)
                f1 = mr.findBest(it.h, it.w, x1, y1, i1);
            int chosenRot = -1, ox=0, oy=0, rw=0, rh=0;
            if(f0 && f1){
                if(y0 < y1 || (y0 == y1 && x0 <= x1)){ chosenRot=0; ox=x0; oy=y0; rw=it.w; rh=it.h; }
                else { chosenRot=1; ox=x1; oy=y1; rw=it.h; rh=it.w; }
            } else if(f0){ chosenRot=0; ox=x0; oy=y0; rw=it.w; rh=it.h; }
            else if(f1){ chosenRot=1; ox=x1; oy=y1; rw=it.h; rh=it.w; }
            if(chosenRot >= 0){
                Placed p; p.typeId=t; p.x=ox; p.y=oy; p.rot=chosenRot;
                res.placements.push_back(p);
                res.totalValue += it.v;
                res.used[t]++; rem[t]--;
                mr.place(ox, oy, rw, rh);
                progress = true;
            }
        }
    }
}

static PackResult greedyFillMaxRects(const vector<int>& order, bool allowRotate, bool tryBothOrient){
    PackResult res;
    res.totalValue = 0;
    res.used.assign(g_items.size(), 0);
    MaxRects mr;
    mr.init(g_bin.W, g_bin.H);
    vector<int> rem(g_items.size());
    for(size_t i = 0; i < g_items.size(); ++i) rem[i] = g_items[i].limit;
    mrGreedyLoop(mr, rem, order, allowRotate, tryBothOrient, res);
    return res;
}

static PackResult fillGaps(const PackResult& base, const vector<int>& order,
                           bool allowRotate, bool tryBothOrient){
    PackResult res = base; // copy placements, totalValue, used
    MaxRects mr;
    mr.init(g_bin.W, g_bin.H);
    for(const Placed& p : base.placements){
        const ItemType& it = g_items[p.typeId];
        int rw = p.rot ? it.h : it.w;
        int rh = p.rot ? it.w : it.h;
        mr.place(p.x, p.y, rw, rh);
    }
    vector<int> rem(g_items.size());
    for(size_t i = 0; i < g_items.size(); ++i) rem[i] = g_items[i].limit - res.used[i];
    mrGreedyLoop(mr, rem, order, allowRotate, tryBothOrient, res);
    return res;
}

static PackResult pruneLow(const PackResult& base, int cnt, const vector<int>& order, bool allowRotate){
    int n = (int)base.placements.size();
    if(cnt <= 0 || n == 0) return base;
    vector<int> ids(n);
    for(int i = 0; i < n; ++i) ids[i] = i;
    sort(ids.begin(), ids.end(), [&](int a, int b){
        int ta = base.placements[a].typeId, tb = base.placements[b].typeId;
        if(g_items[ta].density != g_items[tb].density) return g_items[ta].density < g_items[tb].density;
        return g_items[ta].area > g_items[tb].area;
    });
    vector<char> drop(n, 0);
    for(int i = 0; i < cnt && i < n; ++i) drop[ids[i]] = 1;
    PackResult r; r.totalValue = 0; r.used.assign(g_items.size(), 0);
    for(int i = 0; i < n; ++i) if(!drop[i]){
        const Placed& p = base.placements[i];
        r.placements.push_back(p); r.totalValue += g_items[p.typeId].v; r.used[p.typeId]++;
    }
    return fillGaps(r, order, allowRotate, allowRotate);
}

static PackResult pruneSide(const PackResult& base,int cnt,const vector<int>& order,bool allowRotate,int side){
    int n=(int)base.placements.size(); if(cnt<=0||n==0) return base;
    vector<int> ids(n); for(int i=0;i<n;++i) ids[i]=i;
    sort(ids.begin(),ids.end(),[&](int a,int b){
        const Placed& pa=base.placements[a]; const Placed& pb=base.placements[b];
        const ItemType& ia=g_items[pa.typeId]; const ItemType& ib=g_items[pb.typeId];
        int aw=pa.rot?ia.h:ia.w, ah=pa.rot?ia.w:ia.h, bw=pb.rot?ib.h:ib.w, bh=pb.rot?ib.w:ib.h;
        int ka=(side==0?pa.x:(side==1?pa.y:(side==2?pa.x+aw:pa.y+ah)));
        int kb=(side==0?pb.x:(side==1?pb.y:(side==2?pb.x+bw:pb.y+bh)));
        if(ka!=kb) return side<2?ka<kb:ka>kb;
        if(ia.density!=ib.density) return ia.density<ib.density;
        return ia.area>ib.area;
    });
    vector<char> drop(n,0); for(int i=0;i<cnt&&i<n;++i) drop[ids[i]]=1;
    PackResult r; r.totalValue=0; r.used.assign(g_items.size(),0);
    for(int i=0;i<n;++i) if(!drop[i]){ const Placed&p=base.placements[i]; r.placements.push_back(p); r.totalValue+=g_items[p.typeId].v; r.used[p.typeId]++; }
    return fillGaps(r,order,allowRotate,allowRotate);
}

static PackResult pruneLine(const PackResult& base,int lines,const vector<int>& order,bool allowRotate,int side){
    int n=(int)base.placements.size(); if(lines<=0||n==0) return base;
    vector<int> ks; ks.reserve(n);
    for(const Placed&p:base.placements){
        const ItemType& it=g_items[p.typeId];
        int h=p.rot?it.w:it.h;
        ks.push_back(side==0?p.y:p.y+h);
    }
    sort(ks.begin(),ks.end()); ks.erase(unique(ks.begin(),ks.end()),ks.end());
    if((int)ks.size()<lines) return base;
    int th=side==0?ks[lines-1]:ks[(int)ks.size()-lines];
    vector<char> drop(n,0); int dc=0;
    for(int i=0;i<n;++i){
        const Placed&p=base.placements[i];
        const ItemType& it=g_items[p.typeId];
        int h=p.rot?it.w:it.h;
        int k=side==0?p.y:p.y+h;
        bool d=side==0?k<=th:k>=th;
        if(d){ drop[i]=1; ++dc; }
    }
    if(dc==0||dc>n/4) return base;
    PackResult r; r.totalValue=0; r.used.assign(g_items.size(),0);
    for(int i=0;i<n;++i) if(!drop[i]){ const Placed&p=base.placements[i]; r.placements.push_back(p); r.totalValue+=g_items[p.typeId].v; r.used[p.typeId]++; }
    return fillGaps(r,order,allowRotate,allowRotate);
}


static PackResult choiceMaxRects(bool allowRot, int mode, double tlim){
    PackResult res; res.totalValue = 0; res.used.assign(g_items.size(), 0);
    MaxRects mr; mr.init(g_bin.W, g_bin.H);
    int M = (int)g_items.size();
    vector<int> rem(M);
    for(int i = 0; i < M; ++i) rem[i] = g_items[i].limit;
    while(elapsed() < tlim){
        double bs = -1e100; int bt=-1, br=0, bx=0, by=0, bw=0, bh=0;
        int nf = (int)mr.fx.size();
        for(int t = 0; t < M; ++t){
            if(rem[t] <= 0) continue;
            const ItemType& it = g_items[t];
            for(int rot = 0; rot <= (allowRot ? 1 : 0); ++rot){
                int rw = rot ? it.h : it.w, rh = rot ? it.w : it.h;
                for(int i = 0; i < nf; ++i){
                    if(rw > mr.fw[i] || rh > mr.fh[i]) continue;
                    int lw = mr.fw[i] - rw, lh = mr.fh[i] - rh;
                    int sh = min(lw, lh), lo = max(lw, lh);
                    double sc = it.density * 1000000000.0 - sh * 25000.0 - lo;
                    if(mode == 1) sc = (double)it.v / (1.0 + sh + 0.08 * lo);
                    else if(mode == 2) sc = (double)it.v - (double)sh * 50000.0 - lo * 500.0;
                    else if(mode == 6) sc = it.density * 1000000000.0 - (double)(mr.fw[i] * mr.fh[i] - rw * rh) * 100.0 - sh * 1000.0;
                    else if(mode == 7){
                        int contact = (mr.fx[i] == 0 ? rh : 0) + (mr.fy[i] == 0 ? rw : 0) +
                                      (mr.fx[i] + rw == mr.W ? rh : 0) + (mr.fy[i] + rh == mr.H ? rw : 0);
                        sc = it.density * 1000000000.0 + contact * 1000000.0 - sh * 25000.0 - lo;
                    }
                    int px = mr.fx[i], py = mr.fy[i];
                    if(mode == 1) px += lw;
                    else if(mode == 2 && lw > lh) px += lw;
                    else if(mode == 3) py += lh;
                    else if(mode == 4){ px += lw; py += lh; }
                    else if(mode == 5){ if(lw > lh) px += lw; else py += lh; }
                    sc -= (double)py * 0.01;
                    if(sc > bs){ bs=sc; bt=t; br=rot; bx=px; by=py; bw=rw; bh=rh; }
                }
            }
        }
        if(bt < 0) break;
        Placed p; p.typeId=bt; p.x=bx; p.y=by; p.rot=br;
        res.placements.push_back(p); res.totalValue += g_items[bt].v; res.used[bt]++; rem[bt]--;
        mr.place(bx, by, bw, bh);
    }
    return res;
}

struct ShelfUnit { int typeId; int rot; int count; }; // count = number of shelves of this kind
static PackResult knapsackShelfPlan(bool allowRot){
    int H = g_bin.H, W = g_bin.W;
    int M = (int)g_items.size();
    struct KItem { int w; ll val; int typeId; int rot; int nShelves; int perShelf; bool partial; };
    struct SD { bool valid; int perShelf; int rh; int full; int rem; };
    vector<array<SD,2>> sd(M);
    for(int t = 0; t < M; ++t){
        const ItemType& it = g_items[t];
        for(int rot = 0; rot <= 1; ++rot){
            if(rot == 1 && !allowRot){ sd[t][rot] = {false,0,0,0,0}; continue; }
            int rw = rot ? it.h : it.w;
            int rh = rot ? it.w : it.h;
            if(rw <= 0 || rh <= 0 || rw > W || rh > H){ sd[t][rot] = {false,0,0,0,0}; continue; }
            int perShelf = W / rw;
            if(perShelf <= 0){ sd[t][rot] = {false,0,0,0,0}; continue; }
            int full = it.limit / perShelf;
            int rem = it.limit % perShelf;
            sd[t][rot] = {true, perShelf, rh, full, rem};
        }
    }
    vector<ll> dp(H + 1, -1);
    vector<int> par(H + 1, -1);
    ll bestVal = -1;
    vector<int> bestCombo(M, 0);
    // Instead of 2^M combos (up to 4096 DPs), use greedy hill-climb: start from
    // best-density rotation per item, then try flipping one item at a time until
    // no single flip improves. On the judge this reclaims ~50% of the time budget
    // for downstream phases where actual quality wins live.
    auto evalCombo = [&](const vector<int>& combo, ll* outBest) -> ll {
        vector<KItem> kitems;
        kitems.reserve(M * 12);
        for(int t = 0; t < M; ++t){
            int rot = combo[t];
            if(!sd[t][rot].valid) continue;
            const ItemType& it = g_items[t];
            int perShelf = sd[t][rot].perShelf;
            int rh = sd[t][rot].rh;
            int full = sd[t][rot].full;
            int rem = sd[t][rot].rem;
            int k = full, p = 1;
            while(k > 0){
                int take = min(p, k);
                kitems.push_back({rh * take, (ll)perShelf * take * it.v, t, rot, take, perShelf, false});
                k -= take; p *= 2;
            }
            if(rem > 0){
                kitems.push_back({rh, (ll)rem * it.v, t, rot, 1, rem, true});
            }
        }
        fill(dp.begin(), dp.end(), -1);
        dp[0] = 0;
        for(const KItem& ki : kitems){
            for(int x = H; x >= ki.w; --x){
                if(dp[x - ki.w] >= 0 && dp[x - ki.w] + ki.val > dp[x]){
                    dp[x] = dp[x - ki.w] + ki.val;
                }
            }
        }
        ll bestV = 0;
        for(int x = 0; x <= H; ++x){ if(dp[x] > bestV) bestV = dp[x]; }
        if(outBest) *outBest = bestV;
        return bestV;
    };
    // Seed 0: best density orientation per item.
    vector<int> combo(M, 0);
    if(allowRot){
        for(int t = 0; t < M; ++t){
            const ItemType& it = g_items[t];
            // choose rot maximizing per-shelf density (perShelf / height)
            double d0 = sd[t][0].valid ? (double)sd[t][0].perShelf * it.v / (double)sd[t][0].rh : -1;
            double d1 = sd[t][1].valid ? (double)sd[t][1].perShelf * it.v / (double)sd[t][1].rh : -1;
            combo[t] = (d1 > d0) ? 1 : 0;
        }
    }
    bestVal = evalCombo(combo, nullptr);
    bestCombo = combo;
    if(allowRot){
        // Hill-climb: repeatedly try flipping one item; accept improvements.
        // Use steepest-ascent: try all single flips per pass, pick the best.
        bool improved = true;
        int passes = 0;
        while(improved && passes < 5 && elapsed() < TIME_LIMIT * 0.30){
            improved = false; ++passes;
            int bestFlip = -1; ll bestFlipVal = bestVal;
            for(int t = 0; t < M && elapsed() < TIME_LIMIT * 0.30; ++t){
                if(!sd[t][0].valid || !sd[t][1].valid) continue;
                combo[t] ^= 1;
                ll v = evalCombo(combo, nullptr);
                if(v > bestFlipVal){
                    bestFlipVal = v; bestFlip = t;
                }
                combo[t] ^= 1; // revert for next trial
            }
            if(bestFlip >= 0){
                combo[bestFlip] ^= 1;
                bestVal = bestFlipVal; bestCombo = combo; improved = true;
            }
        }
    }
    vector<KItem> kitems;
    kitems.reserve(M * 12);
    for(int t = 0; t < M; ++t){
        int rot = bestCombo[t];
        if(!sd[t][rot].valid) continue;
        const ItemType& it = g_items[t];
        int perShelf = sd[t][rot].perShelf;
        int rh = sd[t][rot].rh;
        int full = sd[t][rot].full;
        int rem = sd[t][rot].rem;
        int k = full, p = 1;
        while(k > 0){
            int take = min(p, k);
            kitems.push_back({rh * take, (ll)perShelf * take * it.v, t, rot, take, perShelf, false});
            k -= take; p *= 2;
        }
        if(rem > 0){
            kitems.push_back({rh, (ll)rem * it.v, t, rot, 1, rem, true});
        }
    }
    fill(dp.begin(), dp.end(), -1);
    dp[0] = 0;
    fill(par.begin(), par.end(), -1);
    int K = (int)kitems.size();
    vector<vector<ll>> dpBefore(K + 1, vector<ll>(H + 1, -1));
    for(int x = 0; x <= H; ++x) dp[x] = (x == 0) ? 0 : -1;
    for(int idx = 0; idx < K; ++idx){
        for(int x = 0; x <= H; ++x) dpBefore[idx][x] = dp[x];
        const KItem& ki = kitems[idx];
        for(int x = H; x >= ki.w; --x){
            if(dp[x - ki.w] >= 0 && dp[x - ki.w] + ki.val > dp[x]){
                dp[x] = dp[x - ki.w] + ki.val;
            }
        }
    }
    for(int x = 0; x <= H; ++x) dpBefore[K][x] = dp[x];
    ll bestV = 0; int bestX = 0;
    for(int x = 0; x <= H; ++x){ if(dp[x] > bestV){ bestV = dp[x]; bestX = x; } }
    vector<KItem> chosen;
    {
        int x = bestX;
        for(int idx = K - 1; idx >= 0; --idx){
            const KItem& ki = kitems[idx];
            ll before = dpBefore[idx][x];
            ll after = dpBefore[idx + 1][x];
            if(x >= ki.w && after > before && dpBefore[idx][x - ki.w] >= 0 && dpBefore[idx][x - ki.w] + ki.val == after){
                chosen.push_back(ki);
                x -= ki.w;
            }
        }
    }
    PackResult res;
    res.totalValue = 0;
    res.used.assign(g_items.size(), 0);
    int curY = 0;
    sort(chosen.begin(), chosen.end(), [](const KItem& a, const KItem& b){ return a.w > b.w; });
    for(const KItem& ki : chosen){
        const ItemType& it = g_items[ki.typeId];
        int rw = ki.rot ? it.h : it.w;
        int rh = ki.rot ? it.w : it.h;
        int perShelf = ki.perShelf;
        int copiesPerShelf = perShelf;
        int totalShelves = ki.nShelves;
        int canPlace = it.limit - res.used[ki.typeId];
        if(canPlace <= 0) continue;
        for(int s = 0; s < totalShelves && canPlace > 0; ++s){
            if(curY + rh > H) break;
            int curX = 0;
            int thisShelf = min(copiesPerShelf, canPlace);
            for(int c = 0; c < thisShelf; ++c){
                Placed p; p.typeId = ki.typeId; p.x = curX; p.y = curY; p.rot = ki.rot;
                res.placements.push_back(p);
                res.totalValue += it.v;
                res.used[ki.typeId]++;
                curX += rw;
            }
            canPlace -= thisShelf;
            curY += rh;
        }
    }
    return res;
}

static PackResult knapsackColumnPlan(bool allowRot){
    swap(g_bin.W, g_bin.H);
    for(auto& it : g_items) swap(it.w, it.h);
    PackResult r = knapsackShelfPlan(allowRot);
    for(auto& it : g_items) swap(it.w, it.h);
    swap(g_bin.W, g_bin.H);
    for(auto& p : r.placements) swap(p.x, p.y); // rot and type unchanged (see derivation)
    return r;
}

struct MSUnit {
    int typeId, rot;
    int wOr, hOr;    // oriented item size
    int perCol;      // items per full column
    int nItems;      // total items in this unit
    int nCols;       // columns in this unit (width = wOr*nCols)
    ll val;
};

static void buildShelfUnits(int Hs, int W, const vector<int>& rem, bool allowRot,
                            uint32_t flipMask, vector<MSUnit>& units){
    units.clear();
    int M = (int)g_items.size();
    for(int t = 0; t < M; ++t){
        if(rem[t] <= 0) continue;
        const ItemType& it = g_items[t];
        double bestD = -1, secD = -1;
        int bw=0,bh=0,bk=0,brot=-1, sw=0,sh=0,sk=0,srot=-1;
        for(int rot = 0; rot <= (allowRot ? 1 : 0); ++rot){
            int w_or = rot ? it.h : it.w;
            int h_or = rot ? it.w : it.h;
            if(h_or > Hs || w_or > W) continue;
            if(rot == 1 && it.w == it.h) continue;
            int k = Hs / h_or; if(k > rem[t]) k = rem[t];
            double d = (double)k * (double)it.v / (double)w_or;
            if(d > bestD){
                secD=bestD; sw=bw; sh=bh; sk=bk; srot=brot;
                bestD=d; bw=w_or; bh=h_or; bk=k; brot=rot;
            } else if(d > secD){ secD=d; sw=w_or; sh=h_or; sk=k; srot=rot; }
        }
        if(brot < 0) continue;
        if(((flipMask >> t) & 1u) && srot >= 0){ bw=sw; bh=sh; bk=sk; brot=srot; }
        int maxCols = W / bw;
        if(maxCols <= 0) continue;
        int fullCols = rem[t] / bk;
        bool capped = false;
        if(fullCols > maxCols){ fullCols = maxCols; capped = true; }
        int k2 = fullCols, p = 1;
        while(k2 > 0){
            int take = min(p, k2);
            units.push_back({t, brot, bw, bh, bk, bk*take, take, (ll)bk*take*it.v});
            k2 -= take; p *= 2;
        }
        if(!capped){
            int partial = rem[t] - fullCols * bk;
            if(partial > 0 && fullCols < maxCols)
                units.push_back({t, brot, bw, bh, bk, partial, 1, (ll)partial*it.v});
        }
    }
}

static ll msKnapValue(const vector<MSUnit>& units, int W){
    static vector<ll> dp;
    dp.assign(W + 1, 0);
    for(const MSUnit& u : units){
        int uw = u.wOr * u.nCols;
        if(uw > W || uw <= 0) continue;
        ll uv = u.val;
        for(int c = W; c >= uw; --c){
            ll v = dp[c - uw] + uv;
            if(v > dp[c]) dp[c] = v;
        }
    }
    return dp[W];
}

static ll msKnapChoose(const vector<MSUnit>& units, int W, vector<int>& chosen){
    int U = (int)units.size();
    static vector<ll> snap;
    snap.assign((size_t)(U + 1) * (W + 1), 0);
    for(int u = 0; u < U; ++u){
        const ll* prev = &snap[(size_t)u * (W + 1)];
        ll* cur = &snap[(size_t)(u + 1) * (W + 1)];
        int uw = units[u].wOr * units[u].nCols;
        ll uv = units[u].val;
        for(int c = 0; c <= W; ++c){
            ll b = prev[c];
            if(uw > 0 && uw <= W && c >= uw){
                ll v = prev[c - uw] + uv;
                if(v > b) b = v;
            }
            cur[c] = b;
        }
    }
    chosen.clear();
    int c = W;
    for(int u = U - 1; u >= 0; --u){
        const ll* prev = &snap[(size_t)u * (W + 1)];
        const ll* cur = &snap[(size_t)(u + 1) * (W + 1)];
        if(cur[c] != prev[c]){
            chosen.push_back(u);
            c -= units[u].wOr * units[u].nCols;
        }
    }
    return snap[(size_t)U * (W + 1) + W];
}

static double g_msAlpha = 1.0;

static void mixedShelfFillRegion(int x0, int y0, int RW, int RH, vector<int>& rem,
                                 bool allowRot, uint32_t seed, PackResult& res, mt19937& rng){
    int M = (int)g_items.size();
    int y = 0;
    vector<MSUnit> units;
    vector<int> cands, chosen;
    while(y < RH){
        if(elapsed() > TIME_LIMIT) break;
        int Hrem = RH - y;
        cands.clear();
        for(int t = 0; t < M; ++t){
            if(rem[t] <= 0) continue;
            const ItemType& it = g_items[t];
            for(int rot = 0; rot <= (allowRot ? 1 : 0); ++rot){
                int h_or = rot ? it.w : it.h;
                int w_or = rot ? it.h : it.w;
                if(w_or > RW || h_or > Hrem) continue;
                int kmax = Hrem / h_or;
                if(kmax > rem[t]) kmax = rem[t];
                if(kmax > 6) kmax = 6;
                for(int k = 1; k <= kmax; ++k) cands.push_back(k * h_or);
            }
        }
        cands.push_back(Hrem);
        sort(cands.begin(), cands.end());
        cands.erase(unique(cands.begin(), cands.end()), cands.end());
        if(cands.empty()) break;
        uint32_t flipMask = 0;
        if(seed && allowRot){
            flipMask = rng() & rng() & rng(); // ~12.5% of bits
            flipMask &= (M >= 32) ? 0xFFFFFFFFu : ((1u << M) - 1u);
        }
        double bestScore = -1; int bestH = -1;
        for(int Hs : cands){
            buildShelfUnits(Hs, RW, rem, allowRot, flipMask, units);
            if(units.empty()) continue;
            ll v = msKnapValue(units, RW);
            if(v <= 0) continue;
            double score = (g_msAlpha == 1.0) ? (double)v / (double)Hs
                                              : (double)v / pow((double)Hs, g_msAlpha);
            if(seed) score *= 1.0 + 0.05 * (((double)(rng() % 2001) / 1000.0) - 1.0);
            if(score > bestScore){ bestScore = score; bestH = Hs; }
        }
        if(bestH < 0) break;
        buildShelfUnits(bestH, RW, rem, allowRot, flipMask, units);
        msKnapChoose(units, RW, chosen);
        if(chosen.empty()) break;
        int x = 0;
        for(int idx : chosen){
            const MSUnit& u = units[idx];
            int itemsLeft = u.nItems;
            for(int cc = 0; cc < u.nCols && itemsLeft > 0; ++cc){
                int inCol = min(u.perCol, itemsLeft);
                for(int j = 0; j < inCol; ++j){
                    Placed p; p.typeId = u.typeId; p.x = x0 + x; p.y = y0 + y + j * u.hOr; p.rot = u.rot;
                    res.placements.push_back(p);
                    res.totalValue += g_items[u.typeId].v;
                    res.used[u.typeId]++;
                }
                itemsLeft -= inCol;
                x += u.wOr;
            }
            rem[u.typeId] -= u.nItems;
        }
        y += bestH;
    }
}

static void mixedShelfFillRegionT(int x0, int y0, int RW, int RH, vector<int>& rem,
                                  bool allowRot, uint32_t seed, PackResult& res, mt19937& rng){
    size_t before = res.placements.size();
    for(auto& it : g_items) swap(it.w, it.h);
    mixedShelfFillRegion(0, 0, RH, RW, rem, allowRot, seed, res, rng);
    for(auto& it : g_items) swap(it.w, it.h);
    for(size_t i = before; i < res.placements.size(); ++i){
        Placed& p = res.placements[i];
        int nx = x0 + p.y, ny = y0 + p.x;
        p.x = nx; p.y = ny;
    }
}

static PackResult mixedShelfPlan(bool allowRot, uint32_t seed, int stripW, bool transposedMain){
    PackResult res;
    res.totalValue = 0;
    res.used.assign(g_items.size(), 0);
    int M = (int)g_items.size();
    vector<int> rem(M);
    for(int i = 0; i < M; ++i) rem[i] = g_items[i].limit;
    mt19937 rng(seed * 2654435761u + 12345u);
    int W = g_bin.W, H = g_bin.H;
    if(!transposedMain){
        if(stripW >= W) stripW = 0;
        mixedShelfFillRegion(0, 0, W - stripW, H, rem, allowRot, seed, res, rng);
        if(stripW > 0)
            mixedShelfFillRegionT(W - stripW, 0, stripW, H, rem, allowRot, seed, res, rng);
    } else {
        if(stripW >= H) stripW = 0;
        mixedShelfFillRegionT(0, 0, W, H - stripW, rem, allowRot, seed, res, rng);
        if(stripW > 0)
            mixedShelfFillRegion(0, H - stripW, W, stripW, rem, allowRot, seed, res, rng);
    }
    return res;
}

static PackResult splitMixedPlan(bool allowRot, int sw, int mask, uint32_t seed){
    PackResult res; res.totalValue = 0; res.used.assign(g_items.size(), 0);
    int W = g_bin.W, H = g_bin.H, M = (int)g_items.size();
    if(sw <= 0 || sw >= W) return res;
    vector<int> rem(M);
    for(int i = 0; i < M; ++i) rem[i] = g_items[i].limit;
    mt19937 rng(seed * 1103515245u + 97u);
    auto fillR = [&](int x,int w,bool tr){
        if(tr) mixedShelfFillRegionT(x, 0, w, H, rem, allowRot, seed, res, rng);
        else mixedShelfFillRegion(x, 0, w, H, rem, allowRot, seed, res, rng);
    };
    if(mask & 4){ fillR(sw, W - sw, mask & 2); fillR(0, sw, mask & 1); }
    else { fillR(0, sw, mask & 1); fillR(sw, W - sw, mask & 2); }
    return res;
}

static PackResult splitMixedPlanY(bool allowRot, int sh, int mask, uint32_t seed){
    PackResult res; res.totalValue = 0; res.used.assign(g_items.size(), 0);
    int W = g_bin.W, H = g_bin.H, M = (int)g_items.size();
    if(sh <= 0 || sh >= H) return res;
    vector<int> rem(M);
    for(int i = 0; i < M; ++i) rem[i] = g_items[i].limit;
    mt19937 rng(seed * 747796405u + 2891336453u);
    auto fillR = [&](int y,int h,bool tr){
        if(tr) mixedShelfFillRegionT(0, y, W, h, rem, allowRot, seed, res, rng);
        else mixedShelfFillRegion(0, y, W, h, rem, allowRot, seed, res, rng);
    };
    if(mask & 4){ fillR(sh, H - sh, mask & 2); fillR(0, sh, mask & 1); }
    else { fillR(0, sh, mask & 1); fillR(sh, H - sh, mask & 2); }
    return res;
}

static vector<int> g_densOrd; // density-descending type order (filled in main)

static double msFracBound(const vector<int>& rem, double area){
    double v = 0;
    for(int t : g_densOrd){
        if(rem[t] <= 0) continue;
        double a = (double)g_items[t].area * rem[t];
        if(a >= area){ v += area * g_items[t].density; return v; }
        v += (double)g_items[t].v * rem[t];
        area -= a;
    }
    return v;
}

struct MSState {
    vector<int> rem;
    int y;
    PackResult res;
    double rank;
};

static void msCollectCands(int RW, int Hrem, const vector<int>& rem, bool allowRot,
                           vector<int>& cands){
    int M = (int)g_items.size();
    cands.clear();
    for(int t = 0; t < M; ++t){
        if(rem[t] <= 0) continue;
        const ItemType& it = g_items[t];
        for(int rot = 0; rot <= (allowRot ? 1 : 0); ++rot){
            int h_or = rot ? it.w : it.h;
            int w_or = rot ? it.h : it.w;
            if(w_or > RW || h_or > Hrem) continue;
            int kmax = Hrem / h_or;
            if(kmax > rem[t]) kmax = rem[t];
            int cap = allowRot ? 5 : 7;
            if(kmax > cap) kmax = cap;
            for(int k = 1; k <= kmax; ++k) cands.push_back(k * h_or);
        }
    }
    cands.push_back(Hrem);
    sort(cands.begin(), cands.end());
    cands.erase(unique(cands.begin(), cands.end()), cands.end());
}

static void msCommitShelf(MSState& st, int Hs, int RW, bool allowRot,
                          vector<MSUnit>& units, vector<int>& chosen){
    buildShelfUnits(Hs, RW, st.rem, allowRot, 0, units);
    msKnapChoose(units, RW, chosen);
    int x = 0;
    for(int idx : chosen){
        const MSUnit& u = units[idx];
        int itemsLeft = u.nItems;
        for(int cc = 0; cc < u.nCols && itemsLeft > 0; ++cc){
            int inCol = min(u.perCol, itemsLeft);
            for(int j = 0; j < inCol; ++j){
                Placed p; p.typeId = u.typeId; p.x = x; p.y = st.y + j * u.hOr; p.rot = u.rot;
                st.res.placements.push_back(p);
                st.res.totalValue += g_items[u.typeId].v;
                st.res.used[u.typeId]++;
            }
            itemsLeft -= inCol;
            x += u.wOr;
        }
        st.rem[u.typeId] -= u.nItems;
    }
    st.y += Hs;
}

static PackResult beamMixedShelfPlan(bool allowRot, int beamW, int branch,
                                     int branchDepth, double tlimit){
    int W = g_bin.W, H = g_bin.H, M = (int)g_items.size();
    vector<MSState> beam(1);
    beam[0].rem.resize(M);
    for(int i = 0; i < M; ++i) beam[0].rem[i] = g_items[i].limit;
    beam[0].y = 0;
    beam[0].res.totalValue = 0;
    beam[0].res.used.assign(M, 0);
    PackResult best; best.totalValue = -1;
    vector<MSUnit> units; vector<int> cands, chosen;
    int depth = 0;
    while(!beam.empty()){
        int nBranch = (depth < branchDepth) ? branch : 1;
        vector<MSState> next;
        for(MSState& st : beam){
            if(elapsed() > tlimit){
                if(st.res.totalValue > best.totalValue) best = st.res;
                continue;
            }
            int Hrem = H - st.y;
            bool expanded = false;
            if(Hrem > 0){
                msCollectCands(W, Hrem, st.rem, allowRot, cands);
                static vector<pair<double,int>> scored;
                scored.clear();
                for(int Hs : cands){
                    buildShelfUnits(Hs, W, st.rem, allowRot, 0, units);
                    if(units.empty()) continue;
                    ll v = msKnapValue(units, W);
                    if(v <= 0) continue;
                    scored.push_back({(g_msAlpha == 1.0 ? (double)v / (double)Hs : (double)v / pow((double)Hs, g_msAlpha)), Hs});
                }
                sort(scored.begin(), scored.end(),
                     [](const pair<double,int>& a, const pair<double,int>& b){ return a.first > b.first; });
                int taken = 0;
                static vector<int> takenH;
                takenH.clear();
                for(size_t si = 0; si < scored.size() && taken < nBranch; ++si){
                    int Hs = scored[si].second;
                    bool dup = false;
                    for(int th : takenH) if(abs(th - Hs) <= 2){ dup = true; break; }
                    if(dup) continue;
                    takenH.push_back(Hs);
                    MSState child = st;
                    msCommitShelf(child, Hs, W, allowRot, units, chosen);
                    child.rank = (double)child.res.totalValue +
                                 msFracBound(child.rem, (double)(H - child.y) * W);
                    next.push_back(std::move(child));
                    ++taken;
                    expanded = true;
                }
            }
            if(!expanded){
                if(st.res.totalValue > best.totalValue) best = st.res;
            }
        }
        if(next.empty()) break;
        sort(next.begin(), next.end(),
             [](const MSState& a, const MSState& b){ return a.rank > b.rank; });
        vector<MSState> pruned;
        for(auto& s : next){
            bool same = false;
            for(auto& q : pruned)
                if(q.y == s.y && q.res.totalValue == s.res.totalValue){ same = true; break; }
            if(!same) pruned.push_back(std::move(s));
            if((int)pruned.size() >= beamW) break;
        }
        beam = std::move(pruned);
        ++depth;
    }
    if(best.totalValue < 0){ best.totalValue = 0; best.used.assign(M, 0); }
    return best;
}

static PackResult beamMixedShelfPlanT(bool allowRot, int beamW, int branch,
                                      int branchDepth, double tlimit){
    swap(g_bin.W, g_bin.H);
    for(auto& it : g_items) swap(it.w, it.h);
    PackResult r = beamMixedShelfPlan(allowRot, beamW, branch, branchDepth, tlimit);
    for(auto& it : g_items) swap(it.w, it.h);
    swap(g_bin.W, g_bin.H);
    for(auto& p : r.placements) swap(p.x, p.y);
    return r;
}

static vector<int> orderByDensity(){
    vector<int> idx(g_items.size());
    for(size_t i=0;i<g_items.size();++i) idx[i]=i;
    sort(idx.begin(), idx.end(), [](int a, int b){
        if(g_items[a].density != g_items[b].density) return g_items[a].density > g_items[b].density;
        return g_items[a].v > g_items[b].v;
    });
    return idx;
}
static vector<int> orderByHeightDesc(){
    vector<int> idx(g_items.size());
    for(size_t i=0;i<g_items.size();++i) idx[i]=i;
    sort(idx.begin(), idx.end(), [](int a, int b){
        return g_items[a].h > g_items[b].h;
    });
    return idx;
}
static vector<int> orderByWidthDesc(){
    vector<int> idx(g_items.size());
    for(size_t i=0;i<g_items.size();++i) idx[i]=i;
    sort(idx.begin(), idx.end(), [](int a, int b){
        return g_items[a].w > g_items[b].w;
    });
    return idx;
}
static vector<int> orderByAreaDesc(){
    vector<int> idx(g_items.size());
    for(size_t i=0;i<g_items.size();++i) idx[i]=i;
    sort(idx.begin(), idx.end(), [](int a, int b){
        return g_items[a].area > g_items[b].area;
    });
    return idx;
}
static vector<int> orderByValueDesc(){
    vector<int> idx(g_items.size());
    for(size_t i=0;i<g_items.size();++i) idx[i]=i;
    sort(idx.begin(), idx.end(), [](int a, int b){
        return g_items[a].v > g_items[b].v;
    });
    return idx;
}
static vector<int> orderByMinDimAsc(){
    vector<int> idx(g_items.size());
    for(size_t i=0;i<g_items.size();++i) idx[i]=i;
    sort(idx.begin(), idx.end(), [](int a, int b){
        int ma = min(g_items[a].w, g_items[a].h);
        int mb = min(g_items[b].w, g_items[b].h);
        return ma < mb;
    });
    return idx;
}

static void outputResult(const PackResult& res){
    string out;
    out.reserve(res.placements.size() * 40 + 32);
    out += "{\"placements\":[";
    bool first = true;
    for(auto& p : res.placements){
        if(!first) out += ',';
        first = false;
        out += "{\"type\":\"";
        out += g_items[p.typeId].name;
        out += "\",\"x\":";
        out += to_string(p.x);
        out += ",\"y\":";
        out += to_string(p.y);
        out += ",\"rot\":";
        out += to_string(p.rot);
        out += "}";
    }
    out += "]}";
    fputs(out.c_str(), stdout);
    fputc('\n', stdout);
}

int main(){
    T_START = chrono::steady_clock::now();
    string input;
    {
        char buf[1<<16];
        size_t n;
        while((n = fread(buf,1,sizeof(buf),stdin)) > 0) input.append(buf, n);
    }
    JsonParser jp(input);
    jp.expectChar('{');
    for(int iter=0; iter<2; ++iter){
        if(iter>0) jp.expectChar(',');
        string k = jp.readString();
        jp.expectChar(':');
        if(k == "bin"){
            jp.expectChar('{');
            for(int i=0;i<3;++i){
                if(i>0) jp.expectChar(',');
                string bk = jp.readString();
                jp.expectChar(':');
                if(bk=="W") g_bin.W = (int)jp.readLong();
                else if(bk=="H") g_bin.H = (int)jp.readLong();
                else if(bk=="allow_rotate") g_bin.allow_rotate = jp.readBool();
            }
            jp.expectChar('}');
        } else if(k == "items"){
            jp.expectChar('[');
            bool first=true;
            while(true){
                jp.skipWS();
                if(jp.peek()==']'){ jp.getc(); break; }
                if(!first) jp.expectChar(',');
                first=false;
                jp.expectChar('{');
                ItemType it; it.w=it.h=0; it.v=0; it.limit=0;
                for(int f=0; f<5; ++f){
                    if(f>0) jp.expectChar(',');
                    string key = jp.readString();
                    jp.expectChar(':');
                    if(key=="type") it.name = jp.readString();
                    else if(key=="w") it.w = (int)jp.readLong();
                    else if(key=="h") it.h = (int)jp.readLong();
                    else if(key=="v") it.v = jp.readLong();
                    else if(key=="limit") it.limit = (int)jp.readLong();
                }
                jp.expectChar('}');
                it.area = (ll)it.w * it.h;
                it.density = (double)it.v / (double)it.area;
                g_items.push_back(it);
            }
        }
    }
    jp.expectChar('}');

    int M = (int)g_items.size();
    bool allowRot = g_bin.allow_rotate;
    (void)M;

    PackResult best;
    best.totalValue = -1;

    vector<vector<int>> orders = {
        orderByDensity(),
        orderByHeightDesc(),
        orderByWidthDesc(),
        orderByAreaDesc(),
        orderByValueDesc(),
        orderByMinDimAsc()
    };

#ifdef DIAG
    const char* g_label = "";
    auto consider = [&](PackResult&& r){
        if(r.totalValue > best.totalValue){ best = std::move(r); fprintf(stderr,"[%s] newbest=%lld @%.3fs\n", g_label, best.totalValue, elapsed()); }
    };
#else
    auto consider = [&](PackResult&& r){
        if(r.totalValue > best.totalValue) best = std::move(r);
    };
#endif

    #ifdef DIAG
    g_label="knap";
#endif
    consider(knapsackShelfPlan(allowRot));

    vector<int> ordMinDim = orderByMinDimAsc();
    vector<int> ordDens   = orderByDensity();
    vector<int> ordVal    = orderByValueDesc();
    vector<int> ordAreaAsc = ordMinDim; // area-ascending
    sort(ordAreaAsc.begin(), ordAreaAsc.end(), [](int a,int b){ return g_items[a].area < g_items[b].area; });
    vector<vector<int>*> gapOrders = {&ordMinDim, &ordAreaAsc, &ordDens, &ordVal};
    auto polish = [&](PackResult r) -> PackResult {
        ll prev = -1;
        int rounds = 0;
        while(r.totalValue != prev && rounds < 3){
            prev = r.totalValue; ++rounds;
            for(auto* op : gapOrders){
                if(elapsed() > TIME_LIMIT) return r;
                PackResult f = fillGaps(r, *op, allowRot, allowRot);
                if(f.totalValue > r.totalValue) r = std::move(f);
            }
        }
        return r;
    };
    auto gapFill = [&](){
#ifdef DIAG
        g_label="gapfill";
#endif
        ll prev = -1;
        for(int round = 0; round < 5 && best.totalValue != prev; ++round){
            prev = best.totalValue;
            for(auto* op : gapOrders){
                consider(fillGaps(best, *op, allowRot, allowRot));
                if(elapsed() > TIME_LIMIT) return;
            }
        }
    };
    gapFill();
    if(!allowRot){
        consider(polish(knapsackColumnPlan(allowRot)));
        gapFill();
    }

#ifdef DIAG
    g_label="mixed";
#endif
    consider(polish(mixedShelfPlan(allowRot, 0, 0, false)));
    consider(polish(mixedShelfPlan(allowRot, 0, 0, true)));
    {
        double alphas[1] = {0.94};
        for(double a : alphas){
            if(elapsed() > TIME_LIMIT * 0.22) break;
            g_msAlpha = a;
            consider(polish(mixedShelfPlan(allowRot, 0, 0, false)));
            if(elapsed() > TIME_LIMIT * 0.22) break;
            consider(polish(mixedShelfPlan(allowRot, 0, 0, true)));
        }
        g_msAlpha = 1.0;
    }
    // Priority path: tall no-rotation split2 with alpha 0.94 early (protects c11).
    if(!allowRot && g_bin.H * 5 > g_bin.W * 6 && elapsed() < TIME_LIMIT * 0.30){
        double oldA = g_msAlpha; g_msAlpha = 0.94;
        auto cutsE = [&](int L){
            vector<int> v; auto add=[&](int x){ if(x>20&&x<L-20&&find(v.begin(),v.end(),x)==v.end()) v.push_back(x); };
            add(L/3); add(L/2); add((2*L)/3); add(L/4); add((3*L)/4);
            for(int z=0; z<M && z<3; ++z){ const ItemType& it=g_items[ordDens[z]]; int d[2]={it.w,it.h};
                for(int q=0;q<2;++q) for(int k=1;k<=3;++k){ add(d[q]*k); add(L-d[q]*k);} }
            if((int)v.size()>10) v.resize(10); return v;
        };
        int pe=0;
        for(int sh:cutsE(g_bin.H))for(int mask=0;mask<8&&pe<80&&elapsed()<TIME_LIMIT*0.55;++mask,++pe)
            consider(polish(splitMixedPlanY(allowRot,sh,mask,0)));
        pe=0;
        for(int sw:cutsE(g_bin.W))for(int mask=0;mask<8&&pe<80&&elapsed()<TIME_LIMIT*0.75;++mask,++pe)
            consider(polish(splitMixedPlan(allowRot,sw,mask,0)));
        g_msAlpha = oldA;
    }
#ifdef DIAG
    g_label="beam";
#endif
    g_densOrd = orderByDensity();
#ifdef DIAG
    {
        int bw = allowRot ? 14 : 7, br = allowRot ? 5 : 4, bd = allowRot ? 9 : 7;
        PackResult b1 = beamMixedShelfPlan(allowRot, bw, br, bd, allowRot ? TIME_LIMIT * 0.78 : TIME_LIMIT * 0.66);
        fprintf(stderr, "[beamH] raw=%lld @%.3fs\n", b1.totalValue, elapsed());
        consider(polish(std::move(b1)));
        g_msAlpha = allowRot ? 0.935 : 1.0;
        PackResult b2 = beamMixedShelfPlanT(allowRot, bw, br, bd, allowRot ? TIME_LIMIT * 0.9 : TIME_LIMIT * 0.82);
        fprintf(stderr, "[beamT] raw=%lld @%.3fs\n", b2.totalValue, elapsed());
        consider(polish(std::move(b2)));
        g_msAlpha = 1.0;
    }
#else
    if(elapsed() < TIME_LIMIT * 0.55){ g_msAlpha = allowRot ? 0.94 : 1.0;
        consider(polish(beamMixedShelfPlan(allowRot, allowRot ? 14 : 7, allowRot ? 5 : 4, allowRot ? 9 : 7, allowRot ? TIME_LIMIT * 0.78 : TIME_LIMIT * 0.66))); g_msAlpha = 1.0; }
    if(elapsed() < TIME_LIMIT * 0.7){
        g_msAlpha = allowRot ? 0.935 : 1.0;
        consider(polish(beamMixedShelfPlanT(allowRot, allowRot ? 22 : 7, allowRot ? 5 : 4, allowRot ? 9 : 7, allowRot ? TIME_LIMIT * 0.9 : TIME_LIMIT * 0.82)));
        g_msAlpha = 1.0;
    }
    // Focused beam with alpha=0.94 and reduced depth for shorter shelves
    if(allowRot && elapsed() < TIME_LIMIT * 0.80){
        g_msAlpha = 0.94;
        consider(polish(beamMixedShelfPlan(allowRot, 14, 5, 5, TIME_LIMIT * 0.80)));
        g_msAlpha = 1.0;
    }
    // Wide-and-shallow beam for broad exploration
    if(allowRot && elapsed() < TIME_LIMIT * 0.85){
        consider(polish(beamMixedShelfPlan(allowRot, 22, 8, 3, TIME_LIMIT * 0.85)));
    }
#endif
#ifdef DIAG
    g_label="hybrid";
#endif
    {
        vector<int> stripCands;
        for(const ItemType& it : g_items){
            if(it.h >= 2 * it.w && it.w <= g_bin.W / 6) stripCands.push_back(it.w);
            if((allowRot || true) && it.w >= 2 * it.h && it.h <= g_bin.W / 6) stripCands.push_back(it.h);
        }
        sort(stripCands.begin(), stripCands.end());
        stripCands.erase(unique(stripCands.begin(), stripCands.end()), stripCands.end());
        // Add standard fraction strips for more diverse split points
        stripCands.push_back(g_bin.W / 3);
        stripCands.push_back(g_bin.W / 2);
        stripCands.push_back((2 * g_bin.W) / 3);
        sort(stripCands.begin(), stripCands.end());
        stripCands.erase(unique(stripCands.begin(), stripCands.end()), stripCands.end());
        {
            int n = min((int)stripCands.size(), 3);
            vector<int> sums;
            for(int i = 0; i < n; ++i)
                for(int j = i; j < n; ++j) sums.push_back(stripCands[i] + stripCands[j]);
            for(int s : sums) stripCands.push_back(s);
            sort(stripCands.begin(), stripCands.end());
            stripCands.erase(unique(stripCands.begin(), stripCands.end()), stripCands.end());
        }
        if((int)stripCands.size() > (allowRot ? 8 : 10)) stripCands.resize(allowRot ? 8 : 10);
        for(int d : stripCands){
            if(elapsed() > TIME_LIMIT * 0.45) break;
            consider(polish(mixedShelfPlan(allowRot, 0, d, false)));
            if(elapsed() > TIME_LIMIT * 0.45) break;
            consider(polish(mixedShelfPlan(allowRot, 0, d, true)));
        }
    }
#ifdef DIAG
    g_label="split";
#endif
    {
        auto cuts = [&](int L){
            vector<int> v;
            auto add = [&](int x){ if(x > 20 && x < L - 20 && find(v.begin(), v.end(), x) == v.end()) v.push_back(x); };
            add(L/3); add(L/2); add((2*L)/3); add(L/4); add((3*L)/4); add((2*L)/5); add((3*L)/5);
            for(int z = 0; z < M && z < 4; ++z){
                const ItemType& it = g_items[ordDens[z]];
                int d[2] = {it.w, it.h};
                for(int q = 0; q < 2; ++q) for(int k = 1; k <= 3; ++k){ add(d[q] * k); add(L - d[q] * k); }
            }
            if((int)v.size() > 12) v.resize(12);
            return v;
        };
        vector<int> splits = cuts(g_bin.W);
        for(int sw : splits){
            for(int mask = 0; mask < 8 && elapsed() < TIME_LIMIT * 0.56; ++mask)
                consider(polish(splitMixedPlan(allowRot, sw, mask, 0)));
            if(elapsed() > TIME_LIMIT * 0.56) break;
        }
#ifdef DIAG
        g_label="splity";
#endif
        vector<int> ys = cuts(g_bin.H);
        for(int sh : ys){
            for(int mask = 0; mask < 8 && elapsed() < TIME_LIMIT * 0.60; ++mask)
                consider(polish(splitMixedPlanY(allowRot, sh, mask, 0)));
            if(elapsed() > TIME_LIMIT * 0.60) break;
        }
    }
#ifdef DIAG
    g_label="mixedrand";
#endif
    {
        uint32_t s = 1;
        while(elapsed() < TIME_LIMIT * 0.62 && s <= 30){
            consider(polish(mixedShelfPlan(allowRot, s, 0, false)));
            if(elapsed() > TIME_LIMIT * 0.62) break;
            consider(polish(mixedShelfPlan(allowRot, s, 0, true)));
            ++s;
        }
    }
#ifdef DIAG
    g_label="split2";
#endif
    {
        auto cuts = [&](int L, bool ycut){
            vector<int> v;
            auto add = [&](int x){ if(x > 20 && x < L - 20 && find(v.begin(), v.end(), x) == v.end()) v.push_back(x); };
            if(allowRot) for(int z = 0; z < M && z < 3; ++z){
                const ItemType& it = g_items[ordDens[z]];
                int d[2] = {it.w, it.h};
                for(int q = 0; q < 2; ++q) for(int k = 1; k <= 3; ++k){ add(d[q] * k); add(L - d[q] * k); }
            }
            if(ycut && !allowRot && g_bin.H * 5 > g_bin.W * 6) for(int z = 0; z < M && z < 4; ++z){
                const ItemType& it = g_items[ordDens[z]];
                int d[2] = {it.h, it.w};
                for(int q = 0; q < 2; ++q){
                    int km = L / d[q];
                    for(int k = 4; k <= km && k <= 24; ++k){
                        int x = d[q] * k;
                        if(x >= L / 4 && x <= (3 * L) / 4){ add(x); add(L - x); }
                    }
                }
            }
            add(L/3); add(L/2); add((2*L)/3); add(L/4); add((3*L)/4); add((2*L)/5); add((3*L)/5);
            for(int z = 0; !allowRot && z < M && z < 5; ++z){
                const ItemType& it = g_items[ordDens[z]];
                int d[2] = {it.w, it.h};
                for(int q = 0; q < 2; ++q) for(int k = 1; k <= 3; ++k){ add(d[q] * k); add(L - d[q] * k); }
            }
            if((int)v.size() > (ycut && !allowRot && g_bin.H * 5 > g_bin.W * 6 ? 18 : 14)) v.resize(ycut && !allowRot && g_bin.H * 5 > g_bin.W * 6 ? 18 : 14);
            return v;
        };
        vector<int> splits2 = cuts(g_bin.W, false);
        for(uint32_t seed = 1; seed <= 3 && elapsed() < TIME_LIMIT * 0.78; ++seed){
#ifdef DIAG
            g_label="split2";
#endif
            for(int sw : splits2){
                for(int mask = 0; mask < 8 && elapsed() < TIME_LIMIT * 0.78; ++mask)
                    consider(polish(splitMixedPlan(allowRot, sw, mask, seed)));
                if(elapsed() > TIME_LIMIT * 0.78) break;
            }
#ifdef DIAG
            g_label="split2y";
#endif
            vector<int> ys2 = cuts(g_bin.H, true);
            for(int sh : ys2){
                for(int mask = 0; mask < 8 && elapsed() < TIME_LIMIT * 0.78; ++mask)
                    consider(polish(splitMixedPlanY(allowRot, sh, mask, seed)));
                if(elapsed() > TIME_LIMIT * 0.78) break;
            }
        }
    }
#ifdef DIAG
    g_label="choicemr";
#endif
    { int ms[8]={3,6,7,0,4,5,1,2}; for(int ii=0;ii<8&&elapsed()<TIME_LIMIT*0.86;++ii) consider(polish(choiceMaxRects(allowRot,ms[ii],TIME_LIMIT*0.88))); }
#ifdef DIAG
    g_label="prune";
#endif
    if(elapsed() < TIME_LIMIT * 0.90){
        int cs[13] = {1,2,3,4,5,6,8,10,12,16,20,24,32};
        for(int c : cs){
            if(elapsed() > TIME_LIMIT * 0.93) break;
            consider(pruneLow(best, c, ordDens, allowRot));
            if(elapsed() > TIME_LIMIT * 0.93) break;
            consider(pruneLow(best, c, ordMinDim, allowRot));
            if(elapsed() > TIME_LIMIT * 0.93) break;
            consider(pruneLow(best, c, ordVal, allowRot));
            if(elapsed() > TIME_LIMIT * 0.93) break;
            consider(pruneLow(best, c, ordAreaAsc, allowRot));
        }
        int ds[3]={4,8,12};
        for(int c:ds)for(int side=0;side<4&&elapsed()<TIME_LIMIT*0.955;++side)consider(pruneSide(best,c,ordDens,allowRot,side));
        if(!allowRot&&g_bin.H*5>g_bin.W*6)for(int ln=1;ln<=2;++ln)for(int side=0;side<2&&elapsed()<TIME_LIMIT*0.965;++side){
            consider(pruneLine(best,ln,ordDens,allowRot,side));
            if(elapsed()<TIME_LIMIT*0.965) consider(pruneLine(best,ln,ordMinDim,allowRot,side));
            if(elapsed()<TIME_LIMIT*0.965) consider(pruneLine(best,ln,ordVal,allowRot,side));
        }
    }
    #ifdef DIAG
    g_label="maxrects";
#endif
    for(auto& ord : orders){
        consider(greedyFillMaxRects(ord, allowRot, false));
        if(allowRot) consider(greedyFillMaxRects(ord, allowRot, true));
        if(elapsed() > TIME_LIMIT * 0.5) break;
    }

    #ifdef DIAG
    g_label="skyline";
#endif
    for(auto& ord : orders){
        consider(greedyFill(ord, allowRot, false));
        if(allowRot) consider(greedyFill(ord, allowRot, true));
        if(elapsed() > TIME_LIMIT * 0.6) break;
    }

    #ifdef DIAG
    g_label="random";
#endif
    mt19937 rng(987654321u);
    int seed = 0;
    double iterCost = 0.0; // adaptive margin: last iteration's duration
    while(elapsed() + 2.0 * iterCost < TIME_LIMIT - 0.04){
        double t0 = elapsed();
        vector<int> ord = orderByDensity();
        int mode = seed % 8;
        if(mode == 0){
            shuffle(ord.begin(), ord.end(), rng);
        } else if(mode == 1){
            int i = rng() % M, j = rng() % M;
            swap(ord[i], ord[j]);
        } else if(mode == 2){
            sort(ord.begin(), ord.end(), [](int a, int b){
                double da = g_items[a].v / (double)g_items[a].h;
                double db = g_items[b].v / (double)g_items[b].h;
                return da > db;
            });
        } else if(mode == 3){
            sort(ord.begin(), ord.end(), [](int a, int b){
                double da = g_items[a].v / (double)g_items[a].w;
                double db = g_items[b].v / (double)g_items[b].w;
                return da > db;
            });
        } else if(mode == 4){
            sort(ord.begin(), ord.end(), [](int a, int b){
                int ma = min(g_items[a].w, g_items[a].h);
                int mb = min(g_items[b].w, g_items[b].h);
                if(ma != mb) return ma < mb;
                return g_items[a].density > g_items[b].density;
            });
        } else if(mode == 5){
            sort(ord.begin(), ord.end(), [](int a, int b){
                return g_items[a].density > g_items[b].density;
            });
            stable_sort(ord.begin(), ord.end(), [](int a, int b){
                int ma = min(g_items[a].w, g_items[a].h);
                int mb = min(g_items[b].w, g_items[b].h);
                return ma < mb;
            });
        } else if(mode == 6){
            sort(ord.begin(), ord.end(), [](int a, int b){
                return g_items[a].density > g_items[b].density;
            });
        } else {
            shuffle(ord.begin(), ord.end(), rng);
            sort(ord.begin(), ord.end(), [](int a, int b){
                return g_items[a].density > g_items[b].density;
            });
            for(int s=0;s<2;++s){ int i=rng()%M,j=rng()%M; swap(ord[i],ord[j]); }
        }
        bool tryBoth = allowRot ? ((seed & 1) == 0) : false;
        if(mode == 6) tryBoth = true; // force both
        consider(greedyFill(ord, allowRot, tryBoth));
        // Also try MaxRects with this ordering for additional diversity
        if(seed % 4 == 0) consider(greedyFillMaxRects(ord, allowRot, tryBoth));
        iterCost = elapsed() - t0;
        ++seed;
        if(seed > 2000000) break;
    }

    if(elapsed() < TIME_LIMIT - 0.04) gapFill();
    if(best.totalValue < 0){
        best.totalValue = 0;
        best.placements.clear();
    }
    outputResult(best);
    return 0;
}
