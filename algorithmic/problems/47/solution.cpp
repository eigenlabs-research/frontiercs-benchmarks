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

// ---------------- Timer ----------------
static chrono::steady_clock::time_point T_START;
static inline double elapsed(){
    return chrono::duration<double>(chrono::steady_clock::now() - T_START).count();
}
static double TIME_LIMIT = 0.96;

// ---------------- JSON parser (minimal, schema-specific) ----------------
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

// ---------------- Data ----------------
struct ItemType {
    string name;
    int w, h;
    ll v;
    int limit;
    // derived
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

// ---------------- Segment-based skyline packer (best-fit, bottom-left) ----------------
// Segments: list of (x, height) meaning skyline is at 'height' for [x, next_x).
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
    // Find best position for a rect rw x rh using bottom-left best-fit.
    // best = lowest y; tiebreak leftmost x; among candidate segment left edges.
    bool findBest(int rw, int rh, int& outX, int& outY) const {
        if(rw <= 0 || rh <= 0 || rw > W || rh > H) return false;
        int bestX = -1, bestY = INT32_MAX;
        int n = (int)segs.size() - 1; // last index is sentinel
        // For each segment, try left-aligning the rect at seg.x.
        for(int i = 0; i < n; ++i){
            int x = segs[i].x;
            if(x + rw > W) break;
            // The rect spans [x, x+rw). It must rest on the skyline over this interval.
            // Compute max height over [x, x+rw).
            int mx = 0;
            int j = i;
            int curX = x;
            // Walk segments covering [x, x+rw).
            // segs[i] covers [segs[i].x, segs[i+1].x).
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
    // Commit placement at (x,y) size rw x rh: raise skyline.
    void place(int x, int y, int rw, int rh){
        int x2 = x + rw;
        int newH = y + rh;
        // Build new segment list: replace coverage [x, x2) with height newH,
        // keeping segments outside unchanged, merging equal heights.
        vector<Seg> out;
        out.reserve(segs.size() + 4);
        int i = 0;
        // Copy segments fully left of x.
        while(i < (int)segs.size() && segs[i].x < x){
            out.push_back(segs[i]);
            ++i;
        }
        // Now segs[i].x >= x (or end). We may need to clip a segment that started before x.
        // Actually since we copy those with x < x, if segs[i].x == x we just proceed.
        // But a segment with segs[i-1].x < x and segs[i].x > x would have been partially covered;
        // need to re-add the left part. Handle: the last pushed seg has out.back().x < x; its height
        // extends to segs[i].x. We need to terminate it at x.
        if(!out.empty() && out.back().x < x){
            // keep the height but trim end to x by inserting x sentinel? Segments are (x, height) with
            // implied end = next seg x. To trim, we just insert a new seg at x with the new height later.
            // The out.back() height is valid for [out.back().x, x). Good — we leave it, and the next
            // seg we push will be at x.
        }
        // Insert new segment [x, x2) at newH (if same as previous height, merge).
        if(!out.empty() && out.back().height == newH){
            // merge: extend previous by not adding a new x marker; but we must mark x2 boundary.
            // Actually segments are defined by their start x and height; the end is the next seg's x.
            // If previous height == newH, the region [out.back().x, x2) is uniform; we just need
            // to ensure the next seg starts at x2. So skip adding x marker.
        } else {
            out.push_back({x, newH});
        }
        // Skip segments fully inside [x, x2).
        while(i < (int)segs.size() && segs[i].x < x2){
            ++i;
        }
        // Now segs[i].x >= x2 (or end). If segs[i].x > x2, we need to add a boundary at x2 with
        // the height of the segment that was at x2 (i.e., segs[i-1].height, the one we skipped past
        // that covered x2). Find height at x2 from original.
        int hAfter = 0;
        // height at position x2 in original skyline = height of segment containing x2.
        // Find it: the last segment with .x <= x2.
        {
            // binary or linear scan in original
            int lo = 0, hi = (int)segs.size() - 1;
            while(lo < hi){
                int mid = (lo + hi + 1) / 2;
                if(segs[mid].x <= x2) lo = mid; else hi = mid - 1;
            }
            hAfter = segs[lo].height;
        }
        // If the new region's height (newH) equals hAfter, the boundary at x2 disappears (merge).
        // We already pushed new seg at x with newH (or merged with previous). Now push x2 boundary
        // with hAfter only if newH != hAfter.
        if(newH != hAfter){
            out.push_back({x2, hAfter});
        }
        // Copy remaining segments with x > x2 (skip the one exactly at x2 if we added boundary).
        while(i < (int)segs.size()){
            if(segs[i].x > x2){
                out.push_back(segs[i]);
            } else if(segs[i].x == x2){
                // skip; we already added boundary at x2 with hAfter (or merged). But if we merged
                // (newH == hAfter) and didn't add boundary, this seg at x2 with hAfter is the continuation
                // — but its height equals newH, and we want to merge it away too. So skip.
                // If we did add boundary at x2, then this seg is duplicate -> skip.
            }
            ++i;
        }
        // Cleanup: remove trailing sentinel if missing, ensure last seg.x == W.
        // The original had segs.back().x == W. We should preserve that.
        // If out.back().x > W, trim; if out.back().x < W and last height extends to W, add W sentinel.
        // Actually the sentinel seg {W,0} was in original; if we skipped it (because W >= x2 typically
        // and W == x2 case), need to restore.
        if(out.empty() || out.back().x != W){
            // ensure we have the W boundary; the height at W is 0 (sentinel) — but actually the last
            // real segment's height extends to W. The sentinel {W,0} is a marker. Restore it.
            out.push_back({W, 0});
        }
        // Also remove consecutive duplicates (same height) that may have slipped through.
        // Compact.
        vector<Seg> compact;
        compact.reserve(out.size());
        for(auto& s : out){
            if(!compact.empty() && compact.back().height == s.height){
                // merge: skip adding (extends previous)
                continue;
            }
            compact.push_back(s);
        }
        // Ensure first seg at x=0.
        if(compact.empty() || compact[0].x != 0){
            compact.insert(compact.begin(), {0, 0});
        }
        // Ensure last at W.
        if(compact.back().x != W){
            compact.push_back({W, 0});
        }
        segs = std::move(compact);
    }
};

// ---------------- MaxRects packer (free-rectangle list, best-short-side-fit) ----------------
struct MaxRects {
    int W, H;
    // free rectangles stored as x,y,w,h
    vector<int> fx, fy, fw, fh;
    void init(int W_, int H_){
        W = W_; H = H_;
        fx.clear(); fy.clear(); fw.clear(); fh.clear();
        fx.push_back(0); fy.push_back(0); fw.push_back(W_); fh.push_back(H_);
    }
    // Best-short-side-fit: minimize the shorter leftover, tiebreak by longer leftover.
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
        // collect new fragments
        static vector<int> ax, ay, aw, ah;
        ax.clear(); ay.clear(); aw.clear(); ah.clear();
        for(int i = 0; i < n; ){
            int gx1 = fx[i], gy1 = fy[i], gx2 = fx[i] + fw[i], gy2 = fy[i] + fh[i];
            if(px1 < gx2 && px2 > gx1 && py1 < gy2 && py2 > gy1){
                if(px1 > gx1){ ax.push_back(gx1); ay.push_back(gy1); aw.push_back(px1 - gx1); ah.push_back(gy2 - gy1); }
                if(px2 < gx2){ ax.push_back(px2); ay.push_back(gy1); aw.push_back(gx2 - px2); ah.push_back(gy2 - gy1); }
                if(py1 > gy1){ ax.push_back(gx1); ay.push_back(gy1); aw.push_back(gx2 - gx1); ah.push_back(py1 - gy1); }
                if(py2 < gy2){ ax.push_back(gx1); ay.push_back(py2); aw.push_back(gx2 - gx1); ah.push_back(gy2 - py2); }
                // remove i by swap with last
                fx[i] = fx[n-1]; fy[i] = fy[n-1]; fw[i] = fw[n-1]; fh[i] = fh[n-1];
                fx.pop_back(); fy.pop_back(); fw.pop_back(); fh.pop_back();
                --n;
            } else ++i;
        }
        // append new fragments (skip degenerate)
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

// ---------------- Pack result ----------------
struct PackResult {
    vector<Placed> placements;
    ll totalValue;
    vector<int> used;
};

// Greedy fill: try to place items in given order, each up to its limit, until no placement possible.
// tryBothOrient: if allowRotate, try both orientations and pick the one that fits lower (best-fit).
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
            // orientation 0
            if(it.w <= sk.W && it.h <= sk.H){
                sk.findBest(it.w, it.h, x0, y0);
            }
            int chosenRot = -1, ox=0, oy=0, rw=0, rh=0;
            if(tryBothOrient && allowRotate){
                if(it.h <= sk.W && it.w <= sk.H){
                    sk.findBest(it.h, it.w, x1, y1);
                }
                if(x0>=0 && x1>=0){
                    // pick lower y; tiebreak: smaller area-fragmentation — prefer the one with
                    // resulting top closer to neighbors. Simple: lower y, then larger area (fills more).
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

// Core MaxRects greedy loop: places items (order, up to rem[t]) into the given free-rect state,
// appending placements to res and decrementing rem. Tries both orientations if allowed.
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

// Greedy fill using MaxRects (best-short-side-fit), from an empty bin.
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

// Take an existing packing and greedily fill its leftover free space with additional items
// (respecting remaining per-type limits). Carves the occupied rectangles out of the full bin,
// then runs the MaxRects greedy loop on the resulting free space.
static PackResult fillGaps(const PackResult& base, const vector<int>& order,
                           bool allowRotate, bool tryBothOrient){
    PackResult res = base; // copy placements, totalValue, used
    MaxRects mr;
    mr.init(g_bin.W, g_bin.H);
    // Carve out every occupied rectangle so fr becomes the free space.
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

// ---------------- Shelf packer (horizontal shelves) ----------------
// Sort items by some key; for each, choose shelf height = item height (in chosen orientation),
// place item left to right; when shelf full, start new shelf on top.
// This is good for strip-like / uniform height items.
static PackResult shelfPack(const vector<int>& order, bool allowRotate, bool orientRot){
    PackResult res;
    res.totalValue = 0;
    res.used.assign(g_items.size(), 0);
    vector<int> rem(g_items.size());
    for(size_t i=0;i<g_items.size();++i) rem[i] = g_items[i].limit;
    int curY = 0;     // bottom of current shelf
    int shelfH = 0;   // height of current shelf
    int curX = 0;     // next free x in current shelf
    // We iterate in order; items fill current shelf if they fit height-wise (use rotation to match),
    // else close shelf and open a new one.
    bool progress = true;
    while(progress){
        progress = false;
        for(int t : order){
            if(rem[t] <= 0) continue;
            const ItemType& it = g_items[t];
            // Decide orientation: if orientRot and allowRotate, prefer fitting shelf height.
            int rw, rh, rot;
            // Try to place into current shelf.
            // Option A: as-is. Option B: rotated.
            // We want height <= shelfH if shelf open, and width <= W - curX.
            bool placed = false;
            // If current shelf has room:
            if(curY + 1 <= g_bin.H){
                // try orientation 0 in current shelf
                int rw0 = it.w, rh0 = it.h, rot0 = 0;
                int rw1 = it.h, rh1 = it.w, rot1 = 1;
                // choose orientation that fits current shelf height if any
                int candRot = -1, crw=0, crh=0;
                if(curX < g_bin.W){
                    if(shelfH == 0){
                        // open shelf: pick orientation with smaller height (to leave room above)
                        // but only if both fit width; prefer smaller height.
                        bool fit0 = (rw0 <= g_bin.W && rh0 <= g_bin.H - curY);
                        bool fit1 = allowRotate && (rw1 <= g_bin.W && rh1 <= g_bin.H - curY);
                        if(fit0 && fit1){
                            if(rh0 <= rh1){ candRot=0; crw=rw0; crh=rh0; }
                            else { candRot=1; crw=rw1; crh=rh1; }
                        } else if(fit0){ candRot=0; crw=rw0; crh=rh0; }
                        else if(fit1){ candRot=1; crw=rw1; crh=rh1; }
                    } else {
                        // shelf open with height shelfH: prefer orientation fitting height and width
                        bool fit0 = (rh0 <= shelfH && rw0 <= g_bin.W - curX);
                        bool fit1 = allowRotate && (rh1 <= shelfH && rw1 <= g_bin.W - curX);
                        if(fit0 && fit1){
                            // pick larger area (fills more) — both fit, prefer non-rotated for stability
                            candRot=0; crw=rw0; crh=rh0;
                        } else if(fit0){ candRot=0; crw=rw0; crh=rh0; }
                        else if(fit1){ candRot=1; crw=rw1; crh=rh1; }
                    }
                }
                if(candRot >= 0){
                    // place
                    if(shelfH == 0){
                        shelfH = crh;
                    }
                    Placed p; p.typeId=t; p.x=curX; p.y=curY; p.rot=candRot;
                    res.placements.push_back(p);
                    res.totalValue += it.v;
                    res.used[t]++; rem[t]--;
                    curX += crw;
                    progress = true;
                    continue;
                }
            }
            // Can't fit current shelf: close it and open a new one with this item.
            if(shelfH > 0){
                curY += shelfH;
                curX = 0;
                shelfH = 0;
                if(curY >= g_bin.H){ break; }
            }
            // try opening new shelf with this item
            if(curY < g_bin.H){
                int rw0 = it.w, rh0 = it.h, rot0 = 0;
                int rw1 = it.h, rh1 = it.w, rot1 = 1;
                int candRot=-1, crw=0, crh=0;
                bool fit0 = (rw0 <= g_bin.W && rh0 <= g_bin.H - curY);
                bool fit1 = allowRotate && (rw1 <= g_bin.W && rh1 <= g_bin.H - curY);
                if(fit0 && fit1){
                    if(rh0 <= rh1){ candRot=0; crw=rw0; crh=rh0; }
                    else { candRot=1; crw=rw1; crh=rh1; }
                } else if(fit0){ candRot=0; crw=rw0; crh=rh0; }
                else if(fit1){ candRot=1; crw=rw1; crh=rh1; }
                if(candRot >= 0){
                    shelfH = crh;
                    Placed p; p.typeId=t; p.x=curX; p.y=curY; p.rot=candRot;
                    res.placements.push_back(p);
                    res.totalValue += it.v;
                    res.used[t]++; rem[t]--;
                    curX += crw;
                    progress = true;
                }
            }
        }
        // After a full pass, if we made progress but curY reached H, we stop; otherwise loop.
        if(curY >= g_bin.H) break;
    }
    return res;
}

// ---------------- Single-item-per-shelf knapsack planner ----------------
// For each item type, choose an orientation (rot 0 or 1), then a number of shelves of that orientation.
// Each shelf holds perShelf = W // rw copies (last shelf may be partial). Bounded knapsack over height H.
// We enumerate orientation combos (2^M, M<=12 -> 4096) and pick the best knapsack result.
// Time budget: ~0.4s. If M large / time tight, fall back to greedy orientation choice.
struct ShelfUnit { int typeId; int rot; int count; }; // count = number of shelves of this kind
static PackResult knapsackShelfPlan(bool allowRot){
    int H = g_bin.H, W = g_bin.W;
    int M = (int)g_items.size();
    struct KItem { int w; ll val; int typeId; int rot; int nShelves; int perShelf; bool partial; };
    // Precompute per-type per-orientation shelf data.
    // shelfData[t][rot] = {perShelf, rh, fullShelves, rem} or invalid.
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
    // DP scratch (reused). dp[x] = best value, par[x] = kitem index in current combo's list.
    vector<ll> dp(H + 1, -1);
    vector<int> par(H + 1, -1);
    // For reconstruction we need the kitems list of the winning combo; recompute it after picking combo.
    ll bestVal = -1;
    vector<int> bestCombo(M, 0);
    // enumerate combos
    int total = 1;
    for(int t = 0; t < M; ++t) total *= 2;
    // If M == 12, total = 4096. If allowRot false, only 1 combo.
    if(!allowRot) total = 1;
    // Time-guarded enumeration.
    int comboLimit = total;
    // Cap combos by time: each combo knapsack ~ O(kitems * H). Estimate.
        for(int c = 0; c < comboLimit; ++c){
        if((c & 63) == 0 && elapsed() > TIME_LIMIT * 0.5) { comboLimit = c; break; }
        // decode combo (only meaningful if allowRot)
        int code = c;
        vector<int> combo(M, 0);
        if(allowRot){
            for(int t = 0; t < M; ++t){ combo[t] = code & 1; code >>= 1; }
        }
        // build kitems
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
        // Upper-bound pruning: skip DP if sum of kitems' values can't beat best
        ll ub = 0; for(const auto& ki : kitems) ub += ki.val;
        if(ub <= bestVal) continue;
        // 0/1 knapsack
        fill(dp.begin(), dp.end(), -1);
        dp[0] = 0;
        fill(par.begin(), par.end(), -1);
        for(int idx = 0; idx < (int)kitems.size(); ++idx){
            const KItem& ki = kitems[idx];
            for(int x = H; x >= ki.w; --x){
                if(dp[x - ki.w] >= 0 && dp[x - ki.w] + ki.val > dp[x]){
                    dp[x] = dp[x - ki.w] + ki.val;
                    par[x] = idx;
                }
            }
        }
        ll bestV = 0; int bestX = 0;
        for(int x = 0; x <= H; ++x){ if(dp[x] > bestV){ bestV = dp[x]; bestX = x; } }
        if(bestV > bestVal){
            bestVal = bestV;
            bestCombo = combo;
        }
    }
    // Rebuild winning combo's kitems and reconstruct.
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
    // Correct 0/1 knapsack reconstruction: store dp snapshot before each item.
    // dpBefore[i] = dp state before processing item i (i.e., after items 0..i-1).
    // dpBefore[K] = state before item K (none) = state after all items = final dp.
    int K = (int)kitems.size();
    vector<vector<ll>> dpBefore(K + 1, vector<ll>(H + 1, -1));
    for(int x = 0; x <= H; ++x) dp[x] = (x == 0) ? 0 : -1;
    for(int idx = 0; idx < K; ++idx){
        // snapshot before processing item idx
        for(int x = 0; x <= H; ++x) dpBefore[idx][x] = dp[x];
        const KItem& ki = kitems[idx];
        for(int x = H; x >= ki.w; --x){
            if(dp[x - ki.w] >= 0 && dp[x - ki.w] + ki.val > dp[x]){
                dp[x] = dp[x - ki.w] + ki.val;
            }
        }
    }
    // final snapshot = after all items
    for(int x = 0; x <= H; ++x) dpBefore[K][x] = dp[x];
    ll bestV = 0; int bestX = 0;
    for(int x = 0; x <= H; ++x){ if(dp[x] > bestV){ bestV = dp[x]; bestX = x; } }
    // reconstruct: for item idx, before = dpBefore[idx], after = dpBefore[idx+1].
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
    // Build placements.
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
        // For partial shelf items, perShelf field holds 'rem'; nShelves=1.
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

// Column-based variant: transpose the bin and items, run the shelf knapsack (which stacks
// shelves along height), then map placements back. This packs vertical columns instead of
// horizontal shelves — a different layout that can waste less for some aspect ratios.
static PackResult knapsackColumnPlan(bool allowRot){
    swap(g_bin.W, g_bin.H);
    for(auto& it : g_items) swap(it.w, it.h);
    PackResult r = knapsackShelfPlan(allowRot);
    for(auto& it : g_items) swap(it.w, it.h);
    swap(g_bin.W, g_bin.H);
    for(auto& p : r.placements) swap(p.x, p.y); // rot and type unchanged (see derivation)
    return r;
}

// ---------------- Mixed-shelf (3-stage guillotine) planner ----------------
// A shelf of height Hs spanning the full bin width is filled with vertical
// "columns"; each column stacks copies of a single item type (one orientation).
// Column choice = bounded knapsack over the shelf width (exact, binary-encoded).
// Shelves are chosen greedily by value-per-height with exact remaining-limit
// tracking, which lets several item types mix inside one shelf.
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
        // pick orientation with best value density per width (optionally flipped)
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

// Value-only 0/1 knapsack over width (units already binary-encoded).
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

// Knapsack with reconstruction (snapshot rows). Returns chosen unit indices.
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

// score exponent for shelf selection: score = value / Hs^alpha
static double g_msAlpha = 1.0;

// Fill a sub-region [x0,x0+RW) x [y0,y0+RH) with horizontal mixed shelves,
// consuming 'rem' and appending placements/value to res.
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

// Fill a sub-region with vertical mixed shelves (transposed): swap item dims,
// fill the transposed region, then map placements back.
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
        // footprint in transposed space with swapped dims equals real footprint
        // with the same rot flag (see knapsackColumnPlan derivation).
    }
}

// Build one mixed-shelf plan over the whole bin. stripW > 0 reserves a vertical
// band of that width on the right, filled with transposed (vertical) shelves —
// good for tall narrow leftovers that horizontal shelves waste.
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
        // vertical shelves over the whole bin; optional horizontal strip on top
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

// ---------------- Beam search over shelf-height sequences ----------------
// Greedy shelf choice can lock in a bad height partition of H. Branch on the
// first few shelf heights (top-scoring alternatives), keep the most promising
// states (value + fractional bound on the remaining area), finish greedily.
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
            if(kmax > 5) kmax = 5;
            for(int k = 1; k <= kmax; ++k) cands.push_back(k * h_or);
        }
    }
    cands.push_back(Hrem);
    sort(cands.begin(), cands.end());
    cands.erase(unique(cands.begin(), cands.end()), cands.end());
}

// Commit the best fill for shelf height Hs onto state st (placements at y=st.y).
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
                // score candidates
                static vector<pair<double,int>> scored;
                scored.clear();
                for(int Hs : cands){
                    buildShelfUnits(Hs, W, st.rem, allowRot, 0, units);
                    if(units.empty()) continue;
                    ll v = msKnapValue(units, W);
                    if(v <= 0) continue;
                    scored.push_back({(double)v / (double)Hs, Hs});
                }
                sort(scored.begin(), scored.end(),
                     [](const pair<double,int>& a, const pair<double,int>& b){ return a.first > b.first; });
                int taken = 0;
                static vector<int> takenH;
                takenH.clear();
                for(size_t si = 0; si < scored.size() && taken < nBranch; ++si){
                    int Hs = scored[si].second;
                    // require diversity between branched heights
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
        // dedupe identical (y, value) states
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

// ---------------- Orderings ----------------
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

// ---------------- Output ----------------
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

// ---------------- Main ----------------
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

    // Deterministic strategies.
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

    // Knapsack-based single-item-per-shelf plan (strongest in most cases). Run FIRST with ample budget.
    #ifdef DIAG
    g_label="knap";
#endif
    consider(knapsackShelfPlan(allowRot));

    // Orderings that favor small / slender pieces for filling narrow leftover strips.
    vector<int> ordMinDim = orderByMinDimAsc();
    vector<int> ordDens   = orderByDensity();
    vector<int> ordVal    = orderByValueDesc();
    vector<int> ordAreaAsc = ordMinDim; // area-ascending
    sort(ordAreaAsc.begin(), ordAreaAsc.end(), [](int a,int b){ return g_items[a].area < g_items[b].area; });
    vector<vector<int>*> gapOrders = {&ordMinDim, &ordAreaAsc, &ordDens, &ordVal};
    // Polish a candidate result: iterated MaxRects gap-filling on its leftover space.
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
    // Fill gaps left by the knapsack shelf layout (wasted horizontal strips + top strip).
    gapFill();

    // Mixed-shelf (3-stage guillotine) plans: several item types per shelf via an
    // exact width-knapsack over stacked columns; horizontal + transposed; polished.
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
    // Beam search over shelf-height partitions (branch early shelves, greedy tail).
#ifdef DIAG
    g_label="beam";
#endif
    g_densOrd = orderByDensity();
#ifdef DIAG
    {
        PackResult b1 = beamMixedShelfPlan(allowRot, 4, 3, 5, TIME_LIMIT * 0.55);
        fprintf(stderr, "[beamH] raw=%lld @%.3fs\n", b1.totalValue, elapsed());
        consider(polish(std::move(b1)));
        PackResult b2 = beamMixedShelfPlanT(allowRot, 4, 3, 5, TIME_LIMIT * 0.7);
        fprintf(stderr, "[beamT] raw=%lld @%.3fs\n", b2.totalValue, elapsed());
        consider(polish(std::move(b2)));
    }
#else
    if(elapsed() < TIME_LIMIT * 0.55)
        consider(polish(beamMixedShelfPlan(allowRot, 4, 3, 5, TIME_LIMIT * 0.55)));
    if(elapsed() < TIME_LIMIT * 0.7)
        consider(polish(beamMixedShelfPlanT(allowRot, 4, 3, 5, TIME_LIMIT * 0.7)));
#endif
    // Hybrid plans: reserve a narrow band (right side / top) for perpendicular
    // shelves — plugs tall narrow leftovers that one-directional shelves waste.
#ifdef DIAG
    g_label="hybrid";
#endif
    {
        // candidate strip widths from slender item dimensions
        vector<int> stripCands;
        for(const ItemType& it : g_items){
            if(it.h >= 2 * it.w && it.w <= g_bin.W / 6) stripCands.push_back(it.w);
            if((allowRot || true) && it.w >= 2 * it.h && it.h <= g_bin.W / 6) stripCands.push_back(it.h);
        }
        sort(stripCands.begin(), stripCands.end());
        stripCands.erase(unique(stripCands.begin(), stripCands.end()), stripCands.end());
        // pairwise sums of the narrowest few
        {
            int n = min((int)stripCands.size(), 3);
            vector<int> sums;
            for(int i = 0; i < n; ++i)
                for(int j = i; j < n; ++j) sums.push_back(stripCands[i] + stripCands[j]);
            for(int s : sums) stripCands.push_back(s);
            sort(stripCands.begin(), stripCands.end());
            stripCands.erase(unique(stripCands.begin(), stripCands.end()), stripCands.end());
        }
        if((int)stripCands.size() > 6) stripCands.resize(6);
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
        int splits[3] = {g_bin.W/3, g_bin.W/2, (2*g_bin.W)/3};
        for(int sw : splits){
            for(int mask = 0; mask < 8 && elapsed() < TIME_LIMIT * 0.56; ++mask)
                consider(polish(splitMixedPlan(allowRot, sw, mask, 0)));
            if(elapsed() > TIME_LIMIT * 0.56) break;
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

    // Column-based (transposed) knapsack plan + its own gap fill; keep whichever is better.
    if(elapsed() < TIME_LIMIT * 0.55){
#ifdef DIAG
        g_label="column";
#endif
        consider(knapsackColumnPlan(allowRot));
        gapFill();
    }

    #ifdef DIAG
    g_label="maxrects";
#endif
    // MaxRects packer with multiple orderings + orientation variants (fills irregular gaps best).
    for(auto& ord : orders){
        consider(greedyFillMaxRects(ord, allowRot, false));
        if(allowRot) consider(greedyFillMaxRects(ord, allowRot, true));
        if(elapsed() > TIME_LIMIT * 0.5) break;
    }

    #ifdef DIAG
    g_label="skyline";
#endif
    // Skyline packer with multiple orderings + orientation variants (fills irregular gaps).
    for(auto& ord : orders){
        consider(greedyFill(ord, allowRot, false));
        if(allowRot) consider(greedyFill(ord, allowRot, true));
        if(elapsed() > TIME_LIMIT * 0.6) break;
    }

    #ifdef DIAG
    g_label="random";
#endif
    // Randomized multi-start on skyline after extra mixed-shelf probes.
    mt19937 rng(987654321u);
    int seed = 0;
    double iterCost = 0.0; // adaptive margin: last iteration's duration
    while(elapsed() + 1.5 * iterCost < TIME_LIMIT - 0.04){
        double t0 = elapsed();
        vector<int> ord = orderByDensity();
        int mode = seed % 8;
        if(mode == 0){
            shuffle(ord.begin(), ord.end(), rng);
        } else if(mode == 1){
            // density with random swap
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
            // density order, but move strip items (large min dim) to end
            sort(ord.begin(), ord.end(), [](int a, int b){
                return g_items[a].density > g_items[b].density;
            });
            stable_sort(ord.begin(), ord.end(), [](int a, int b){
                int ma = min(g_items[a].w, g_items[a].h);
                int mb = min(g_items[b].w, g_items[b].h);
                // slender (small min) first
                return ma < mb;
            });
        } else if(mode == 6){
            // random rotation per type chosen, density order
            sort(ord.begin(), ord.end(), [](int a, int b){
                return g_items[a].density > g_items[b].density;
            });
        } else {
            // shuffle and re-sort by density partially
            shuffle(ord.begin(), ord.end(), rng);
            sort(ord.begin(), ord.end(), [](int a, int b){
                return g_items[a].density > g_items[b].density;
            });
            // small perturbation
            for(int s=0;s<2;++s){ int i=rng()%M,j=rng()%M; swap(ord[i],ord[j]); }
        }
        bool tryBoth = allowRot ? ((seed & 1) == 0) : false;
        if(mode == 6) tryBoth = true; // force both
        consider(greedyFill(ord, allowRot, tryBoth));
        iterCost = elapsed() - t0;
        ++seed;
        if(seed > 2000000) break;
    }

    // Final gap-fill on the champion (in case the multistart produced a new best).
    if(elapsed() < TIME_LIMIT - 0.04) gapFill();

    if(best.totalValue < 0){
        best.totalValue = 0;
        best.placements.clear();
    }
    outputResult(best);
    return 0;
}
