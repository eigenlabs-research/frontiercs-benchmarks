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
static double TIME_LIMIT = 0.72;

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

// ---------------- Orderings ----------------
// For a chosen shelf-height sequence, fill each shelf with all items whose (oriented) height
// <= shelf height, highest density first, respecting limits. Shelf heights are chosen greedily
// as we go: at each new shelf, pick the shelf height that maximizes value-per-area of the
// items that can fit in it (considering remaining counts and width-packing within shelf).
// Simpler: try shelf height = each item's height (in both orientations if rotation); pick best overall.
struct ShelfFiller {
    int W, H;
    bool allowRot;
    // Place items into shelves. For each shelf, given a target height Hshelf, iterate items
    // (in 'order') whose oriented height <= Hshelf and oriented width fits remaining shelf width,
    // place them, highest density first. Returns placements and used counts.
    PackResult fillWithShelfHeights(const vector<int>& order, const vector<int>& shelfHeights){
        PackResult res;
        res.totalValue = 0;
        res.used.assign(g_items.size(), 0);
        vector<int> rem(g_items.size());
        for(size_t i=0;i<g_items.size();++i) rem[i] = g_items[i].limit;
        int curY = 0;
        int si = 0;
        while(curY < g_bin.H && si < (int)shelfHeights.size()){
            int sh = shelfHeights[si++];
            if(sh <= 0 || sh > g_bin.H - curY) continue;
            int curX = 0;
            // Repeatedly fill this shelf: scan order, place any item whose oriented height <= sh
            // and width <= W - curX. Prefer higher density (order is already density-sorted).
            bool prog = true;
            while(prog){
                prog = false;
                for(int t : order){
                    if(rem[t] <= 0) continue;
                    const ItemType& it = g_items[t];
                    // try orientation 0
                    int rw0=it.w, rh0=it.h;
                    int rw1=it.h, rh1=it.w;
                    int candRot=-1, crw=0, crh=0;
                    bool fit0 = (rh0 <= sh && rw0 <= g_bin.W - curX);
                    bool fit1 = allowRot && (rh1 <= sh && rw1 <= g_bin.W - curX);
                    if(fit0 && fit1){
                        // pick higher density orientation (same density since area same & v same),
                        // so pick the one with smaller height (less waste) then larger width to fill.
                        if(rh0 <= rh1){ candRot=0; crw=rw0; crh=rh0; }
                        else { candRot=1; crw=rw1; crh=rh1; }
                    } else if(fit0){ candRot=0; crw=rw0; crh=rh0; }
                    else if(fit1){ candRot=1; crw=rw1; crh=rh1; }
                    if(candRot >= 0){
                        Placed p; p.typeId=t; p.x=curX; p.y=curY; p.rot=candRot;
                        res.placements.push_back(p);
                        res.totalValue += it.v;
                        res.used[t]++; rem[t]--;
                        curX += crw;
                        prog = true;
                    }
                }
            }
            curY += sh;
        }
        return res;
    }
};

// Generate a shelf-height plan: greedily pick next shelf height as the height (over orientations)
// of the highest-density item that still has remaining copies and fits, considering what fits in it.
// This produces a list of shelf heights.
static vector<int> greedyShelfHeights(const vector<int>& order, bool allowRot){
    vector<int> rem(g_items.size());
    for(size_t i=0;i<g_items.size();++i) rem[i] = g_items[i].limit;
    vector<int> heights;
    int curY = 0;
    while(curY < g_bin.H){
        // Pick shelf height: the oriented height of the first (in order) item with rem>0 that fits.
        int sh = -1;
        for(int t : order){
            if(rem[t] <= 0) continue;
            const ItemType& it = g_items[t];
            // choose oriented height that fits and is smallest (to leave room)
            int h0 = it.h, h1 = it.w;
            int cand = -1;
            // need at least one copy to fit width-wise
            if(it.w <= g_bin.W && h0 <= g_bin.H - curY) cand = h0;
            if(allowRot && it.h <= g_bin.W && h1 <= g_bin.H - curY){
                if(cand < 0 || h1 < cand) cand = h1;
            }
            if(cand > 0){ sh = cand; break; }
        }
        if(sh < 0) break;
        heights.push_back(sh);
        // simulate filling this shelf to update rem (so next shelf height choice reflects leftovers)
        int curX = 0;
        bool prog = true;
        while(prog){
            prog = false;
            for(int t : order){
                if(rem[t] <= 0) continue;
                const ItemType& it = g_items[t];
                int rw0=it.w, rh0=it.h, rw1=it.h, rh1=it.w;
                int candRot=-1, crw=0, crh=0;
                bool fit0 = (rh0 <= sh && rw0 <= g_bin.W - curX);
                bool fit1 = allowRot && (rh1 <= sh && rw1 <= g_bin.W - curX);
                if(fit0 && fit1){ if(rh0<=rh1){candRot=0;crw=rw0;crh=rh0;} else {candRot=1;crw=rw1;crh=rh1;} }
                else if(fit0){candRot=0;crw=rw0;crh=rh0;}
                else if(fit1){candRot=1;crw=rw1;crh=rh1;}
                if(candRot>=0){
                    rem[t]--; curX += crw; prog = true;
                }
            }
        }
         curY += sh;
     }
     return heights;
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
        if((c & 63) == 0 && elapsed() > TIME_LIMIT * 0.6) { comboLimit = c; break; }
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

    // Shelf packer with multiple orderings.
    ShelfFiller sf;
    sf.W = g_bin.W; sf.H = g_bin.H; sf.allowRot = allowRot;
    for(auto& ord : orders){
        // greedy shelf heights based on this order
        vector<int> heights = greedyShelfHeights(ord, allowRot);
        if(!heights.empty()){
            consider(sf.fillWithShelfHeights(ord, heights));
        }
        if(elapsed() > TIME_LIMIT * 0.55) break;
    }
    // Also try: shelf heights = sorted unique item heights (descending), fill each.
    {
        vector<int> allH;
        for(auto& it : g_items){
            allH.push_back(it.h);
            if(allowRot) allH.push_back(it.w);
        }
        sort(allH.begin(), allH.end());
        allH.erase(unique(allH.begin(), allH.end()), allH.end());
        sort(allH.begin(), allH.end(), greater<int>());
        // For each candidate first-shelf-height, build a plan: pick heights greedily from allH
        // that fit remaining vertical space, fill each shelf with density order.
        vector<int> densOrd = orderByDensity();
        // Try a few plans: allH as-is (descending); also greedy from each starting height.
        consider(sf.fillWithShelfHeights(densOrd, allH));
        // Greedy plan from density order
        vector<int> gh = greedyShelfHeights(densOrd, allowRot);
        if(!gh.empty()) consider(sf.fillWithShelfHeights(densOrd, gh));
        // Try ascending too (small shelves first can fill gaps)
        sort(allH.begin(), allH.end());
        consider(sf.fillWithShelfHeights(densOrd, allH));
        // Try tiling with a single shelf height = each candidate (best for uniform items).
        for(int h : allH){
            if(h <= 0) continue;
            int n = g_bin.H / h;
            if(n <= 0) continue;
            vector<int> plan;
            for(int k=0;k<n;++k) plan.push_back(h);
            // leftover strip
            int remH = g_bin.H - n*h;
            if(remH > 0) plan.push_back(remH);
            consider(sf.fillWithShelfHeights(densOrd, plan));
            if(elapsed() > TIME_LIMIT * 0.55) break;
        }
        // Try plans starting with each candidate height, then greedy.
        for(int h0 : allH){
            if(h0 <= 0 || h0 > g_bin.H) continue;
            // build plan: first shelf h0, then greedy from density order with that rem
            vector<int> rem(g_items.size());
            for(size_t i=0;i<g_items.size();++i) rem[i] = g_items[i].limit;
            vector<int> plan;
            int curY = 0;
            // first shelf h0
            plan.push_back(h0);
            // simulate fill of first shelf to update rem
            {
                int curX=0; bool prog=true;
                while(prog){ prog=false;
                    for(int t : densOrd){
                        if(rem[t]<=0) continue;
                        const ItemType& it=g_items[t];
                        int rw0=it.w,rh0=it.h,rw1=it.h,rh1=it.w;
                        int cr=-1,cw=0,ch=0;
                        bool f0=(rh0<=h0 && rw0<=g_bin.W-curX);
                        bool f1=allowRot&&(rh1<=h0 && rw1<=g_bin.W-curX);
                        if(f0&&f1){if(rh0<=rh1){cr=0;cw=rw0;ch=rh0;}else{cr=1;cw=rw1;ch=rh1;}}
                        else if(f0){cr=0;cw=rw0;ch=rh0;}
                        else if(f1){cr=1;cw=rw1;ch=rh1;}
                        if(cr>=0){rem[t]--;curX+=cw;prog=true;}
                    }
                }
            }
            curY += h0;
            // greedy remainder
            while(curY < g_bin.H){
                int sh=-1;
                for(int t : densOrd){
                    if(rem[t]<=0) continue;
                    const ItemType& it=g_items[t];
                    int h0b=it.h, h1b=it.w;
                    int cand=-1;
                    if(it.w<=g_bin.W && h0b<=g_bin.H-curY) cand=h0b;
                    if(allowRot && it.h<=g_bin.W && h1b<=g_bin.H-curY){ if(cand<0||h1b<cand) cand=h1b; }
                    if(cand>0){ sh=cand; break; }
                }
                if(sh<0) break;
                plan.push_back(sh);
                int curX=0; bool prog=true;
                while(prog){ prog=false;
                    for(int t : densOrd){
                        if(rem[t]<=0) continue;
                        const ItemType& it=g_items[t];
                        int rw0=it.w,rh0=it.h,rw1=it.h,rh1=it.w;
                        int cr=-1,cw=0,ch=0;
                        bool f0=(rh0<=sh && rw0<=g_bin.W-curX);
                        bool f1=allowRot&&(rh1<=sh && rw1<=g_bin.W-curX);
                        if(f0&&f1){if(rh0<=rh1){cr=0;cw=rw0;ch=rh0;}else{cr=1;cw=rw1;ch=rh1;}}
                        else if(f0){cr=0;cw=rw0;ch=rh0;}
                        else if(f1){cr=1;cw=rw1;ch=rh1;}
                        if(cr>=0){rem[t]--;curX+=cw;prog=true;}
                    }
                }
                curY += sh;
            }
            consider(sf.fillWithShelfHeights(densOrd, plan));
            if(elapsed() > TIME_LIMIT * 0.55) break;
        }
    }

    #ifdef DIAG
    g_label="random";
#endif
    // Randomized multi-start on skyline.
    mt19937 rng(987654321u);
    int seed = 0;
    while(elapsed() < TIME_LIMIT){
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
        ++seed;
        if(seed > 2000000) break;
    }

    // Final gap-fill on the champion (in case the multistart produced a new best).
    gapFill();

    if(best.totalValue < 0){
        best.totalValue = 0;
        best.placements.clear();
    }
    outputResult(best);
    return 0;
}
