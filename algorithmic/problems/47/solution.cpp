#include <bits/stdc++.h>
using namespace std;

struct Item { string name; int w,h,lim; long long v; };
struct Pl { int id,x,y,rot; };
struct Rect { int x,y,w,h; };

static string getString(const string& s, const string& k) {
    regex r("\\\"" + k + "\\\"\\s*:\\s*\\\"([^\\\"]*)\\\""); smatch m;
    return regex_search(s,m,r) ? m[1].str() : "";
}
static long long getInt(const string& s, const string& k) {
    regex r("\\\"" + k + "\\\"\\s*:\\s*(-?[0-9]+)"); smatch m;
    return regex_search(s,m,r) ? stoll(m[1].str()) : 0;
}
static string esc(const string& s) {
    string r; for(char c:s) { if(c=='\\' || c=='\"') r+='\\'; r+=c; } return r;
}
static bool contains(const Rect&a,const Rect&b) {
    return b.x>=a.x && b.y>=a.y && b.x+b.w<=a.x+a.w && b.y+b.h<=a.y+a.h;
}
static bool intersects(const Rect&a,int x,int y,int w,int h) {
    return x < a.x+a.w && x+w > a.x && y < a.y+a.h && y+h > a.y;
}

/* Maximal empty rectangles retain overlapping descriptions of empty space.  Thus a
   selected rectangle is always safe, while the alternate descriptions let later
   items use L-shaped leftovers that a shelf representation loses. */
static vector<Pl> pack(const vector<Item>& a, int W, int H, bool canRot,
                       const vector<int>& order, int policy) {
    int n=a.size(), cap=0; for(auto&q:a) cap+=q.lim;
    cap=min(cap, 4000);
    vector<int> used(n); vector<Rect> freeR(1,{0,0,W,H}); vector<Pl> ans;
    for(int step=0; step<cap; ++step) {
        int bi=-1,bf=-1,br=0,bw=0,bh=0;
        long long bm1=LLONG_MAX,bm2=LLONG_MAX;
        // Type priority is deliberate: it protects scarce high-density pieces.
        // Within that type, choose the empty rectangle with the least fragmentation.
        for(int id:order) {
            if(used[id]>=a[id].lim) continue;
            bool found=false; int ci=-1,cr=0,cw=0,ch=0; long long cm1=LLONG_MAX,cm2=LLONG_MAX;
            for(int f=0; f<(int)freeR.size(); ++f) for(int r=0;r<=(canRot?1:0);++r) {
                int w=r?a[id].h:a[id].w, h=r?a[id].w:a[id].h;
                const Rect& z=freeR[f]; if(w>z.w || h>z.h) continue;
                long long shortSide=min(z.w-w,z.h-h), longSide=max(z.w-w,z.h-h);
                long long area=1LL*z.w*z.h-1LL*w*h;
                long long m1,m2;
                if(policy==0) m1=shortSide,m2=longSide;
                else if(policy==1) m1=area,m2=shortSide;
                else m1=longSide,m2=area;
                if(!found || m1<cm1 || (m1==cm1 && (m2<cm2 || (m2==cm2 && (z.y<freeR[ci].y || (z.y==freeR[ci].y && z.x<freeR[ci].x)))))) {
                    found=true; ci=f;cr=r;cw=w;ch=h;cm1=m1;cm2=m2;
                }
            }
            if(found) { bi=id;bf=ci;br=cr;bw=cw;bh=ch;bm1=cm1;bm2=cm2; break; }
        }
        if(bi<0) break;
        Rect place{freeR[bf].x,freeR[bf].y,bw,bh};
        ans.push_back({bi,place.x,place.y,br}); ++used[bi];
        vector<Rect> next, add; next.reserve(freeR.size()+4); add.reserve(16);
        for(const Rect& z:freeR) {
            if(!intersects(z,place.x,place.y,place.w,place.h)) { next.push_back(z); continue; }
            if(place.x>z.x) add.push_back({z.x,z.y,place.x-z.x,z.h});
            if(place.x+place.w<z.x+z.w) add.push_back({place.x+place.w,z.y,z.x+z.w-(place.x+place.w),z.h});
            if(place.y>z.y) add.push_back({z.x,z.y,z.w,place.y-z.y});
            if(place.y+place.h<z.y+z.h) add.push_back({z.x,place.y+place.h,z.w,z.y+z.h-(place.y+place.h)});
        }
        // The old uncut rectangles were already maximal.  Consequently only a
        // newly generated residual needs containment comparisons; this avoids the
        // quadratic all-pairs prune commonly responsible for MaxRects timeouts.
        for(const Rect& g:add) {
            if(g.w<=0 || g.h<=0) continue;
            bool covered=false;
            for(const Rect& z:next) if(contains(z,g)) { covered=true; break; }
            if(covered) continue;
            for(int i=0;i<(int)next.size();) {
                if(contains(g,next[i])) next[i]=next.back(),next.pop_back();
                else ++i;
            }
            next.push_back(g);
        }
        freeR.swap(next);
    }
    return ans;
}

int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string s((istreambuf_iterator<char>(cin)),{});
    int W=(int)getInt(s,"W"), H=(int)getInt(s,"H");
    smatch m; regex rr("\\\"allow_rotate\\\"\\s*:\\s*(true|false)");
    bool rot=regex_search(s,m,rr) && m[1].str()=="true";
    vector<Item>a; regex obj("\\{[^{}]*\\}");
    for(sregex_iterator it(s.begin(),s.end(),obj),en;it!=en;++it) {
        string z=it->str(), name=getString(z,"type"); if(name.empty()) continue;
        Item q{name,(int)getInt(z,"w"),(int)getInt(z,"h"),(int)getInt(z,"limit"),getInt(z,"v")};
        if(q.w>0 && q.h>0 && q.lim>0) a.push_back(q);
    }
    int n=a.size(); vector<int> base(n); iota(base.begin(),base.end(),0);
    vector<vector<int>> orders;
    auto add=[&](auto cmp) { auto x=base; stable_sort(x.begin(),x.end(),cmp); orders.push_back(x); };
    add([&](int i,int j) { __int128 l=(__int128)a[i].v*a[j].w*a[j].h, r=(__int128)a[j].v*a[i].w*a[i].h; return l!=r ? l>r : a[i].v>a[j].v; });
    add([&](int i,int j) { return a[i].v!=a[j].v ? a[i].v>a[j].v : 1LL*a[i].w*a[i].h<1LL*a[j].w*a[j].h; });
    add([&](int i,int j) { long double di=(long double)a[i].v/(a[i].w*a[i].h), dj=(long double)a[j].v/(a[j].w*a[j].h); long double si=di*(1+.18L*min(a[i].w,a[i].h)/(long double)max(a[i].w,a[i].h)), sj=dj*(1+.18L*min(a[j].w,a[j].h)/(long double)max(a[j].w,a[j].h)); return si>sj; });
    add([&](int i,int j) { return 1LL*a[i].w*a[i].h!=1LL*a[j].w*a[j].h ? 1LL*a[i].w*a[i].h<1LL*a[j].w*a[j].h : a[i].v>a[j].v; });
    vector<Pl> best; long long bestV=-1;
    for(int k=0;k<(int)orders.size();++k) for(int p=0;p<3;++p) {
        // Six focused restarts are more useful under a one-second limit than a
        // randomized long tail, and are deterministic for reproducible scoring.
        if((k==2 || k==3) && p==2) continue;
        auto cur=pack(a,W,H,rot,orders[k],p); long long v=0; for(auto&q:cur)v+=a[q.id].v;
        if(v>bestV) bestV=v,best.swap(cur);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i) { if(i) cout<<','; const Pl&q=best[i];
        cout << "{\"type\":\"" << esc(a[q.id].name) << "\",\"x\":" << q.x << ",\"y\":" << q.y << ",\"rot\":" << q.rot << '}'; }
    cout << "]}\n";
}
