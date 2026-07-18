#include <bits/stdc++.h>
using namespace std;

struct Cell { int x, y; };
struct Ori {
    vector<Cell> c;
    int w=0, h=0;
    int r=0, f=0;
    int minx=0, miny=0;
};
struct Piece {
    vector<pair<int,int>> orig;
    vector<Ori> ori;
    int k=0;
};
struct Placement { int x=0, y=0, r=0, f=0, minx=0, miny=0; int shelf=0; };
struct Solution { long long area=LLONG_MAX; int W=0,H=0; vector<Placement> p; };

static pair<int,int> transformCell(int x,int y,int r,int f){
    if(f) x = -x;
    switch(r&3){
        case 0: return {x,y};
        case 1: return {y,-x};
        case 2: return {-x,-y};
        default: return {-y,x};
    }
}

static vector<Ori> makeOrientations(const vector<pair<int,int>>& cells){
    vector<Ori> res;
    set<vector<pair<int,int>>> seen;
    for(int f=0; f<2; ++f) for(int r=0; r<4; ++r){
        vector<pair<int,int>> v;
        int mnx=INT_MAX,mny=INT_MAX,mxx=INT_MIN,mxy=INT_MIN;
        for(auto [x,y]: cells){
            auto [tx,ty]=transformCell(x,y,r,f);
            v.push_back({tx,ty});
            mnx=min(mnx,tx); mny=min(mny,ty); mxx=max(mxx,tx); mxy=max(mxy,ty);
        }
        for(auto &p:v){ p.first-=mnx; p.second-=mny; }
        sort(v.begin(), v.end());
        v.erase(unique(v.begin(), v.end()), v.end());
        if(seen.insert(v).second){
            Ori o; o.r=r; o.f=f; o.minx=mnx; o.miny=mny; o.w=mxx-mnx+1; o.h=mxy-mny+1;
            for(auto [x,y]: v) o.c.push_back({x,y});
            res.push_back(o);
        }
    }
    sort(res.begin(), res.end(), [](const Ori&a,const Ori&b){
        if(a.h!=b.h) return a.h<b.h;
        if(a.w!=b.w) return a.w<b.w;
        return a.r*2+a.f < b.r*2+b.f;
    });
    return res;
}

struct Shelf { int used=0, h=0; vector<int> ids; };
struct Info { int id,k,minArea,maxDim,bestH,bestW; };

static Solution packShelves(const vector<Piece>& pieces, const vector<int>& order, int W){
    int n=(int)pieces.size();
    vector<Shelf> shelves;
    Solution sol; sol.W=W; sol.p.assign(n, {});

    for(int id: order){
        const Piece &pc=pieces[id];
        int bestS=-1,bestOi=-1;
        long long bestCost=LLONG_MAX;
        int bestWaste=INT_MAX;

        int firstShelf = max(0, (int)shelves.size() - 256);
        for(int s=firstShelf; s<(int)shelves.size(); ++s){
            const Shelf &sh=shelves[s];
            for(int oi=0; oi<(int)pc.ori.size(); ++oi){
                const Ori &o=pc.ori[oi];
                if(o.w > W - sh.used) continue;
                int nh=max(sh.h,o.h);
                long long cost=1LL*(nh-sh.h)*W;
                int waste=(W - (sh.used+o.w)) + max(0, sh.h-o.h);
                if(cost<bestCost || (cost==bestCost && waste<bestWaste)){
                    bestCost=cost; bestWaste=waste; bestS=s; bestOi=oi;
                }
            }
        }

        // Start a new shelf if no existing shelf fits, or if fitting would raise a shelf
        // much more than the height of a compact new shelf.  This keeps the method fast
        // and avoids the zero-score timeout of a cell-by-cell bottom-left search.
        int newOi=-1;
        for(int oi=0; oi<(int)pc.ori.size(); ++oi){
            const Ori&o=pc.ori[oi];
            if(o.w>W) continue;
            if(newOi<0 || o.h < pc.ori[newOi].h || (o.h==pc.ori[newOi].h && o.w<pc.ori[newOi].w)) newOi=oi;
        }
        if(newOi<0){ sol.area=LLONG_MAX; return sol; }
        int newH=pc.ori[newOi].h;
        if(bestS<0 || bestCost > 1LL*newH*W/2){
            Shelf sh; sh.used=0; sh.h=0; shelves.push_back(sh); bestS=(int)shelves.size()-1; bestOi=newOi;
        }

        const Ori &o=pc.ori[bestOi];
        Shelf &sh=shelves[bestS];
        sol.p[id] = {sh.used, 0, o.r, o.f, o.minx, o.miny, bestS};
        sh.used += o.w;
        sh.h = max(sh.h, o.h);
        sh.ids.push_back(id);
    }

    int y=0;
    for(int s=0; s<(int)shelves.size(); ++s){
        for(int id: shelves[s].ids) sol.p[id].y = y;
        y += shelves[s].h;
    }
    sol.H=max(1,y);
    sol.area=1LL*sol.W*sol.H;
    return sol;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<Piece> pieces(n);
    long long total=0; int globalMaxDim=1;
    for(int i=0;i<n;++i){
        int k; cin>>k; pieces[i].k=k; total+=k;
        pieces[i].orig.resize(k);
        for(int j=0;j<k;++j) cin>>pieces[i].orig[j].first>>pieces[i].orig[j].second;
        pieces[i].ori=makeOrientations(pieces[i].orig);
        for(auto &o:pieces[i].ori) globalMaxDim=max(globalMaxDim, max(o.w,o.h));
    }

    vector<Info> info; info.reserve(n);
    for(int i=0;i<n;++i){
        int ma=INT_MAX, md=0, bh=INT_MAX, bw=INT_MAX;
        for(auto&o:pieces[i].ori){ ma=min(ma,o.w*o.h); md=max(md,max(o.w,o.h)); if(o.h<bh || (o.h==bh&&o.w<bw)){bh=o.h; bw=o.w;} }
        info.push_back({i,pieces[i].k,ma,md,bh,bw});
    }
    vector<int> order(n); iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a,int b){
        const auto&A=info[a], &B=info[b];
        if(A.bestH!=B.bestH) return A.bestH>B.bestH;
        if(A.bestW!=B.bestW) return A.bestW>B.bestW;
        if(A.minArea!=B.minArea) return A.minArea>B.minArea;
        if(A.k!=B.k) return A.k>B.k;
        return a<b;
    });

    int root=max(1,(int)llround(sqrt((long double)total)));
    vector<int> widths;
    auto addW=[&](int w){ if(w>=globalMaxDim && w<=4000000) widths.push_back(w); };
    for(double m: {0.75,0.90,1.00,1.10,1.25,1.45,1.70,2.05,2.50,3.20}) addW((int)llround(root*m));
    for(int d=-12; d<=12; d+=4) addW(root+d);
    for(int w: {globalMaxDim, globalMaxDim+1, 16,24,32,48,64,96,128,192,256,384}) addW(w);
    sort(widths.begin(), widths.end()); widths.erase(unique(widths.begin(), widths.end()), widths.end());

    Solution best;
    for(int W: widths){
        Solution s=packShelves(pieces, order, W);
        if(s.area<best.area || (s.area==best.area && (s.H<best.H || (s.H==best.H && s.W<best.W)))) best=move(s);
    }

    if(best.area==LLONG_MAX){
        best.W=1; best.H=(int)total; best.p.assign(n,{}); int y=0;
        for(int i=0;i<n;++i){ const Ori&o=pieces[i].ori[0]; best.W=max(best.W,o.w); best.p[i]={0,y,o.r,o.f,o.minx,o.miny,0}; y+=o.h; }
        best.H=max(1,y); best.area=1LL*best.W*best.H;
    }

    cout<<best.W<<' '<<best.H<<'\n';
    for(int i=0;i<n;++i){
        const auto&p=best.p[i];
        cout << (p.x - p.minx) << ' ' << (p.y - p.miny) << ' ' << p.r << ' ' << p.f << '\n';
    }
    return 0;
}
