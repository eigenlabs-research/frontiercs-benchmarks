#include <bits/stdc++.h>
using namespace std;

struct Variant {
    int w=0,h=0,r=0,f=0;
    int minx=0,miny=0;              // transformed coords were normalized by subtracting these
    vector<pair<int,int>> cells;     // normalized
    vector<unsigned short> rowmask;  // h rows, bits in [0,w)
};
struct Piece { vector<pair<int,int>> orig; vector<Variant> vars; int k; int id; };
struct Place { int x=0,y=0,r=0,f=0,minx=0,miny=0,w=0,h=0; };

static pair<int,int> trcell(int x,int y,int r,int f){
    if(f) x=-x;
    switch(r&3){
        case 0: return {x,y};
        case 1: return {y,-x};
        case 2: return {-x,-y};
        default: return {-y,x};
    }
}

static vector<Variant> make_variants(const vector<pair<int,int>>& cells){
    vector<Variant> res;
    set<vector<pair<int,int>>> seen;
    for(int f=0; f<2; ++f) for(int r=0; r<4; ++r){
        vector<pair<int,int>> v;
        int mnx=INT_MAX,mny=INT_MAX,mxx=INT_MIN,mxy=INT_MIN;
        for(auto [x,y]: cells){
            auto p=trcell(x,y,r,f);
            v.push_back(p); mnx=min(mnx,p.first); mny=min(mny,p.second);
            mxx=max(mxx,p.first); mxy=max(mxy,p.second);
        }
        for(auto &p:v){ p.first-=mnx; p.second-=mny; }
        sort(v.begin(), v.end());
        if(!seen.insert(v).second) continue;
        Variant a; a.w=mxx-mnx+1; a.h=mxy-mny+1; a.r=r; a.f=f; a.minx=mnx; a.miny=mny; a.cells=v;
        a.rowmask.assign(a.h,0);
        for(auto [x,y]:v) a.rowmask[y] |= (unsigned short)(1u<<x);
        res.push_back(a);
    }
    sort(res.begin(), res.end(), [](const Variant& a,const Variant& b){
        if(a.h!=b.h) return a.h<b.h;
        if(a.w!=b.w) return a.w<b.w;
        return a.cells<b.cells;
    });
    return res;
}

struct Packer {
    int W, blocks, H=0;
    vector<vector<unsigned long long>> rows;
    Packer(int w=1):W(w),blocks((w+63)>>6){}
    void ensure(int h){ while((int)rows.size()<h) rows.push_back(vector<unsigned long long>(blocks,0)); }
    inline bool fits(const Variant& v,int x,int y){
        if(x<0 || x+v.w>W || y<0) return false;
        int b=x>>6, off=x&63;
        for(int ry=0; ry<v.h; ++ry){
            unsigned long long m=v.rowmask[ry];
            if(!m) continue;
            if(y+ry >= (int)rows.size()) continue;
            unsigned long long lo = m << off;
            if(rows[y+ry][b] & lo) return false;
            if(off && b+1<blocks){
                unsigned long long hi = m >> (64-off);
                if(hi && (rows[y+ry][b+1] & hi)) return false;
            }
        }
        return true;
    }
    inline void put(const Variant& v,int x,int y){
        ensure(y+v.h);
        int b=x>>6, off=x&63;
        for(int ry=0; ry<v.h; ++ry){
            unsigned long long m=v.rowmask[ry]; if(!m) continue;
            rows[y+ry][b] |= (m << off);
            if(off && b+1<blocks) rows[y+ry][b+1] |= (m >> (64-off));
        }
        H=max(H,y+v.h);
    }
};

static vector<Place> pack_exact(const vector<Piece>& pieces, const vector<int>& order, int W, int scanLimit){
    Packer pk(W); vector<Place> ans(pieces.size());
    for(int idx: order){
        const Piece& pc=pieces[idx];
        bool done=false;
        int maxY = pk.H + 12;
        if(scanLimit>0) maxY = min(maxY, max(0, pk.H-scanLimit) + scanLimit + 12);
        int startY = (scanLimit>0 && pk.H>scanLimit) ? pk.H-scanLimit : 0;
        // For modest instances, search the whole built rectangle; for large ones, only the active top band.
        for(int y=startY; !done && y<=maxY; ++y){
            for(int x=0; !done && x<W; ++x){
                for(const Variant& v: pc.vars){
                    if(x+v.w<=W && pk.fits(v,x,y)){
                        pk.put(v,x,y); ans[pc.id]={x,y,v.r,v.f,v.minx,v.miny,v.w,v.h}; done=true; break;
                    }
                }
            }
        }
        if(!done){ // guaranteed skyline-like fallback above everything
            int bestY=INT_MAX,bestX=0; const Variant* best=nullptr;
            for(const Variant& v: pc.vars) if(v.w<=W){
                for(int x=0; x+v.w<=W; ++x){ if(pk.H < bestY){ bestY=pk.H; bestX=x; best=&v; break; } }
            }
            if(!best){ best=&pc.vars[0]; }
            pk.put(*best,bestX,bestY); ans[pc.id]={bestX,bestY,best->r,best->f,best->minx,best->miny,best->w,best->h};
        }
    }
    return ans;
}

static pair<vector<Place>,int> pack_skyline(const vector<Piece>& pieces, const vector<int>& order, int W){
    Packer pk(W); vector<int> colH(W,0); vector<Place> ans(pieces.size());
    for(int idx: order){
        const Piece& pc=pieces[idx];
        int bestTop=INT_MAX,bestY=INT_MAX,bestX=0; const Variant* best=nullptr;
        for(const Variant& v: pc.vars) if(v.w<=W){
            for(int x=0; x+v.w<=W; ++x){
                int y=0;
                for(auto [cx,cy]: v.cells) y=max(y, colH[x+cx]-cy);
                int top=y+v.h;
                if(top<bestTop || (top==bestTop && (y<bestY || (y==bestY && x<bestX)))){
                    bestTop=top; bestY=y; bestX=x; best=&v;
                }
            }
        }
        if(!best){ best=&pc.vars[0]; bestX=0; bestY=pk.H; }
        pk.put(*best,bestX,bestY);
        for(auto [cx,cy]: best->cells) colH[bestX+cx]=max(colH[bestX+cx], bestY+cy+1);
        ans[pc.id]={bestX,bestY,best->r,best->f,best->minx,best->miny,best->w,best->h};
    }
    return {ans, pk.H};
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<Piece> pieces(n); long long S=0; int maxMinW=1;
    for(int i=0;i<n;++i){
        int k; cin>>k; pieces[i].k=k; pieces[i].id=i; pieces[i].orig.resize(k); S+=k;
        for(auto &p: pieces[i].orig) cin>>p.first>>p.second;
        pieces[i].vars=make_variants(pieces[i].orig);
        int mw=1000000; for(auto &v:pieces[i].vars) mw=min(mw,v.w); maxMinW=max(maxMinW,mw);
    }
    vector<int> order(n); iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a,int b){
        const Piece &A=pieces[a], &B=pieces[b];
        int ab=100, bb=100, ad=0, bd=0;
        for(auto &v:A.vars){ ab=min(ab,v.w*v.h); ad=max(ad,max(v.w,v.h)); }
        for(auto &v:B.vars){ bb=min(bb,v.w*v.h); bd=max(bd,max(v.w,v.h)); }
        if(A.k!=B.k) return A.k>B.k;
        if(ab!=bb) return ab>bb;
        if(ad!=bd) return ad>bd;
        return a<b;
    });

    double rt=sqrt((double)S);
    vector<int> cand;
    double factors[] = {0.75,0.88,1.0,1.12,1.28,1.45,1.7,2.05};
    for(double f:factors) cand.push_back(max(maxMinW,(int)llround(rt*f)));
    for(int extra: {maxMinW, (int)rt, (int)rt+5, (int)(rt*1.5)+3}) cand.push_back(max(maxMinW,extra));
    sort(cand.begin(), cand.end()); cand.erase(unique(cand.begin(), cand.end()), cand.end());

    vector<Place> bestAns; long long bestArea=LLONG_MAX; int bestW=0,bestH=0;
    auto consider=[&](int W, vector<Place>& ans){
        int H=0; for(auto &p:ans) H=max(H,p.y+p.h);
        long long area=1LL*W*H;
        if(area<bestArea || (area==bestArea && (H<bestH || (H==bestH && W<bestW)))){ bestArea=area; bestW=W; bestH=H; bestAns=ans; }
    };

    auto t0 = chrono::steady_clock::now();
    auto elapsed = [&](){ return chrono::duration<double>(chrono::steady_clock::now()-t0).count(); };
    // Dense bottom-left search for small cases; top-band search for larger cases.
    int tries=0;
    for(int W: cand){
        if(tries>0 && elapsed()>1.55) break;
        int scan = (S<=6000 ? 0 : 48);
        vector<Place> ans = pack_exact(pieces, order, W, scan);
        consider(W, ans);
        if(++tries >= (S>50000 ? 5 : 8)) break;
    }
    // A few fast skyline alternatives, often better on very large instances.
    for(int W: cand){
        if(elapsed()>1.85) break;
        auto pr=pack_skyline(pieces, order, W);
        vector<Place> ans=std::move(pr.first);
        consider(W, ans);
    }

    cout << bestW << ' ' << bestH << '\n';
    for(int i=0;i<n;++i){
        auto &p=bestAns[i];
        cout << (p.x - p.minx) << ' ' << (p.y - p.miny) << ' ' << p.r << ' ' << p.f << '\n';
    }
    return 0;
}
