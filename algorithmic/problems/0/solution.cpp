#include <bits/stdc++.h>
using namespace std;

struct Orient {
    int R, F;
    int minx, miny, w, h;
    vector<pair<int,int>> cells;      // normalized
    vector<unsigned long long> rows;  // row bit masks, w <= 10 normally
};
struct Piece { int k; vector<pair<int,int>> orig; vector<Orient> ors; };
struct AnsPlace { long long X=0, Y=0; int R=0, F=0; };

static pair<int,int> rotcw(int x,int y,int r){
    switch(r&3){
        case 0: return {x,y};
        case 1: return {y,-x};
        case 2: return {-x,-y};
        default: return {-y,x};
    }
}

static vector<Orient> make_orients(const vector<pair<int,int>>& c){
    vector<Orient> res;
    set<vector<pair<int,int>>> seen;
    for(int F=0; F<=1; ++F) for(int R=0; R<4; ++R){
        vector<pair<int,int>> t;
        int mnx=INT_MAX,mny=INT_MAX,mxx=INT_MIN,mxy=INT_MIN;
        for(auto [x,y]: c){
            if(F) x=-x;
            auto p=rotcw(x,y,R);
            t.push_back(p);
            mnx=min(mnx,p.first); mny=min(mny,p.second);
            mxx=max(mxx,p.first); mxy=max(mxy,p.second);
        }
        for(auto &p:t){ p.first-=mnx; p.second-=mny; }
        sort(t.begin(), t.end());
        if(!seen.insert(t).second) continue;
        Orient o; o.R=R; o.F=F; o.minx=mnx; o.miny=mny;
        o.w=mxx-mnx+1; o.h=mxy-mny+1; o.cells=t;
        o.rows.assign(o.h,0);
        for(auto [x,y]:t) o.rows[y] |= (1ULL<<x);
        res.push_back(o);
    }
    sort(res.begin(), res.end(), [](const Orient&a,const Orient&b){
        int aa=a.w*a.h, bb=b.w*b.h;
        if(aa!=bb) return aa<bb;
        if(a.h!=b.h) return a.h<b.h;
        return a.w<b.w;
    });
    return res;
}

struct Packer {
    int W, chunks, maxH=0;
    vector<vector<unsigned long long>> grid;
    vector<pair<int,int>> cand;
    unordered_set<unsigned long long> seen;
    vector<AnsPlace> place;

    Packer(int w,int n): W(w), chunks((w+63)>>6), place(n) {
        grid.reserve(4096);
        seen.reserve(1<<18);
        addCand(0,0);
    }
    void ensureRows(int h){ while((int)grid.size()<h) grid.push_back(vector<unsigned long long>(chunks,0)); }
    void addCand(int x,int y){
        if(x<0 || x>=W || y<0) return;
        unsigned long long key=((unsigned long long)(unsigned int)y<<32) | (unsigned int)x;
        if(seen.insert(key).second) cand.push_back({x,y});
    }
    bool fits(const Orient& o,int x,int y){
        if(x<0 || y<0 || x+o.w>W) return false;
        ensureRows(y+o.h);
        int ch=x>>6, sh=x&63;
        for(int r=0;r<o.h;++r){
            unsigned long long m=o.rows[r]; if(!m) continue;
            unsigned long long lo = (sh? (m<<sh) : m);
            if(grid[y+r][ch] & lo) return false;
            if(sh && ch+1<chunks){
                unsigned long long hi = (m>>(64-sh));
                if(hi && (grid[y+r][ch+1]&hi)) return false;
            }
        }
        return true;
    }
    int contact(const Orient& o,int x,int y){
        int s=0;
        for(auto [cx,cy]: o.cells){
            int gx=x+cx, gy=y+cy;
            if(gx==0 || gx==W-1 || gy==0) ++s;
            static const int dx[4]={1,-1,0,0};
            static const int dy[4]={0,0,1,-1};
            for(int d=0; d<4; ++d){
                int nx=gx+dx[d], ny=gy+dy[d];
                if(nx<0||nx>=W||ny<0) continue;
                if(ny<(int)grid.size() && (grid[ny][nx>>6]&(1ULL<<(nx&63)))) ++s;
            }
        }
        return s;
    }
    void put(const Orient& o,int x,int y,int idx){
        ensureRows(y+o.h);
        int ch=x>>6, sh=x&63;
        for(int r=0;r<o.h;++r){
            unsigned long long m=o.rows[r]; if(!m) continue;
            grid[y+r][ch] |= (sh? (m<<sh) : m);
            if(sh && ch+1<chunks) grid[y+r][ch+1] |= (m>>(64-sh));
        }
        maxH=max(maxH,y+o.h);
        place[idx] = {x - o.minx, y - o.miny, o.R, o.F};
        for(auto [cx,cy]: o.cells){
            int gx=x+cx, gy=y+cy;
            addCand(gx+1,gy); addCand(gx-1,gy); addCand(gx,gy+1); addCand(gx,gy-1);
        }
        addCand(0,maxH);
    }
};

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<Piece> pieces(n);
    long long total=0; int maxMinW=1;
    for(int i=0;i<n;++i){
        cin>>pieces[i].k; total+=pieces[i].k; pieces[i].orig.resize(pieces[i].k);
        for(auto &p: pieces[i].orig) cin>>p.first>>p.second;
        pieces[i].ors=make_orients(pieces[i].orig);
        int mw=1000000; for(auto &o:pieces[i].ors) mw=min(mw,o.w); maxMinW=max(maxMinW,mw);
    }
    vector<int> ord(n); iota(ord.begin(),ord.end(),0);
    sort(ord.begin(), ord.end(), [&](int a,int b){
        const auto &pa=pieces[a], &pb=pieces[b];
        int besta=1000,bestb=1000, maxha=0,maxhb=0;
        for(auto&o:pa.ors){besta=min(besta,o.w*o.h); maxha=max(maxha,max(o.w,o.h));}
        for(auto&o:pb.ors){bestb=min(bestb,o.w*o.h); maxhb=max(maxhb,max(o.w,o.h));}
        if(pa.k!=pb.k) return pa.k>pb.k;
        if(besta!=bestb) return besta>bestb;
        return maxha>maxhb;
    });

    int root=max(1,(int)(sqrt((double)total)+0.5));
    vector<int> widths;
    double facs[] = {0.75,0.9,1.05,1.2,1.4,1.65};
    for(double f:facs){ int w=max(maxMinW,(int)llround(root*f)); widths.push_back(w); }
    widths.push_back(max(maxMinW, root)); widths.push_back(max(maxMinW, (int)ceil(sqrt((double)total*1.35))));
    sort(widths.begin(),widths.end()); widths.erase(unique(widths.begin(),widths.end()),widths.end());

    long long bestArea=LLONG_MAX; int bestW=0,bestH=0; vector<AnsPlace> bestPlace;
    auto start=chrono::steady_clock::now();
    int trials=0;
    for(int W: widths){
        ++trials;
        Packer pk(W,n);
        int step=0;
        for(int id: ord){
            if((step++ % 80)==0) sort(pk.cand.begin(), pk.cand.end(), [](auto&a,auto&b){return a.second==b.second?a.first<b.first:a.second<b.second;});
            long long bestScore=LLONG_MAX; int bx=0,by=pk.maxH, bo=0, bestCont=-1;
            auto consider = [&](int cx,int cy){
                for(int oi=0; oi<(int)pieces[id].ors.size(); ++oi){
                    const Orient &o=pieces[id].ors[oi];
                    if(cx+o.w>W) continue;
                    int nh=max(pk.maxH, cy+o.h);
                    if((long long)nh*W > bestArea) continue;
                    if(!pk.fits(o,cx,cy)) continue;
                    int cont=pk.contact(o,cx,cy);
                    long long sc=(long long)nh*100000000LL + (long long)cy*10000LL + cx*10LL - cont*1500LL;
                    if(sc<bestScore){ bestScore=sc; bx=cx; by=cy; bo=oi; bestCont=cont; }
                }
            };
            int L=(n>5000?700:1300);
            int m=pk.cand.size();
            for(int i=0;i<m && i<L;i++) consider(pk.cand[i].first, pk.cand[i].second);
            for(int i=max(0,m-L); i<m; ++i) consider(pk.cand[i].first, pk.cand[i].second);
            consider(0,pk.maxH);
            for(int x=0; x<W; x+=max(1,W/16)) consider(x,pk.maxH);
            if(bestScore==LLONG_MAX){
                // Always-valid fallback: a fresh row, using the narrowest orientation.
                int oi=0; for(int j=1;j<(int)pieces[id].ors.size();++j) if(pieces[id].ors[j].w < pieces[id].ors[oi].w) oi=j;
                bx=0; by=pk.maxH; bo=oi;
            }
            pk.put(pieces[id].ors[bo], bx, by, id);
        }
        long long area=(long long)W*pk.maxH;
        if(area<bestArea || (area==bestArea && (pk.maxH<bestH || (pk.maxH==bestH && W<bestW)))){
            bestArea=area; bestW=W; bestH=pk.maxH; bestPlace.swap(pk.place);
        }
        auto now=chrono::steady_clock::now();
        double elapsed=chrono::duration<double>(now-start).count();
        if(elapsed>1.65 && trials>=2) break;
    }
    cout << bestW << ' ' << bestH << '\n';
    for(int i=0;i<n;++i) cout << bestPlace[i].X << ' ' << bestPlace[i].Y << ' ' << bestPlace[i].R << ' ' << bestPlace[i].F << '\n';
    return 0;
}
