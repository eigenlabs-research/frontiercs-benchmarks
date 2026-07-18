#include <bits/stdc++.h>
using namespace std;

struct Cell { int x, y; };
struct Orient {
    vector<Cell> c;
    int w=0, h=0, r=0, f=0;
    int minx=0, miny=0;
};
struct Piece {
    vector<Cell> orig;
    vector<Orient> os;
    int k=0, id=0;
};
struct Place { long long X=0, Y=0; int R=0, F=0; };
struct Cand { long long area=LLONG_MAX, W=0, H=0; vector<Place> p; };

static inline pair<int,int> trans(int x,int y,int r,int f){
    if(f) x = -x;
    switch(r&3){
        case 0: return {x,y};
        case 1: return {y,-x};
        case 2: return {-x,-y};
        default: return {-y,x};
    }
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if(!(cin>>n)) return 0;
    vector<Piece> pieces(n);
    long long total=0;
    for(int i=0;i<n;i++){
        int k; cin>>k; pieces[i].k=k; pieces[i].id=i; total+=k;
        pieces[i].orig.resize(k);
        for(int j=0;j<k;j++) cin>>pieces[i].orig[j].x>>pieces[i].orig[j].y;
        set<vector<pair<int,int>>> seen;
        for(int f=0; f<2; ++f) for(int r=0; r<4; ++r){
            vector<pair<int,int>> v;
            int mnx=INT_MAX,mny=INT_MAX,mxx=INT_MIN,mxy=INT_MIN;
            for(auto cc: pieces[i].orig){
                auto q=trans(cc.x,cc.y,r,f);
                v.push_back(q);
                mnx=min(mnx,q.first); mny=min(mny,q.second);
                mxx=max(mxx,q.first); mxy=max(mxy,q.second);
            }
            sort(v.begin(), v.end());
            if(!seen.insert(v).second) continue;
            Orient o; o.r=r; o.f=f; o.minx=mnx; o.miny=mny; o.w=mxx-mnx+1; o.h=mxy-mny+1;
            for(auto q: v) o.c.push_back({q.first-mnx,q.second-mny});
            pieces[i].os.push_back(o);
        }
        sort(pieces[i].os.begin(), pieces[i].os.end(), [](const Orient&a,const Orient&b){
            if(a.w*a.h!=b.w*b.h) return a.w*a.h < b.w*b.h;
            if(max(a.w,a.h)!=max(b.w,b.h)) return max(a.w,a.h)<max(b.w,b.h);
            return a.h<b.h;
        });
    }

    vector<int> ord(n);
    iota(ord.begin(), ord.end(), 0);
    sort(ord.begin(), ord.end(), [&](int a,int b){
        auto bestMetric=[&](const Piece& p){
            int ba=1000000, bm=0, bh=0;
            for(auto &o:p.os){ ba=min(ba,o.w*o.h); bm=max(bm,max(o.w,o.h)); bh=max(bh,o.h); }
            return tuple<int,int,int,int>(p.k, bm, ba, bh);
        };
        return bestMetric(pieces[a]) > bestMetric(pieces[b]);
    });

    int minW=1;
    for(auto &p: pieces){
        int mw=1000000;
        for(auto &o:p.os) mw=min(mw,o.w);
        minW=max(minW,mw);
    }
    double root = sqrt((double)max(1LL,total));
    vector<int> widths;
    double mults[] = {1.00,1.12,0.90,1.28,0.78,1.50,0.65,1.80,2.20};
    for(double m: mults) widths.push_back(max(minW, (int)llround(root*m)));
    for(int w: {16,24,32,48,64,96,128,192,256,384,512}) if(w>=minW) widths.push_back(w);
    sort(widths.begin(), widths.end());
    widths.erase(unique(widths.begin(), widths.end()), widths.end());
    stable_sort(widths.begin(), widths.end(), [&](int a,int b){ return abs(a-root)<abs(b-root); });

    auto start = chrono::steady_clock::now();
    auto elapsed = [&](){ return chrono::duration<double>(chrono::steady_clock::now()-start).count(); };
    Cand best;

    auto packWidth = [&](int W)->Cand{
        Cand res; res.W=W; res.H=0; res.p.assign(n, {});
        int words=(W+63)>>6;
        vector<vector<unsigned long long>> rows;
        vector<int> colH(W,0);
        long long maxX=0,maxY=0;
        auto occupied = [&](int x,int y)->bool{
            if(y<0 || x<0 || x>=W) return true;
            if(y>=(int)rows.size()) return false;
            return (rows[y][x>>6]>>(x&63))&1ULL;
        };
        auto canPlace = [&](const Orient&o,int px,int py)->bool{
            if(px<0 || px+o.w>W || py<0) return false;
            for(auto cc:o.c) if(occupied(px+cc.x, py+cc.y)) return false;
            return true;
        };
        for(int idx=0; idx<n; ++idx){
            if((idx&127)==0 && elapsed()>1.86){ res.area=LLONG_MAX; return res; }
            Piece &pc = pieces[ord[idx]];
            long long bestNH=LLONG_MAX; int bestX=0,bestY=0,bestContact=-1; const Orient* bestO=nullptr;
            auto contactScore = [&](const Orient&o,int px,int py)->int{
                int s=0;
                for(auto cc:o.c){
                    int gx=px+cc.x, gy=py+cc.y;
                    const int dx[4]={1,-1,0,0}, dy[4]={0,0,1,-1};
                    for(int d=0; d<4; ++d){
                        int nx=gx+dx[d], ny=gy+dy[d];
                        if(nx<0 || nx>=W || ny<0) ++s;
                        else if(occupied(nx,ny)) ++s;
                    }
                }
                return s;
            };
            for(const auto &o: pc.os){
                if(o.w>W) continue;
                for(int x=0; x+o.w<=W; ++x){
                    int y=0;
                    for(auto cc:o.c) y=max(y, colH[x+cc.x]-cc.y);
                    if(y<0) y=0;
                    // Usually skyline y is feasible; if holes/overhangs collide, lift a little.
                    while(!canPlace(o,x,y)){
                        ++y;
                        if(max<long long>(maxY, y+o.h)>bestNH) break;
                    }
                    long long nh=max<long long>(maxY, y+o.h);
                    if(nh>bestNH || !canPlace(o,x,y)) continue;
                    int cs=contactScore(o,x,y);
                    if(nh < bestNH || cs > bestContact || (cs==bestContact && (y<bestY || (y==bestY && x<bestX)))){
                        bestNH=nh; bestX=x; bestY=y; bestO=&o; bestContact=cs;
                    }
                }
            }
            if(!bestO){ res.area=LLONG_MAX; return res; }
            if((int)rows.size() < bestY + bestO->h) rows.resize(bestY+bestO->h, vector<unsigned long long>(words,0));
            for(auto cc: bestO->c){
                int gx=bestX+cc.x, gy=bestY+cc.y;
                rows[gy][gx>>6] |= 1ULL<<(gx&63);
                colH[gx]=max(colH[gx], gy+1);
                maxX=max<long long>(maxX,gx+1); maxY=max<long long>(maxY,gy+1);
            }
            Place pl;
            pl.X = (long long)bestX - bestO->minx;
            pl.Y = (long long)bestY - bestO->miny;
            pl.R = bestO->r; pl.F = bestO->f;
            res.p[pc.id]=pl;
        }
        res.W=maxX; res.H=maxY; res.area=maxX*maxY;
        return res;
    };

    for(int W: widths){
        if(elapsed()>1.82) break;
        Cand c=packWidth(W);
        if(c.area<best.area || (c.area==best.area && (c.H<best.H || (c.H==best.H && c.W<best.W)))) best=std::move(c);
    }

    // Extremely safe fallback: put normalized pieces one after another in a single row.
    if(best.area==LLONG_MAX){
        best.p.assign(n,{}); long long x=0,H=1;
        for(int i=0;i<n;i++){
            const Orient &o=pieces[i].os[0];
            best.p[i]={(long long)x-o.minx, (long long)0-o.miny, o.r, o.f};
            x+=o.w; H=max<long long>(H,o.h);
        }
        best.W=x; best.H=H; best.area=x*H;
    }
    cout << best.W << ' ' << best.H << '\n';
    for(int i=0;i<n;i++) cout << best.p[i].X << ' ' << best.p[i].Y << ' ' << best.p[i].R << ' ' << best.p[i].F << '\n';
    return 0;
}
