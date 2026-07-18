#include <bits/stdc++.h>
using namespace std;

struct Orient {
    int R, F;
    int w, h;
    long long minx, miny;
    vector<pair<int,int>> cells; // normalized cells
};
struct Piece {
    int id, k;
    vector<pair<long long,long long>> raw;
    vector<Orient> ors;
};
struct Ans { long long X=0,Y=0; int R=0,F=0; };

static pair<long long,long long> rotcw(long long x,long long y,int r){
    switch(r&3){
        case 0: return {x,y};
        case 1: return {y,-x};
        case 2: return {-x,-y};
        default: return {-y,x};
    }
}

static vector<Orient> gen_orients(const vector<pair<long long,long long>>& raw){
    vector<Orient> res;
    set<vector<pair<int,int>>> seen;
    for(int F=0; F<=1; ++F){
        for(int R=0; R<4; ++R){
            vector<pair<long long,long long>> t;
            long long mnx=LLONG_MAX,mny=LLONG_MAX,mxx=LLONG_MIN,mxy=LLONG_MIN;
            for(auto [x,y]: raw){
                if(F) x=-x;
                auto [qx,qy]=rotcw(x,y,R);
                t.push_back({qx,qy});
                mnx=min(mnx,qx); mny=min(mny,qy); mxx=max(mxx,qx); mxy=max(mxy,qy);
            }
            vector<pair<int,int>> c;
            for(auto [qx,qy]: t) c.push_back({(int)(qx-mnx),(int)(qy-mny)});
            sort(c.begin(), c.end());
            if(seen.insert(c).second){
                Orient o; o.R=R; o.F=F; o.minx=mnx; o.miny=mny; o.w=(int)(mxx-mnx+1); o.h=(int)(mxy-mny+1); o.cells=c;
                res.push_back(o);
            }
        }
    }
    sort(res.begin(), res.end(), [](const Orient& a,const Orient& b){
        int aa=a.w*a.h, bb=b.w*b.h;
        if(aa!=bb) return aa<bb;
        if(max(a.w,a.h)!=max(b.w,b.h)) return max(a.w,a.h)<max(b.w,b.h);
        return a.h<b.h;
    });
    return res;
}

struct PackResult { long long W,H,area; vector<Ans> ans; };

static PackResult pack_width(const vector<Piece>& pieces, const vector<int>& order, int W){
    int n=pieces.size();
    vector<Ans> ans(n);
    vector<int> colH(W,0);
    vector<vector<unsigned char>> grid;
    grid.reserve(1024);
    long long curH=0;

    auto ensureH = [&](int h){
        while((int)grid.size()<h) grid.push_back(vector<unsigned char>(W,0));
    };
    auto fits = [&](const Orient& o,int x,int y)->bool{
        if(x<0 || x+o.w>W || y<0) return false;
        if(y+o.h <= (int)grid.size()){
            for(auto [dx,dy]: o.cells) if(grid[y+dy][x+dx]) return false;
        } else {
            for(auto [dx,dy]: o.cells) if(y+dy < (int)grid.size() && grid[y+dy][x+dx]) return false;
        }
        return true;
    };
    auto contact = [&](const Orient& o,int x,int y)->int{
        int s=0;
        for(auto [dx,dy]: o.cells){
            int gx=x+dx, gy=y+dy;
            if(gy==0) s++;
            if(gx==0) s++;
            if(gx+1==W) s++;
            static int dd[4][2]={{1,0},{-1,0},{0,1},{0,-1}};
            for(auto &d: dd){
                int nx=gx+d[0], ny=gy+d[1];
                if(nx>=0 && nx<W && ny>=0 && ny<(int)grid.size() && grid[ny][nx]) s++;
            }
        }
        return s;
    };

    for(int idx: order){
        const Piece& p=pieces[idx];
        bool found=false;
        const Orient* bestO=nullptr;
        int bestX=0,bestY=0,bestContact=-1;
        long long bestH=LLONG_MAX;
        int bestWaste=INT_MAX;
        for(const Orient& o: p.ors){
            if(o.w>W) continue;
            for(int x=0; x+o.w<=W; ++x){
                int y=0;
                for(auto [dx,dy]: o.cells) y=max(y, colH[x+dx]-dy);
                if(y<0) y=0;
                while(!fits(o,x,y)) ++y;
                long long nh=max(curH, (long long)y+o.h);
                int waste=0;
                for(int xx=x; xx<x+o.w; ++xx) waste += max(0, y - colH[xx]);
                int con=0;
                if(nh < bestH || (nh==bestH && (y<bestY || (y==bestY && (waste<bestWaste || (waste==bestWaste && (con=contact(o,x,y))>bestContact)))))){
                    if(con==0) con=contact(o,x,y);
                    found=true; bestO=&o; bestX=x; bestY=y; bestH=nh; bestWaste=waste; bestContact=con;
                }
            }
        }
        if(!found){ // should not happen if W is sane; append with first orientation by enlarging height conceptually
            bestO=&p.ors[0]; bestX=0; bestY=(int)curH; bestH=curH+bestO->h;
        }
        ensureH(bestY+bestO->h);
        for(auto [dx,dy]: bestO->cells){
            int gx=bestX+dx, gy=bestY+dy;
            grid[gy][gx]=1;
            colH[gx]=max(colH[gx], gy+1);
        }
        curH=max(curH, (long long)bestY+bestO->h);
        ans[p.id] = { (long long)bestX - bestO->minx, (long long)bestY - bestO->miny, bestO->R, bestO->F };
    }
    return {(long long)W, curH, (long long)W*curH, move(ans)};
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<Piece> pieces(n);
    long long total=0;
    for(int i=0;i<n;++i){
        int k; cin>>k; total+=k;
        pieces[i].id=i; pieces[i].k=k; pieces[i].raw.resize(k);
        for(int j=0;j<k;++j) cin>>pieces[i].raw[j].first>>pieces[i].raw[j].second;
        pieces[i].ors=gen_orients(pieces[i].raw);
    }
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    vector<int> hard(n,0), minw(n,1000000), minh(n,1000000), maxbox(n,0);
    int minRequiredW=1;
    for(int i=0;i<n;++i){
        for(auto &o: pieces[i].ors){ minw[i]=min(minw[i],o.w); minh[i]=min(minh[i],o.h); maxbox[i]=max(maxbox[i], o.w*o.h); }
        minRequiredW=max(minRequiredW, minw[i]);
        hard[i]=maxbox[i]*20 - pieces[i].k*10 + max(minw[i],minh[i]);
    }
    sort(order.begin(), order.end(), [&](int a,int b){
        if(pieces[a].k!=pieces[b].k) return pieces[a].k>pieces[b].k;
        if(hard[a]!=hard[b]) return hard[a]>hard[b];
        return a<b;
    });

    double root=sqrt((double)max(1LL,total));
    vector<int> widths;
    auto addW=[&](int w){
        w=max(w,minRequiredW);
        if(w<1) w=1;
        widths.push_back(w);
    };
    double factors[] = {0.55,0.70,0.85,1.00,1.12,1.25,1.42,1.62,1.85,2.15,2.55,3.10};
    for(double f: factors) addW((int)llround(root*f));
    addW((int)ceil(root)); addW(minRequiredW); addW(max(minRequiredW, (int)root/2));
    sort(widths.begin(), widths.end()); widths.erase(unique(widths.begin(), widths.end()), widths.end());

    auto start=chrono::steady_clock::now();
    bool have=false; PackResult best;
    for(int wi=0; wi<(int)widths.size(); ++wi){
        auto now=chrono::steady_clock::now();
        if(have && chrono::duration<double>(now-start).count()>1.75) break;
        PackResult r=pack_width(pieces, order, widths[wi]);
        if(!have || r.area<best.area || (r.area==best.area && (r.H<best.H || (r.H==best.H && r.W<best.W)))){
            best=move(r); have=true;
        }
    }
    if(!have) best=pack_width(pieces, order, minRequiredW);
    cout << best.W << ' ' << best.H << '\n';
    for(int i=0;i<n;++i) cout << best.ans[i].X << ' ' << best.ans[i].Y << ' ' << best.ans[i].R << ' ' << best.ans[i].F << '\n';
    return 0;
}
