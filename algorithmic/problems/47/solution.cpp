#include <bits/stdc++.h>
using namespace std;

struct JsonParser{
    string s; size_t p=0;
    JsonParser(const string& t):s(t){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    bool eat(char c){ ws(); if(p<s.size()&&s[p]==c){p++; return true;} return false; }
    void expect(char c){ if(!eat(c)) exit(0); }
    string str(){ ws(); expect('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\' && p<s.size()){ char e=s[p++]; if(e=='"'||e=='\\'||e=='/') r+=e; else if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); int sign=1; if(p<s.size()&&(s[p]=='-'||s[p]=='+')){ if(s[p]=='-') sign=-1; p++; } long long v=0; while(p<s.size()&&isdigit((unsigned char)s[p])) v=v*10+(s[p++]-'0'); return sign*v; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} if(s.compare(p,5,"false")==0){p+=5; return false;} return false; }
};

struct Item{ string id; int w,h; long long v; int lim; };
struct Rect{ int x,y,w,h; };
struct Pl{ int idx,x,y,rot,w,h; };
int W,H; bool allowRot; vector<Item> items;

static inline long long area(const Rect& r){ return 1LL*r.w*r.h; }
static inline bool inter(const Rect&a,const Rect&b){ return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }
static inline bool contains(const Rect&a,const Rect&b){ return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h; }

struct Packer{
    vector<Rect> freeR;
    vector<Pl> out;
    vector<int> used;
    long long val=0;
    int placeMode=0;
    Packer(int pm=0):used(items.size(),0),placeMode(pm){ freeR.push_back({0,0,W,H}); }
    
    bool bestSpot(int w,int h,int &bx,int &by,long long &k1,long long &k2,long long &k3){
        bool ok=false; k1=k2=k3=LLONG_MAX; bx=by=0;
        for(auto &r: freeR) if(w<=r.w && h<=r.h){
            long long waste=1LL*r.w*r.h-1LL*w*h;
            long long ss=min(r.w-w,r.h-h), ls=max(r.w-w,r.h-h);
            long long a,b,c;
            if(placeMode==0){ a=ss; b=ls; c=waste; }          // best short side fit
            else if(placeMode==1){ a=waste; b=ss; c=ls; }     // best area fit
            else if(placeMode==2){ a=r.y; b=r.x; c=ss; }      // bottom-left
            else { a=max(r.w-w,r.h-h); b=waste; c=r.y; }      // best long side fit
            if(!ok || tuple<long long,long long,long long,int,int>(a,b,c,r.y,r.x) < tuple<long long,long long,long long,int,int>(k1,k2,k3,by,bx)){
                ok=true; k1=a; k2=b; k3=c; bx=r.x; by=r.y;
            }
        }
        return ok;
    }
    void prune(){
        vector<Rect> v;
        for(auto r: freeR) if(r.w>0&&r.h>0) v.push_back(r);
        vector<char> rem(v.size(),0);
        for(size_t i=0;i<v.size();++i) if(!rem[i]) for(size_t j=0;j<v.size();++j) if(i!=j&&!rem[j]){
            if(contains(v[i],v[j])) rem[j]=1;
        }
        freeR.clear();
        for(size_t i=0;i<v.size();++i) if(!rem[i]) freeR.push_back(v[i]);
        sort(freeR.begin(), freeR.end(), [](const Rect&a,const Rect&b){
            if(a.y!=b.y) return a.y<b.y; if(a.x!=b.x) return a.x<b.x; return area(a)>area(b);
        });
        if(freeR.size()>700) freeR.resize(700);
    }
    void place(int idx,int x,int y,int rot,int w,int h){
        Rect p{x,y,w,h}; vector<Rect> nf;
        for(auto r: freeR){
            if(!inter(r,p)){ nf.push_back(r); continue; }
            int rx2=r.x+r.w, ry2=r.y+r.h, px2=p.x+p.w, py2=p.y+p.h;
            if(p.x>r.x) nf.push_back({r.x,r.y,p.x-r.x,r.h});
            if(px2<rx2) nf.push_back({px2,r.y,rx2-px2,r.h});
            if(p.y>r.y) nf.push_back({r.x,r.y,r.w,p.y-r.y});
            if(py2<ry2) nf.push_back({r.x,py2,r.w,ry2-py2});
        }
        freeR.swap(nf); prune();
        used[idx]++; val+=items[idx].v; out.push_back({idx,x,y,rot,w,h});
    }
};

vector<pair<int,int>> orientations(int i){
    vector<pair<int,int>> o; o.push_back({items[i].w,items[i].h});
    if(allowRot && items[i].w!=items[i].h) o.push_back({items[i].h,items[i].w});
    return o;
}

Packer dynamicPack(int scoreMode,int pm){
    Packer pk(pm); int n=items.size();
    const int hardSteps=20000;
    for(int step=0; step<hardSteps; ++step){
        int bi=-1, brot=0, bx=0, by=0, bw=0,bh=0; double best=-1e300; long long bk1=0,bk2=0,bk3=0;
        for(int i=0;i<n;i++) if(pk.used[i]<items[i].lim && items[i].v>0){
            auto os=orientations(i);
            for(int oi=0; oi<(int)os.size(); ++oi){
                int w=os[oi].first,h=os[oi].second; long long k1,k2,k3; int x,y;
                if(!pk.bestSpot(w,h,x,y,k1,k2,k3)) continue;
                double ar=double(w)*h, den=double(items[i].v)/ar;
                double sc;
                if(scoreMode==0) sc=den*1e9 + items[i].v*1e-3 - k1*10 - k2*1e-3;
                else if(scoreMode==1) sc=items[i].v - k1*1000.0 - k2*0.01;
                else if(scoreMode==2) sc=den*7e8 + double(items[i].v)*0.15 - k3*den*0.02;
                else if(scoreMode==3) sc=double(items[i].v) / (ar + max(0LL,k3)*0.18) * 1e9;
                else if(scoreMode==4) sc=den*5e8 + sqrt((double)items[i].v)*1000 - (w+h)*10 - k1;
                else sc=items[i].v*0.7 + den*2e8 - abs(w-h)*100.0 - k3*1e-4;
                // Prefer placements lower/left and candidates consuming scarce awkward free rectangles.
                sc -= y*1e-5 + x*1e-6;
                if(sc>best){ best=sc; bi=i; brot=oi; bx=x; by=y; bw=w; bh=h; bk1=k1; bk2=k2; bk3=k3; }
            }
        }
        if(bi<0) break;
        pk.place(bi,bx,by,brot,bw,bh);
    }
    return pk;
}

Packer shelfPack(int mode, bool vertical){
    Packer pk(0); int n=items.size();
    struct Cand{int id,w,h,rot;};
    vector<int> sizes;
    vector<Cand> all;
    for(int i=0;i<n;i++){
        auto os=orientations(i);
        for(int r=0;r<(int)os.size();++r){
            int w=os[r].first,h=os[r].second;
            if(w<=W && h<=H){ all.push_back({i,w,h,r}); sizes.push_back(vertical?w:h); }
        }
    }
    sort(sizes.begin(), sizes.end()); sizes.erase(unique(sizes.begin(), sizes.end()), sizes.end());
    auto makeStrip = [&](int sz, vector<int> cnt, bool commit){
        vector<Pl> row; long long value=0; int pos=0, limitLen=vertical?H:W;
        while(pos<limitLen){
            int best=-1; double bsc=-1e300;
            for(int ci=0;ci<(int)all.size();++ci){
                auto &c=all[ci]; if(cnt[c.id]>=items[c.id].lim) continue;
                int along=vertical?c.h:c.w, thick=vertical?c.w:c.h;
                if(thick>sz || pos+along>limitLen) continue;
                double waste = (double)(sz-thick)*along;
                double sc;
                if(mode==0) sc=(double)items[c.id].v/along;
                else if(mode==1) sc=(double)items[c.id].v/(c.w*(double)c.h);
                else if(mode==2) sc=(double)items[c.id].v/(along + 0.35*(sz-thick));
                else sc=(double)items[c.id].v - waste*0.01;
                sc += (items[c.id].lim-cnt[c.id])*1e-9;
                if(sc>bsc){ bsc=sc; best=ci; }
            }
            if(best<0) break;
            auto c=all[best];
            if(vertical) row.push_back({c.id,0,pos,c.rot,c.w,c.h});
            else row.push_back({c.id,pos,0,c.rot,c.w,c.h});
            cnt[c.id]++; value+=items[c.id].v; pos += vertical?c.h:c.w;
        }
        if(commit){ for(auto &p: row) pk.place(p.idx,p.x,p.y,p.rot,p.w,p.h); }
        return pair<long long,vector<Pl>>(value,row);
    };
    int offset=0, maxThick=vertical?W:H;
    while(offset<maxThick){
        int bestSz=0; long long bestVal=0; vector<Pl> bestRow;
        for(int sz: sizes) if(offset+sz<=maxThick){
            auto sim = makeStrip(sz, pk.used, false);
            if(sim.first==0) continue;
            double lhs=(double)sim.first/(mode==3?1.0:sz), rhs=(double)bestVal/(mode==3?1.0:max(1,bestSz));
            if(bestSz==0 || lhs>rhs){ bestSz=sz; bestVal=sim.first; bestRow=std::move(sim.second); }
        }
        if(bestSz==0) break;
        for(auto p: bestRow){ if(vertical) p.x=offset; else p.y=offset; pk.place(p.idx,p.x,p.y,p.rot,p.w,p.h); }
        offset += bestSz;
    }
    return pk;
}

Packer orderPack(int orderMode,int pm){
    Packer pk(pm); int n=items.size(); vector<int> ord(n); iota(ord.begin(),ord.end(),0);
    sort(ord.begin(),ord.end(),[&](int a,int b){
        double aa=1.0*items[a].w*items[a].h, bb=1.0*items[b].w*items[b].h;
        if(orderMode==0) return items[a].v/aa > items[b].v/bb;
        if(orderMode==1) return items[a].v > items[b].v;
        if(orderMode==2) return aa > bb;
        if(orderMode==3) return max(items[a].w,items[a].h) > max(items[b].w,items[b].h);
        return min(items[a].w,items[a].h) < min(items[b].w,items[b].h);
    });
    bool progress=true; int loops=0;
    while(progress && loops++<3){
        progress=false;
        for(int id: ord){
            while(pk.used[id]<items[id].lim){
                int bestOi=-1,bx=0,by=0,bw=0,bh=0; long long bestTuple1=LLONG_MAX,bestTuple2=LLONG_MAX,bestTuple3=LLONG_MAX;
                auto os=orientations(id);
                for(int oi=0; oi<(int)os.size(); ++oi){
                    long long k1,k2,k3; int x,y,w=os[oi].first,h=os[oi].second;
                    if(pk.bestSpot(w,h,x,y,k1,k2,k3) && tuple<long long,long long,long long>(k1,k2,k3)<tuple<long long,long long,long long>(bestTuple1,bestTuple2,bestTuple3)){
                        bestOi=oi; bx=x; by=y; bw=w; bh=h; bestTuple1=k1; bestTuple2=k2; bestTuple3=k3;
                    }
                }
                if(bestOi<0) break;
                pk.place(id,bx,by,bestOi,bw,bh); progress=true;
            }
        }
    }
    return pk;
}

void parseInput(){
    string s((istreambuf_iterator<char>(cin)), {}); JsonParser jp(s); jp.expect('{');
    for(int top=0; top<2; ++top){ if(top) jp.expect(','); string key=jp.str(); jp.expect(':');
        if(key=="bin"){
            jp.expect('{'); for(int i=0;i<3;i++){ if(i) jp.expect(','); string k=jp.str(); jp.expect(':'); if(k=="W") W=jp.num(); else if(k=="H") H=jp.num(); else if(k=="allow_rotate") allowRot=jp.boolean(); } jp.expect('}');
        }else if(key=="items"){
            jp.expect('['); bool first=true; while(true){ jp.ws(); if(jp.eat(']')) break; if(!first) jp.expect(','); first=false; jp.expect('{'); Item it; for(int f=0; f<5; ++f){ if(f) jp.expect(','); string k=jp.str(); jp.expect(':'); if(k=="type") it.id=jp.str(); else if(k=="w") it.w=jp.num(); else if(k=="h") it.h=jp.num(); else if(k=="v") it.v=jp.num(); else if(k=="limit") it.lim=jp.num(); } jp.expect('}'); items.push_back(it); }
        }
    }
}

string esc(const string& s){ string r; for(char c:s){ if(c=='"'||c=='\\') r+='\\',r+=c; else r+=c; } return r; }

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    parseInput();
    Packer best(0);
    for(int mode=0; mode<4; ++mode){
        auto p=shelfPack(mode,false); if(p.val>best.val) best=std::move(p);
        auto q=shelfPack(mode,true); if(q.val>best.val) best=std::move(q);
    }
    for(int pm=0; pm<4; ++pm){
        for(int sm=0; sm<6; ++sm){ auto p=dynamicPack(sm,pm); if(p.val>best.val) best=std::move(p); }
        for(int om=0; om<5; ++om){ auto p=orderPack(om,pm); if(p.val>best.val) best=std::move(p); }
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.out.size();++i){
        if(i) cout << ',';
        const auto &p=best.out[i];
        cout << "{\"type\":\"" << esc(items[p.idx].id) << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}";
    }
    cout << "]}\n";
    return 0;
}
