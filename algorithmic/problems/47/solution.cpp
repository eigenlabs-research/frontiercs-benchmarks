#include <bits/stdc++.h>
using namespace std;

struct JsonParser{
    string s; size_t p=0;
    JsonParser(const string& ss):s(ss){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    void expect(char c){ ws(); if(p>=s.size()||s[p]!=c) exit(0); p++; }
    string str(){ ws(); expect('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\' && p<s.size()){ char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); int sign=1; if(s[p]=='-'){sign=-1;p++;} long long v=0; while(p<s.size() && isdigit((unsigned char)s[p])) v=v*10+(s[p++]-'0'); return sign*v; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} p+=5; return false; }
};

struct Item{ string type; int w,h; long long v; int limit; };
struct Rect{ int x,y,w,h; };
struct Placement{ int id,x,y,rot,w,h; };
static inline bool inter(const Rect&a,const Rect&b){ return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h; }
static inline bool contains(const Rect&a,const Rect&b){ return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h; }

int W,H; bool allowRot; vector<Item> items;
chrono::steady_clock::time_point gStart;

void parseInput(){
    string input((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    JsonParser jp(input);
    jp.expect('{');
    for(int top=0; top<2; ++top){
        if(top) jp.expect(',');
        string k=jp.str(); jp.expect(':');
        if(k=="bin"){
            jp.expect('{');
            for(int i=0;i<3;i++){
                if(i) jp.expect(','); string bk=jp.str(); jp.expect(':');
                if(bk=="W") W=(int)jp.num(); else if(bk=="H") H=(int)jp.num(); else if(bk=="allow_rotate") allowRot=jp.boolean();
            }
            jp.expect('}');
        }else if(k=="items"){
            jp.expect('['); jp.ws();
            while(jp.p<jp.s.size() && jp.s[jp.p]!=']'){
                if(!items.empty()) jp.expect(',');
                jp.expect('{'); Item it; int seen=0;
                for(int f=0; f<5; ++f){
                    if(f) jp.expect(','); string key=jp.str(); jp.expect(':');
                    if(key=="type") it.type=jp.str(); else if(key=="w") it.w=(int)jp.num(); else if(key=="h") it.h=(int)jp.num(); else if(key=="v") it.v=jp.num(); else if(key=="limit") it.limit=(int)jp.num(); seen++;
                }
                jp.expect('}'); items.push_back(it); jp.ws();
            }
            jp.expect(']');
        }
    }
}

void prune(vector<Rect>& fr){
    // remove non-positive, contained duplicates; keep a moderate number of useful free rectangles
    vector<Rect> a; a.reserve(fr.size());
    for(auto &r:fr) if(r.w>0 && r.h>0) a.push_back(r);
    fr.swap(a);
    int n=fr.size(); vector<char> del(n,0);
    for(int i=0;i<n;i++) if(!del[i]){
        for(int j=0;j<n;j++) if(i!=j && !del[j]){
            if(contains(fr[i],fr[j])) del[j]=1;
        }
    }
    a.clear();
    for(int i=0;i<n;i++) if(!del[i]) a.push_back(fr[i]);
    fr.swap(a);
    if(fr.size()>3500){
        sort(fr.begin(), fr.end(), [](const Rect&a,const Rect&b){
            long long aa=1LL*a.w*a.h, bb=1LL*b.w*b.h;
            if(aa!=bb) return aa>bb; return min(a.w,a.h)>min(b.w,b.h);
        });
        fr.resize(3500);
    }
}

vector<Placement> packVariant(int variant){
    int M=items.size(); vector<int> rem(M); for(int i=0;i<M;i++) rem[i]=items[i].limit;
    vector<Rect> freeR(1, {0,0,W,H}); vector<Placement> out; out.reserve(4096);
    long long binArea=1LL*W*H;
    int stagnant=0;
    while(true){
        if((out.size()&63)==0){
            auto ms=chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-gStart).count();
            if(ms>900) break;
        }
        long double best=-1e100; int bi=-1,bf=-1,brot=0,bw=0,bh=0;
        long long bestWaste=0; int bestY=0,bestX=0;
        for(int fi=0; fi<(int)freeR.size(); ++fi){
            const Rect&r=freeR[fi];
            for(int i=0;i<M;i++) if(rem[i]>0){
                for(int ro=0; ro<=(allowRot?1:0); ++ro){
                    int iw= ro?items[i].h:items[i].w, ih= ro?items[i].w:items[i].h;
                    if(ro && iw==items[i].w && ih==items[i].h) continue;
                    if(iw>r.w || ih>r.h) continue;
                    long double area=(long double)iw*ih;
                    long double dens=(long double)items[i].v/area;
                    long long waste=1LL*r.w*r.h-(long long)iw*ih;
                    int dw=r.w-iw, dh=r.h-ih;
                    int exact=(dw==0)+(dh==0);
                    long double sc;
                    if(variant==0) sc=dens*1e12L + exact*1e9L - min(dw,dh)*1000.0L - waste*0.01L;
                    else if(variant==1) sc=dens*1e12L - waste*10.0L + exact*1e8L;
                    else if(variant==2) sc=(long double)items[i].v*1e6L + dens*1e8L - waste*0.001L;
                    else if(variant==3) sc=dens*1e12L + (long double)area*100.0L - max(dw,dh)*1000.0L;
                    else if(variant==4) sc=dens*1e12L - (long double)(dw*dw + dh*dh)*100.0L + exact*1e9L;
                    else sc=(long double)items[i].v*1e5L - waste*1.0L + dens*1e9L;
                    // deterministic spatial tie breakers vary slightly by variant
                    if(variant%2==0) sc -= r.y*0.1L + r.x*0.0001L; else sc -= r.x*0.1L + r.y*0.0001L;
                    if(sc>best || (fabsl(sc-best)<1e-9 && (r.y<bestY || (r.y==bestY && r.x<bestX)))){
                        best=sc; bi=i; bf=fi; brot=ro; bw=iw; bh=ih; bestWaste=waste; bestY=r.y; bestX=r.x;
                    }
                }
            }
        }
        if(bi<0) break;
        Rect placed{freeR[bf].x, freeR[bf].y, bw, bh};
        out.push_back({bi,placed.x,placed.y,brot,bw,bh}); rem[bi]--;
        vector<Rect> nf; nf.reserve(freeR.size()+8);
        for(const Rect& r: freeR){
            if(!inter(r, placed)){ nf.push_back(r); continue; }
            int rx2=r.x+r.w, ry2=r.y+r.h, px2=placed.x+placed.w, py2=placed.y+placed.h;
            if(placed.x>r.x) nf.push_back({r.x,r.y,placed.x-r.x,r.h});
            if(px2<rx2) nf.push_back({px2,r.y,rx2-px2,r.h});
            if(placed.y>r.y) nf.push_back({r.x,r.y,r.w,placed.y-r.y});
            if(py2<ry2) nf.push_back({r.x,py2,r.w,ry2-py2});
        }
        freeR.swap(nf);
        prune(freeR);
        if(out.size()>=25000) break;
        if(freeR.empty()) break;
    }
    return out;
}

long long valueOf(const vector<Placement>& p){ long long s=0; for(auto &q:p) s+=items[q.id].v; return s; }

vector<Placement> shelfPack(vector<int> order, bool useRot){
    vector<int> rem(items.size()); for(int i=0;i<(int)items.size();++i) rem[i]=items[i].limit;
    vector<Placement> out; int x=0,y=0,rowH=0;
    for(int id: order){
        while(rem[id]>0){
            int bestRot=0, bw=items[id].w, bh=items[id].h;
            if(useRot && allowRot && items[id].h<=W && items[id].w<=H){
                // Prefer the orientation that fits the current shelf, otherwise the shorter height.
                bool f0=(x+items[id].w<=W && y+max(rowH,items[id].h)<=H);
                bool f1=(x+items[id].h<=W && y+max(rowH,items[id].w)<=H);
                if((!f0 && f1) || (f1 && items[id].w<items[id].h)){ bestRot=1; bw=items[id].h; bh=items[id].w; }
            }
            if(bw>W || bh>H) break;
            if(x+bw>W){ y+=rowH; x=0; rowH=0; }
            if(y+bh>H) break;
            out.push_back({id,x,y,bestRot,bw,bh}); x+=bw; rowH=max(rowH,bh); rem[id]--;
            if(out.size()>=25000) return out;
        }
    }
    return out;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    gStart=chrono::steady_clock::now();
    parseInput();
    vector<Placement> best; long long bestVal=-1;
    vector<int> ord(items.size()); iota(ord.begin(), ord.end(), 0);
    auto tryCand = [&](vector<Placement> cand){ long long val=valueOf(cand); if(val>bestVal){ bestVal=val; best.swap(cand); } };
    tryCand(shelfPack(ord,false));
    sort(ord.begin(), ord.end(), [](int a,int b){ return (long double)items[a].v/(items[a].w*items[a].h) > (long double)items[b].v/(items[b].w*items[b].h); });
    tryCand(shelfPack(ord,allowRot));
    int variants=6;
    for(int v=0; v<variants; ++v){
        auto ms=chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-gStart).count();
        if(ms>920 && bestVal>=0) break;
        auto cand=packVariant(v); tryCand(cand);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){
        if(i) cout << ',';
        const auto&p=best[i];
        cout << "{\"type\":\"";
        for(char c: items[p.id].type){ if(c=='"'||c=='\\') cout << '\\'; cout << c; }
        cout << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}";
    }
    cout << "]}\n";
    return 0;
}
