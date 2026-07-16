#include <bits/stdc++.h>
using namespace std;

struct Item{string type; int w,h; long long v; int limit;};
struct Placement{int id,x,y,rot,w,h;};
struct Rect{int x,y,w,h;};
static chrono::steady_clock::time_point T0;
static inline bool timeUp(){ return chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-T0).count() > 880; }

struct Parser{
    string s; size_t p=0;
    Parser(const string& t):s(t){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) p++; }
    string str(){ ws(); ch('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\' && p<s.size()){ char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); long long sign=1,val=0; if(p<s.size()&&(s[p]=='-'||s[p]=='+')){ if(s[p]=='-') sign=-1; p++; } while(p<s.size() && isdigit((unsigned char)s[p])) val=val*10+(s[p++]-'0'); return sign*val; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} if(s.compare(p,5,"false")==0){p+=5; return false;} return false; }
};

static inline bool inter(const Rect&a,const Rect&b){return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h;}
static inline bool contains(const Rect&a,const Rect&b){return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h;}

void splitFree(vector<Rect>& freeR, const Rect& used){
    vector<Rect> nr; nr.reserve(freeR.size()+8);
    for(auto fr: freeR){
        if(!inter(fr,used)){ nr.push_back(fr); continue; }
        int fx2=fr.x+fr.w, fy2=fr.y+fr.h, ux2=used.x+used.w, uy2=used.y+used.h;
        if(used.x>fr.x) nr.push_back({fr.x,fr.y,used.x-fr.x,fr.h});
        if(ux2<fx2) nr.push_back({ux2,fr.y,fx2-ux2,fr.h});
        if(used.y>fr.y) nr.push_back({fr.x,fr.y,fr.w,used.y-fr.y});
        if(uy2<fy2) nr.push_back({fr.x,uy2,fr.w,fy2-uy2});
    }
    freeR.clear();
    for(auto r:nr) if(r.w>0 && r.h>0) freeR.push_back(r);
    // Remove contained rectangles and exact duplicates. Small problem sizes make O(n^2) acceptable.
    vector<char> dead(freeR.size(),0);
    for(size_t i=0;i<freeR.size();++i) if(!dead[i]){
        for(size_t j=i+1;j<freeR.size();++j) if(!dead[j]){
            if(contains(freeR[i],freeR[j])) dead[j]=1;
            else if(contains(freeR[j],freeR[i])) {dead[i]=1; break;}
        }
    }
    vector<Rect> pr; pr.reserve(freeR.size());
    for(size_t i=0;i<freeR.size();++i) if(!dead[i]) pr.push_back(freeR[i]);
    if(pr.size()>450){
        sort(pr.begin(), pr.end(), [](const Rect&a,const Rect&b){
            long long aa=1LL*a.w*a.h, bb=1LL*b.w*b.h;
            if(aa!=bb) return aa>bb; if(a.y!=b.y) return a.y<b.y; return a.x<b.x;
        });
        pr.resize(450);
    }
    freeR.swap(pr);
}

long double evalScore(int mode, const Item& it, int rw, int rh, const Rect& fr, int step){
    long double area=1.0L*rw*rh, dens=it.v/area;
    int remW=fr.w-rw, remH=fr.h-rh; long double waste=1.0L*fr.w*fr.h-area;
    switch(mode){
        case 0: return dens*1e12L + it.v - (min(remW,remH))*1e-3L; // density first
        case 1: return it.v*1.0L - waste*0.02L;                 // high absolute value
        case 2: return it.v/area*1e9L - (long double)min(remW,remH)*1000 - max(remW,remH); // tight fit
        case 3: return (long double)it.v*sqrt(area) - waste;     // balanced
        case 4: return dens*1e10L + (rw>rh?1:0)*100 + it.v*1e-3L;
        case 5: return dens*1e10L + (rh>rw?1:0)*100 + it.v*1e-3L;
        default: return dens*1e11L + it.v*0.1L - waste*0.001L;
    }
}

vector<Placement> runMaxRects(const vector<Item>& items, int W, int H, bool allowRot, int mode){
    int n=items.size(); vector<int> left(n); for(int i=0;i<n;i++) left[i]=items[i].limit;
    vector<Rect> freeR={{0,0,W,H}}; vector<Placement> sol; sol.reserve(5000);
    int guard=0;
    while(++guard<=12000){
        if((guard&63)==0 && timeUp()) break;
        long double best=-1e100L; int bi=-1, brot=0, bx=0, by=0, bw=0,bh=0;
        for(int f=0; f<(int)freeR.size(); ++f){ Rect fr=freeR[f];
            for(int i=0;i<n;i++) if(left[i]>0 && items[i].v>0){
                for(int rot=0; rot<=(allowRot?1:0); ++rot){
                    int rw=rot?items[i].h:items[i].w, rh=rot?items[i].w:items[i].h;
                    if(rw>fr.w || rh>fr.h) continue;
                    long double sc=evalScore(mode,items[i],rw,rh,fr,(int)sol.size());
                    // deterministic tie-breakers: lower y/x, then smaller waste.
                    sc -= (long double)(fr.y*0.000001L + fr.x*0.000000001L);
                    if(sc>best){best=sc; bi=i; brot=rot; bx=fr.x; by=fr.y; bw=rw; bh=rh;}
                }
            }
        }
        if(bi<0) break;
        sol.push_back({bi,bx,by,brot,bw,bh}); left[bi]--; splitFree(freeR,{bx,by,bw,bh});
    }
    return sol;
}

// Simple shelf baseline variants: sort item copies by a key, fill rows left-to-right.
vector<Placement> runShelves(vector<Item> items, int W, int H, bool allowRot, int mode){
    int n=items.size(); vector<int> left(n); for(int i=0;i<n;i++) left[i]=items[i].limit;
    vector<Placement> sol; int y=0;
    while(y<H){
        int bestH=0; vector<tuple<long double,int,int,int,int>> cand;
        for(int i=0;i<n;i++) if(left[i]>0) for(int rot=0; rot<=(allowRot?1:0); ++rot){
            int rw=rot?items[i].h:items[i].w, rh=rot?items[i].w:items[i].h;
            if(rw<=W && y+rh<=H){
                long double key;
                if(mode==0) key=items[i].v/(long double)(rw*rh);
                else if(mode==1) key=items[i].v/(long double)rh;
                else key=items[i].v/(long double)rw;
                cand.emplace_back(-key,i,rot,rw,rh);
            }
        }
        if(cand.empty()) break;
        sort(cand.begin(), cand.end());
        int seed=get<1>(cand[0]), srot=get<2>(cand[0]);
        bestH=srot?items[seed].w:items[seed].h;
        int x=0;
        bool placed=true;
        while(placed){ placed=false; int bi=-1,br=0,bw=0,bh=0; long double bk=-1e100;
            for(int i=0;i<n;i++) if(left[i]>0) for(int rot=0; rot<=(allowRot?1:0); ++rot){
                int rw=rot?items[i].h:items[i].w, rh=rot?items[i].w:items[i].h;
                if(rh<=bestH && x+rw<=W && y+rh<=H){
                    long double k=(mode==0?items[i].v/(long double)(rw*rh): mode==1?items[i].v/(long double)rw:items[i].v/(long double)rh);
                    if(k>bk){bk=k; bi=i;br=rot;bw=rw;bh=rh;}
                }
            }
            if(bi>=0){ sol.push_back({bi,x,y,br,bw,bh}); left[bi]--; x+=bw; placed=true; }
        }
        y += bestH;
    }
    return sol;
}

long long valueOf(const vector<Placement>& sol, const vector<Item>& items){ long long v=0; for(auto&p:sol) v+=items[p.id].v; return v; }
string escJson(const string& s){ string r; for(char c:s){ if(c=='"'||c=='\\'){ r+='\\'; r+=c; } else if(c=='\n') r+="\\n"; else if(c=='\r') r+="\\r"; else if(c=='\t') r+="\\t"; else r+=c; } return r; }

int main(){
    T0=chrono::steady_clock::now();
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    Parser ps(input); int W=0,H=0; bool allowRot=false; vector<Item> items;
    ps.ch('{');
    for(int top=0; top<2; ++top){ if(top) ps.ch(','); string k=ps.str(); ps.ch(':');
        if(k=="bin"){
            ps.ch('{'); for(int i=0;i<3;i++){ if(i) ps.ch(','); string bk=ps.str(); ps.ch(':'); if(bk=="W") W=ps.num(); else if(bk=="H") H=ps.num(); else if(bk=="allow_rotate") allowRot=ps.boolean(); } ps.ch('}');
        }else if(k=="items"){
            ps.ch('['); ps.ws(); bool first=true; while(ps.p<ps.s.size() && ps.s[ps.p]!=']'){
                if(!first) ps.ch(','); first=false; ps.ch('{'); Item it; for(int f=0; f<5; ++f){ if(f) ps.ch(','); string key=ps.str(); ps.ch(':'); if(key=="type") it.type=ps.str(); else if(key=="w") it.w=ps.num(); else if(key=="h") it.h=ps.num(); else if(key=="v") it.v=ps.num(); else if(key=="limit") it.limit=ps.num(); } ps.ch('}'); items.push_back(it); ps.ws(); }
            ps.ch(']');
        }
    }
    vector<Placement> best; long long bestV=-1;
    for(int m=0;m<6 && !timeUp();m++){ auto s=runMaxRects(items,W,H,allowRot,m); long long v=valueOf(s,items); if(v>bestV){bestV=v; best=s;} }
    for(int m=0;m<3 && !timeUp();m++){ auto s=runShelves(items,W,H,allowRot,m); long long v=valueOf(s,items); if(v>bestV){bestV=v; best=s;} }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){ auto&p=best[i]; if(i) cout << ','; cout << "{\"type\":\"" << escJson(items[p.id].type) << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}"; }
    cout << "]}\n";
    return 0;
}
