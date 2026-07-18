#include <bits/stdc++.h>
using namespace std;

struct Json {
    string s; size_t p=0;
    Json(string t):s(move(t)){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) p++; }
    string str(){ ws(); ch('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\' && p<s.size()){ char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); int sg=1; if(p<s.size() && s[p]=='-') sg=-1,p++; long long v=0; while(p<s.size() && isdigit((unsigned char)s[p])) v=v*10+(s[p++]-'0'); return sg*v; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} p+=5; return false; }
};

struct Item{ string id; int w,h,lim; long long v; };
struct Rect{ int x,y,w,h; };
struct Place{ int t,x,y,rot,w,h; };

static inline bool inter(const Rect&a,const Rect&b){ return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h; }
static inline bool contains(const Rect&a,const Rect&b){ return b.x>=a.x && b.y>=a.y && b.x+b.w<=a.x+a.w && b.y+b.h<=a.y+a.h; }

struct Solver{
    int W,H; bool allow; vector<Item> it; int n;
    struct Result{ vector<Place> p; long long val=0; } best;

    void splitFree(vector<Rect>& freeR, const Rect& u){
        vector<Rect> nf; nf.reserve(freeR.size()+4);
        for(auto f: freeR){
            if(!inter(f,u)){ nf.push_back(f); continue; }
            if(u.x > f.x) nf.push_back({f.x,f.y,u.x-f.x,f.h});
            if(u.x+u.w < f.x+f.w) nf.push_back({u.x+u.w,f.y,f.x+f.w-(u.x+u.w),f.h});
            if(u.y > f.y) nf.push_back({f.x,f.y,f.w,u.y-f.y});
            if(u.y+u.h < f.y+f.h) nf.push_back({f.x,u.y+u.h,f.w,f.y+f.h-(u.y+u.h)});
        }
        vector<Rect> good;
        good.reserve(nf.size());
        for(auto r:nf) if(r.w>0 && r.h>0) good.push_back(r);
        vector<char> dead(good.size(),0);
        for(size_t i=0;i<good.size();++i) if(!dead[i]){
            for(size_t j=0;j<good.size();++j) if(i!=j && !dead[j]){
                if(contains(good[i],good[j])) dead[j]=1;
            }
        }
        freeR.clear();
        for(size_t i=0;i<good.size();++i) if(!dead[i]) freeR.push_back(good[i]);
    }

    Result runHeuristic(int mode, unsigned seed){
        vector<int> rem(n); for(int i=0;i<n;i++) rem[i]=it[i].lim;
        vector<Rect> freeR{{0,0,W,H}};
        Result res;
        mt19937 rng(seed);
        int totalLim=0; for(auto &a:it) totalLim+=a.lim;
        for(int step=0; step<totalLim; ++step){
            int bt=-1, br=0, bfi=-1, bw=0,bh=0; long double bs=-1e100; int bx=0,by=0;
            for(int fi=0; fi<(int)freeR.size(); ++fi){
                Rect f=freeR[fi];
                for(int t=0;t<n;t++) if(rem[t]>0){
                    for(int rr=0; rr<(allow && it[t].w!=it[t].h ? 2:1); ++rr){
                        int w = rr?it[t].h:it[t].w, h = rr?it[t].w:it[t].h;
                        if(w>f.w || h>f.h) continue;
                        long double area=(long double)w*h, dens=(long double)it[t].v/area;
                        long double waste=(long double)(f.w*f.h - w*h);
                        int shortside=min(f.w-w, f.h-h), longside=max(f.w-w, f.h-h);
                        long double sc;
                        if(mode==0) sc = dens*1000000.0L + it[t].v*0.001L - shortside*0.1L - f.y*0.01L;
                        else if(mode==1) sc = it[t].v*1.0L + dens*50000.0L - waste*0.02L;
                        else if(mode==2) sc = dens*700000.0L - shortside*20.0L - longside*0.05L + it[t].v*0.0005L;
                        else if(mode==3) sc = (long double)it[t].v * (1.0L + 0.15L*(long double)area/(W*(long double)H)) - waste*0.01L;
                        else if(mode==4) sc = dens*400000.0L + it[t].v*0.002L - (f.x+f.y)*0.05L - abs((f.w-w)-(f.h-h))*0.2L;
                        else sc = dens*600000.0L + it[t].v*0.001L - shortside*5.0L + (long double)(rng()%1000)*0.001L;
                        // Prefer bottom-left and then smaller free rectangle for deterministic tie-breaking.
                        sc -= (long double)f.y*0.0001L + (long double)f.x*0.00005L + (long double)f.w*f.h*1e-9L;
                        if(sc>bs){ bs=sc; bt=t; br=rr; bfi=fi; bw=w; bh=h; bx=f.x; by=f.y; }
                    }
                }
            }
            if(bt<0) break;
            res.p.push_back({bt,bx,by,br,bw,bh}); res.val += it[bt].v; rem[bt]--;
            splitFree(freeR, {bx,by,bw,bh});
            if(freeR.size()>5000){
                sort(freeR.begin(), freeR.end(), [](const Rect&a,const Rect&b){ return (long long)a.w*a.h > (long long)b.w*b.h; });
                freeR.resize(5000);
            }
        }
        return res;
    }

    Result shelf(int mode){
        struct Ori{int t,w,h,r; long double key;}; vector<Ori> o;
        for(int t=0;t<n;t++) for(int r=0;r<(allow&&it[t].w!=it[t].h?2:1);r++){
            int w=r?it[t].h:it[t].w, h=r?it[t].w:it[t].h; if(w<=W&&h<=H){
                long double d=(long double)it[t].v/(w*(long double)h);
                long double k = mode==0? -d : mode==1? -(long double)it[t].v : mode==2? -((long double)it[t].v/h) : -((long double)it[t].v/w);
                o.push_back({t,w,h,r,k});
            }
        }
        sort(o.begin(),o.end(),[](const Ori&a,const Ori&b){ if(a.key!=b.key) return a.key<b.key; return a.h>b.h; });
        vector<int> rem(n); for(int i=0;i<n;i++) rem[i]=it[i].lim;
        Result res; int y=0;
        while(y<H){
            bool any=false; int shelfH=0, x=0;
            for(auto &a:o) if(rem[a.t]>0 && a.h<=H-y){ shelfH=a.h; any=true; break; }
            if(!any) break;
            bool placed=true;
            while(placed){
                placed=false; int bi=-1;
                for(int i=0;i<(int)o.size();++i){ auto&a=o[i]; if(rem[a.t]>0 && a.h<=shelfH && x+a.w<=W){ bi=i; break; } }
                if(bi>=0){ auto&a=o[bi]; res.p.push_back({a.t,x,y,a.r,a.w,a.h}); res.val+=it[a.t].v; rem[a.t]--; x+=a.w; placed=true; }
            }
            y += shelfH;
        }
        return res;
    }

    Result adaptiveShelf(int mode){
        struct Ori{int t,w,h,r;}; vector<Ori> o; vector<int> heights;
        for(int t=0;t<n;t++) for(int r=0;r<(allow&&it[t].w!=it[t].h?2:1);r++){
            int w=r?it[t].h:it[t].w, h=r?it[t].w:it[t].h;
            if(w<=W && h<=H){ o.push_back({t,w,h,r}); heights.push_back(h); }
        }
        sort(heights.begin(), heights.end()); heights.erase(unique(heights.begin(), heights.end()), heights.end());
        vector<int> rem(n); for(int i=0;i<n;i++) rem[i]=it[i].lim;
        Result res; int y=0;
        auto buildRow = [&](int sh, const vector<int>& baseRem){
            vector<int> tr=baseRem; Result row; int x=0;
            while(true){
                int bi=-1; long double bs=-1e100;
                for(int i=0;i<(int)o.size();++i){ auto &a=o[i]; if(tr[a.t]<=0 || a.h>sh || x+a.w>W) continue;
                    long double dens=(long double)it[a.t].v/(a.w*(long double)a.h);
                    long double sc;
                    if(mode==0) sc=(long double)it[a.t].v/a.w;
                    else if(mode==1) sc=dens*10000.0L + (long double)it[a.t].v/a.w;
                    else sc=(long double)it[a.t].v - (long double)(sh-a.h)*a.w*dens*0.35L;
                    sc -= (long double)(W-x-a.w)*1e-7L;
                    if(sc>bs){ bs=sc; bi=i; }
                }
                if(bi<0) break;
                auto &a=o[bi]; row.p.push_back({a.t,x,0,a.r,a.w,a.h}); row.val+=it[a.t].v; tr[a.t]--; x+=a.w;
            }
            return row;
        };
        while(y<H){
            int bestH=-1; Result bestRow; long double bestScore=-1e100;
            for(int sh: heights) if(sh<=H-y){
                bool possible=false; for(auto &a:o) if(a.h<=sh && rem[a.t]>0){ possible=true; break; }
                if(!possible) continue;
                Result row=buildRow(sh, rem); if(row.p.empty()) continue;
                long double score;
                if(mode==0) score=(long double)row.val/sh;
                else if(mode==1) score=(long double)row.val/(sh+0.20L*(H-y-sh));
                else score=(long double)row.val - (long double)sh*W*0.0005L;
                if(score>bestScore){ bestScore=score; bestH=sh; bestRow=move(row); }
            }
            if(bestH<0) break;
            for(auto p: bestRow.p){ p.y=y; res.val+=it[p.t].v; rem[p.t]--; res.p.push_back(p); }
            y += bestH;
        }
        return res;
    }

    void solve(){
        best = Result();
        for(int m=0;m<7;m++){ auto r=runHeuristic(m, 12345u+1009u*m); if(r.val>best.val) best=move(r); }
        for(int m=0;m<4;m++){ auto r=shelf(m); if(r.val>best.val) best=move(r); }
        for(int m=0;m<3;m++){ auto r=adaptiveShelf(m); if(r.val>best.val) best=move(r); }
    }
};

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    Json js(input); Solver sol; vector<Item> items;
    js.ch('{');
    for(int top=0; top<2; ++top){
        if(top) js.ch(','); string k=js.str(); js.ch(':');
        if(k=="bin"){
            js.ch('{'); for(int i=0;i<3;i++){ if(i) js.ch(','); string bk=js.str(); js.ch(':'); if(bk=="W") sol.W=js.num(); else if(bk=="H") sol.H=js.num(); else sol.allow=js.boolean(); } js.ch('}');
        }else if(k=="items"){
            js.ch('['); js.ws(); bool first=true; while(js.p<input.size() && input[js.p]!=']'){
                if(!first) js.ch(','); first=false; js.ch('{'); Item a; for(int f=0; f<5; ++f){ if(f) js.ch(','); string ik=js.str(); js.ch(':'); if(ik=="type") a.id=js.str(); else if(ik=="w") a.w=js.num(); else if(ik=="h") a.h=js.num(); else if(ik=="v") a.v=js.num(); else if(ik=="limit") a.lim=js.num(); } js.ch('}'); items.push_back(a); js.ws(); }
            js.ch(']');
        }
    }
    sol.it=items; sol.n=items.size(); sol.solve();
    cout << "{\"placements\":[";
    for(size_t i=0;i<sol.best.p.size();++i){ auto &p=sol.best.p[i]; if(i) cout << ','; cout << "{\"type\":\"" << sol.it[p.t].id << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}"; }
    cout << "]}\n";
    return 0;
}
