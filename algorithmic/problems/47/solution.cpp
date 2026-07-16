#include <bits/stdc++.h>
using namespace std;

struct JsonParser {
    string s; size_t p=0;
    JsonParser(const string& in):s(in){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) p++; }
    string str(){ ws(); ch('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\' && p<s.size()){ char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); int sg=1; if(p<s.size() && s[p]=='-'){sg=-1;p++;} long long v=0; while(p<s.size() && isdigit((unsigned char)s[p])) v=v*10+(s[p++]-'0'); return sg*v; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} p+=5; return false; }
};

struct Item { string type; long long w,h,v,lim; double dens; };
struct Pl { int idx; long long x,y,w,h; int rot; };
struct Rect { long long x,y,w,h; };

static bool intersect(const Rect& a, const Rect& b){
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}
static bool containsR(const Rect& a, const Rect& b){
    return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h;
}

struct Solver {
    long long W,H; bool allow; vector<Item> it;
    chrono::steady_clock::time_point start;
    double limitSec = 0.88;
    bool timeout() const { return chrono::duration<double>(chrono::steady_clock::now()-start).count() > limitSec; }

    vector<Pl> runMode(int mode){
        int n=it.size();
        vector<long long> rem(n); for(int i=0;i<n;i++) rem[i]=it[i].lim;
        vector<Rect> freeR{{0,0,W,H}};
        vector<Pl> out; out.reserve(4096);
        long long totalLimit=0; for(auto x: rem) totalLimit += min<long long>(x, 1000000);
        int steps=0;
        while(!freeR.empty() && steps < totalLimit && steps < 20000 && !timeout()){
            struct Cand{ bool ok=false; int fi=-1, idx=-1, rot=0; long long w=0,h=0; double a=0,b=0,c=0; } best;
            for(int fi=0; fi<(int)freeR.size() && !timeout(); ++fi){
                Rect fr=freeR[fi];
                for(int i=0;i<n;i++) if(rem[i]>0 && it[i].v>0){
                    for(int rr=0; rr<(allow?2:1); ++rr){
                        long long ww = rr?it[i].h:it[i].w, hh = rr?it[i].w:it[i].h;
                        if(ww>fr.w || hh>fr.h) continue;
                        long long waste = fr.w*fr.h - ww*hh;
                        long long dw=fr.w-ww, dh=fr.h-hh;
                        double den = it[i].dens, val = (double)it[i].v, area=(double)ww*hh;
                        double A=0,B=0,C=0;
                        if(mode==0){ A=den; B=-min(dw,dh); C=-waste; }
                        else if(mode==1){ A=val; B=den; C=-waste; }
                        else if(mode==2){ A=den; B=-waste; C=area; }
                        else if(mode==3){ A=val/sqrt(max(1.0,area)); B=den; C=-max(dw,dh); }
                        else if(mode==4){ A=den; B=-(double)(dw*dw+dh*dh); C=val; }
                        else if(mode==5){ A=den*(1.0 + 0.07*sin((i+1)*12.9898 + (fi+1)*78.233)); B=val; C=-waste; }
                        else if(mode==6){ A=(double)it[i].v / max(1LL, min(ww,hh)); B=den; C=-waste; }
                        else { A=area; B=den; C=val; }
                        // Prefer lower-left coordinates for reproducible compact packings.
                        double D = -(double)fr.y*1e-9 - (double)fr.x*1e-12;
                        bool better = !best.ok || A>best.a+1e-12 || (fabs(A-best.a)<=1e-12 && (B>best.b+1e-9 || (fabs(B-best.b)<=1e-9 && (C>best.c+1e-6 || (fabs(C-best.c)<=1e-6 && D>0)))));
                        if(better){ best={true,fi,i,rr,ww,hh,A,B,C}; }
                    }
                }
            }
            if(!best.ok) break;
            Rect used{freeR[best.fi].x, freeR[best.fi].y, best.w, best.h};
            out.push_back({best.idx, used.x, used.y, used.w, used.h, best.rot});
            rem[best.idx]--; steps++;

            vector<Rect> nf; nf.reserve(freeR.size()+4);
            for(const Rect& r: freeR){
                if(!intersect(r, used)){ nf.push_back(r); continue; }
                long long rx2=r.x+r.w, ry2=r.y+r.h, ux2=used.x+used.w, uy2=used.y+used.h;
                if(used.x > r.x) nf.push_back({r.x,r.y,used.x-r.x,r.h});
                if(ux2 < rx2) nf.push_back({ux2,r.y,rx2-ux2,r.h});
                if(used.y > r.y) nf.push_back({r.x,r.y,r.w,used.y-r.y});
                if(uy2 < ry2) nf.push_back({r.x,uy2,r.w,ry2-uy2});
            }
            vector<Rect> clean;
            clean.reserve(nf.size());
            for(auto r: nf) if(r.w>0 && r.h>0){
                bool dup=false; for(auto &q: clean) if(q.x==r.x&&q.y==r.y&&q.w==r.w&&q.h==r.h){dup=true;break;}
                if(!dup) clean.push_back(r);
            }
            vector<char> bad(clean.size(),0);
            for(int i=0;i<(int)clean.size();++i) if(!bad[i])
                for(int j=0;j<(int)clean.size();++j) if(i!=j && !bad[j])
                    if(containsR(clean[i], clean[j])) bad[j]=1;
            freeR.clear();
            for(int i=0;i<(int)clean.size();++i) if(!bad[i]) freeR.push_back(clean[i]);
            if((int)freeR.size()>2500){
                sort(freeR.begin(), freeR.end(), [](const Rect&a,const Rect&b){
                    long long aa=a.w*a.h, bb=b.w*b.h; if(aa!=bb) return aa>bb; if(a.y!=b.y) return a.y<b.y; return a.x<b.x;
                });
                freeR.resize(2500);
            }
        }
        return out;
    }
};

static string esc(const string& s){ string r; for(char c:s){ if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    JsonParser jp(input);
    long long W=0,H=0; bool allow=false; vector<Item> items;
    jp.ch('{');
    bool firstTop=true;
    while(true){ jp.ws(); if(jp.p>=jp.s.size() || jp.s[jp.p]=='}'){ if(jp.p<jp.s.size()) jp.p++; break; } if(!firstTop) jp.ch(','); firstTop=false; string key=jp.str(); jp.ch(':');
        if(key=="bin"){
            jp.ch('{'); bool first=true; while(true){ jp.ws(); if(jp.s[jp.p]=='}'){jp.p++; break;} if(!first) jp.ch(','); first=false; string k=jp.str(); jp.ch(':'); if(k=="W") W=jp.num(); else if(k=="H") H=jp.num(); else if(k=="allow_rotate") allow=jp.boolean(); }
        } else if(key=="items"){
            jp.ch('['); bool firstArr=true; while(true){ jp.ws(); if(jp.s[jp.p]==']'){jp.p++; break;} if(!firstArr) jp.ch(','); firstArr=false; Item it; jp.ch('{'); bool first=true; while(true){ jp.ws(); if(jp.s[jp.p]=='}'){jp.p++; break;} if(!first) jp.ch(','); first=false; string k=jp.str(); jp.ch(':'); if(k=="type") it.type=jp.str(); else if(k=="w") it.w=jp.num(); else if(k=="h") it.h=jp.num(); else if(k=="v") it.v=jp.num(); else if(k=="limit") it.lim=jp.num(); }
                it.dens = (double)it.v / max(1.0, (double)it.w*it.h); items.push_back(it); }
        }
    }
    Solver sol; sol.W=W; sol.H=H; sol.allow=allow; sol.it=items; sol.start=chrono::steady_clock::now();
    vector<Pl> best; long long bestVal=-1;
    for(int mode=0; mode<8 && !sol.timeout(); ++mode){
        vector<Pl> cur=sol.runMode(mode); long long val=0; for(auto&p:cur) val += items[p.idx].v;
        if(val>bestVal){ bestVal=val; best.swap(cur); }
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){
        if(i) cout << ',';
        auto &p=best[i];
        cout << "{\"type\":\"" << esc(items[p.idx].type) << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}";
    }
    cout << "]}\n";
    return 0;
}
