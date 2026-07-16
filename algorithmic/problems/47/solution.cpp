#include <bits/stdc++.h>
using namespace std;

struct JsonParser {
    string s; size_t p=0;
    JsonParser(const string& t):s(t){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) p++; }
    string str(){ ws(); ch('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\' && p<s.size()){ char e=s[p++]; if(e=='"'||e=='\\'||e=='/') r+=e; else if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); int sign=1; if(p<s.size()&&(s[p]=='-'||s[p]=='+')){ if(s[p]=='-') sign=-1; p++; } long long v=0; while(p<s.size() && isdigit((unsigned char)s[p])) v=v*10+(s[p++]-'0'); return sign*v; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} if(s.compare(p,5,"false")==0){p+=5; return false;} return false; }
};

struct Item { string id; int w,h; long long v; int lim; };
struct Rect { int x,y,w,h; };
struct Place { int idx,x,y,rot,w,h; };

static inline bool inter(const Rect&a,const Rect&b){ return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h; }
static inline bool contains(const Rect&a,const Rect&b){ return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h; }

struct Solver {
    int W,H; bool rotAllowed; vector<Item> it;
    chrono::steady_clock::time_point deadline;

    vector<Place> run(int mode){
        vector<int> rem; for(auto &x:it) rem.push_back(x.lim);
        vector<Rect> freeR{{0,0,W,H}};
        vector<Place> out;
        long long totalLim=0; for(int x:rem) totalLim+=x;
        while(totalLim>0){
            if(chrono::steady_clock::now() > deadline) break;
            bool found=false; int bi=-1,bf=-1,brot=0,bw=0,bh=0; double best=-1e300; long long bv=0;
            for(int i=0;i<(int)it.size();++i) if(rem[i]>0 && it[i].v>=0){
                for(int ro=0; ro<=(rotAllowed?1:0); ++ro){
                    int w = ro?it[i].h:it[i].w, h = ro?it[i].w:it[i].h;
                    if(w>W || h>H) continue;
                    for(int f=0; f<(int)freeR.size(); ++f){
                        const Rect &r=freeR[f]; if(w>r.w || h>r.h) continue;
                        int dw=r.w-w, dh=r.h-h;
                        double area=(double)w*h, dens=it[i].v/area;
                        double waste=(double)r.w*r.h-area;
                        double shortFit=min(dw,dh), longFit=max(dw,dh);
                        double score;
                        if(mode==0) score=dens*1000000000.0 - shortFit*2000.0 - longFit - waste*0.0001;
                        else if(mode==1) score=(double)it[i].v - waste*0.02 - shortFit*100.0 - longFit*0.1;
                        else if(mode==2) score=dens*100000000.0 + (double)it[i].v*0.05 - shortFit*5000.0 - longFit*5.0;
                        else if(mode==3) score=(double)it[i].v/(sqrt(area)) - shortFit*50.0 - longFit*0.5;
                        else if(mode==4) score=(double)it[i].v + dens*10000000.0 - (dw==0||dh==0? -5000000.0:0) - waste*0.005;
                        else score=dens*1000000000.0 + (double)min(rem[i],5)*1000.0 - waste*0.001;
                        // deterministic tie breakers: prefer higher value, lower y then x
                        if(!found || score>best+1e-9 || (fabs(score-best)<=1e-9 && (it[i].v>bv || (it[i].v==bv && (freeR[f].y<freeR[bf].y || (freeR[f].y==freeR[bf].y && freeR[f].x<freeR[bf].x)))))){
                            found=true; best=score; bi=i; bf=f; brot=ro; bw=w; bh=h; bv=it[i].v;
                        }
                    }
                }
            }
            if(!found) break;
            Rect used{freeR[bf].x, freeR[bf].y, bw, bh};
            out.push_back({bi,used.x,used.y,brot,bw,bh}); rem[bi]--; totalLim--;
            vector<Rect> nf;
            nf.reserve(freeR.size()+4);
            for(const Rect &r: freeR){
                if(!inter(r,used)){ nf.push_back(r); continue; }
                int rx2=r.x+r.w, ry2=r.y+r.h, ux2=used.x+used.w, uy2=used.y+used.h;
                if(used.x>r.x) nf.push_back({r.x,r.y,used.x-r.x,r.h});
                if(ux2<rx2) nf.push_back({ux2,r.y,rx2-ux2,r.h});
                if(used.y>r.y) nf.push_back({r.x,r.y,r.w,used.y-r.y});
                if(uy2<ry2) nf.push_back({r.x,uy2,r.w,ry2-uy2});
            }
            // discard degenerate and contained free rectangles
            vector<Rect> clean;
            for(auto &r:nf){
                if(r.w<=0 || r.h<=0) continue;
                bool dup=false;
                for(auto &q:clean) if(r.x==q.x&&r.y==q.y&&r.w==q.w&&r.h==q.h){ dup=true; break; }
                if(!dup) clean.push_back(r);
            }
            if(clean.size()>900){
                sort(clean.begin(), clean.end(), [](const Rect&a,const Rect&b){ return (long long)a.w*a.h > (long long)b.w*b.h; });
                clean.resize(900);
            }
            vector<char> bad(clean.size(),0);
            for(int a=0;a<(int)clean.size();++a) if(!bad[a])
                for(int b=0;b<(int)clean.size();++b) if(a!=b && !bad[b] && contains(clean[a],clean[b])) bad[b]=1;
            freeR.clear();
            for(int a=0;a<(int)clean.size();++a) if(!bad[a]) freeR.push_back(clean[a]);
            sort(freeR.begin(), freeR.end(), [](const Rect&a,const Rect&b){
                if(a.y!=b.y) return a.y<b.y; if(a.x!=b.x) return a.x<b.x; return (long long)a.w*a.h>(long long)b.w*b.h;
            });
            if(freeR.size()>1200) freeR.resize(1200);
        }
        return out;
    }
};

long long valOf(const vector<Place>&p,const vector<Item>&it){ long long s=0; for(auto &x:p) s+=it[x.idx].v; return s; }

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    JsonParser jp(input);
    int W=0,H=0; bool allow=false; vector<Item> items;
    jp.ch('{');
    for(int top=0; top<2; ++top){ if(top) jp.ch(','); string key=jp.str(); jp.ch(':');
        if(key=="bin"){
            jp.ch('{'); for(int i=0;i<3;i++){ if(i) jp.ch(','); string k=jp.str(); jp.ch(':'); if(k=="W") W=jp.num(); else if(k=="H") H=jp.num(); else if(k=="allow_rotate") allow=jp.boolean(); } jp.ch('}');
        } else if(key=="items"){
            jp.ch('['); jp.ws(); bool first=true; while(jp.p<input.size() && input[jp.p]!=']'){
                if(!first) jp.ch(','); first=false; jp.ch('{'); Item it; for(int f=0;f<5;f++){ if(f) jp.ch(','); string k=jp.str(); jp.ch(':'); if(k=="type") it.id=jp.str(); else if(k=="w") it.w=jp.num(); else if(k=="h") it.h=jp.num(); else if(k=="v") it.v=jp.num(); else if(k=="limit") it.lim=jp.num(); } jp.ch('}'); items.push_back(it); jp.ws(); }
            jp.ch(']');
        }
    }
    Solver sol{W,H,allow,items, chrono::steady_clock::now() + chrono::milliseconds(850)};
    vector<Place> best; long long bestv=-1;
    for(int m=0;m<6;m++){
        if(chrono::steady_clock::now() > sol.deadline) break;
        vector<Place> p=sol.run(m); long long v=valOf(p,items);
        if(v>bestv){ bestv=v; best.swap(p); }
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){
        if(i) cout << ',';
        const auto &p=best[i];
        cout << "{\"type\":\"";
        for(char c: items[p.idx].id){ if(c=='"'||c=='\\') cout << '\\'; cout << c; }
        cout << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}";
    }
    cout << "]}\n";
    return 0;
}
