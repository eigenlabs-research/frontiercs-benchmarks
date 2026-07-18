#include <bits/stdc++.h>
using namespace std;

struct JsonParser {
    string s; size_t p=0;
    JsonParser(const string& in):s(in){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) p++; }
    string str(){ ws(); ch('"'); string r; while(p<s.size()){ char c=s[p++]; if(c=='"') break; if(c=='\\' && p<s.size()){ char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else r+=e; } else r+=c; } return r; }
    long long num(){ ws(); int sign=1; if(p<s.size()&&(s[p]=='-'||s[p]=='+')){ if(s[p]=='-') sign=-1; p++; } long long v=0; while(p<s.size()&&isdigit((unsigned char)s[p])) v=v*10+(s[p++]-'0'); return sign*v; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4; return true;} if(s.compare(p,5,"false")==0){p+=5; return false;} return false; }
};

struct Item { string type; int w,h; long long v; int lim; };
struct Placement { int id,x,y,rot,w,h; };
struct Rect { int x,y,w,h; };

static bool inter(const Rect&a,const Rect&b){ return a.x < b.x+b.w && a.x+a.w > b.x && a.y < b.y+b.h && a.y+a.h > b.y; }
static bool contains(const Rect&a,const Rect&b){ return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h; }

struct Solver {
    int W,H,n; bool allowRot; vector<Item> items;
    vector<Placement> best;
    long long bestVal=0;

    void prune(vector<Rect>& fr){
        vector<Rect> v;
        for(auto &r: fr) if(r.w>0 && r.h>0) v.push_back(r);
        vector<char> rem(v.size(),0);
        for(size_t i=0;i<v.size();++i) if(!rem[i]){
            for(size_t j=0;j<v.size();++j) if(i!=j && !rem[j]){
                if(contains(v[i],v[j])) rem[j]=1;
                else if(contains(v[j],v[i])) { rem[i]=1; break; }
            }
        }
        fr.clear();
        for(size_t i=0;i<v.size();++i) if(!rem[i]) fr.push_back(v[i]);
        sort(fr.begin(), fr.end(), [](const Rect&a,const Rect&b){
            if(a.y!=b.y) return a.y<b.y; if(a.x!=b.x) return a.x<b.x; return (long long)a.w*a.h > (long long)b.w*b.h;
        });
        if(fr.size()>900) fr.resize(900);
    }

    void placeSplit(vector<Rect>& fr, const Rect& used){
        vector<Rect> nf;
        nf.reserve(fr.size()+8);
        for(auto r: fr){
            if(!inter(r,used)){ nf.push_back(r); continue; }
            if(used.x > r.x) nf.push_back({r.x,r.y,used.x-r.x,r.h});
            if(used.x+used.w < r.x+r.w) nf.push_back({used.x+used.w,r.y,r.x+r.w-(used.x+used.w),r.h});
            if(used.y > r.y) nf.push_back({r.x,r.y,r.w,used.y-r.y});
            if(used.y+used.h < r.y+r.h) nf.push_back({r.x,used.y+used.h,r.w,r.y+r.h-(used.y+used.h)});
        }
        fr.swap(nf); prune(fr);
    }

    double key(int mode, const Item& it, int w, int h, const Rect& r, int step){
        double area = (double)w*h, dens = it.v / area;
        double waste = (double)r.w*r.h - area;
        double ss = min(r.w-w, r.h-h), ls = max(r.w-w, r.h-h);
        if(mode==0) return dens*1e9 - waste - ss*1000;
        if(mode==1) return (double)it.v*1000.0 - waste*0.01;
        if(mode==2) return (double)it.v*sqrt(dens) - waste*0.1;
        if(mode==3) return -ss*1e9 - ls*1000 + dens;
        if(mode==4) return -waste*1000 + (double)it.v;
        if(mode==5) return dens*1e8 + (double)max(w,h)*1000 - waste;
        return (double)it.v/ (1.0 + waste) + dens*1000;
    }

    void runMode(int mode){
        vector<int> rem(n); for(int i=0;i<n;i++) rem[i]=items[i].lim;
        vector<Rect> fr; fr.push_back({0,0,W,H});
        vector<Placement> out; long long val=0;
        int totalLim=0; for(auto &it:items) totalLim += max(0,it.lim);
        for(int step=0; step<totalLim; ++step){
            double bk=-1e300; int bi=-1, brot=0, bw=0,bh=0; Rect br{0,0,0,0};
            for(int i=0;i<n;i++) if(rem[i]>0 && items[i].v>=0){
                for(int ro=0; ro<(allowRot?2:1); ++ro){
                    int w = ro?items[i].h:items[i].w, h = ro?items[i].w:items[i].h;
                    if(w>W || h>H) continue;
                    for(const auto& r: fr) if(w<=r.w && h<=r.h){
                        double k = key(mode, items[i], w,h,r,step);
                        // Prefer lower/left coordinates and exact edge matches for stability.
                        k += (r.w==w)*5000.0 + (r.h==h)*5000.0 - r.y*1e-6 - r.x*1e-7;
                        if(k>bk){ bk=k; bi=i; brot=ro; bw=w; bh=h; br=r; }
                    }
                }
            }
            if(bi<0) break;
            out.push_back({bi,br.x,br.y,brot,bw,bh}); val += items[bi].v; rem[bi]--;
            placeSplit(fr, {br.x,br.y,bw,bh});
        }
        if(val>bestVal){ bestVal=val; best.swap(out); }
    }

    void runShelves(){
        // Independent deterministic shelf fill: try each type/orientation as shelf height, fill rows greedily by density.
        for(int seed=0; seed<n; ++seed) for(int srot=0; srot<(allowRot?2:1); ++srot){
            int sh = srot?items[seed].w:items[seed].h; if(sh<=0 || sh>H) continue;
            vector<int> rem(n); for(int i=0;i<n;i++) rem[i]=items[i].lim;
            vector<Placement> out; long long val=0; int y=0;
            while(y+sh<=H){
                int x=0, rowh=sh;
                while(true){
                    int bi=-1, bro=0, bw=0,bh=0; double bk=-1e300;
                    for(int i=0;i<n;i++) if(rem[i]>0){
                        for(int ro=0; ro<(allowRot?2:1); ++ro){
                            int w=ro?items[i].h:items[i].w, h=ro?items[i].w:items[i].h;
                            if(h<=rowh && x+w<=W){ double k=items[i].v/(double)(w*h)*1e9 + items[i].v*0.001 + (h==rowh)*1000; if(k>bk){bk=k;bi=i;bro=ro;bw=w;bh=h;} }
                        }
                    }
                    if(bi<0) break;
                    out.push_back({bi,x,y,bro,bw,bh}); val+=items[bi].v; rem[bi]--; x+=bw;
                }
                y += rowh;
            }
            if(val>bestVal){ bestVal=val; best.swap(out); }
        }
    }
};

static string escJson(const string& s){
    string r;
    for(char c: s){
        if(c=='"' || c=='\\') { r+='\\'; r+=c; }
        else if(c=='\n') r += "\\n";
        else if(c=='\r') r += "\\r";
        else if(c=='\t') r += "\\t";
        else r += c;
    }
    return r;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    JsonParser jp(input);
    int W=0,H=0; bool allow=false; vector<Item> items;
    jp.ch('{');
    for(int top=0; top<2; ++top){
        if(top) jp.ch(','); string k=jp.str(); jp.ch(':');
        if(k=="bin"){
            jp.ch('{'); for(int i=0;i<3;i++){ if(i) jp.ch(','); string bk=jp.str(); jp.ch(':'); if(bk=="W") W=jp.num(); else if(bk=="H") H=jp.num(); else if(bk=="allow_rotate") allow=jp.boolean(); } jp.ch('}');
        } else if(k=="items"){
            jp.ch('['); jp.ws(); bool first=true;
            while(jp.p<jp.s.size() && jp.s[jp.p]!=']'){
                if(!first) jp.ch(','); first=false; jp.ch('{'); Item it; for(int f=0; f<5; ++f){ if(f) jp.ch(','); string kk=jp.str(); jp.ch(':'); if(kk=="type") it.type=jp.str(); else if(kk=="w") it.w=jp.num(); else if(kk=="h") it.h=jp.num(); else if(kk=="v") it.v=jp.num(); else if(kk=="limit") it.lim=jp.num(); } jp.ch('}'); items.push_back(it); jp.ws();
            } jp.ch(']');
        }
    }
    Solver sol; sol.W=W; sol.H=H; sol.allowRot=allow; sol.items=items; sol.n=items.size();
    for(int m=0;m<7;m++) sol.runMode(m);
    sol.runShelves();
    cout << "{\"placements\":[";
    for(size_t i=0;i<sol.best.size();++i){
        if(i) cout << ',';
        auto &p=sol.best[i];
        cout << "{\"type\":\"" << escJson(sol.items[p.id].type) << "\",\"x\":" << p.x << ",\"y\":" << p.y << ",\"rot\":" << p.rot << "}";
    }
    cout << "]}\n";
    return 0;
}
