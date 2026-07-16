#include <bits/stdc++.h>
using namespace std;

struct Item { string name; long long w,h,v,lim; };
struct Rect { long long x,y,w,h; };
struct Pl { int t,rot; long long x,y,w,h; };

// The checker accepts the usual one-character escapes and treats a backslash-u escape literally.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); string r; if(p>=s.size()||s[p++]!='"') return r;
        while(p<s.size() && s[p]!='"') { char c=s[p++];
            if(c=='\\' && p<s.size()) { char e=s[p++];
                if(e=='"'||e=='\\'||e=='/') r+=e;
                else if(e=='b') r+='\b'; else if(e=='f') r+='\f'; else if(e=='n') r+='\n';
                else if(e=='r') r+='\r'; else if(e=='t') r+='\t'; else r+=e;
            } else r+=c;
        }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); bool neg=false; if(p<s.size()&&(s[p]=='-'||s[p]=='+')) neg=s[p++]=='-'; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+(s[p++]-'0'); return neg?-x:x; }
    bool boolean(){ ws(); bool r=p<s.size()&&s[p]=='t'; p+=r?4:5; return r; }
};
static bool inside(const Rect&a,const Rect&b){ return a.x>=b.x&&a.y>=b.y&&a.x+a.w<=b.x+b.w&&a.y+a.h<=b.y+b.h; }
static bool hit(const Rect&a,const Rect&b){ return a.x<b.x+b.w&&b.x<a.x+a.w&&a.y<b.y+b.h&&b.y<a.y+a.h; }

static void jsonString(const string& s){
    cout << '"';
    for(unsigned char c:s) {
        if(c=='"') cout << "\\\""; else if(c=='\\') cout << "\\\\";
        else if(c=='\b') cout << "\\b"; else if(c=='\f') cout << "\\f";
        else if(c=='\n') cout << "\\n"; else if(c=='\r') cout << "\\r";
        else if(c=='\t') cout << "\\t"; else cout << c;
    }
    cout << '"';
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Json j; long long W=0,H=0; bool allow=false; vector<Item> a;
    j.ch('{');
    for(int top=0;top<2;top++){
        if(top) j.ch(','); string key=j.str(); j.ch(':');
        if(key=="bin"){
            j.ch('{'); for(int z=0;z<3;z++){ if(z)j.ch(','); string k=j.str(); j.ch(':');
                if(k=="W") W=j.num(); else if(k=="H") H=j.num(); else allow=j.boolean(); }
            j.ch('}');
        } else {
            j.ch('['); bool first=true;
            while(true){ j.ws(); if(j.p>=j.s.size()||j.s[j.p]==']'){j.ch(']');break;} if(!first)j.ch(','); first=false;
                j.ch('{'); Item it; for(int z=0;z<5;z++){if(z)j.ch(','); string k=j.str();j.ch(':');
                    if(k=="type")it.name=j.str(); else if(k=="w")it.w=j.num(); else if(k=="h")it.h=j.num(); else if(k=="v")it.v=j.num(); else it.lim=j.num(); }
                j.ch('}'); a.push_back(it);
            }
        }
    }
    int n=a.size(); long double maxD=1,maxV=1, binA=(long double)W*H;
    for(auto &q:a){ maxD=max(maxD,(long double)q.v/(q.w*q.h)); maxV=max(maxV,(long double)q.v); }
    vector<Pl> answer; long long best=-1;
    // Profiles are an ablation of the selection rule, not different geometry engines.
    const long double prof[5][5]={{1,.10,.02,.18,.02},{1,.35,.00,.35,.01},{.75,.75,.02,.10,.00},{1.2,.02,.08,.42,.04},{.9,.5,.10,.05,.07}};
    for(int run=0;run<5;run++){
        vector<Rect> freeR; freeR.push_back({0,0,W,H}); vector<int> used(n); vector<Pl> out; long long val=0;
        // The supplied ranges make 1200 enough for normal near-optimal layouts; it is also a hard time guard.
        for(int step=0;step<1200 && !freeR.empty();step++){
            long double bs=-1e100L; int bi=-1,bf=-1,br=0; long long bw=0,bh=0;
            for(int fi=0;fi<(int)freeR.size();fi++) for(int t=0;t<n;t++) if(used[t]<a[t].lim) for(int ro=0;ro<= (allow?1:0);ro++){
                if(ro && a[t].w==a[t].h) continue;
                long long w=ro?a[t].h:a[t].w, h=ro?a[t].w:a[t].h; Rect f=freeR[fi];
                if(w>f.w||h>f.h) continue;
                long double ar=(long double)w*h, den=((long double)a[t].v/ar)/maxD;
                long double vv=(long double)a[t].v/maxV, sz=ar/binA;
                long double waste=((long double)f.w*f.h-ar)/binA, tall=(long double)(f.y+h)/H;
                long double sc=prof[run][0]*den+prof[run][1]*vv+prof[run][2]*sz-prof[run][3]*waste-prof[run][4]*tall;
                // deterministic tie diversity changes which equally ranked cavity is consumed.
                sc -= 1e-12L*((run&1)?(f.x+3*f.y):(3*f.x+f.y));
                if(sc>bs){bs=sc;bi=t;bf=fi;br=ro;bw=w;bh=h;}
            }
            if(bi<0) break;
            Rect q{freeR[bf].x,freeR[bf].y,bw,bh}; out.push_back({bi,br,q.x,q.y,bw,bh}); used[bi]++; val+=a[bi].v;
            vector<Rect> nxt; nxt.reserve(freeR.size()*2+4);
            for(Rect f:freeR) {
                if(!hit(f,q)){nxt.push_back(f);continue;}
                if(q.x>f.x) nxt.push_back({f.x,f.y,q.x-f.x,f.h});
                if(q.x+q.w<f.x+f.w) nxt.push_back({q.x+q.w,f.y,f.x+f.w-(q.x+q.w),f.h});
                if(q.y>f.y) nxt.push_back({f.x,f.y,f.w,q.y-f.y});
                if(q.y+q.h<f.y+f.h) nxt.push_back({f.x,q.y+q.h,f.w,f.y+f.h-(q.y+q.h)});
            }
            if(nxt.size()>300) { sort(nxt.begin(),nxt.end(),[](const Rect&x,const Rect&y){return x.w*x.h>y.w*y.h;}); nxt.resize(300); }
            vector<Rect> clean; clean.reserve(nxt.size());
            for(int x=0;x<(int)nxt.size();x++) if(nxt[x].w>0&&nxt[x].h>0){ bool covered=false; for(int y=0;y<(int)nxt.size();y++) if(x!=y && inside(nxt[x],nxt[y])) {covered=true;break;} if(!covered) clean.push_back(nxt[x]); }
            // A bounded frontier is essential under the one-second limit; discarding empty
            // rectangles can only lose options, never invalidate a placement already made.
            if(clean.size()>240) { sort(clean.begin(),clean.end(),[](const Rect&x,const Rect&y){return x.w*x.h>y.w*y.h;}); clean.resize(240); }
            freeR.swap(clean);
        }
        if(val>best){best=val;answer.swap(out);}
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<answer.size();i++){ if(i)cout<<','; cout<<"{\"type\":"; jsonString(a[answer[i].t].name); cout<<",\"x\":"<<answer[i].x<<",\"y\":"<<answer[i].y<<",\"rot\":"<<answer[i].rot<<'}'; }
    cout << "]}";
}
