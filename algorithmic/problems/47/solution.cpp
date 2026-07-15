#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h,lim; long long v; };
struct Rect { int x,y,w,h; };
struct Pl { int t,x,y,rot,w,h; };

// Small permissive JSON reader.  The instance format contains only strings, integers and bools.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); string r; if(p>=s.size() || s[p++]!='"') return r;
        while(p<s.size() && s[p]!='"') {
            char c=s[p++];
            if(c=='\\' && p<s.size()) { char e=s[p++];
                if(e=='b') r+='\b'; else if(e=='f') r+='\f'; else if(e=='n') r+='\n';
                else if(e=='r') r+='\r'; else if(e=='t') r+='\t'; else r+=e;
            } else r+=c;
        }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); int sign=1; if(p<s.size()&&(s[p]=='-'||s[p]=='+')) sign=s[p++]=='-'?-1:1; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sign*x; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4;return true;} p+=5;return false; }
};
static bool inter(const Rect&a,const Rect&b) {
    return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h;
}
static bool contains(const Rect&a,const Rect&b) {
    return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h;
}

struct Result { long long value=0; vector<Pl> p; };
int W,H; bool rotateOK; vector<Item> it;

// A bounded maximal-rectangles construction.  Keeping a bounded collection is important on
// instances with many tiny pieces: it preserves the useful large holes without quadratic blowup.
Result pack(int variant) {
    vector<Rect> fr(1, {0,0,W,H});
    vector<int> used(it.size()); Result ans;
    const long double binArea=(long double)W*H;
    for(;;) {
        int bi=-1,bf=-1,br=0; long double bestRank=-1e100L, bestA=1e100L,bestB=1e100L;
        int bestY=INT_MAX,bestX=INT_MAX;
        for(int i=0;i<(int)it.size();++i) if(used[i]<it[i].lim) {
            long double area=(long double)it[i].w*it[i].h;
            long double den=(long double)it[i].v/area;
            long double rank=den;
            if(variant==3) rank=den*(1.0L+0.20L*area/binArea);
            if(variant==4) rank=den/(1.0L+0.16L*area/binArea);
            if(variant==5) rank=(long double)it[i].v; // occasionally a scarce large prize is worth reserving space for
            for(int r=0;r<=(rotateOK && it[i].w!=it[i].h);++r) {
                int ww=r?it[i].h:it[i].w, hh=r?it[i].w:it[i].h;
                for(int f=0;f<(int)fr.size();++f) if(ww<=fr[f].w && hh<=fr[f].h) {
                    int dx=fr[f].w-ww, dy=fr[f].h-hh;
                    long double a,b;
                    int policy=variant%3;
                    if(policy==0) { a=min(dx,dy); b=max(dx,dy); }
                    else if(policy==1) { a=(long double)fr[f].w*fr[f].h-(long double)ww*hh; b=min(dx,dy); }
                    else { a=fr[f].y; b=fr[f].x; } // bottom-left alternative
                    bool take=false;
                    if(rank>bestRank+1e-12L) take=true;
                    else if(fabsl(rank-bestRank)<1e-12L) {
                        if(a<bestA-1e-12L || (fabsl(a-bestA)<1e-12L && (b<bestB-1e-12L ||
                           (fabsl(b-bestB)<1e-12L && (fr[f].y<bestY || (fr[f].y==bestY && fr[f].x<bestX)))))) take=true;
                    }
                    if(take) bi=i,bf=f,br=r,bestRank=rank,bestA=a,bestB=b,bestY=fr[f].y,bestX=fr[f].x;
                }
            }
        }
        if(bi<0) break;
        Rect q{fr[bf].x,fr[bf].y,br?it[bi].h:it[bi].w,br?it[bi].w:it[bi].h};
        ans.p.push_back({bi,q.x,q.y,br,q.w,q.h}); ans.value+=it[bi].v; ++used[bi];
        vector<Rect> nf; nf.reserve(fr.size()*2+4);
        for(const Rect& a:fr) {
            if(!inter(a,q)) { nf.push_back(a); continue; }
            int ar=a.x+a.w, at=a.y+a.h, qr=q.x+q.w, qt=q.y+q.h;
            if(q.x>a.x) nf.push_back({a.x,a.y,q.x-a.x,a.h});
            if(qr<ar) nf.push_back({qr,a.y,ar-qr,a.h});
            if(q.y>a.y) nf.push_back({a.x,a.y,a.w,q.y-a.y});
            if(qt<at) nf.push_back({a.x,qt,a.w,at-qt});
        }
        // First discard small fragments, then do the (much cheaper) maximality cleanup.
        if(nf.size()>110) {
            nth_element(nf.begin(),nf.begin()+110,nf.end(),[](const Rect&a,const Rect&b){return (long long)a.w*a.h>(long long)b.w*b.h;});
            nf.resize(110);
        }
        vector<char> dead(nf.size());
        for(int a=0;a<(int)nf.size();++a) for(int b=0;b<(int)nf.size();++b) if(a!=b && contains(nf[b],nf[a])) { dead[a]=1; break; }
        fr.clear(); for(int a=0;a<(int)nf.size();++a) if(!dead[a] && nf[a].w>0 && nf[a].h>0) fr.push_back(nf[a]);
        if(fr.size()>80) { sort(fr.begin(),fr.end(),[](const Rect&a,const Rect&b){return (long long)a.w*a.h>(long long)b.w*b.h;}); fr.resize(80); }
        if(ans.p.size()>=100000) break;
    }
    return ans;
}

static string esc(const string& x) { string r; for(unsigned char c:x) { if(c=='"'||c=='\\') r+='\\',r+=c; else if(c=='\n') r+="\\n"; else if(c=='\r') r+="\\r"; else if(c=='\t') r+="\\t"; else if(c=='\b') r+="\\b"; else if(c=='\f') r+="\\f"; else r+=c; } return r; }
int main(){
    Json j; j.ch('{');
    for(int z=0;z<2;z++) { if(z) j.ch(','); string key=j.str(); j.ch(':');
        if(key=="bin") { j.ch('{'); for(int k=0;k<3;k++){ if(k)j.ch(','); string q=j.str();j.ch(':'); if(q=="W")W=j.num(); else if(q=="H")H=j.num(); else rotateOK=j.boolean(); } j.ch('}'); }
        else { j.ch('['); bool first=true; while(true){ j.ws(); if(j.p<j.s.size()&&j.s[j.p]==']'){++j.p;break;} if(!first)j.ch(',');first=false; j.ch('{'); Item a; for(int k=0;k<5;k++){if(k)j.ch(',');string q=j.str();j.ch(':'); if(q=="type")a.id=j.str(); else if(q=="w")a.w=j.num(); else if(q=="h")a.h=j.num(); else if(q=="v")a.v=j.num(); else a.lim=j.num();}j.ch('}');it.push_back(a); } }
    }
    Result best; for(int v=0;v<6;v++){ Result r=pack(v); if(r.value>best.value) best=move(r); }
    cout << "{\"placements\":[";
    for(size_t k=0;k<best.p.size();++k){ if(k)cout<<','; const Pl&p=best.p[k]; cout<<"{\"type\":\""<<esc(it[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
