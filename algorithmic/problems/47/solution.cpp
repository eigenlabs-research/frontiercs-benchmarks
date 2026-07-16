#include <bits/stdc++.h>
using namespace std;

// Small JSON reader: the input is deliberately parsed by names, rather than relying on
// the illustrative field ordering in the statement.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); string r; if(p>=s.size() || s[p++]!='"') return r;
        while(p<s.size() && s[p]!='"') { char c=s[p++]; if(c=='\\' && p<s.size()) { char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else if(e=='r') r+='\r'; else if(e=='b') r+='\b'; else if(e=='f') r+='\f'; else r+=e; } else r+=c; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); long long z=1,r=0; if(p<s.size() && s[p]=='-') z=-1,++p; while(p<s.size() && isdigit((unsigned char)s[p])) r=r*10+s[p++]-'0'; return z*r; }
    bool boolean(){ ws(); bool v=p<s.size()&&s[p]=='t'; p += v?4:5; return v; }
};
struct Item { string name; long long w,h,v,lim; };
struct R { long long x,y,w,h; };
struct P { int t,rot; long long x,y,w,h; };

// Emit item identifiers as JSON strings; identifiers are arbitrary JSON strings, not
// necessarily the alphanumeric examples used in the statement.
static void emitString(const string& s) {
    cout << '"';
    for(unsigned char c : s) {
        if(c=='"' || c=='\\') cout << '\\' << (char)c;
        else if(c=='\n') cout << "\\n";
        else if(c=='\r') cout << "\\r";
        else if(c=='\t') cout << "\\t";
        else if(c=='\b') cout << "\\b";
        else if(c=='\f') cout << "\\f";
        else if(c<0x20) { static const char hex[]="0123456789abcdef"; cout << "\\u00" << hex[c>>4] << hex[c&15]; }
        else cout << (char)c;
    }
    cout << '"';
}

static bool intersects(const R&a,const R&b) {
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}
static void normalize(vector<R>& a) {
    vector<R> b; b.reserve(a.size());
    for(auto q:a) if(q.w>0 && q.h>0) b.push_back(q);
    // Do this before the quadratic containment pass.  Discarding an empty region can
    // only reduce a heuristic candidate, never invalidate its remaining placements.
    if(b.size()>280) {
        nth_element(b.begin(), b.begin()+220, b.end(), [](const R&u,const R&v){return u.w*u.h>v.w*v.h;});
        b.resize(220);
    }
    vector<char> dead(b.size());
    for(size_t i=0;i<b.size();++i) if(!dead[i])
        for(size_t j=0;j<b.size();++j) if(i!=j && !dead[j] &&
           b[i].x>=b[j].x && b[i].y>=b[j].y && b[i].x+b[i].w<=b[j].x+b[j].w && b[i].y+b[i].h<=b[j].y+b[j].h) { dead[i]=1; break; }
    a.clear(); for(size_t i=0;i<b.size();++i) if(!dead[i]) a.push_back(b[i]);
    // A pathological sequence can create many equivalent maximal spaces.  Keeping the
    // largest ones is safe (each is empty) and keeps the heuristic within the time limit.
    if(a.size()>220) {
        nth_element(a.begin(), a.begin()+180, a.end(), [](const R&u,const R&v){return u.w*u.h>v.w*v.h;});
        a.resize(180);
    }
}
static void occupy(vector<R>& f, R q) {
    vector<R> out; out.reserve(f.size()*2+4);
    for(const R& z:f) {
        if(!intersects(z,q)) { out.push_back(z); continue; }
        if(q.x>z.x) out.push_back({z.x,z.y,q.x-z.x,z.h});
        if(q.x+q.w<z.x+z.w) out.push_back({q.x+q.w,z.y,z.x+z.w-(q.x+q.w),z.h});
        if(q.y>z.y) out.push_back({z.x,z.y,z.w,q.y-z.y});
        if(q.y+q.h<z.y+z.h) out.push_back({z.x,q.y+q.h,z.w,z.y+z.h-(q.y+q.h)});
    }
    f.swap(out); normalize(f);
}

static vector<P> solveOne(long long W,long long H,bool rot,const vector<Item>& it, long double alpha,long double beta,int flavor) {
    vector<R> freeR={{0,0,W,H}}; vector<long long> used(it.size()); vector<P> ans;
    long long maxSteps=0; for(auto &z:it) maxSteps+=z.lim;
    // A bin cannot contain more than WH/min-area useful copies; this also protects malformed huge limits.
    long long mn=LLONG_MAX; for(auto &z:it) if(z.v>0) mn=min(mn,z.w*z.h);
    if(mn!=LLONG_MAX && W>0 && H>0) maxSteps=min(maxSteps, W*H/mn+1);
    for(long long step=0; step<maxSteps; ++step) {
        int bi=-1,bj=-1,br=0; long long bw=0,bh=0; long double best=-1e100L;
        for(int i=0;i<(int)it.size();++i) if(used[i]<it[i].lim && it[i].v>0) {
            for(int r=0;r<=(rot && it[i].w!=it[i].h);++r) {
                long long w=r?it[i].h:it[i].w, h=r?it[i].w:it[i].h;
                long double base=log((long double)it[i].v+1.0L)-alpha*log((long double)w*h);
                for(int j=0;j<(int)freeR.size();++j) if(w<=freeR[j].w && h<=freeR[j].h) {
                    const R& z=freeR[j]; long double unused=((long double)z.w*z.h-(long double)w*h)/((long double)z.w*z.h);
                    long double shape=(long double)(z.w-w+z.h-h)/(W+H+1);
                    long double sc=base-beta*unused-0.08L*shape;
                    // Different tie policies make equal-density grids start from different axes.
                    if(flavor==1) sc-=1e-7L*(z.x+z.y);
                    if(flavor==2) sc-=1e-7L*(z.y*W+z.x);
                    if(sc>best) best=sc,bi=i,bj=j,br=r,bw=w,bh=h;
                }
            }
        }
        if(bi<0) break;
        R q{freeR[bj].x,freeR[bj].y,bw,bh};
        ans.push_back({bi,br,q.x,q.y,bw,bh}); ++used[bi]; occupy(freeR,q);
    }
    return ans;
}
int main(){
    Json j; j.ws(); j.ch('{'); long long W=0,H=0; bool allow=false; vector<Item> items;
    for(bool first=true;;first=false) {
        j.ws(); if(j.p>=j.s.size() || j.s[j.p]=='}'){j.ch('}');break;} if(!first) j.ch(',');
        string key=j.str(); j.ch(':');
        if(key=="bin") {
            j.ch('{'); for(bool f=true;;f=false){ j.ws(); if(j.s[j.p]=='}'){j.ch('}');break;} if(!f)j.ch(','); string k=j.str();j.ch(':'); if(k=="W")W=j.num(); else if(k=="H")H=j.num(); else if(k=="allow_rotate")allow=j.boolean(); }
        } else if(key=="items") {
            j.ch('['); for(bool f=true;;f=false){ j.ws(); if(j.s[j.p]==']'){j.ch(']');break;} if(!f)j.ch(','); j.ch('{'); Item a{}; for(bool g=true;;g=false){j.ws();if(j.s[j.p]=='}'){j.ch('}');break;}if(!g)j.ch(',');string k=j.str();j.ch(':');if(k=="type")a.name=j.str();else if(k=="w")a.w=j.num();else if(k=="h")a.h=j.num();else if(k=="v")a.v=j.num();else if(k=="limit")a.lim=j.num();} items.push_back(a); }
        }
    }
    vector<pair<long double,long double>> modes={{0.0L,0.10L},{0.45L,0.12L},{0.75L,0.10L},{1.0L,0.08L},{1.2L,0.06L},{1.45L,0.04L},{0.9L,0.35L},{1.15L,0.30L},{0.6L,0.45L}};
    vector<P> best; __int128 bestVal=-1;
    for(int k=0;k<(int)modes.size();++k) { auto q=solveOne(W,H,allow,items,modes[k].first,modes[k].second,k%3); __int128 val=0;for(auto&p:q)val+=items[p.t].v; if(val>bestVal)bestVal=val,best.swap(q); }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){ if(i)cout<<','; auto&p=best[i]; cout<<"{\"type\":"; emitString(items[p.t].name); cout<<",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}";
}
