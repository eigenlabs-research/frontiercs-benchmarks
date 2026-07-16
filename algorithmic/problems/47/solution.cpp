#include <bits/stdc++.h>
using namespace std;

// Minimal JSON reader for the (small, fixed-schema) input.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char take(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){
        ws(); string r; if(p>=s.size()||s[p++]!='"') return r;
        while(p<s.size() && s[p]!='"') { char c=s[p++]; if(c=='\\' && p<s.size()) { char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else r+=e; } else r+=c; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); long long z=1,x=0; if(s[p]=='-') z=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return z*x; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};
struct Item { string id; int w,h; long long v; int lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };
int W,H; bool canrot; vector<Item> a;

static void skipValue(Json& j) { // input has no unknown fields, retained for robust key order
    j.ws(); if(j.p>=j.s.size()) return;
    if(j.s[j.p]=='"') { j.str(); return; }
    if(j.s[j.p]=='{' || j.s[j.p]=='[') { char o=j.s[j.p++], c=o=='{'?'}':']'; int d=1; while(j.p<j.s.size()&&d){ if(j.s[j.p]=='"'){j.str();continue;} if(j.s[j.p++]==o)d++; else if(j.s[j.p-1]==c)d--; } return; }
    while(j.p<j.s.size() && string(",}]").find(j.s[j.p])==string::npos) ++j.p;
}
void parse(){
    Json j; j.take();
    while(true){ j.ws(); if(j.take()=='}') break; --j.p; string key=j.str(); j.take();
        if(key=="bin") {
            j.take(); while(true){ j.ws(); if(j.take()=='}') break; --j.p; string k=j.str(); j.take();
                if(k=="W") W=(int)j.num(); else if(k=="H") H=(int)j.num(); else if(k=="allow_rotate") canrot=j.boolean(); else skipValue(j);
                j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p;
            }
        } else if(key=="items") {
            j.take(); while(true){ j.ws(); if(j.take()==']') break; --j.p; j.take(); Item q;
                while(true){ j.ws(); if(j.take()=='}') break; --j.p; string k=j.str(); j.take();
                    if(k=="type")q.id=j.str(); else if(k=="w")q.w=(int)j.num(); else if(k=="h")q.h=(int)j.num(); else if(k=="v")q.v=j.num(); else if(k=="limit")q.lim=(int)j.num(); else skipValue(j);
                    j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p;
                } a.push_back(q); j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p;
            }
        } else skipValue(j);
        j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p;
    }
}
// Remove all maximal free rectangles containing another one.
void prune(vector<R>& f){
    vector<char> dead(f.size());
    for(int i=0;i<(int)f.size();++i) for(int k=0;k<(int)f.size();++k) if(i!=k && !dead[i])
        if(f[k].x<=f[i].x && f[k].y<=f[i].y && f[k].x+f[k].w>=f[i].x+f[i].w && f[k].y+f[k].h>=f[i].y+f[i].h) dead[i]=1;
    vector<R> g; g.reserve(f.size()); for(int i=0;i<(int)f.size();++i) if(!dead[i]&&f[i].w>0&&f[i].h>0) g.push_back(f[i]); f.swap(g);
}
void occupy(vector<R>& f, R u){
    vector<R> g; g.reserve(f.size()*2);
    for(R r:f){
        int lx=max(r.x,u.x), rx=min(r.x+r.w,u.x+u.w), by=max(r.y,u.y), ty=min(r.y+r.h,u.y+u.h);
        if(lx>=rx || by>=ty) { g.push_back(r); continue; }
        if(u.x>r.x) g.push_back({r.x,r.y,u.x-r.x,r.h});
        if(u.x+u.w<r.x+r.w) g.push_back({u.x+u.w,r.y,r.x+r.w-u.x-u.w,r.h});
        if(u.y>r.y) g.push_back({r.x,r.y,r.w,u.y-r.y});
        if(u.y+u.h<r.y+r.h) g.push_back({r.x,u.y+u.h,r.w,r.y+r.h-u.y-u.h});
    }
    f.swap(g); prune(f);
}
long long value(const vector<P>& p){ long long z=0; for(auto q:p) z+=a[q.t].v; return z; }
bool valid(const vector<P>& p){
    vector<int> cnt(a.size());
    for(int i=0;i<(int)p.size();++i){ auto q=p[i]; if(q.t<0||q.t>=(int)a.size()||q.x<0||q.y<0||q.x+q.w>W||q.y+q.h>H||q.rot&& !canrot||++cnt[q.t]>a[q.t].lim) return false;
        for(int k=0;k<i;k++){auto z=p[k]; if(!(q.x+q.w<=z.x||z.x+z.w<=q.x||q.y+q.h<=z.y||z.y+z.h<=q.y))return false;}}
    return true;
}
vector<P> pack(int mode){
    vector<R> f={{0,0,W,H}}; vector<int> used(a.size()); vector<P> out;
    // Static priorities give a deliberate supply-first ordering; modes vary density/value and shape bias.
    vector<int> ord(a.size()); iota(ord.begin(),ord.end(),0);
    auto pri=[&](int i){ long double ar=(long double)a[i].w*a[i].h, d=a[i].v/ar;
        if(mode/4==0) return d;
        if(mode/4==1) return d*sqrt(ar);
        if(mode/4==2) return (long double)a[i].v;
        return d/sqrt(ar);
    };
    sort(ord.begin(),ord.end(),[&](int x,int y){long double px=pri(x),py=pri(y); if(px!=py)return px>py; return a[x].v>a[y].v;});
    while(!f.empty()){
        int bt=-1,bfi=-1,br=0; long long b1=LLONG_MAX,b2=LLONG_MAX,b3=LLONG_MAX;
        // Highest-priority remaining type that fits wins; choose its best maximal rectangle.
        for(int t:ord) if(used[t]<a[t].lim){
            bool any=false; long long x1=LLONG_MAX,x2=LLONG_MAX,x3=LLONG_MAX; int xi=-1,xr=0;
            for(int k=0;k<(int)f.size();++k) for(int rot=0;rot<=canrot;rot++){
                int w=rot?a[t].h:a[t].w,h=rot?a[t].w:a[t].h; if(w>f[k].w||h>f[k].h)continue;
                long long dw=f[k].w-w, dh=f[k].h-h;
                long long q1,q2,q3;
                if(mode%4==0) q1=min(dw,dh),q2=max(dw,dh),q3=dw*dh;
                else if(mode%4==1) q1=dw*dh,q2=min(dw,dh),q3=max(dw,dh);
                else if(mode%4==2) q1=max(dw,dh),q2=min(dw,dh),q3=dw*dh;
                else q1=dw+dh,q2=dw*dh,q3=min(dw,dh);
                if(!any||make_tuple(q1,q2,q3,k,rot)<make_tuple(x1,x2,x3,xi,xr)) any=true,x1=q1,x2=q2,x3=q3,xi=k,xr=rot;
            }
            if(any){bt=t;bfi=xi;br=xr;b1=x1;b2=x2;b3=x3;break;}
        }
        if(bt<0)break;
        int w=br?a[bt].h:a[bt].w,h=br?a[bt].w:a[bt].h; R u={f[bfi].x,f[bfi].y,w,h};
        out.push_back({bt,u.x,u.y,br,w,h}); ++used[bt]; occupy(f,u);
        // Guard against pathological free-list growth; usual instances remain far below this.
        if(f.size()>5000) prune(f);
    }
    return out;
}
string esc(const string& x){ string r; for(char c:x){if(c=='"'||c=='\\')r+='\\';r+=c;}return r; }
int main(){
    parse(); vector<P> best;
    // Eight deterministic starts keep the quadratic free-rectangle cleanup bounded for 1s cases.
    for(int m=0;m<8;m++){ vector<P> q=pack(m); if(valid(q)&&value(q)>value(best)) best.swap(q); }
    cout<<"{\"placements\":[";
    for(int i=0;i<(int)best.size();++i){if(i)cout<<','; auto q=best[i]; cout<<"{\"type\":\""<<esc(a[q.t].id)<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.rot<<'}';}
    cout<<"]}\n";
}
