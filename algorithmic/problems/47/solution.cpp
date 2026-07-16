#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h,lim; long long v; };
struct R { int x,y,w,h; };
struct P { int t,x,y; bool rot; };
struct Result { long long value=0; vector<P> p; };

// Small JSON reader: input is restricted to the documented object format.
struct Json {
    string s; size_t i=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(i<s.size() && isspace((unsigned char)s[i])) ++i; }
    void ch(char c){ ws(); if(i<s.size() && s[i]==c) ++i; }
    string str(){ ws(); string r; if(i<s.size()&&s[i]=='"') ++i; while(i<s.size()&&s[i]!='"'){ if(s[i]=='\\'&&i+1<s.size()) ++i; r+=s[i++]; } if(i<s.size()) ++i; return r; }
    long long num(){ ws(); int z=1; if(s[i]=='-') z=-1,++i; long long x=0; while(i<s.size()&&isdigit((unsigned char)s[i])) x=x*10+s[i++]-'0'; return z*x; }
    bool boolean(){ ws(); bool v=s.compare(i,4,"true")==0; i+=v?4:5; return v; }
};

static bool intersects(const R&a,const R&b){ return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }
static bool contains(const R&a,const R&b){ return b.x>=a.x&&b.y>=a.y&&b.x+b.w<=a.x+a.w&&b.y+b.h<=a.y+a.h; }

struct MaxRects {
    int W,H; vector<R> free;
    MaxRects(int W,int H):W(W),H(H){ free.push_back({0,0,W,H}); }
    bool put(int w,int h, int &ox,int &oy) {
        int best=-1, bs=INT_MAX, bl=INT_MAX, by=INT_MAX, bx=INT_MAX;
        for(int k=0;k<(int)free.size();++k) if(w<=free[k].w && h<=free[k].h){
            int a=min(free[k].w-w,free[k].h-h), b=max(free[k].w-w,free[k].h-h);
            if(make_tuple(a,b,free[k].y,free[k].x)<make_tuple(bs,bl,by,bx)) best=k,bs=a,bl=b,by=free[k].y,bx=free[k].x;
        }
        if(best<0) return false;
        R u{free[best].x,free[best].y,w,h}; ox=u.x; oy=u.y;
        vector<R> nf; nf.reserve(free.size()+4);
        for(const R& f:free){
            if(!intersects(f,u)){ nf.push_back(f); continue; }
            if(u.x>f.x) nf.push_back({f.x,f.y,u.x-f.x,f.h});
            if(u.x+u.w<f.x+f.w) nf.push_back({u.x+u.w,f.y,f.x+f.w-(u.x+u.w),f.h});
            if(u.y>f.y) nf.push_back({f.x,f.y,f.w,u.y-f.y});
            if(u.y+u.h<f.y+f.h) nf.push_back({f.x,u.y+u.h,f.w,f.y+f.h-(u.y+u.h)});
        }
        // Remove rectangles subsumed by another free rectangle.
        vector<char> dead(nf.size());
        for(int a=0;a<(int)nf.size();++a) if(!dead[a]) for(int b=0;b<(int)nf.size();++b)
            if(a!=b && !dead[b] && contains(nf[a],nf[b])) dead[b]=1;
        free.clear(); free.reserve(nf.size());
        for(int k=0;k<(int)nf.size();++k) if(!dead[k] && nf[k].w>0 && nf[k].h>0) free.push_back(nf[k]);
        return true;
    }
};

int main(){
    Json j; int W=0,H=0; bool allow=false; vector<Item> a;
    j.ch('{');
    for(int top=0;top<2;top++){
        string key=j.str(); j.ch(':');
        if(key=="bin"){
            j.ch('{'); for(int q=0;q<3;q++){ string k=j.str(); j.ch(':'); if(k=="W") W=(int)j.num(); else if(k=="H") H=(int)j.num(); else allow=j.boolean(); if(q<2)j.ch(','); } j.ch('}');
        } else {
            j.ch('['); j.ws();
            while(j.i<j.s.size()&&j.s[j.i]!=']'){
                j.ch('{'); Item z; for(int q=0;q<5;q++){ string k=j.str(); j.ch(':'); if(k=="type")z.id=j.str(); else if(k=="w")z.w=j.num(); else if(k=="h")z.h=j.num(); else if(k=="v")z.v=j.num(); else z.lim=j.num(); if(q<4)j.ch(','); } j.ch('}'); a.push_back(z); j.ws(); if(j.s[j.i]==',')j.ch(','); j.ws();
            } j.ch(']');
        }
        if(top==0)j.ch(',');
    }
    Result ans;
    // Different exponents trade raw value against value density, avoiding dependence on sample-specific ordering.
    const long double exps[]={1.0L,.80L,1.20L,.60L,1.45L};
    for(int run=0;run<10;run++){
        long double e=exps[run%5]; vector<int> ord(a.size()); iota(ord.begin(),ord.end(),0);
        sort(ord.begin(),ord.end(),[&](int x,int y){
            long double ax=(long double)a[x].v/powl((long double)a[x].w*a[x].h,e);
            long double ay=(long double)a[y].v/powl((long double)a[y].w*a[y].h,e);
            if(fabsl(ax-ay)>1e-18L) return ax>ay;
            return (run<5 ? a[x].v>a[y].v : a[x].w*a[x].h<a[y].w*a[y].h);
        });
        MaxRects mr(W,H); Result cur;
        for(int id:ord){
            for(int c=0;c<a[id].lim;c++){
                bool firstRot=allow && (run&1); int x,y;
                int fw=firstRot?a[id].h:a[id].w, fh=firstRot?a[id].w:a[id].h;
                if(mr.put(fw,fh,x,y)) { cur.p.push_back({id,x,y,firstRot}); cur.value+=a[id].v; continue; }
                if(allow && a[id].w!=a[id].h && mr.put(firstRot?a[id].w:a[id].h,firstRot?a[id].h:a[id].w,x,y)) {
                    cur.p.push_back({id,x,y,!firstRot}); cur.value+=a[id].v; continue;
                }
                break;
            }
        }
        if(cur.value>ans.value) ans=move(cur);
    }
    auto esc=[](const string&x){ string r; for(char c:x){ if(c=='"'||c=='\\')r+='\\'; r+=c; } return r; };
    cout<<"{\"placements\":[";
    for(size_t k=0;k<ans.p.size();++k){ auto q=ans.p[k]; if(k)cout<<','; cout<<"{\"type\":\""<<esc(a[q.t].id)<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<(q.rot?1:0)<<"}"; }
    cout<<"]}";
}
