#include <bits/stdc++.h>
using namespace std;

struct Item { string name; int w,h; long long v; int lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot; };

// The input is deliberately parsed here rather than depending on a JSON library, which
// is not part of the standard contest toolchain.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); string r; if(p>=s.size() || s[p++]!='"') return r;
        while(p<s.size() && s[p]!='"') { char c=s[p++]; if(c=='\\' && p<s.size()) { char e=s[p++]; if(e=='n') r+='\n'; else if(e=='t') r+='\t'; else r+=e; } else r+=c; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); int sg=1; if(s[p]=='-') sg=-1,++p; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sg*x; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0) {p+=4;return true;} p+=5;return false; }
};

static bool intersects(const R&a,const R&b) {
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}

struct Pack {
    int W,H; const vector<Item>& a; bool canrot; vector<R> free; vector<P> out; vector<int> used;
    Pack(int W,int H,const vector<Item>&a,bool r):W(W),H(H),a(a),canrot(r),free{{0,0,W,H}},used(a.size()){}
    void split(R u) {
        vector<R> nf;
        for(const R& f:free) {
            if(!intersects(f,u)) { nf.push_back(f); continue; }
            // These four rectangles cover the old free rectangle except for its overlap
            // with u.  They may overlap each other; containment pruning below is the
            // standard maximal-rectangles representation.
            if(u.x>f.x) nf.push_back({f.x,f.y,u.x-f.x,f.h});
            if(u.x+u.w<f.x+f.w) nf.push_back({u.x+u.w,f.y,f.x+f.w-u.x-u.w,f.h});
            if(u.y>f.y) nf.push_back({f.x,f.y,f.w,u.y-f.y});
            if(u.y+u.h<f.y+f.h) nf.push_back({f.x,u.y+u.h,f.w,f.y+f.h-u.y-u.h});
        }
        vector<R> q; q.reserve(nf.size());
        for(int i=0;i<(int)nf.size();++i) {
            if(nf[i].w<=0||nf[i].h<=0) continue;
            bool inside=false;
            for(int j=0;j<(int)nf.size();++j) if(i!=j && nf[j].x<=nf[i].x && nf[j].y<=nf[i].y && nf[j].x+nf[j].w>=nf[i].x+nf[i].w && nf[j].y+nf[j].h>=nf[i].y+nf[i].h) { inside=true; break; }
            if(!inside) q.push_back(nf[i]);
        }
        free.swap(q);
    }
    bool put(int t, int style) {
        int bi=-1, br=0, bw=0,bh=0; long long best1=LLONG_MAX,best2=LLONG_MAX; int by=INT_MAX,bx=INT_MAX;
        for(int r=0;r<= (canrot && a[t].w!=a[t].h);++r) {
            int w=r?a[t].h:a[t].w, h=r?a[t].w:a[t].h;
            for(int i=0;i<(int)free.size();++i) if(w<=free[i].w&&h<=free[i].h) {
                long long s1,s2;
                if(style==0) { s1=min(free[i].w-w,free[i].h-h); s2=max(free[i].w-w,free[i].h-h); }
                else if(style==1) { s1=(long long)free[i].w*free[i].h-(long long)w*h; s2=min(free[i].w-w,free[i].h-h); }
                else { s1=free[i].y; s2=free[i].x; }
                // Stable bottom/left ties make the layouts compact and reproducible.
                if(s1<best1 || (s1==best1&&(s2<best2 || (s2==best2&&(free[i].y<by || (free[i].y==by&&free[i].x<bx)))))) {
                    best1=s1;best2=s2;bi=i;br=r;bw=w;bh=h;by=free[i].y;bx=free[i].x;
                }
            }
        }
        if(bi<0) return false;
        R z{bx,by,bw,bh}; split(z); out.push_back({t,bx,by,br}); ++used[t]; return true;
    }
    vector<P> run(const vector<int>& ord,int style) {
        for(int t:ord) while(used[t]<a[t].lim && put(t,style));
        return out;
    }
};

static string esc(const string& x) { string r; for(char c:x) { if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }
int main(){
    Json j; j.ch('{');
    int W=0,H=0; bool rotate=false; vector<Item> a;
    for(int top=0;top<2;++top) {
        if(top) j.ch(','); string key=j.str(); j.ch(':');
        if(key=="bin") {
            j.ch('{'); for(int z=0;z<3;++z) { if(z)j.ch(','); string k=j.str();j.ch(':'); if(k=="W")W=j.num(); else if(k=="H")H=j.num(); else rotate=j.boolean(); } j.ch('}');
        } else {
            j.ch('['); j.ws(); bool first=true;
            while(j.p<j.s.size()&&j.s[j.p]!=']') { if(!first)j.ch(','); first=false; j.ch('{'); Item it; for(int z=0;z<5;++z) {if(z)j.ch(',');string k=j.str();j.ch(':'); if(k=="type")it.name=j.str(); else if(k=="w")it.w=j.num(); else if(k=="h")it.h=j.num(); else if(k=="v")it.v=j.num(); else it.lim=j.num();} j.ch('}');a.push_back(it);j.ws(); } j.ch(']');
        }
    }
    vector<P> answer; long long best=-1;
    vector<int> base(a.size()); iota(base.begin(),base.end(),0);
    // Several deterministic orderings are cheap and protect against density ties and
    // against a valuable, low-density item being fragmented out by tiny rectangles.
    for(int trial=0;trial<14;++trial) {
        vector<int> o=base;
        sort(o.begin(),o.end(),[&](int x,int y){
            long double dx=(long double)a[x].v/(a[x].w*a[x].h), dy=(long double)a[y].v/(a[y].w*a[y].h);
            long double sx=dx,sy=dy;
            if(trial==1) sx=a[x].v,sy=a[y].v;
            else if(trial==2) sx=dx*sqrt((long double)a[x].w*a[x].h),sy=dy*sqrt((long double)a[y].w*a[y].h);
            else if(trial>=3) { unsigned hx=(unsigned)(x*1103515245u+trial*12345u), hy=(unsigned)(y*1103515245u+trial*12345u); sx=dx*(0.78L+(hx%1000)/1000.0L*.44L); sy=dy*(0.78L+(hy%1000)/1000.0L*.44L); }
            if(fabsl(sx-sy)>1e-15L) return sx>sy; return x<y;
        });
        Pack pk(W,H,a,rotate); auto got=pk.run(o,trial%3); long long val=0; for(auto&p:got) val+=a[p.t].v;
        if(val>best) best=val,answer=move(got);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<answer.size();++i) { if(i)cout<<','; auto&p=answer[i]; cout<<"{\"type\":\""<<esc(a[p.t].name)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
