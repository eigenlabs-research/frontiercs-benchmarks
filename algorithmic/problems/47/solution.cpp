#include <bits/stdc++.h>
using namespace std;

// A small permissive JSON reader; input is the fixed object described in the task.
struct JS {
    string s; size_t p=0;
    JS(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ws(); return p<s.size()?s[p]:0;}
    void eat(char c){ws(); if(p<s.size()&&s[p]==c) ++p;}
    string str(){ eat('"'); string r; while(p<s.size()&&s[p]!='"') { char c=s[p++]; if(c=='\\'&&p<s.size()) { char q=s[p++]; if(q=='n') r+='\n'; else if(q=='t') r+='\t'; else r+=q; } else r+=c; } eat('"'); return r; }
    long long num(){ ws(); int sg=1; if(s[p]=='-') sg=-1,++p; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return x*sg; }
    bool boolean(){ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r;}
};
struct It { string name; int w,h; long long v, lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };
struct Pack { vector<P> a; long long val=0; };
int W,H,N; bool allowRot; vector<It> it;

static bool contains(const R&a,const R&b) {
    return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h;
}
static bool intersects(const R&a,const R&b) {
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}

// Maximal empty rectangles: splitting every intersected empty rectangle is conservative,
// and every subsequently selected rectangle is therefore disjoint from all placed items.
struct MR {
    vector<R> f;
    vector<P> out;
    vector<int> used;
    long long value=0;
    MR(): f(1,{0,0,W,H}),used(N,0) {}
    void put(int t,int fi,int rot) {
        int w=rot?it[t].h:it[t].w, h=rot?it[t].w:it[t].h;
        R q{f[fi].x,f[fi].y,w,h};
        vector<R> nf; nf.reserve(f.size()+8);
        for(const R& z:f) {
            if(!intersects(z,q)) { nf.push_back(z); continue; }
            if(q.x>z.x) nf.push_back({z.x,z.y,q.x-z.x,z.h});
            if(q.x+q.w<z.x+z.w) nf.push_back({q.x+q.w,z.y,z.x+z.w-(q.x+q.w),z.h});
            if(q.y>z.y) nf.push_back({z.x,z.y,z.w,q.y-z.y});
            if(q.y+q.h<z.y+z.h) nf.push_back({z.x,q.y+q.h,z.w,z.y+z.h-(q.y+q.h)});
        }
        // Delete contained maximal rectangles.  Keeping overlapping (but empty) maxima is intended.
        vector<char> bad(nf.size());
        for(int i=0;i<(int)nf.size();++i) if(nf[i].w>0&&nf[i].h>0)
            for(int j=0;j<(int)nf.size();++j) if(i!=j && contains(nf[j],nf[i]) &&
                (nf[j].x!=nf[i].x || nf[j].y!=nf[i].y || nf[j].w!=nf[i].w || nf[j].h!=nf[i].h || j<i)) { bad[i]=1; break; }
        f.clear(); f.reserve(nf.size());
        for(int i=0;i<(int)nf.size();++i) if(!bad[i]&&nf[i].w>0&&nf[i].h>0) f.push_back(nf[i]);
        out.push_back({t,q.x,q.y,rot,w,h}); ++used[t]; value+=it[t].v;
    }
    // rule: 0 short-side, 1 area waste, 2 bottom-left. Returns an admissible placement.
    bool addType(int t,int rule) {
        if(used[t]>=it[t].lim) return false;
        int bi=-1, br=0; long long bk1=LLONG_MAX,bk2=LLONG_MAX,bk3=LLONG_MAX;
        for(int i=0;i<(int)f.size();++i) for(int r=0;r<= (allowRot && it[t].w!=it[t].h);++r) {
            int w=r?it[t].h:it[t].w,h=r?it[t].w:it[t].h;
            if(w>f[i].w||h>f[i].h) continue;
            long long a=f[i].w-w,b=f[i].h-h, k1,k2,k3;
            if(rule==0) k1=min(a,b),k2=max(a,b),k3=(long long)f[i].w*f[i].h-w*h;
            else if(rule==1) k1=(long long)f[i].w*f[i].h-w*h,k2=min(a,b),k3=max(a,b);
            else k1=f[i].y,k2=f[i].x,k3=(long long)f[i].w*f[i].h-w*h;
            if(bi<0 || tie(k1,k2,k3,i,r) < tie(bk1,bk2,bk3,bi,br)) bi=i,br=r,bk1=k1,bk2=k2,bk3=k3;
        }
        if(bi<0) return false; put(t,bi,br); return true;
    }
    static tuple<long long,long long,long long,int,int> tie(long long a,long long b,long long c,int d,int e){return {a,b,c,d,e};}
};

Pack ordered(vector<int> ord,int rule) {
    MR m;
    for(int t:ord) while(m.addType(t,rule));
    return {move(m.out),m.value};
}
// Interleaved ablation: choose the currently fitting item globally.  A small fit term
// distinguishes this from exhausting a type, while density remains the dominant objective.
Pack interleaved(int rule) {
    MR m;
    const long double binA=(long double)W*H;
    while(true) {
        int bt=-1; long double bs=-1e100L; long long bv=-1;
        for(int t=0;t<N;++t) if(m.used[t]<it[t].lim) {
            // Test feasibility and estimate the best free-rectangle fit without modifying m.
            long long waste=LLONG_MAX; bool ok=false;
            for(const R& z:m.f) for(int r=0;r<=(allowRot&&it[t].w!=it[t].h);++r) { int w=r?it[t].h:it[t].w,h=r?it[t].w:it[t].h; if(w<=z.w&&h<=z.h) ok=true,waste=min(waste,(long long)z.w*z.h-(long long)w*h); }
            if(!ok) continue;
            long double den=(long double)it[t].v/((long double)it[t].w*it[t].h);
            long double fit=1.0L-(long double)waste/binA;
            long double score=den;
            if(rule==1) score=den*(0.92L+0.08L*fit);
            if(rule==2) score=den*(0.80L+0.20L*fit)+(long double)it[t].v/binA*0.02L;
            if(score>bs+1e-18L || (fabsl(score-bs)<1e-18L && it[t].v>bv)) bs=score,bv=it[t].v,bt=t;
        }
        if(bt<0 || !m.addType(bt,rule==2?1:rule)) break;
    }
    return {move(m.out),m.value};
}
static string esc(const string& s) { string r; for(unsigned char c:s) { if(c=='"'||c=='\\') r+='\\',r+=c; else if(c=='\n') r+="\\n"; else if(c=='\r') r+="\\r"; else if(c=='\t') r+="\\t"; else r+=c; } return r; }
int main(){
    JS j; j.eat('{');
    // Keys may legally arrive in either order.
    for(int top=0;top<2;++top){ if(top)j.eat(','); string key=j.str();j.eat(':');
        if(key=="bin") { j.eat('{'); for(int z=0;z<3;++z){if(z)j.eat(',');string k=j.str();j.eat(':');if(k=="W")W=j.num();else if(k=="H")H=j.num();else allowRot=j.boolean();}j.eat('}'); }
        else { j.eat('['); bool first=true; while(j.ch()!=']'){if(!first)j.eat(',');first=false;j.eat('{');It a;for(int z=0;z<5;++z){if(z)j.eat(',');string k=j.str();j.eat(':');if(k=="type")a.name=j.str();else if(k=="w")a.w=j.num();else if(k=="h")a.h=j.num();else if(k=="v")a.v=j.num();else a.lim=j.num();}j.eat('}');it.push_back(a);}j.eat(']'); }
    } j.eat('}'); N=it.size();
    vector<int> base(N); iota(base.begin(),base.end(),0); Pack best;
    for(int mode=0;mode<4;++mode){ vector<int> o=base; sort(o.begin(),o.end(),[&](int a,int b){
        long double da=(long double)it[a].v/(it[a].w*it[a].h), db=(long double)it[b].v/(it[b].w*it[b].h);
        if(mode==0 && da!=db)return da>db;
        if(mode==1 && it[a].v!=it[b].v)return it[a].v>it[b].v;
        if(mode==2 && it[a].w*it[a].h!=it[b].w*it[b].h)return it[a].w*it[a].h<it[b].w*it[b].h;
        int ma=max(it[a].w,it[a].h),mb=max(it[b].w,it[b].h); if(ma!=mb)return ma<mb; return a<b;
    }); for(int rule=0;rule<3;++rule){Pack q=ordered(o,rule);if(q.val>best.val)best=move(q);} }
    for(int rule=0;rule<3;++rule){Pack q=interleaved(rule);if(q.val>best.val)best=move(q);}
    cout<<"{\"placements\":[";
    for(size_t i=0;i<best.a.size();++i){auto&p=best.a[i];if(i)cout<<',';cout<<"{\"type\":\""<<esc(it[p.t].name)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}
    cout<<"]}";
}
