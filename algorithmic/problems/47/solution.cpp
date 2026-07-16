#include <bits/stdc++.h>
using namespace std;

struct J {
    string s; size_t p=0;
    J(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); string r; if(p<s.size()&&s[p]=='"') ++p;
        while(p<s.size() && s[p]!='"') { if(s[p]=='\\' && p+1<s.size()) { ++p; char c=s[p++]; if(c=='n')r+='\n'; else if(c=='t')r+='\t'; else r+=c; } else r+=s[p++]; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); long long z=0,sg=1; if(s[p]=='-')sg=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p]))z=z*10+s[p++]-'0'; return z*sg; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};
struct It { string id; int w,h; long long v; int lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot; };
struct Sol { vector<P> p; long long val=0; };
int W,H,n; bool canrot; vector<It>a;

static void prune(vector<R>& f) {
    vector<R> q; q.reserve(f.size());
    for(int i=0;i<(int)f.size();++i) {
        if(f[i].w<=0||f[i].h<=0) continue;
        bool inside=false;
        for(int j=0;j<(int)f.size();++j) if(i!=j && f[j].x<=f[i].x && f[j].y<=f[i].y && f[j].x+f[j].w>=f[i].x+f[i].w && f[j].y+f[j].h>=f[i].y+f[i].h) { inside=true; break; }
        if(!inside) q.push_back(f[i]);
    }
    f.swap(q);
}
static void occupy(vector<R>& f, R u) {
    vector<R> q; q.reserve(f.size()*2);
    for(R r:f) {
        if(u.x>=r.x+r.w || u.x+u.w<=r.x || u.y>=r.y+r.h || u.y+u.h<=r.y) { q.push_back(r); continue; }
        if(u.x>r.x) q.push_back({r.x,r.y,u.x-r.x,r.h});
        if(u.x+u.w<r.x+r.w) q.push_back({u.x+u.w,r.y,r.x+r.w-(u.x+u.w),r.h});
        if(u.y>r.y) q.push_back({r.x,r.y,r.w,u.y-r.y});
        if(u.y+u.h<r.y+r.h) q.push_back({r.x,u.y+u.h,r.w,r.y+r.h-(u.y+u.h)});
    }
    f.swap(q); prune(f);
}

Sol build(vector<int> ord, int mode, int rotmode) {
    vector<R> free{{0,0,W,H}}; vector<int> used(n); Sol z;
    for(int t:ord) for(int rep=0;rep<a[t].lim;rep++) {
        int bi=-1, br=0; long long best=LLONG_MAX;
        for(int i=0;i<(int)free.size();++i) for(int ro=0;ro<= (canrot && a[t].w!=a[t].h);++ro) {
            if(rotmode==1 && canrot && a[t].w!=a[t].h) ro=1;
            int w=ro?a[t].h:a[t].w, h=ro?a[t].w:a[t].h;
            if(w>free[i].w||h>free[i].h) { if(rotmode==1) break; else continue; }
            long long dx=free[i].w-w, dy=free[i].h-h, key;
            if(mode==0) key=min(dx,dy)*10000000LL+max(dx,dy)*1000LL+free[i].y;
            else if(mode==1) key=(long long)free[i].y*10000000LL+free[i].x*1000LL+min(dx,dy);
            else key=dx*dy*1000LL+min(dx,dy);
            if(key<best) best=key,bi=i,br=ro;
            if(rotmode==1) break;
        }
        if(bi<0) break;
        int w=br?a[t].h:a[t].w,h=br?a[t].w:a[t].h;
        R u{free[bi].x,free[bi].y,w,h}; occupy(free,u);
        z.p.push_back({t,u.x,u.y,br}); z.val+=a[t].v;
    }
    return z;
}
int main(){
    J j; j.ch('{');
    // Input key order is not semantically significant.
    for(int top=0;top<2;++top){ string k=j.str(); j.ch(':');
        if(k=="bin") { j.ch('{'); for(int q=0;q<3;++q){ string x=j.str();j.ch(':'); if(x=="W")W=(int)j.num(); else if(x=="H")H=(int)j.num(); else canrot=j.boolean(); if(q<2)j.ch(','); } j.ch('}'); }
        else { j.ch('['); j.ws(); while(j.p<j.s.size()&&j.s[j.p]!=']') { j.ch('{'); It e; for(int q=0;q<5;++q){ string x=j.str();j.ch(':'); if(x=="type")e.id=j.str(); else if(x=="w")e.w=j.num(); else if(x=="h")e.h=j.num(); else if(x=="v")e.v=j.num(); else e.lim=j.num(); if(q<4)j.ch(','); }j.ch('}');a.push_back(e);j.ws();if(j.s[j.p]==',')j.ch(','); }j.ch(']'); }
        if(top==0)j.ch(',');
    }
    n=a.size(); vector<int> base(n); iota(base.begin(),base.end(),0); Sol ans;
    // Several orderings are alternate tie-breaking views of the same MaxRects heuristic.
    for(int variant=0;variant<12;++variant) {
        vector<int> o=base;
        sort(o.begin(),o.end(),[&](int x,int y){
            long double dx=(long double)a[x].v/(a[x].w*a[x].h), dy=(long double)a[y].v/(a[y].w*a[y].h);
            if(variant%4==0) { if(dx!=dy)return dx>dy; }
            else if(variant%4==1) { if(a[x].v!=a[y].v)return a[x].v>a[y].v; }
            else if(variant%4==2) { if(a[x].w*a[x].h!=a[y].w*a[y].h)return a[x].w*a[x].h<a[y].w*a[y].h; }
            else { long double sx=dx*(1.0L+0.13L*((x*37+variant*11)%7)), sy=dy*(1.0L+0.13L*((y*37+variant*11)%7)); if(sx!=sy)return sx>sy; }
            int ax=max(a[x].w,a[x].h), ay=max(a[y].w,a[y].h); return (variant&1)?ax<ay:ax>ay;
        });
        Sol cur=build(o,(variant/4)%3,(variant%3==2)?1:0);
        if(cur.val>ans.val) ans=move(cur);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<ans.p.size();++i){ if(i)cout<<','; auto p=ans.p[i]; cout<<"{\"type\":\""<<a[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}";
}
