#include <bits/stdc++.h>
using namespace std;
struct Item { string id; long long w,h,v,lim; };
struct Placement { int t,rot; long long x,y,w,h; };

// The input is a small JSON object; this deliberately accepts field order changes.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p]:0; }
    void eat(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string r; if(p>=s.size()||s[p++]!='"') return r; while(p<s.size()&&s[p]!='"'){ char c=s[p++]; if(c=='\\'&&p<s.size()){char e=s[p++]; r+=(e=='n'?'\n':e=='t'?'\t':e);} else r+=c; } if(p<s.size())++p; return r; }
    long long num(){ ws(); long long z=1,x=0; if(p<s.size()&&s[p]=='-')z=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p]))x=x*10+s[p++]-'0'; return z*x; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4;return true;} p+=5;return false; }
};
struct Rect { long long x,y,w,h; };
static bool intersects(const Rect&a,const Rect&b){ return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }
static bool contains(const Rect&a,const Rect&b){ return a.x<=b.x&&a.y<=b.y&&a.x+a.w>=b.x+b.w&&a.y+a.h>=b.y+b.h; }

// Maximal-free-rectangles packing.  Every rectangle in free is known empty; it is
// harmless for free rectangles to overlap, and splitting every intersected one
// before the next choice preserves that invariant.
vector<Placement> build(const vector<Item>& it,long long W,long long H,bool canrot,int order,int fit,clock_t deadline){
    int n=it.size(); vector<int> q(n); iota(q.begin(),q.end(),0);
    sort(q.begin(),q.end(),[&](int a,int b){
        __int128 da=(__int128)it[a].v*it[b].w*it[b].h, db=(__int128)it[b].v*it[a].w*it[a].h;
        if(order==1 && it[a].v!=it[b].v) return it[a].v>it[b].v;
        if(order==2 && it[a].w*it[a].h!=it[b].w*it[b].h) return it[a].w*it[a].h>it[b].w*it[b].h;
        if(da!=db) return da>db;
        return it[a].v>it[b].v;
    });
    vector<Rect> free{{0,0,W,H}}; vector<Placement> out; long long work=0;
    for(int t:q) for(long long copy=0;copy<it[t].lim;++copy){
        if(((++work)&31)==0 && clock()>deadline) return out;
        int bi=-1,br=0; long long bw=0,bh=0; __int128 best1=-1,best2=-1;
        for(int r=0;r<=(canrot?1:0);++r){
            long long w=r?it[t].h:it[t].w, h=r?it[t].w:it[t].h;
            for(int i=0;i<(int)free.size();++i) if(w<=free[i].w&&h<=free[i].h){
                long long dx=free[i].w-w, dy=free[i].h-h;
                __int128 a1,a2;
                if(fit==0) a1=min(dx,dy),a2=max(dx,dy);          // short-side fit
                else if(fit==1) a1=(__int128)dx*free[i].h+(__int128)dy*free[i].w,a2=min(dx,dy); // area waste
                else a1=free[i].y+h,a2=free[i].x+w;              // bottom-left
                if(bi<0||a1<best1||(a1==best1&&(a2<best2||(a2==best2&&free[i].y<free[bi].y))))
                    bi=i,br=r,bw=w,bh=h,best1=a1,best2=a2;
            }
        }
        if(bi<0) break;
        Rect p{free[bi].x,free[bi].y,bw,bh}; out.push_back({t,br,p.x,p.y,p.w,p.h});
        vector<Rect> next; next.reserve(free.size()+4);
        for(const Rect&r:free){
            if(!intersects(r,p)){ next.push_back(r); continue; }
            if(p.x>r.x) next.push_back({r.x,r.y,p.x-r.x,r.h});
            if(p.x+p.w<r.x+r.w) next.push_back({p.x+p.w,r.y,r.x+r.w-(p.x+p.w),r.h});
            if(p.y>r.y) next.push_back({r.x,r.y,r.w,p.y-r.y});
            if(p.y+p.h<r.y+r.h) next.push_back({r.x,p.y+p.h,r.w,r.y+r.h-(p.y+p.h)});
        }
        // Standard maximal-rectangle pruning.  It keeps the list manageable
        // without changing the union of available space.
        vector<char> dead(next.size());
        for(int i=0;i<(int)next.size();++i) if(!dead[i]) for(int k=0;k<(int)next.size();++k)
            if(i!=k&&!dead[k]&&contains(next[k],next[i])) { dead[i]=1; break; }
        free.clear(); for(int i=0;i<(int)next.size();++i) if(!dead[i]&&next[i].w&&next[i].h) free.push_back(next[i]);
    }
    return out;
}
int main(){
    Json j; long long W=0,H=0; bool allow=false; vector<Item> items; j.eat('{');
    while(j.ch()=='"'){
        string key=j.str(); j.eat(':');
        if(key=="bin"){ j.eat('{'); while(j.ch()=='"'){ string k=j.str();j.eat(':'); if(k=="W")W=j.num();else if(k=="H")H=j.num();else if(k=="allow_rotate")allow=j.boolean(); if(j.ch()==',')j.eat(',');else break;} j.eat('}'); }
        else if(key=="items"){ j.eat('['); while(j.ch()=='{'){ j.eat('{'); Item a; while(j.ch()=='"'){string k=j.str();j.eat(':');if(k=="type")a.id=j.str();else if(k=="w")a.w=j.num();else if(k=="h")a.h=j.num();else if(k=="v")a.v=j.num();else if(k=="limit")a.lim=j.num();if(j.ch()==',')j.eat(',');else break;}j.eat('}');items.push_back(a);if(j.ch()==',')j.eat(',');else break;}j.eat(']'); }
        if(j.ch()==',')j.eat(',');else break;
    }
    vector<Placement> ans; __int128 best=-1; clock_t deadline=clock()+CLOCKS_PER_SEC*55/100;
    for(int order=0;order<3&&clock()<deadline;++order) for(int fit=0;fit<3&&clock()<deadline;++fit){
        auto p=build(items,W,H,allow,order,fit,deadline); __int128 value=0; for(auto&a:p)value+=items[a.t].v;
        if(value>best) best=value,ans.swap(p);
    }
    cout<<"{\"placements\":["; for(size_t i=0;i<ans.size();++i){if(i)cout<<',';auto&p=ans[i];cout<<"{\"type\":\""<<items[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';} cout<<"]}";
}
