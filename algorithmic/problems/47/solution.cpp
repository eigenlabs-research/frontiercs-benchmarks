#include <bits/stdc++.h>
using namespace std;
struct Item { string id; int w,h,lim; long long v; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot; };
struct Ans { long long value=0; vector<P> p; };

static string input;
static long long numberAfter(const string& s, const string& key) {
    size_t q=s.find("\""+key+"\""); if(q==string::npos) return 0;
    q=s.find(':',q); ++q; while(q<s.size() && isspace((unsigned char)s[q])) ++q;
    long long z=0; while(q<s.size() && isdigit((unsigned char)s[q])) z=z*10+s[q++]-'0'; return z;
}
static string stringAfter(const string& s, const string& key) {
    size_t q=s.find("\""+key+"\""); q=s.find(':',q); q=s.find('"',q); size_t e=s.find('"',q+1); return s.substr(q+1,e-q-1);
}
static bool intersects(const R&a,const R&b) { return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }
static bool contains(const R&a,const R&b) { return b.x>=a.x&&b.y>=a.y&&b.x+b.w<=a.x+a.w&&b.y+b.h<=a.y+a.h; }

// Split every intersected maximal empty rectangle, then remove contained remnants.
static void insertRect(vector<R>& f, R u) {
    vector<R> n; n.reserve(f.size()*2+4);
    for(const R&r:f) {
        if(!intersects(r,u)) { n.push_back(r); continue; }
        if(u.x>r.x) n.push_back({r.x,r.y,u.x-r.x,r.h});
        if(u.x+u.w<r.x+r.w) n.push_back({u.x+u.w,r.y,r.x+r.w-u.x-u.w,r.h});
        if(u.y>r.y) n.push_back({r.x,r.y,r.w,u.y-r.y});
        if(u.y+u.h<r.y+r.h) n.push_back({r.x,u.y+u.h,r.w,r.y+r.h-u.y-u.h});
    }
    vector<char> bad(n.size());
    for(size_t i=0;i<n.size();++i) if(n[i].w<=0||n[i].h<=0) bad[i]=1;
    for(size_t i=0;i<n.size();++i) if(!bad[i]) for(size_t j=0;j<n.size();++j)
        if(i!=j&&!bad[j]&&contains(n[j],n[i]) && (n[j].x!=n[i].x||n[j].y!=n[i].y||n[j].w!=n[i].w||n[j].h!=n[i].h||j<i)) { bad[i]=1; break; }
    f.clear(); f.reserve(n.size()); for(size_t i=0;i<n.size();++i) if(!bad[i]) f.push_back(n[i]);
}
struct Fit { int fi=-1,rot=0; long long a=LLONG_MAX,b=LLONG_MAX,c=LLONG_MAX; };
static Fit bestFit(const vector<R>& f, const Item& it, bool rotate, int style) {
    Fit z;
    for(int k=0;k<(rotate?2:1);++k) {
        int w=k?it.h:it.w, h=k?it.w:it.h;
        for(int i=0;i<(int)f.size();++i) if(w<=f[i].w&&h<=f[i].h) {
            long long dx=f[i].w-w, dy=f[i].h-h;
            long long A,B,C;
            if(style==0) A=min(dx,dy),B=max(dx,dy),C=dx*dy; // short side fit
            else if(style==1) A=dx*dy,B=min(dx,dy),C=max(dx,dy); // area fit
            else A=f[i].y+h,B=f[i].x+w,C=min(dx,dy); // bottom-left
            if(make_tuple(A,B,C,k)<make_tuple(z.a,z.b,z.c,z.rot)) z={i,k,A,B,C};
        }
    }
    return z;
}
static Ans pack(vector<Item>& it,int W,int H,bool canrot,int orderMode,int fitStyle) {
    int m=it.size(); vector<int> ord(m); iota(ord.begin(),ord.end(),0);
    auto density=[&](int i){return (long double)it[i].v/(it[i].w*1.0L*it[i].h);};
    sort(ord.begin(),ord.end(),[&](int a,int b){
        long double ka,kb;
        if(orderMode==0) ka=density(a),kb=density(b);
        else if(orderMode==1) ka=it[a].v,kb=it[b].v;
        else if(orderMode==2) ka=density(a)*sqrt((long double)it[a].w*it[a].h),kb=density(b)*sqrt((long double)it[b].w*it[b].h);
        else ka=(long double)it[a].v/max(it[a].w,it[a].h),kb=(long double)it[b].v/max(it[b].w,it[b].h);
        if(ka!=kb) return ka>kb; return it[a].v>it[b].v;
    });
    vector<R> freeR={{0,0,W,H}}; Ans ans;
    // A type is exhausted before moving on: this intentionally tests the density-first packing model.
    for(int t:ord) for(int cnt=0;cnt<it[t].lim;++cnt) {
        Fit q=bestFit(freeR,it[t],canrot,fitStyle); if(q.fi<0) break;
        int w=q.rot?it[t].h:it[t].w, h=q.rot?it[t].w:it[t].h;
        R u={freeR[q.fi].x,freeR[q.fi].y,w,h};
        ans.p.push_back({t,u.x,u.y,q.rot}); ans.value+=it[t].v;
        insertRect(freeR,u);
        if(freeR.empty()) return ans;
    }
    return ans;
}
static string esc(const string&s) { string r; for(char c:s) { if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    input.assign(istreambuf_iterator<char>(cin), istreambuf_iterator<char>());
    int W=(int)numberAfter(input,"W"), H=(int)numberAfter(input,"H");
    size_t b=input.find("\"allow_rotate\""); bool rot=input.find("true",b)!=string::npos && input.find("true",b)<input.find('}',b);
    vector<Item> it; size_t p=input.find("\"items\"");
    while((p=input.find('{',p))!=string::npos) {
        size_t e=input.find('}',p); if(e==string::npos) break; string o=input.substr(p,e-p+1);
        if(o.find("\"type\"")!=string::npos) it.push_back({stringAfter(o,"type"),(int)numberAfter(o,"w"),(int)numberAfter(o,"h"),(int)numberAfter(o,"limit"),numberAfter(o,"v")});
        p=e+1;
    }
    Ans best;
    // Independent packings make choice robust to unknown hidden aspect-ratio distributions.
    for(int order=0;order<4;++order) for(int fit=0;fit<3;++fit) {
        Ans a=pack(it,W,H,rot,order,fit); if(a.value>best.value) best=move(a);
    }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<best.p.size();++i) { auto&q=best.p[i]; if(i) cout<<','; cout<<"{\"type\":\""<<esc(it[q.t].id)<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.rot<<'}'; }
    cout<<"]}\n";
}
