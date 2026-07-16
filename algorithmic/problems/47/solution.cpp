#include <bits/stdc++.h>
using namespace std;

struct Item { string id; long long w,h,v,lim; };
struct Pl { int t; long long x,y,w,h; int r; };
struct FR { long long x,y,w,h; };

// Small JSON reader for the specified input grammar.  Object member order is deliberately ignored.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char get(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){
        ws(); string r; if(p>=s.size() || s[p++]!='"') return r;
        while(p<s.size() && s[p]!='"') { if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); int sg=1; if(s[p]=='-'){sg=-1;++p;} long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sg*x; }
    bool boolean(){ ws(); bool z=p<s.size()&&s[p]=='t'; while(p<s.size()&&isalpha((unsigned char)s[p]))++p; return z; }
};

struct Solver {
    long long W=0,H=0; bool rot=false; vector<Item> a;
    vector<Pl> best; long long bestV=-1;

    static bool inter(const FR& f, const Pl& q) {
        return f.x < q.x+q.w && q.x < f.x+f.w && f.y < q.y+q.h && q.y < f.y+f.h;
    }
    void trim(vector<FR>& f) {
        vector<FR> g; g.reserve(f.size());
        for(auto z:f) if(z.w>0 && z.h>0) g.push_back(z);
        // Retaining only the largest empty rectangles bounds cleanup cost; discarded rectangles
        // are merely unused feasible space.
        if(g.size()>90) {
            nth_element(g.begin(), g.begin()+90, g.end(), [](const FR&x,const FR&y){return x.w*x.h>y.w*y.h;});
            g.resize(90);
        }
        vector<char> gone(g.size());
        for(int i=0;i<(int)g.size();++i) if(!gone[i]) for(int j=0;j<(int)g.size();++j) if(i!=j && !gone[i]) {
            if(g[j].x<=g[i].x && g[j].y<=g[i].y && g[j].x+g[j].w>=g[i].x+g[i].w && g[j].y+g[j].h>=g[i].y+g[i].h && (g[j].x!=g[i].x || g[j].y!=g[i].y || g[j].w!=g[i].w || g[j].h!=g[i].h || j<i)) gone[i]=1;
        }
        f.clear(); for(int i=0;i<(int)g.size();++i) if(!gone[i]) f.push_back(g[i]);
        // This cap only discards unused space, never creates an invalid placement.  It bounds
        // the quadratic maximal-rectangle cleanup on instances with many tiny items.
        if(f.size()>90) {
            nth_element(f.begin(), f.begin()+90, f.end(), [](const FR&x,const FR&y){return x.w*x.h>y.w*y.h;});
            f.resize(90);
        }
    }
    void insert(vector<FR>& f, const Pl& q) {
        vector<FR> n; n.reserve(f.size()*3);
        for(const FR& z:f) {
            if(!inter(z,q)) { n.push_back(z); continue; }
            long long rx=z.x+z.w, ty=z.y+z.h, qx=q.x+q.w, qy=q.y+q.h;
            if(q.x>z.x) n.push_back({z.x,z.y,q.x-z.x,z.h});
            if(qx<rx) n.push_back({qx,z.y,rx-qx,z.h});
            if(q.y>z.y) n.push_back({z.x,z.y,z.w,q.y-z.y});
            if(qy<ty) n.push_back({z.x,qy,z.w,ty-qy});
        }
        f.swap(n); trim(f);
    }
    long double priority(const Item& z, int mode) {
        long double area=(long double)z.w*z.h, d=z.v/area;
        if(mode==0) return d;
        if(mode==1) return z.v/sqrt(area);           // absolute value gets more weight
        if(mode==2) return z.v/powl(area,0.75L);     // middle ground between value and density
        if(mode==3) return d/sqrt(area);             // strongly favors compact gap fillers
        long double aspect=max(z.w,z.h)/(long double)min(z.w,z.h);
        return d * (1.0L + 0.08L*min((long double)5,aspect)); // gives useful strips a chance
    }
    void run(int mode) {
        vector<FR> freeR(1,{0,0,W,H}); vector<long long> used(a.size()); vector<Pl> out;
        long long val=0;
        // Global choice of both type and maximal empty rectangle, rather than an item-order shelf.
        while(!freeR.empty() && out.size()<2000) {
            int bi=-1,bf=-1,br=0; long double bs=-1;
            for(int i=0;i<(int)a.size();++i) if(used[i]<a[i].lim) {
                for(int r=0;r<=(rot?1:0);++r) {
                    long long iw=r?a[i].h:a[i].w, ih=r?a[i].w:a[i].h;
                    if(r && a[i].w==a[i].h) continue;
                    long double p=priority(a[i],mode);
                    for(int j=0;j<(int)freeR.size();++j) {
                        const FR& z=freeR[j]; if(iw>z.w||ih>z.h) continue;
                        long double waste=(long double)(z.w*z.h-iw*ih)/(z.w*z.h);
                        long double shape=(long double)min(z.w-iw,z.h-ih)/(max(z.w, z.h)+1);
                        // Fit influences positions without allowing low-density fillers to dominate.
                        long double sc=p*(1.0L-0.035L*waste-0.010L*shape);
                        if(sc>bs) bs=sc,bi=i,bf=j,br=r;
                    }
                }
            }
            if(bi<0) break;
            long long iw=br?a[bi].h:a[bi].w, ih=br?a[bi].w:a[bi].h;
            Pl q{bi,freeR[bf].x,freeR[bf].y,iw,ih,br};
            out.push_back(q); used[bi]++; val+=a[bi].v; insert(freeR,q);
        }
        if(val>bestV) bestV=val,best.swap(out);
    }
};

int main(){
    Json j; Solver S; j.get();
    for(int top=0;top<2;++top){
        string key=j.str(); j.get();
        if(key=="bin"){
            j.get(); for(int k=0;k<3;++k){ string q=j.str(); j.get(); if(q=="W")S.W=j.num(); else if(q=="H")S.H=j.num(); else S.rot=j.boolean(); if(k<2)j.get(); } j.get();
        } else {
            j.get(); j.ws(); bool first=true;
            while(j.p<j.s.size() && j.s[j.p]!=']'){
                if(!first) j.get(); first=false; j.get(); Item z;
                for(int k=0;k<5;++k){ string q=j.str(); j.get(); if(q=="type")z.id=j.str(); else if(q=="w")z.w=j.num(); else if(q=="h")z.h=j.num(); else if(q=="v")z.v=j.num(); else z.lim=j.num(); if(k<4)j.get(); } j.get(); S.a.push_back(z); j.ws();
            } j.get();
        }
        if(top==0) j.get();
    }
    for(int mode=0;mode<4;++mode) S.run(mode);
    auto emitString=[](const string& x){ for(char c:x) { if(c=='"' || c=='\\') cout << '\\'; cout << c; } };
    cout << "{\"placements\":[";
    for(size_t k=0;k<S.best.size();++k){ auto &p=S.best[k]; if(k) cout<<','; cout<<"{\"type\":\""; emitString(S.a[p.t].id); cout<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<'}'; }
    cout << "]}";
}
