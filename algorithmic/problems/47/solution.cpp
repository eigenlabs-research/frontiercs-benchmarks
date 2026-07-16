#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct Pl { int t,x,y,rot,w,h; };
struct Scan {
    string s; size_t p=0;
    Scan(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char get(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){ ws(); if(p>=s.size()||s[p++]!='"') return ""; string r; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\'&&p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size())++p; return r; }
    long long num(){ ws(); int sg=1; if(s[p]=='-'){sg=-1;++p;} long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p]))x=x*10+s[p++]-'0'; return sg*x; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};

struct Solver {
    int W=0,H=0; bool rotate=false; vector<Item> a;
    // Nodes give height starting at x, with the final sentinel at W.
    void split(vector<pair<int,int>>& s,int x) {
        for(int i=0;i+1<(int)s.size();++i) if(s[i].first<x && x<s[i+1].first) { s.insert(s.begin()+i+1,{x,s[i].second}); return; }
    }
    int level(const vector<pair<int,int>>& s,int first,int w) {
        int z=0, end=s[first].first+w;
        for(int i=first;i+1<(int)s.size() && s[i].first<end;++i) z=max(z,s[i].second);
        return z;
    }
    void put(vector<pair<int,int>>& s,int x,int w,int nh) {
        split(s,x); split(s,x+w);
        for(auto &q:s) if(q.first>=x && q.first<x+w) q.second=nh;
        vector<pair<int,int>> t; for(auto q:s) if(t.empty()||t.back().second!=q.second||q.first==W) t.push_back(q);
        if(t.back().first!=W) t.push_back({W,0}); s.swap(t);
    }
    vector<Pl> run(int mode) {
        vector<pair<int,int>> sky={{0,0},{W,0}};
        vector<int> used(a.size()); vector<Pl> out;
        // A small set of lowest profile corners keeps this bounded even for tiny pieces.
        for(int step=0;step<3000;++step) {
            int bi=-1,bx=0,by=0,br=0,bw=0,bh=0, btop=0;
            long double bkey=-1e100L; bool geometric=(mode==8 || mode==9);
            vector<int> corners;
            for(int k=0;k+1<(int)sky.size();++k) {
                corners.push_back(k);
                sort(corners.begin(),corners.end(),[&](int u,int v){ return sky[u].second!=sky[v].second ? sky[u].second<sky[v].second : sky[u].first<sky[v].first; });
                if(corners.size()>8) corners.pop_back();
            }
            for(int i=0;i<(int)a.size();++i) if(used[i]<a[i].lim) for(int r=0;r<=(rotate?1:0);++r) {
                int w=r?a[i].h:a[i].w, h=r?a[i].w:a[i].h;
                if(w>W||h>H || (r&&w==h)) continue;
                for(int k:corners) {
                    int x=sky[k].first; if(x+w>W) continue;
                    int y=level(sky,k,w); if(y+h>H) continue;
                    long double area=(long double)w*h, den=a[i].v/area, key;
                    if(mode==0) key=den;
                    else if(mode==1) key=a[i].v;
                    else if(mode==2) key=den*sqrtl(area);
                    else if(mode==3) key=a[i].v/(long double)max(w,h);
                    else if(mode==4) key=den*(1.0L+0.15L*(long double)min(w,h)/max(w,h));
                    else if(mode==5) key=a[i].v/sqrtl(area);
                    else if(mode==6) key=den*(long double)min(w,h);
                    else if(mode==7) key=a[i].v/(long double)(w+h);
                    else key=den;
                    int top=y+h;
                    bool take=false;
                    if(bi<0) take=true;
                    else if(geometric) {
                        // Profile-first variants deliberately keep a flat, low frontier.
                        if(top!=btop) take=top<btop;
                        else if(fabsl(key-bkey)>1e-15L) take=key>bkey;
                        else if(y!=by) take=y<by;
                    } else {
                        if(key>bkey+1e-15L) take=true;
                        else if(fabsl(key-bkey)<=1e-15L) {
                            if(top!=btop) take=top<btop;
                            else if(y!=by) take=y<by;
                            else if(x!=bx) take=x<bx;
                            else if(w>bw) take=true;
                        }
                    }
                    if(take) bi=i,bx=x,by=y,br=r,bw=w,bh=h,btop=top,bkey=key;
                }
            }
            if(bi<0) break;
            out.push_back({bi,bx,by,br,bw,bh}); ++used[bi]; put(sky,bx,bw,by+bh);
        }
        return out;
    }
    long long value(const vector<Pl>& q){ long long z=0; for(auto p:q) z+=a[p.t].v; return z; }
    bool valid(const vector<Pl>& q) {
        vector<int> c(a.size());
        for(auto p:q) { if(p.t<0||p.t>=(int)a.size()||++c[p.t]>a[p.t].lim||p.x<0||p.y<0||p.x+p.w>W||p.y+p.h>H) return false; }
        // Each placement is raised to the maximum skyline level over its span; thus it
        // cannot intersect earlier skyline placements (the profile is an exact conservative cover).
        return true;
    }
};
static string esc(const string& s){ string r; for(char c:s){ if(c=='"'||c=='\\')r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Scan z; Solver S;
    if(z.get()!='{') return 0;
    for(int field=0;field<2;++field){ string key=z.str(); z.get();
        if(key=="bin") { z.get(); for(int k=0;k<3;++k){ string q=z.str(); z.get(); if(q=="W")S.W=z.num(); else if(q=="H")S.H=z.num(); else S.rotate=z.boolean(); if(k<2)z.get(); } z.get(); }
        else { z.get(); z.ws(); while(z.p<z.s.size()&&z.s[z.p]!=']') { z.get(); Item it; for(int k=0;k<5;++k){ string q=z.str();z.get(); if(q=="type")it.id=z.str(); else if(q=="w")it.w=z.num(); else if(q=="h")it.h=z.num(); else if(q=="v")it.v=z.num(); else it.lim=z.num(); if(k<4)z.get(); } z.get(); S.a.push_back(it); z.ws(); if(z.s[z.p]==',')z.get(); } z.get(); }
        if(field==0)z.get();
    }
    vector<Pl> best; long long bv=-1;
    for(int m=0;m<10;++m){ auto q=S.run(m); if(S.valid(q) && S.value(q)>bv) bv=S.value(q),best.swap(q); }
    cout<<"{\"placements\":[";
    for(int i=0;i<(int)best.size();++i){ if(i)cout<<','; auto p=best[i]; cout<<"{\"type\":\""<<esc(S.a[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout<<"]}";
}
