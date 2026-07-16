#include <bits/stdc++.h>
using namespace std;

struct Item { string name; int w,h; long long v; int lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };

struct Parser {
    string s; size_t p=0;
    Parser(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char get(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){ ws(); if(p<s.size() && s[p]=='"') ++p; string r; while(p<s.size() && s[p]!='"') { if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size()) ++p; return r; }
    long long num(){ ws(); int z=1; if(s[p]=='-'){z=-1;++p;} long long r=0; while(p<s.size()&&isdigit((unsigned char)s[p])) r=r*10+s[p++]-'0'; return z*r; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};
static bool intersects(const R&a,const R&b) { return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }
static bool contains(const R&a,const R&b) { return b.x>=a.x&&b.y>=a.y&&b.x+b.w<=a.x+a.w&&b.y+b.h<=a.y+a.h; }

struct Result { vector<P> p; long long val=0; };

Result pack(const vector<Item>& it, int W, int H, bool allow, int mode) {
    vector<R> free{{0,0,W,H}};
    vector<int> used(it.size()); Result ans;
    // mode chooses how strongly a placement is rewarded for consuming its current free rectangle.
    const int penalties[] = {0, 1, 4, 10, 20};
    int pen=penalties[mode%5], q=20;
    while(true) {
        int bi=-1,bf=-1,br=0; long long bn=0,bd=1; int bestA=INT_MAX, bestShort=INT_MAX;
        for(int i=0;i<(int)it.size();++i) if(used[i]<it[i].lim) {
            for(int rot=0;rot<=(allow?1:0);++rot) {
                int w=rot?it[i].h:it[i].w, h=rot?it[i].w:it[i].h;
                for(int f=0;f<(int)free.size();++f) if(w<=free[f].w && h<=free[f].h) {
                    long long area=1LL*w*h, waste=1LL*free[f].w*free[f].h-area;
                    // maximize v / (q*area + pen*waste); pen=0 is pure density
                    long long den=q*area+pen*waste;
                    int shortside=min(free[f].w-w,free[f].h-h);
                    int longa=max(free[f].w-w,free[f].h-h);
                    bool take=false;
                    if(bi<0 || (__int128)it[i].v*bd > (__int128)bn*den) take=true;
                    else if((__int128)it[i].v*bd == (__int128)bn*den) {
                        // Different deterministic tie breaking produces distinct useful layouts.
                        if(mode==1) take = waste < bestA || (waste==bestA && shortside<bestShort);
                        else if(mode==2) take = shortside < bestShort || (shortside==bestShort && waste<bestA);
                        else if(mode==3) take = waste < bestA || (waste==bestA && longa<bestShort);
                        else take = waste < bestA || (waste==bestA && i<bi);
                    }
                    if(take) { bi=i; bf=f; br=rot; bn=it[i].v; bd=den; bestA=(int)min<long long>(waste,INT_MAX); bestShort=shortside; }
                }
            }
        }
        if(bi<0) break;
        int w=br?it[bi].h:it[bi].w, h=br?it[bi].w:it[bi].h;
        R placed{free[bf].x,free[bf].y,w,h};
        ans.p.push_back({bi,placed.x,placed.y,br,w,h}); ans.val+=it[bi].v; ++used[bi];
        vector<R> nf;
        for(R a:free) {
            if(!intersects(a,placed)) { nf.push_back(a); continue; }
            if(placed.x>a.x) nf.push_back({a.x,a.y,placed.x-a.x,a.h});
            if(placed.x+placed.w<a.x+a.w) nf.push_back({placed.x+placed.w,a.y,a.x+a.w-(placed.x+placed.w),a.h});
            if(placed.y>a.y) nf.push_back({a.x,a.y,a.w,placed.y-a.y});
            if(placed.y+placed.h<a.y+a.h) nf.push_back({a.x,placed.y+placed.h,a.w,a.y+a.h-(placed.y+placed.h)});
        }
        vector<char> dead(nf.size());
        for(int a=0;a<(int)nf.size();++a) for(int b=0;b<(int)nf.size();++b) if(a!=b && contains(nf[b],nf[a])) { dead[a]=1; break; }
        free.clear(); for(int a=0;a<(int)nf.size();++a) if(!dead[a]&&nf[a].w&&nf[a].h) free.push_back(nf[a]);
    }
    return ans;
}

int main(){
    Parser z; z.get(); int W=0,H=0; bool allow=false; vector<Item> items;
    for(int root=0;root<2;++root) {
        string key=z.str(); z.get();
        if(key=="bin") {
            z.get(); for(int j=0;j<3;++j){ string k=z.str(); z.get(); if(k=="W") W=z.num(); else if(k=="H") H=z.num(); else allow=z.boolean(); if(j<2)z.get(); } z.get();
        } else {
            z.get(); z.ws(); while(z.p<z.s.size() && z.s[z.p]!=']') {
                z.get(); Item a; for(int j=0;j<5;++j) { string k=z.str(); z.get(); if(k=="type")a.name=z.str(); else if(k=="w")a.w=z.num(); else if(k=="h")a.h=z.num(); else if(k=="v")a.v=z.num(); else a.lim=z.num(); if(j<4)z.get(); } z.get(); items.push_back(a); z.ws(); if(z.s[z.p]==',')z.get(); z.ws();
            } z.get();
        }
        z.ws(); if(root==0) z.get();
    }
    Result best;
    for(int m=0;m<3;++m) { Result r=pack(items,W,H,allow,m); if(r.val>best.val) best=move(r); }
    cout << "{\"placements\":[";
    for(size_t k=0;k<best.p.size();++k){ if(k) cout<<','; auto &a=best.p[k]; cout<<"{\"type\":\""<<items[a.t].name<<"\",\"x\":"<<a.x<<",\"y\":"<<a.y<<",\"rot\":"<<a.rot<<'}'; }
    cout << "]}\n";
}
