#include <bits/stdc++.h>
using namespace std;

struct Item { string id; long long w,h,v,lim; };
struct Pl { int t; long long x,y,w,h; int r; };
struct Seg { long long l,r,h; };

// Small JSON reader for the fixed, object-only input schema.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); string z; if(p>=s.size() || s[p++]!='\"') return z;
        while(p<s.size() && s[p]!='\"') { char c=s[p++]; if(c=='\\' && p<s.size()) { char e=s[p++]; if(e=='n') z+='\n'; else if(e=='t') z+='\t'; else z+=e; } else z+=c; }
        if(p<s.size()) ++p; return z;
    }
    long long num(){ ws(); long long sign=1,x=0; if(s[p]=='-') sign=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sign*x; }
    bool boolean(){ ws(); bool b=s.compare(p,4,"true")==0; p+=b?4:5; return b; }
};

static long long W,H; static bool rotateOK; static vector<Item> a;

static void putSky(vector<Seg>& q, long long x, long long w, long long top) {
    vector<Seg> z; z.reserve(q.size()+3);
    long long e=x+w;
    for(auto s:q) {
        if(s.r<=x || s.l>=e) z.push_back(s);
        else { if(s.l<x) z.push_back({s.l,x,s.h}); if(s.r>e) z.push_back({e,s.r,s.h}); }
    }
    z.push_back({x,e,top});
    sort(z.begin(),z.end(),[](const Seg&A,const Seg&B){return A.l<B.l;});
    q.clear();
    for(auto s:z) {
        if(!q.empty() && q.back().r==s.l && q.back().h==s.h) q.back().r=s.r;
        else q.push_back(s);
    }
}

static vector<Pl> pack(int mode) {
    vector<int> ord(a.size()); iota(ord.begin(),ord.end(),0);
    sort(ord.begin(),ord.end(),[&](int i,int j){
        const Item &A=a[i],&B=a[j];
        if(mode==2 && A.v!=B.v) return A.v>B.v;
        __int128 l=(__int128)A.v*B.w*B.h, r=(__int128)B.v*A.w*A.h;
        if(l!=r) return l>r;
        if(mode==0) { if(A.h!=B.h) return A.h>B.h; if(A.w!=B.w) return A.w>B.w; }
        else if(mode==1) { if(A.w!=B.w) return A.w>B.w; if(A.h!=B.h) return A.h>B.h; }
        else { if(A.w*A.h!=B.w*B.h) return A.w*A.h<B.w*B.h; }
        return A.id<B.id;
    });
    vector<Seg> sky={{0,W,0}}; vector<long long> used(a.size()); vector<Pl> ans;
    for(int id:ord) {
        for(long long copy=0; copy<a[id].lim; ++copy) {
            long long bx=-1,by=0,bw=0,bh=0; int br=0;
            // Every skyline left edge is a sufficient candidate x for bottom-left placement.
            for(int rr=0;rr<(rotateOK?2:1);++rr) {
                long long w=rr?a[id].h:a[id].w, h=rr?a[id].w:a[id].h;
                if(w>W || h>H) continue;
                for(const auto &start:sky) {
                    long long x=start.l; if(x+w>W) continue;
                    long long y=0;
                    for(const auto &s:sky) { if(s.l>=x+w) break; if(s.r>x) y=max(y,s.h); }
                    if(y+h>H) continue;
                    // Bottom-left, with a lower resulting skyline as a useful tie-break.
                    if(bx<0 || y<by || (y==by && x<bx) || (y==by && x==bx && y+h<by+bh) ||
                       (y==by && x==bx && y+h==by+bh && rr<br)) {
                        bx=x; by=y; bw=w; bh=h; br=rr;
                    }
                }
            }
            if(bx<0) break; // Skyline heights only rise, so this type cannot become feasible later.
            putSky(sky,bx,bw,by+bh);
            ans.push_back({id,bx,by,bw,bh,br}); ++used[id];
        }
    }
    return ans;
}

static long long value(const vector<Pl>& p){ long long z=0; for(auto x:p) z+=a[x.t].v; return z; }
static string escaped(const string& s) { string r; for(char c:s) { if(c=='\"'||c=='\\') r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Json j; j.ch('{');
    // Keys may be ordered arbitrarily, although official inputs use bin then items.
    for(int root=0;root<2;++root){
        if(root) j.ch(','); string key=j.str(); j.ch(':');
        if(key=="bin") {
            j.ch('{'); for(int k=0;k<3;++k){ if(k)j.ch(','); string n=j.str();j.ch(':'); if(n=="W")W=j.num(); else if(n=="H")H=j.num(); else rotateOK=j.boolean(); } j.ch('}');
        } else {
            j.ch('['); bool first=true; while(true){ j.ws(); if(j.p<j.s.size()&&j.s[j.p]==']'){++j.p;break;} if(!first)j.ch(','); first=false; j.ch('{'); Item it; for(int k=0;k<5;++k){if(k)j.ch(',');string n=j.str();j.ch(':'); if(n=="type")it.id=j.str(); else if(n=="w")it.w=j.num(); else if(n=="h")it.h=j.num(); else if(n=="v")it.v=j.num(); else it.lim=j.num();}j.ch('}');a.push_back(it); }
        }
    }
    vector<Pl> best;
    for(int m=0;m<3;++m){ auto cur=pack(m); if(value(cur)>value(best)) best.swap(cur); }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){ if(i)cout<<','; const auto&p=best[i]; cout<<"{\"type\":\""<<escaped(a[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<'}'; }
    cout << "]}";
}
