#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct P { int t,x,y,r,w,h; };
struct Seg { int l,r,y; };

// Small JSON reader: the input grammar only contains objects, arrays, strings, numbers and booleans.
struct JS {
    string s; size_t p=0;
    JS(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p]:'\0'; }
    void need(char c){ ws(); if(p<s.size()&&s[p]==c) ++p; }
    string str(){ need('"'); string r; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\' && p+1<s.size()){ ++p; char c=s[p++]; if(c=='n')r+='\n'; else if(c=='t')r+='\t'; else r+=c; } else r+=s[p++]; } need('"'); return r; }
    long long num(){ ws(); int sg=1; if(s[p]=='-')sg=-1,++p; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p]))x=x*10+s[p++]-'0'; return x*sg; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4;return true;} p+=5;return false; }
};

struct Pack {
    int W,H; bool rot; const vector<Item>& a; vector<Seg> sky; vector<P> out; long long val=0;
    Pack(int W,int H,bool r,const vector<Item>& a):W(W),H(H),rot(r),a(a){sky.push_back({0,W,0});}
    bool position(int w,int h,int &bx,int &by){
        int bestY=INT_MAX,bestX=INT_MAX;
        for(const auto &q:sky){
            int x=q.l; if(x+w>W) continue;
            int y=0;
            for(const auto &z:sky) if(z.r>x && z.l<x+w) y=max(y,z.y);
            if(y+h<=H && (y<bestY || (y==bestY&&x<bestX))) bestY=y,bestX=x;
        }
        if(bestY==INT_MAX) return false; bx=bestX;by=bestY;return true;
    }
    void addSky(int x,int w,int top){
        vector<Seg> q;
        for(auto z:sky){
            if(z.r<=x || z.l>=x+w) q.push_back(z);
            else { if(z.l<x) q.push_back({z.l,x,z.y}); if(z.r>x+w) q.push_back({x+w,z.r,z.y}); }
        }
        q.push_back({x,x+w,top}); sort(q.begin(),q.end(),[](Seg A,Seg B){return A.l<B.l;});
        sky.clear(); for(auto z:q){ if(!sky.empty()&&sky.back().r==z.l&&sky.back().y==z.y) sky.back().r=z.r; else sky.push_back(z); }
    }
    bool put(int t){
        const Item &it=a[t]; int bx=0,by=0,bw=0,bh=0,br=0; bool ok=false;
        auto test=[&](int w,int h,int r){ int x,y; if(!position(w,h,x,y)) return; if(!ok || y<by || (y==by&&x<bx) || (y==by&&x==bx&&y+h<by+bh)){ok=true;bx=x;by=y;bw=w;bh=h;br=r;} };
        test(it.w,it.h,0); if(rot && it.w!=it.h) test(it.h,it.w,1);
        if(!ok)return false; addSky(bx,bw,by+bh); out.push_back({t,bx,by,br,bw,bh}); val+=it.v; return true;
    }
};

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    JS j; vector<Item> a; int W=0,H=0; bool rotate=false;
    j.need('{');
    while(j.ch()!='}' && j.p<j.s.size()){
        string k=j.str(); j.need(':');
        if(k=="bin"){
            j.need('{'); while(j.ch()!='}') { string z=j.str();j.need(':'); if(z=="W")W=(int)j.num(); else if(z=="H")H=(int)j.num(); else rotate=j.boolean(); if(j.ch()==',')j.need(','); } j.need('}');
        } else if(k=="items") {
            j.need('['); while(j.ch()!=']') { Item it; j.need('{'); while(j.ch()!='}') { string z=j.str();j.need(':'); if(z=="type")it.id=j.str(); else if(z=="w")it.w=j.num(); else if(z=="h")it.h=j.num(); else if(z=="v")it.v=j.num(); else it.lim=j.num(); if(j.ch()==',')j.need(','); } j.need('}'); a.push_back(it); if(j.ch()==',')j.need(','); } j.need(']');
        }
        if(j.ch()==',')j.need(',');
    }
    vector<P> best; long long bestv=-1; int n=a.size();
    vector<double> den(n), vv(n), ar(n); double md=0,mv=0,ma=0;
    for(int i=0;i<n;i++){ den[i]=(double)a[i].v/((double)a[i].w*a[i].h); vv[i]=a[i].v; ar[i]=(double)a[i].w*a[i].h; md=max(md,den[i]);mv=max(mv,vv[i]);ma=max(ma,ar[i]); }
    // A compact deterministic portfolio of priorities, then both batch and interleaved copies.
    vector<array<double,3>> weights={{ {1,0,0},{0,1,0},{.8,.15,.05},{.7,.3,-.1},{.55,.1,.35},{.4,.55,.05},{.9,.05,-.25},{.3,.3,.4},{.65,.5,-.3},{.5,.15,.2} }};
    for(auto q:weights) for(int mode=0;mode<2;mode++){
        vector<int> ord(n); iota(ord.begin(),ord.end(),0);
        sort(ord.begin(),ord.end(),[&](int x,int y){ double sx=q[0]*den[x]/md+q[1]*vv[x]/mv+q[2]*ar[x]/ma, sy=q[0]*den[y]/md+q[1]*vv[y]/mv+q[2]*ar[y]/ma; if(fabs(sx-sy)>1e-12)return sx>sy; return a[x].v>a[y].v; });
        Pack p(W,H,rotate,a); vector<int> used(n); bool progress=true;
        if(!mode) { for(int t:ord) while(used[t]<a[t].lim && p.put(t)) ++used[t]; }
        else { while(progress){ progress=false; for(int t:ord) if(used[t]<a[t].lim && p.put(t)){++used[t];progress=true;} } }
        if(p.val>bestv){bestv=p.val;best=move(p.out);}
    }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<best.size();++i){ if(i)cout<<','; const auto &p=best[i]; cout<<"{\"type\":\""<<a[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<"}"; }
    cout<<"]}\n";
}
