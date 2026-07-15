#include <bits/stdc++.h>
using namespace std;

struct Item { string id; long long w,h,v,lim; };
struct Node { long long l,r,y; };
struct Pl { int t; long long x,y; int rot; };
struct Parser {
    string s; size_t p=0;
    Parser(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string z; if(p>=s.size()||s[p++]!='"') return z; while(p<s.size()&&s[p]!='"'){ char c=s[p++]; if(c=='\\'&&p<s.size()){ char e=s[p++]; if(e=='n')z+='\n'; else if(e=='t')z+='\t'; else z+=e; } else z+=c; } if(p<s.size())++p; return z; }
    long long num(){ ws(); int sg=1; if(s[p]=='-'){sg=-1;++p;} long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sg*x; }
    bool boolean(){ ws(); bool q=s.compare(p,4,"true")==0; p+=q?4:5; return q; }
};
static long long W,H; static bool rotateOK; static vector<Item> a;

static long long heightAt(const vector<Node>& q,long long x,long long w){
    long long z=0, e=x+w;
    for(const auto &n:q) if(n.r>x && n.l<e) z=max(z,n.y);
    return z;
}
static void put(vector<Node>& q,long long x,long long w,long long ny){
    long long e=x+w; vector<Node> b; b.reserve(q.size()+3);
    for(auto n:q){
        if(n.r<=x || n.l>=e) b.push_back(n);
        else { if(n.l<x) b.push_back({n.l,x,n.y}); if(n.r>e) b.push_back({e,n.r,n.y}); }
    }
    b.push_back({x,e,ny}); sort(b.begin(),b.end(),[](Node A,Node B){return A.l<B.l;});
    q.clear(); for(auto n:b){ if(n.l>=n.r) continue; if(!q.empty()&&q.back().r==n.l&&q.back().y==n.y)q.back().r=n.r; else q.push_back(n); }
}
// Return a bottom-left position. Aligning either side to a skyline edge covers useful corners.
static bool bestPos(const vector<Node>& q,long long w,long long h,long long &bx,long long &by){
    if(w>W||h>H) return false; vector<long long> xs; xs.reserve(q.size()*2+1);
    xs.push_back(0); xs.push_back(W-w);
    for(auto n:q){ if(n.l<=W-w) xs.push_back(n.l); if(n.r>=w) xs.push_back(n.r-w); }
    sort(xs.begin(),xs.end()); xs.erase(unique(xs.begin(),xs.end()),xs.end()); bool ok=false;
    for(long long x:xs) if(x>=0&&x+w<=W){ long long y=heightAt(q,x,w); if(y+h<=H && (!ok||y<by||(y==by&&x<bx))) bx=x,by=y,ok=true; }
    return ok;
}
static vector<Pl> pack(int mode, int orientationMode){
    vector<Node> sky{{0,W,0}}; vector<long long> used(a.size()); vector<Pl> out;
    // Static priorities give deliberately different useful packings; geometry breaks ties.
    auto pri=[&](int i)->long double{
        long double ar=(long double)a[i].w*a[i].h;
        if(mode==0) return a[i].v/ar;                 // value density
        if(mode==1) return a[i].v;                    // expensive items
        if(mode==2) return a[i].v/sqrt(ar);           // compact high value
        if(mode==3) return a[i].v/(max(a[i].w,a[i].h));
        if(mode==4) return a[i].v/(min(a[i].w,a[i].h));
        return a[i].v/(ar+50.0L*mode);                // density variants
    };
    while(true){
        int take=-1,tr=0; long long tx=0,ty=0; long double bp=-1;
        for(int i=0;i<(int)a.size();++i) if(used[i]<a[i].lim){
            for(int r=0;r<=(rotateOK&&a[i].w!=a[i].h);++r){
                long long ww=r?a[i].h:a[i].w, hh=r?a[i].w:a[i].h, x,y;
                if(!bestPos(sky,ww,hh,x,y)) continue;
                long double z=pri(i);
                // Favor a low skyline only when priorities are essentially tied.
                bool better = take<0 || z>bp+1e-15L || (fabsl(z-bp)<1e-15L && (y<ty || (y==ty && (orientationMode==0 ? hh<a[take].h : ww<a[take].w))));
                if(better) take=i,tr=r,tx=x,ty=y,bp=z;
            }
        }
        if(take<0) break;
        long long ww=tr?a[take].h:a[take].w, hh=tr?a[take].w:a[take].h;
        put(sky,tx,ww,ty+hh); used[take]++; out.push_back({take,tx,ty,tr});
        // A pathological input can have huge limits; output itself cannot usefully be enormous.
        if(out.size()>200000) break;
    }
    return out;
}
// Unlike a skyline, this maintains every maximal empty rectangle, so it can use
// pockets below a previously raised skyline segment.
struct Free { long long x,y,w,h; };
static bool contains(const Free& A,const Free& B){ return A.x<=B.x&&A.y<=B.y&&A.x+A.w>=B.x+B.w&&A.y+A.h>=B.y+B.h; }
static vector<Pl> maxrectPack(int mode){
    vector<Free> fr{{0,0,W,H}}; vector<long long> used(a.size()); vector<Pl> out;
    auto priority=[&](int i)->long double{
        long double ar=(long double)a[i].w*a[i].h;
        if(mode==0) return a[i].v/ar;
        if(mode==1) return a[i].v/sqrt(ar);
        if(mode==2) return a[i].v;
        if(mode==3) return a[i].v/(long double)max(a[i].w,a[i].h);
        return a[i].v/(ar+200.0L);
    };
    // The cap avoids pathological quadratic free-list work.  Normal hard-set
    // packings are much smaller, and skyline candidates remain a fallback.
    while(out.size()<3500 && !fr.empty()){
        int ti=-1, rr=0, fi=-1; long double bp=-1; long long best1=LLONG_MAX,best2=LLONG_MAX;
        for(int i=0;i<(int)a.size();++i) if(used[i]<a[i].lim) for(int r=0;r<=(rotateOK&&a[i].w!=a[i].h);++r){
            long long w=r?a[i].h:a[i].w, h=r?a[i].w:a[i].h;
            for(int j=0;j<(int)fr.size();++j) if(w<=fr[j].w&&h<=fr[j].h){
                long long s=min(fr[j].w-w,fr[j].h-h), l=max(fr[j].w-w,fr[j].h-h);
                long double q=priority(i);
                if(ti<0 || q>bp+1e-15L || (fabsl(q-bp)<1e-15L && (s<best1 || (s==best1&&l<best2))))
                    ti=i,rr=r,fi=j,bp=q,best1=s,best2=l;
            }
        }
        if(ti<0) break;
        Free u=fr[fi]; long long w=rr?a[ti].h:a[ti].w, h=rr?a[ti].w:a[ti].h;
        long long px=u.x,py=u.y; out.push_back({ti,px,py,rr}); used[ti]++;
        vector<Free> nf; nf.reserve(fr.size()+4);
        for(const Free& q:fr){
            long long qx2=q.x+q.w,qy2=q.y+q.h, px2=px+w,py2=py+h;
            if(px2<=q.x||px>=qx2||py2<=q.y||py>=qy2) { nf.push_back(q); continue; }
            if(px>q.x) nf.push_back({q.x,q.y,px-q.x,q.h});
            if(px2<qx2) nf.push_back({px2,q.y,qx2-px2,q.h});
            if(py>q.y) nf.push_back({q.x,q.y,q.w,py-q.y});
            if(py2<qy2) nf.push_back({q.x,py2,q.w,qy2-py2});
        }
        fr.clear();
        for(int i=0;i<(int)nf.size();++i) if(nf[i].w>0&&nf[i].h>0){
            bool dead=false;
            for(int j=0;j<(int)nf.size();++j) if(i!=j&&contains(nf[j],nf[i])) { dead=true; break; }
            if(!dead) fr.push_back(nf[i]);
        }
        // Fragmentation beyond this size is not competitive with the skyline
        // branch and makes the exhaustive containment pruning needlessly slow.
        if(fr.size()>2500) break;
    }
    return out;
}
static long long value(const vector<Pl>& p){ long long z=0; for(auto x:p) z+=a[x.t].v; return z; }
static string esc(const string& x){ string z; for(char c:x){ if(c=='"'||c=='\\')z+='\\'; z+=c; } return z; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr); Parser p; p.ch('{');
    for(int top=0;top<2;top++){
        string k=p.str(); p.ch(':');
        if(k=="bin"){
            p.ch('{'); for(int j=0;j<3;j++){ string q=p.str(); p.ch(':'); if(q=="W")W=p.num(); else if(q=="H")H=p.num(); else rotateOK=p.boolean(); if(j<2)p.ch(','); } p.ch('}');
        } else {
            p.ch('['); p.ws(); while(p.p<p.s.size()&&p.s[p.p]!=']'){
                p.ch('{'); Item it; for(int j=0;j<5;j++){ string q=p.str();p.ch(':'); if(q=="type")it.id=p.str(); else if(q=="w")it.w=p.num(); else if(q=="h")it.h=p.num(); else if(q=="v")it.v=p.num(); else it.lim=p.num(); if(j<4)p.ch(','); } p.ch('}'); a.push_back(it); p.ws(); if(p.s[p.p]==',')p.ch(','); p.ws();
            } p.ch(']');
        } if(top==0)p.ch(',');
    }
    vector<Pl> ans; long long bv=-1;
    // Retain the inexpensive skyline family, then explicitly test the
    // free-rectangle representation against it.
    for(int m=0;m<7;m++) for(int o=0;o<(rotateOK?2:1);o++){ auto q=pack(m,o); long long v=value(q); if(v>bv) bv=v,ans.swap(q); }
    for(int m=0;m<5;m++){ auto q=maxrectPack(m); long long v=value(q); if(v>bv) bv=v,ans.swap(q); }
    cout << "{\"placements\":[";
    for(size_t i=0;i<ans.size();++i){ if(i)cout<<','; auto z=ans[i]; cout<<"{\"type\":\""<<esc(a[z.t].id)<<"\",\"x\":"<<z.x<<",\"y\":"<<z.y<<",\"rot\":"<<z.rot<<'}'; }
    cout << "]}";
}
