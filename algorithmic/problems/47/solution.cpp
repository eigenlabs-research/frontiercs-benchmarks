#include <bits/stdc++.h>
using namespace std;

struct Item { string id; long long w,h,v,lim; };
struct Tok {
    string s; size_t p=0;
    Tok(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){
        ws(); if(p>=s.size()) return "";
        if(s[p]!='"'){ size_t q=p; while(p<s.size() && (isalnum((unsigned char)s[p])||s[p]=='-'||s[p]=='.')) ++p; return s.substr(q,p-q); }
        ++p; string r;
        while(p<s.size() && s[p]!='"') { if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ return stoll(str()); }
    void colon(){ ch(); }
};
struct Seg { long long l,r,y; };
struct Pl { int t; long long x,y; int rot; };
struct Ans { long long val=0; vector<Pl> p; };

int W,H; bool canrot; vector<Item> a;
static chrono::steady_clock::time_point began;
static const int BUDGET=930;
bool expired(){ return chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-began).count()>=BUDGET; }

// A contour (skyline) is a compact representation of the already occupied lower envelope.
// Raising an interval is conservative: it can never create an overlap.
void raiseSky(vector<Seg>& q,long long x,long long w,long long top){
    long long z=x+w; vector<Seg> n; n.reserve(q.size()+2);
    for(auto s:q){
        if(s.r<=x || s.l>=z) n.push_back(s);
        else {
            if(s.l<x) n.push_back({s.l,x,s.y});
            if(s.r>z) n.push_back({z,s.r,s.y});
        }
    }
    n.push_back({x,z,top});
    sort(n.begin(),n.end(),[](const Seg&A,const Seg&B){return A.l<B.l;});
    q.clear();
    for(auto s:n) {
        if(!q.empty() && q.back().r==s.l && q.back().y==s.y) q.back().r=s.r;
        else q.push_back(s);
    }
}

Ans build(int trial){
    vector<Seg> sky(1,{0,W,0});
    vector<long long> used(a.size()); Ans out;
    uint64_t seed=0x9e3779b97f4a7c15ULL ^ (uint64_t)(trial+17)*0xbf58476d1ce4e5b9ULL;
    auto rnd=[&](){ seed^=seed<<7; seed^=seed>>9; return (double)(seed&0xffff)/65535.0; };
    // Vary the exponent: alpha=1 is density, lower values deliberately favor large profit items.
    double alpha = 0.48 + (trial%17)*0.075;
    vector<double> bias(a.size());
    for(size_t i=0;i<a.size();++i) bias[i]=(rnd()-.5)*(trial%3==0?.55:.20);
    for(int step=0; !expired(); ++step){
        int bt=-1, br=0; long long bx=0,by=0, bw=0,bh=0;
        double best=-1e100;
        for(int i=0;i<(int)a.size();++i) if(used[i]<a[i].lim) {
            for(int r=0;r<=1;r++) {
                if(r && (!canrot || a[i].w==a[i].h)) continue;
                long long w=r?a[i].h:a[i].w, h=r?a[i].w:a[i].h;
                if(w>W||h>H) continue;
                long long lx=-1,ly=LLONG_MAX, ltop=LLONG_MAX;
                // Candidate left edges are contour boundaries plus right alignment to the bin.
                // At a contour boundary the intersecting segments begin at k, so this is
                // linear in the contour rather than a quadratic free-list scan.
                for(int k=0;k<(int)sky.size();++k){
                    long long x=sky[k].l;
                    if(x+w<=W) {
                        long long y=0;
                        for(int j=k;j<(int)sky.size() && sky[j].l<x+w;++j) y=max(y,sky[j].y);
                        if(y+h<=H && (y<ly || (y==ly && y+h<ltop) || (y==ly && y+h==ltop && x<lx))) lx=x,ly=y,ltop=y+h;
                    }
                }
                // Right alignment is the only useful non-boundary candidate.
                long long x=W-w, y=0;
                if(x>=0) {
                    for(auto s:sky) if(s.l<x+w && s.r>x) y=max(y,s.y);
                    if(y+h<=H && (y<ly || (y==ly && y+h<ltop) || (y==ly && y+h==ltop && x<lx))) lx=x,ly=y,ltop=y+h;
                }
                if(lx<0) continue;
                double pr=log((double)a[i].v+1.0)-alpha*log((double)w*h+1.0)+bias[i];
                // Slightly prefer placements that do not grow a very tall skyline prematurely.
                pr -= 0.035*(double)(ly+h)/max(1,H);
                if(pr>best+1e-12 || (fabs(pr-best)<1e-12 && a[i].v> (bt<0?0:a[bt].v))){
                    best=pr; bt=i;br=r;bx=lx;by=ly;bw=w;bh=h;
                }
            }
        }
        if(bt<0) break;
        out.p.push_back({bt,bx,by,br}); out.val+=a[bt].v; ++used[bt];
        raiseSky(sky,bx,bw,by+bh);
        if(out.p.size()>200000) break;
    }
    return out;
}

int main(){
    Tok t; if(t.s.empty()) return 0; t.ch();
    while(true){
        string key=t.str(); t.colon();
        if(key=="bin"){
            t.ch(); for(int z=0;z<3;z++){ string k=t.str();t.colon();string v=t.str(); if(k=="W")W=stoi(v); else if(k=="H")H=stoi(v); else canrot=(v=="true"); if(z<2)t.ch(); } t.ch();
        } else if(key=="items"){
            t.ch(); t.ws(); if(t.p<t.s.size()&&t.s[t.p]==']') ++t.p;
            else while(true){
                t.ch(); Item it; for(int z=0;z<5;z++){ string k=t.str();t.colon(); string v=t.str(); if(k=="type")it.id=v; else if(k=="w")it.w=stoll(v); else if(k=="h")it.h=stoll(v); else if(k=="v")it.v=stoll(v); else if(k=="limit")it.lim=stoll(v); if(z<4)t.ch(); } t.ch(); a.push_back(it); t.ws(); char c=t.ch(); if(c==']')break; }
        }
        t.ws(); if(t.p>=t.s.size()||t.ch()=='}') break;
    }
    began=chrono::steady_clock::now(); Ans best;
    for(int trial=0; !expired() && trial<200; ++trial){ Ans cur=build(trial); if(cur.val>best.val) best=move(cur); }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<best.p.size();++i){ auto z=best.p[i]; if(i)cout<<','; cout<<"{\"type\":\""<<a[z.t].id<<"\",\"x\":"<<z.x<<",\"y\":"<<z.y<<",\"rot\":"<<z.rot<<'}'; }
    cout<<"]}";
}
