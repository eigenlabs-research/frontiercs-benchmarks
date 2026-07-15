#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct Pl { int t,x,y,r,w,h; };
struct Rect { int x,y,w,h; };

// Tiny JSON reader: input strings in this task have only ordinary identifier keys/values.
struct Parser {
    string s; size_t p=0;
    Parser(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); ch('"'); string r; while(p<s.size() && s[p]!='"') { if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } ch('"'); return r; }
    long long num(){ ws(); int z=1; if(s[p]=='-')z=-1,++p; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p]))x=x*10+s[p++]-'0'; return x*z; }
    bool boolean(){ ws(); bool r=s[p]=='t'; p+=r?4:5; return r; }
};

static bool intersects(const Rect&a,const Rect&b){
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}

struct Packer {
    int W,H; const vector<Item>& a; bool rot;
    vector<Rect> fr; vector<Pl> out; vector<int> used;
    Packer(int W,int H,const vector<Item>&a,bool r):W(W),H(H),a(a),rot(r),used(a.size()) { fr.push_back({0,0,W,H}); }
    void prune(){
        for(int i=0;i<(int)fr.size();++i) for(int j=i+1;j<(int)fr.size();) {
            Rect A=fr[i],B=fr[j];
            if(A.x>=B.x && A.y>=B.y && A.x+A.w<=B.x+B.w && A.y+A.h<=B.y+B.h) { fr.erase(fr.begin()+i); --i; break; }
            if(B.x>=A.x && B.y>=A.y && B.x+B.w<=A.x+A.w && B.y+B.h<=A.y+A.h) fr.erase(fr.begin()+j);
            else ++j;
        }
    }
    void put(int t,int x,int y,int w,int h,int r){
        Rect q{x,y,w,h}; vector<Rect> nf;
        for(Rect f:fr){
            if(!intersects(f,q)){ nf.push_back(f); continue; }
            // The four residual rectangles are intentional: this is the MaxRects split.
            if(q.x>f.x) nf.push_back({f.x,f.y,q.x-f.x,f.h});
            if(q.x+q.w<f.x+f.w) nf.push_back({q.x+q.w,f.y,f.x+f.w-(q.x+q.w),f.h});
            if(q.y>f.y) nf.push_back({f.x,f.y,f.w,q.y-f.y});
            if(q.y+q.h<f.y+f.h) nf.push_back({f.x,q.y+q.h,f.w,f.y+f.h-(q.y+q.h)});
        }
        fr.swap(nf); prune(); out.push_back({t,x,y,r,w,h}); ++used[t];
    }
    // Best short-side fit, with a bottom-left tie break, produces compact packings.
    bool one(int t, int mode){
        if(used[t]>=a[t].lim) return false;
        int bx=0,by=0,bw=0,bh=0,br=0; long long bs=LLONG_MAX, bl=LLONG_MAX;
        for(const Rect& f:fr) for(int rr=0;rr<(rot&&a[t].w!=a[t].h?2:1);++rr){
            int w=rr?a[t].h:a[t].w, h=rr?a[t].w:a[t].h;
            if(w>f.w||h>f.h) continue;
            long long s=min(f.w-w,f.h-h), l=max(f.w-w,f.h-h);
            // alternate trials occasionally prefer area fit, useful for large valuable pieces
            if(mode==1) s=(long long)(f.w-w)*(f.h-h), l=min(f.w-w,f.h-h);
            if(s<bs || (s==bs && (l<bl || (l==bl && (f.y<by || (f.y==by&&f.x<bx))))))
                bs=s,bl=l,bx=f.x,by=f.y,bw=w,bh=h,br=rr;
        }
        if(bs==LLONG_MAX) return false;
        put(t,bx,by,bw,bh,br); return true;
    }
    long long value() const { long long z=0; for(auto&p:out) z+=a[p.t].v; return z; }
};

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Parser q; q.ch('{'); string k=q.str(); q.ch(':'); q.ch('{');
    int W=0,H=0; bool allow=false;
    for(int z=0;z<3;z++){ if(z)q.ch(','); string x=q.str();q.ch(':'); if(x=="W")W=q.num(); else if(x=="H")H=q.num(); else allow=q.boolean(); }
    q.ch('}');q.ch(',');q.str();q.ch(':');q.ch('[');
    vector<Item>a;
    q.ws(); bool first=true;
    while(q.p<q.s.size() && q.s[q.p]!=']'){
        if(!first)q.ch(','); first=false; q.ch('{'); Item e;
        for(int z=0;z<5;z++){ if(z)q.ch(','); string x=q.str();q.ch(':'); if(x=="type")e.id=q.str(); else if(x=="w")e.w=q.num(); else if(x=="h")e.h=q.num(); else if(x=="v")e.v=q.num(); else e.lim=q.num(); }
        q.ch('}'); q.ws(); a.push_back(e);
    }
    vector<Pl> best; long long bestv=-1;
    vector<int> base(a.size()); iota(base.begin(),base.end(),0);
    mt19937 rng(712367);
    auto start=chrono::steady_clock::now(); int trial=0;
    // Several ordering perturbations matter more than a single elaborate placement rule.
    while(trial<55 && chrono::duration<double>(chrono::steady_clock::now()-start).count()<0.72){
        vector<int> ord=base;
        vector<double> key(a.size());
        for(int i=0;i<(int)a.size();i++){
            double d=(double)a[i].v/((double)a[i].w*a[i].h);
            if(trial==0) key[i]=log(max(d,1e-30));
            else { double u=(double)(rng()%1000000)/1000000.0-.5; key[i]=log(max(d,1e-30))+u*(trial<20?.18:.70); }
        }
        sort(ord.begin(),ord.end(),[&](int x,int y){return key[x]>key[y];});
        Packer p(W,H,a,allow);
        // Exhaust a type before moving to the next one; later types then fill the
        // residual strips left by the more profitable types.
        for(int id:ord) while(p.one(id,(trial%7==0)?1:0)) {}
        long long v=p.value(); if(v>bestv){bestv=v;best=p.out;}
        ++trial;
    }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<best.size();++i){ if(i)cout<<','; const Pl&p=best[i];
        cout<<"{\"type\":\""<<a[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<'}'; }
    cout<<"]}\n";
}
