#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct R { int x,y,w,h; };
struct Put { int t,x,y,rot; };

// Small JSON reader: input contains only the primitive JSON values used by this task.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){
        ws(); ++p; string r;
        while(p<s.size() && s[p]!='"') {
            if(s[p]=='\\' && p+1<s.size()) { ++p; char c=s[p++];
                if(c=='n') r+='\n'; else if(c=='t') r+='\t'; else r+=c;
            } else r+=s[p++];
        }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); int sg=1; if(s[p]=='-') sg=-1,++p; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sg*x; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0) {p+=4;return true;} p+=5;return false; }
};

static bool contains(const R&a,const R&b) {
    return b.x>=a.x && b.y>=a.y && b.x+b.w<=a.x+a.w && b.y+b.h<=a.y+a.h;
}
struct Packer {
    int W,H; bool rotate; const vector<Item>& it;
    vector<R> fr; vector<Put> out; long long value=0;
    Packer(int W,int H,bool ro,const vector<Item>& a):W(W),H(H),rotate(ro),it(a) { fr.push_back({0,0,W,H}); }
    void carve(R q) {
        vector<R> nf; nf.reserve(fr.size()*2);
        for(const R& f:fr) {
            int ix=max(f.x,q.x), iy=max(f.y,q.y), ax=min(f.x+f.w,q.x+q.w), ay=min(f.y+f.h,q.y+q.h);
            if(ix>=ax || iy>=ay) { nf.push_back(f); continue; }
            if(q.y>f.y) nf.push_back({f.x,f.y,f.w,q.y-f.y});
            if(q.y+q.h<f.y+f.h) nf.push_back({f.x,q.y+q.h,f.w,f.y+f.h-q.y-q.h});
            if(q.x>f.x) nf.push_back({f.x,f.y,q.x-f.x,f.h});
            if(q.x+q.w<f.x+f.w) nf.push_back({q.x+q.w,f.y,f.x+f.w-q.x-q.w,f.h});
        }
        // Keeping only maximal empty rectangles controls the otherwise rapid growth.
        if(nf.size()>420) {
            nth_element(nf.begin(),nf.begin()+420,nf.end(),[](const R&a,const R&b){return 1LL*a.w*a.h>1LL*b.w*b.h;});
            nf.resize(420);
        }
        vector<char> dead(nf.size());
        for(size_t i=0;i<nf.size();++i) if(!dead[i])
            for(size_t j=0;j<nf.size();++j) if(i!=j && !dead[i] && contains(nf[j],nf[i]) && (!contains(nf[i],nf[j]) || j<i)) dead[i]=1;
        fr.clear(); fr.reserve(nf.size());
        for(size_t i=0;i<nf.size();++i) if(!dead[i] && nf[i].w>0 && nf[i].h>0) fr.push_back(nf[i]);
    }
    bool add(int t,int style) {
        int best=-1,brot=0; long long ba=LLONG_MAX; int bs=INT_MAX, bl=INT_MAX;
        for(int r=0;r<=(rotate && it[t].w!=it[t].h);++r) {
            int w=r?it[t].h:it[t].w, h=r?it[t].w:it[t].h;
            for(int k=0;k<(int)fr.size();++k) if(w<=fr[k].w && h<=fr[k].h) {
                int a=fr[k].w-w, b=fr[k].h-h, sh=min(a,b), lo=max(a,b); long long ar=1LL*a*b;
                bool take=false;
                if(best<0) take=true;
                else if(style==0) take=make_tuple(sh,lo,ar)<make_tuple(bs,bl,ba);
                else if(style==1) take=make_tuple(ar,sh,lo)<make_tuple(ba,bs,bl);
                else take=make_tuple(lo,sh,ar)<make_tuple(bl,bs,ba);
                if(take) best=k,brot=r,ba=ar,bs=sh,bl=lo;
            }
        }
        if(best<0) return false;
        int w=brot?it[t].h:it[t].w, h=brot?it[t].w:it[t].h;
        R q={fr[best].x,fr[best].y,w,h};
        out.push_back({t,q.x,q.y,brot}); value+=it[t].v; carve(q); return true;
    }
};

static string esc(const string& x) { string r; for(char c:x) { if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Json j; int W=0,H=0; bool rot=false; vector<Item> a;
    j.ch(); // root object
    while(true) {
        j.ws(); if(j.p>=j.s.size()||j.s[j.p]=='}'){j.ch();break;} string key=j.str(); j.ch();
        if(key=="bin") {
            j.ch(); while(true){ j.ws(); if(j.s[j.p]=='}'){j.ch();break;} string k=j.str();j.ch();
                if(k=="W") W=(int)j.num(); else if(k=="H") H=(int)j.num(); else if(k=="allow_rotate") rot=j.boolean();
                j.ws(); if(j.s[j.p]==',') j.ch();
            }
        } else { // items
            j.ch(); while(true){ j.ws(); if(j.s[j.p]==']'){j.ch();break;} j.ch(); Item z;
                while(true){j.ws(); if(j.s[j.p]=='}'){j.ch();break;} string k=j.str();j.ch();
                    if(k=="type")z.id=j.str(); else if(k=="w")z.w=j.num(); else if(k=="h")z.h=j.num(); else if(k=="v")z.v=j.num(); else if(k=="limit")z.lim=j.num();
                    j.ws();if(j.s[j.p]==',')j.ch();
                } a.push_back(z); j.ws();if(j.s[j.p]==',')j.ch();
            }
        } j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',')j.ch();
    }
    vector<Put> ans; long long best=-1; mt19937 rng(712367);
    vector<double> alphas={0.25,0.55,0.8,1.0,1.25,1.55};
    clock_t start=clock(); int run=0;
    while(run<42 && (double)(clock()-start)/CLOCKS_PER_SEC<0.72) {
        double alpha=alphas[run%alphas.size()]; int style=(run/alphas.size())%3;
        vector<int> ord(a.size()); iota(ord.begin(),ord.end(),0);
        vector<double> noise(a.size()); for(double&x:noise) x=(double)(rng()%100000)/100000.0;
        sort(ord.begin(),ord.end(),[&](int x,int y){
            double sx=log((double)a[x].v+1)-alpha*log((double)a[x].w*a[x].h)+(run<6?0:0.32*noise[x]);
            double sy=log((double)a[y].v+1)-alpha*log((double)a[y].w*a[y].h)+(run<6?0:0.32*noise[y]);
            return sx>sy;
        });
        Packer p(W,H,rot,a);
        for(int t:ord) for(int k=0;k<a[t].lim;k++) if(!p.add(t,style)) break;
        if(p.value>best) best=p.value,ans.swap(p.out);
        ++run;
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<ans.size();++i) { if(i) cout<<','; const Put&q=ans[i];
        cout<<"{\"type\":\""<<esc(a[q.t].id)<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.rot<<'}';
    }
    cout << "]}\n";
}
