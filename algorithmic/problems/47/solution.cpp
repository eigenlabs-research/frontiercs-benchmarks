#include <bits/stdc++.h>
using namespace std;
struct Item { string id; int w,h; long long v,lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };

// Small JSON reader: input has only the primitive values specified by the task.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string r; if(p<s.size()) ++p; while(p<s.size() && s[p]!='"') { if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size()) ++p; return r; }
    long long num(){ ws(); long long z=0, q=1; if(s[p]=='-') q=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p])) z=z*10+s[p++]-'0'; return z*q; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};
static bool contains(const R&a,const R&b){ return a.x<=b.x&&a.y<=b.y&&a.x+a.w>=b.x+b.w&&a.y+a.h>=b.y+b.h; }
static bool inter(const R&a,const R&b){ return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }

struct Pack {
    int W,H; vector<Item> &it; bool rotok; vector<R> fr; vector<P> out; vector<int> used;
    Pack(int W,int H,vector<Item>&a,bool r):W(W),H(H),it(a),rotok(r),fr{{0,0,W,H}},used(a.size()){}
    bool bestSpot(int t,int &bi,int &bro,int &bw,int &bh) {
        long long best=LLONG_MAX; bi=-1; bro=0;
        for(int i=0;i<(int)fr.size();++i) for(int ro=0;ro<= (rotok && it[t].w!=it[t].h);++ro){
            int w=ro?it[t].h:it[t].w, h=ro?it[t].w:it[t].h;
            if(w>fr[i].w||h>fr[i].h) continue;
            long long a=(long long)fr[i].w*fr[i].h-(long long)w*h;
            long long sh=min(fr[i].w-w,fr[i].h-h), lg=max(fr[i].w-w,fr[i].h-h);
            // area fit avoids wasting scarce small holes; short side breaks ties.
            long long key=a*1000000LL+sh*1000+lg;
            if(key<best){best=key;bi=i;bro=ro;bw=w;bh=h;}
        }
        return bi>=0;
    }
    void put(int t,int fi,int ro,int w,int h){
        R q{fr[fi].x,fr[fi].y,w,h}; out.push_back({t,q.x,q.y,ro,w,h}); ++used[t];
        vector<R> add; vector<char> dead(fr.size());
        for(int i=0;i<(int)fr.size();++i) if(inter(fr[i],q)){
            R f=fr[i]; dead[i]=1;
            if(q.x>f.x) add.push_back({f.x,f.y,q.x-f.x,f.h});
            if(q.x+q.w<f.x+f.w) add.push_back({q.x+q.w,f.y,f.x+f.w-q.x-q.w,f.h});
            if(q.y>f.y) add.push_back({f.x,f.y,f.w,q.y-f.y});
            if(q.y+q.h<f.y+f.h) add.push_back({f.x,q.y+q.h,f.w,f.y+f.h-q.y-q.h});
        }
        vector<R> nf; nf.reserve(fr.size()+add.size());
        for(int i=0;i<(int)fr.size();++i) if(!dead[i]) nf.push_back(fr[i]);
        for(R a:add) if(a.w>0&&a.h>0){
            bool bad=false;
            for(auto &b:nf) if(contains(b,a)){bad=true;break;}
            if(bad) continue;
            for(int i=(int)nf.size()-1;i>=0;--i) if(contains(a,nf[i])) nf.erase(nf.begin()+i);
            nf.push_back(a);
        }
        fr.swap(nf);
    }
    long long run(vector<int> ord){
        // A full type is consumed before the next one.  This deliberately tests the
        // density-first assumption against several independent priority orders.
        for(int t:ord) while(used[t]<it[t].lim){ int f,ro,w,h; if(!bestSpot(t,f,ro,w,h)) break; put(t,f,ro,w,h); }
        long long z=0; for(auto&p:out) z+=it[p.t].v; return z;
    }
};
static void emitString(const string& s){
    cout << '"';
    for(unsigned char c:s){ if(c=='"'||c=='\\') cout<<'\\'<<(char)c; else if(c=='\n') cout<<"\\n"; else if(c=='\r') cout<<"\\r"; else if(c=='\t') cout<<"\\t"; else cout<<(char)c; }
    cout << '"';
}
int main(){
    Json j; int W=0,H=0; bool ar=false; vector<Item>a;
    j.ch('{'); for(int top=0;top<2;++top){ if(top)j.ch(','); string k=j.str();j.ch(':');
        if(k=="bin"){ j.ch('{'); for(int z=0;z<3;++z){if(z)j.ch(',');string q=j.str();j.ch(':');if(q=="W")W=j.num();else if(q=="H")H=j.num();else ar=j.boolean();}j.ch('}'); }
        else { j.ch('['); bool first=true; while(true){j.ws();if(j.p<j.s.size()&&j.s[j.p]==']'){++j.p;break;}if(!first)j.ch(',');first=false;j.ch('{');Item x;for(int z=0;z<5;++z){if(z)j.ch(',');string q=j.str();j.ch(':');if(q=="type")x.id=j.str();else if(q=="w")x.w=j.num();else if(q=="h")x.h=j.num();else if(q=="v")x.v=j.num();else x.lim=j.num();}j.ch('}');a.push_back(x);} }
    }
    int n=a.size(); vector<vector<int>> orders;
    vector<int> base(n); iota(base.begin(),base.end(),0);
    auto add=[&](auto cmp){auto v=base;sort(v.begin(),v.end(),cmp);orders.push_back(v);};
    add([&](int x,int y){return (__int128)a[x].v*a[y].w*a[y].h>(__int128)a[y].v*a[x].w*a[x].h;});
    add([&](int x,int y){return a[x].v>a[y].v;});
    add([&](int x,int y){return (__int128)a[x].v*a[y].w*a[y].h*min(a[y].w,a[y].h)>(__int128)a[y].v*a[x].w*a[x].h*min(a[x].w,a[x].h);});
    add([&](int x,int y){return (__int128)a[x].v*max(a[y].w,a[y].h)>(__int128)a[y].v*max(a[x].w,a[x].h);});
    // Tie alternatives matter because all item fields are often deliberately close.
    for(int shift=1;shift<min(n,3);++shift){auto v=orders[0]; rotate(v.begin(),v.begin()+shift,v.end());orders.push_back(v);}
    long long bv=-1; vector<P> ans;
    for(auto &o:orders){ Pack q(W,H,a,ar); long long z=q.run(o); if(z>bv){bv=z;ans=q.out;} }
    cout<<"{\"placements\":[";
    for(int i=0;i<(int)ans.size();++i){if(i)cout<<',';auto&p=ans[i];cout<<"{\"type\":";emitString(a[p.t].id);cout<<",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}
    cout<<"]}";
}
