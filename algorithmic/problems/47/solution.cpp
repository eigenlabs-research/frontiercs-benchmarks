#include <bits/stdc++.h>
using namespace std;

struct Item { string id; long long w,h,v,lim; };
struct Pl { int t; long long x,y,w,h; int r; };
struct FR { long long x,y,w,h; };

// Small JSON reader: the instance format contains only objects, arrays, strings, integers, booleans.
struct JS {
    string s; size_t p=0;
    JS(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string r; if(p>=s.size()||s[p++]!='\"') return r; while(p<s.size()&&s[p]!='\"'){ if(s[p]=='\\'&&p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size())++p; return r; }
    long long num(){ ws(); long long z=0,sg=1; if(s[p]=='-')sg=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p]))z=z*10+s[p++]-'0'; return z*sg; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};

static bool intersects(const FR&a,const Pl&b){
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}
static bool contains(const FR&a,const FR&b){ return a.x<=b.x&&a.y<=b.y&&a.x+a.w>=b.x+b.w&&a.y+a.h>=b.y+b.h; }
static void splitFree(vector<FR>& f,const Pl& q){
    vector<FR> n; n.reserve(f.size()+4);
    for(const FR&a:f){
        if(!intersects(a,q)){ n.push_back(a); continue; }
        if(q.x>a.x) n.push_back({a.x,a.y,q.x-a.x,a.h});
        if(q.x+q.w<a.x+a.w) n.push_back({q.x+q.w,a.y,a.x+a.w-q.x-q.w,a.h});
        if(q.y>a.y) n.push_back({a.x,a.y,a.w,q.y-a.y});
        if(q.y+q.h<a.y+a.h) n.push_back({a.x,q.y+q.h,a.w,a.y+a.h-q.y-q.h});
    }
    vector<char> bad(n.size());
    for(size_t i=0;i<n.size();++i) for(size_t j=0;j<n.size();++j) if(i!=j && contains(n[j],n[i])) {bad[i]=1;break;}
    f.clear(); f.reserve(n.size());
    for(size_t i=0;i<n.size();++i) if(!bad[i]&&n[i].w>0&&n[i].h>0) f.push_back(n[i]);
}

struct Result { long long value=0; vector<Pl> p; };
static uint64_t mix64(uint64_t x){ x+=0x9e3779b97f4a7c15ULL; x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL; x=(x^(x>>27))*0x94d049bb133111ebULL; return x^(x>>31); }
static bool pref(initializer_list<long long> a, initializer_list<long long> b){
    auto x=a.begin(), y=b.begin();
    for(;x!=a.end();++x,++y) if(*x!=*y) return *x<*y;
    return false;
}

static Result pack(const vector<Item>& it,long long W,long long H,bool rot,int run){
    int m=it.size(); vector<int> seq;
    // A copy cannot be useful beyond this geometric bound. This also bounds pathological inputs.
    for(int i=0;i<m;i++){
        long long cap=min(it[i].lim, max(0LL, W*H/max(1LL,it[i].w*it[i].h)));
        cap=min(cap,3000LL);
        while(cap--) seq.push_back(i);
    }
    uint64_t seed=1469598103934665603ULL ^ (uint64_t)run*11995408973635179863ULL;
    for(auto&a:it) for(unsigned char c:a.id) seed=(seed^c)*1099511628211ULL;
    vector<long double> key(m);
    for(int i=0;i<m;i++){
        long double ar=(long double)it[i].w*it[i].h;
        long double den=(long double)it[i].v/ar;
        long double val=it[i].v;
        long double sh=min(it[i].w,it[i].h), lg=max(it[i].w,it[i].h);
        long double noise=((mix64(seed+i*12347)>>11)*(1.0L/(1ULL<<53))-.5L);
        int z=run%6;
        if(z==0) key[i]=den;
        else if(z==1) key[i]=val;
        else if(z==2) key[i]=den*(1.0L+0.28L*noise);
        else if(z==3) key[i]=sqrt(max((long double)0,val))*den*(1.0L+0.18L*noise);
        else if(z==4) key[i]=den/(1.0L+0.001L*lg)+den*0.22L*noise;
        else key[i]=den*(1.0L+0.55L*noise)+val*0.00001L;
        // In every ordering, prohibitably large items get a modest preference only as a tie breaker.
        key[i]+=1e-12L*(sh+0.01L*lg);
    }
    stable_sort(seq.begin(),seq.end(),[&](int a,int b){return key[a]>key[b];});
    vector<FR> freeR={{0,0,W,H}}; Result ans;
    for(int id:seq){
        int best=-1, br=0; long long ba=LLONG_MAX, bs=LLONG_MAX, bl=LLONG_MAX, by=LLONG_MAX,bx=LLONG_MAX;
        for(int k=0;k<(int)freeR.size();k++) for(int rr=0;rr<=int(rot);rr++){
            long long w=rr?it[id].h:it[id].w, h=rr?it[id].w:it[id].h; auto&a=freeR[k];
            if(w>a.w||h>a.h) continue;
            long long dx=a.w-w,dy=a.h-h, ss=min(dx,dy), ll=max(dx,dy), area=a.w*a.h-w*h;
            bool take=false; int mode=(run/6)%3;
            if(mode==0) take=pref({ss,ll,area,a.y,a.x},{bs,bl,ba,by,bx});
            else if(mode==1) take=pref({area,ss,ll,a.y,a.x},{ba,bs,bl,by,bx});
            else take=pref({a.y,a.x,ss,ll,area},{by,bx,bs,bl,ba});
            if(take){best=k;br=rr;ba=area;bs=ss;bl=ll;by=a.y;bx=a.x;}
        }
        if(best<0) continue;
        long long w=br?it[id].h:it[id].w,h=br?it[id].w:it[id].h;
        Pl q{id,freeR[best].x,freeR[best].y,w,h,br}; ans.p.push_back(q); ans.value+=it[id].v; splitFree(freeR,q);
    }
    return ans;
}
static string esc(const string&s){ string r; for(char c:s){if(c=='\"'||c=='\\')r+='\\';r+=c;}return r; }
int main(){
    JS j; j.ch('{'); long long W=0,H=0; bool allow=false; vector<Item> it;
    for(int top=0;top<2;top++){
        if(top)j.ch(','); string k=j.str(); j.ch(':');
        if(k=="bin"){
            j.ch('{'); for(int z=0;z<3;z++){if(z)j.ch(',');string q=j.str();j.ch(':');if(q=="W")W=j.num();else if(q=="H")H=j.num();else allow=j.boolean();} j.ch('}');
        }else{
            j.ch('['); j.ws(); bool first=true;
            while(j.p<j.s.size()&&j.s[j.p]!=']'){
                if(!first)j.ch(','); first=false; j.ch('{'); Item a;
                for(int z=0;z<5;z++){if(z)j.ch(',');string q=j.str();j.ch(':');if(q=="type")a.id=j.str();else if(q=="w")a.w=j.num();else if(q=="h")a.h=j.num();else if(q=="v")a.v=j.num();else a.lim=j.num();} j.ch('}');it.push_back(a);j.ws();
            } j.ch(']');
        }
    }
    Result best; for(int r=0;r<18;r++){ Result q=pack(it,W,H,allow,r); if(q.value>best.value)best=move(q); }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<best.p.size();i++){auto&q=best.p[i];if(i)cout<<',';cout<<"{\"type\":\""<<esc(it[q.t].id)<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.r<<'}';}
    cout<<"]}";
}
