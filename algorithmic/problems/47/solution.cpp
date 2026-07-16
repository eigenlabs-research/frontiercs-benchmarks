#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v, lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };

// Small JSON reader: input has the fixed, simple schema from the statement.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string r; if(p<s.size()&&s[p]=='"') ++p; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size()) ++p; return r; }
    long long num(){ ws(); int sg=1; if(s[p]=='-'){sg=-1;++p;} long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sg*x; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};
static bool fit(const R&a,int w,int h){ return w<=a.w && h<=a.h; }
static bool intersects(const R&a,const R&b){ return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }
static bool contains(const R&a,const R&b){ return a.x<=b.x&&a.y<=b.y&&a.x+a.w>=b.x+b.w&&a.y+a.h>=b.y+b.h; }

struct Pack {
    vector<R> f;
    vector<P> out;
    long long val=0;
    vector<Item> *it;
    int W,H; bool rotok; int style;
    Pack(vector<Item>& a,int W_,int H_,bool ro,int st):f{{0,0,W_,H_}},it(&a),W(W_),H(H_),rotok(ro),style(st){}
    void add(R q){
        if(q.w<=0||q.h<=0) return;
        for(const R& z:f) if(contains(z,q)) return;
        size_t k=0;
        for(size_t i=0;i<f.size();++i) if(!contains(q,f[i])) f[k++]=f[i];
        f.resize(k); f.push_back(q);
    }
    void occupy(R b){
        vector<R> old; old.swap(f);
        for(const R&a:old){
            if(!intersects(a,b)){ add(a); continue; }
            if(b.x>a.x) add({a.x,a.y,b.x-a.x,a.h});
            if(b.x+b.w<a.x+a.w) add({b.x+b.w,a.y,a.x+a.w-(b.x+b.w),a.h});
            if(b.y>a.y) add({a.x,a.y,a.w,b.y-a.y});
            if(b.y+b.h<a.y+a.h) add({a.x,b.y+b.h,a.w,a.y+a.h-(b.y+b.h)});
        }
    }
    bool put(int t){
        const Item&a=(*it)[t];
        int bi=-1, br=0, bw=0,bh=0;
        long long best1=LLONG_MAX,best2=LLONG_MAX,best3=LLONG_MAX;
        for(int i=0;i<(int)f.size();++i) for(int rr=0;rr<=(rotok&&a.w!=a.h);++rr){
            int w=rr?a.h:a.w, h=rr?a.w:a.h; if(!fit(f[i],w,h)) continue;
            long long dx=f[i].w-w, dy=f[i].h-h;
            long long area=(long long)f[i].w*f[i].h-(long long)w*h;
            long long sh=min(dx,dy), lg=max(dx,dy);
            long long q1,q2,q3;
            if(style%3==0) q1=area,q2=sh,q3=lg;       // best area fit
            else if(style%3==1) q1=sh,q2=lg,q3=area;  // best short side fit
            else q1=lg,q2=sh,q3=area;                 // best long side fit
            if(q1<best1 || (q1==best1&&(q2<best2||(q2==best2&&q3<best3)))){
                best1=q1;best2=q2;best3=q3;bi=i;br=rr;bw=w;bh=h;
            }
        }
        if(bi<0) return false;
        R b{f[bi].x,f[bi].y,bw,bh};
        occupy(b);
        out.push_back({t,b.x,b.y,br,bw,bh}); val+=a.v;
        return f.size()<=7000;
    }
};

static string esc(const string& x){ string r; for(char c:x){ if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Json j; j.ch('{'); int W=0,H=0; bool allow=false; vector<Item> a;
    for(int top=0;top<2;++top){ if(top)j.ch(','); string key=j.str(); j.ch(':');
        if(key=="bin"){
            j.ch('{'); for(int z=0;z<3;++z){ if(z)j.ch(','); string k=j.str();j.ch(':'); if(k=="W")W=j.num(); else if(k=="H")H=j.num(); else allow=j.boolean(); } j.ch('}');
        } else {
            j.ch('['); bool first=true; while(true){ j.ws(); if(j.p<j.s.size()&&j.s[j.p]==']'){++j.p;break;} if(!first)j.ch(','); first=false; j.ch('{'); Item q; for(int z=0;z<5;++z){if(z)j.ch(',');string k=j.str();j.ch(':'); if(k=="type")q.id=j.str(); else if(k=="w")q.w=j.num(); else if(k=="h")q.h=j.num(); else if(k=="v")q.v=j.num(); else q.lim=j.num();} j.ch('}'); a.push_back(q); }
        }
    }
    vector<P> answer; long long best=-1; int n=a.size();
    // Several deterministic priorities only differ in tie/shape bias; all retain profit density as the main signal.
    for(int mode=0;mode<9;++mode){
        vector<int> ord(n); iota(ord.begin(),ord.end(),0);
        sort(ord.begin(),ord.end(),[&](int x,int y){
            long double dx=(long double)a[x].v/(a[x].w*a[x].h), dy=(long double)a[y].v/(a[y].w*a[y].h);
            long double sx=dx, sy=dy;
            int mx=max(a[x].w,a[x].h), my=max(a[y].w,a[y].h);
            int mn=min(a[x].w,a[x].h), nny=min(a[y].w,a[y].h);
            if(mode>=3&&mode<6){ sx*=1.0L+(long double)(mode-2)*mn/(mx+1)/50; sy*=1.0L+(long double)(mode-2)*nny/(my+1)/50; }
            if(mode>=6&&mode<9){ sx*=1.0L+(long double)(mode-5)*mx/(W+H)/30; sy*=1.0L+(long double)(mode-5)*my/(W+H)/30; }
            if(fabsl(sx-sy)>1e-18L) return sx>sy;
            if(mode&1) return a[x].v>a[y].v;
            if(mode&2) return a[x].w*a[x].h>a[y].w*a[y].h;
            return x<y;
        });
        Pack pk(a,W,H,allow,mode);
        for(int t:ord) for(long long c=0;c<a[t].lim;++c) if(!pk.put(t)) break;
        if(pk.val>best){best=pk.val;answer.swap(pk.out);}
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<answer.size();++i){ if(i)cout<<','; const P&p=answer[i]; cout<<"{\"type\":\""<<esc(a[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
