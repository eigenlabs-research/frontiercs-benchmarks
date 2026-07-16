#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot; };

struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); ch('"'); string r; while(p<s.size() && s[p]!='"') { char c=s[p++]; if(c=='\\' && p<s.size()) { char q=s[p++]; if(q=='n') r+='\n'; else if(q=='t') r+='\t'; else r+=q; } else r+=c; } ch('"'); return r; }
    long long num(){ ws(); long long z=1,n=0; if(s[p]=='-') z=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p])) n=n*10+s[p++]-'0'; return z*n; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4;return true;} p+=5;return false; }
};

static bool intersect(const R&a,const R&b) { return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h; }
static bool contained(const R&a,const R&b) { return a.x>=b.x && a.y>=b.y && a.x+a.w<=b.x+b.w && a.y+a.h<=b.y+b.h; }

// Standard maximal-free-rectangle update; it deliberately permits overlapping free regions,
// which lets later placements exploit L-shaped residual areas safely.
static void put(vector<R>& f, R u) {
    vector<R> q; q.reserve(f.size()*2+4);
    for(const R&r:f) {
        if(!intersect(r,u)) { q.push_back(r); continue; }
        if(u.x>r.x) q.push_back({r.x,r.y,u.x-r.x,r.h});
        if(u.x+u.w<r.x+r.w) q.push_back({u.x+u.w,r.y,r.x+r.w-(u.x+u.w),r.h});
        if(u.y>r.y) q.push_back({r.x,r.y,r.w,u.y-r.y});
        if(u.y+u.h<r.y+r.h) q.push_back({r.x,u.y+u.h,r.w,r.y+r.h-(u.y+u.h)});
    }
    vector<R> z; z.reserve(q.size());
    for(size_t i=0;i<q.size();++i) if(q[i].w>0&&q[i].h>0) {
        bool inside=false;
        for(size_t j=0;j<q.size();++j) if(i!=j && contained(q[i],q[j]) &&
            ((long long)q[i].w*q[i].h < (long long)q[j].w*q[j].h ||
             ((long long)q[i].w*q[i].h == (long long)q[j].w*q[j].h && j<i))) { inside=true; break; }
        if(!inside) z.push_back(q[i]);
    }
    f.swap(z);
}

static pair<long long,vector<P>> pack(int W,int H,bool rotation,const vector<Item>& a, const vector<double>& rank, double fitBias) {
    vector<R> free={{0,0,W,H}}; vector<int> used(a.size()); vector<P> ans;
    // Ranking values have arbitrary scales (density, raw value, etc.).  Normalize them so
    // that the geometric term has the same meaning in every ordering.
    double scale=1;
    for(double z:rank) scale=max(scale,z);
    long long val=0;
    // Availability is already bounded by the input contract.  Do not truncate a feasible,
    // profitable packing before all available copies have been considered.
    int outputGuard=0;
    for(const Item& it:a) outputGuard+=it.lim;
    while(!free.empty() && (int)ans.size()<outputGuard) {
        int bt=-1, bf=-1, br=0, bw=0,bh=0;
        double best=-1e300; int bestshort=INT_MAX,bestlong=INT_MAX;
        for(int t=0;t<(int)a.size();++t) if(used[t]<a[t].lim && a[t].v>0) {
            for(int o=0;o<=(rotation && a[t].w!=a[t].h);++o) {
                int w=o?a[t].h:a[t].w, h=o?a[t].w:a[t].h;
                for(int k=0;k<(int)free.size();++k) if(w<=free[k].w && h<=free[k].h) {
                    int ss=min(free[k].w-w,free[k].h-h), ls=max(free[k].w-w,free[k].h-h);
                    // A rectangle which nearly fills a currently maximal free rectangle is
                    // less likely to turn it into a thin, unusable L-shaped remainder.  This
                    // is deliberately a bounded secondary signal: density/value still leads.
                    double fill=(double)((long long)w*h)/((long long)free[k].w*free[k].h);
                    double sc=rank[t]/scale + fitBias*fill;
                    if(sc>best+1e-10 || (fabs(sc-best)<1e-10 && (ss<bestshort || (ss==bestshort&&ls<bestlong))))
                        best=sc,bestshort=ss,bestlong=ls,bt=t,bf=k,br=o,bw=w,bh=h;
                }
            }
        }
        if(bt<0) break;
        R u={free[bf].x,free[bf].y,bw,bh};
        ans.push_back({bt,u.x,u.y,br}); val+=a[bt].v; ++used[bt]; put(free,u);
    }
    return {val,ans};
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Json j; int W=0,H=0; bool rot=false; vector<Item>a;
    j.ch('{');
    while(true){ j.ws(); if(j.p>=j.s.size()||j.s[j.p]=='}'){j.ch('}');break;} string key=j.str();j.ch(':');
        if(key=="bin") { j.ch('{'); while(true){j.ws();if(j.s[j.p]=='}'){j.ch('}');break;} string k=j.str();j.ch(':'); if(k=="W")W=(int)j.num(); else if(k=="H")H=(int)j.num(); else rot=j.boolean(); j.ws();if(j.s[j.p]==',')j.ch(',');} }
        else if(key=="items") { j.ch('['); while(true){j.ws();if(j.s[j.p]==']'){j.ch(']');break;} Item x; j.ch('{'); while(true){j.ws();if(j.s[j.p]=='}'){j.ch('}');break;} string k=j.str();j.ch(':'); if(k=="type")x.id=j.str();else if(k=="w")x.w=(int)j.num();else if(k=="h")x.h=(int)j.num();else if(k=="v")x.v=j.num();else x.lim=(int)j.num();j.ws();if(j.s[j.p]==',')j.ch(',');} a.push_back(x);j.ws();if(j.s[j.p]==',')j.ch(',');} }
        j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',')j.ch(',');
    }
    vector<P> best; long long bestv=-1;
    // The zero-bias density pass is retained as a control.  The other passes test whether
    // fit-aware choices repair the fragmentation caused by an absolute item ordering.
    for(int mode=0;mode<7;++mode){ vector<double> r(a.size()); double bias=0;
        for(int i=0;i<(int)a.size();++i){ double area=(double)a[i].w*a[i].h, d=a[i].v/area;
            if(mode==0) r[i]=d;
            else if(mode==1) r[i]=d, bias=.10;
            else if(mode==2) r[i]=d*pow(area,0.12), bias=.14;
            else if(mode==3) r[i]=d/pow(area,0.12), bias=.18;
            else if(mode==4) r[i]=(double)a[i].v, bias=.12;
            else if(mode==5) r[i]=d*pow((double)min(a[i].w,a[i].h),0.22), bias=.16;
            else { unsigned z=(unsigned)(i*1103515245u+mode*12345u); r[i]=log(d+1.0)+(double)(z%1000)/1000.0*.28; bias=.22; }
        }
        auto got=pack(W,H,rot,a,r,bias); if(got.first>bestv) bestv=got.first,best=move(got.second);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i){ if(i)cout<<','; const P&p=best[i];
        cout<<"{\"type\":\""; for(char c:a[p.t].id){if(c=='\\' || c=='\"')cout<<'\\';cout<<c;} cout<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout<<"]}\n";
}
