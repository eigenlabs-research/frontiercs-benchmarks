#include <bits/stdc++.h>
using namespace std;
struct Item { string name; long long w,h,v,lim; };
struct R { long long x,y,w,h; };
struct P { int t,rot; long long x,y,w,h; };

// Tiny JSON reader: the instance format is deliberately simple, but keys need not be ordered.
struct JS {
    string s; size_t p=0;
    JS(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p]:'\0'; }
    void eat(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string r; if(p<s.size()&&s[p]=='"') ++p; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size()) ++p; return r; }
    long long num(){ ws(); int sg=1; if(s[p]=='-'){sg=-1;++p;} long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return sg*x; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};

struct Sol { vector<P> p; long long val=0; };
static bool contains(const R&a,const R&b){ return b.x>=a.x&&b.y>=a.y&&b.x+b.w<=a.x+a.w&&b.y+b.h<=a.y+a.h; }
static bool intersect(const R&a,const R&b){ return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }

Sol pack(const vector<Item>& a,long long W,long long H,bool allow,const vector<int>& order,int style,int rotbias){
    vector<R> fr(1,{0,0,W,H}); vector<int> used(a.size()); Sol ans;
    for(int id:order){
        for(int rep=0;rep<a[id].lim && ans.p.size()<7000;rep++){
            int bi=-1, br=0; long long best1=LLONG_MAX,best2=LLONG_MAX,best3=LLONG_MAX;
            for(int i=0;i<(int)fr.size();i++) for(int ro=0;ro<=(allow?1:0);ro++){
                long long w=ro?a[id].h:a[id].w, h=ro?a[id].w:a[id].h;
                if(w>fr[i].w||h>fr[i].h) continue;
                long long dx=fr[i].w-w, dy=fr[i].h-h, q1,q2,q3;
                if(style==0){ q1=dx*dy; q2=min(dx,dy); q3=fr[i].y; } // best area fit
                else if(style==1){ q1=min(dx,dy); q2=max(dx,dy); q3=fr[i].y; } // short-side fit
                else { q1=fr[i].y; q2=fr[i].x; q3=dx*dy; } // bottom-left
                // Orientation tie breaks are intentional: both landscape and portrait layouts are tried.
                long long rb=(ro==rotbias?0:1);
                if(q1<best1 || (q1==best1 && (q2<best2 || (q2==best2 && (q3<best3 || (q3==best3 && rb<(br==rotbias?0:1))))))){
                    best1=q1; best2=q2; best3=q3; bi=i; br=ro;
                }
            }
            if(bi<0) break;
            long long w=br?a[id].h:a[id].w, h=br?a[id].w:a[id].h;
            R take{fr[bi].x,fr[bi].y,w,h}; ans.p.push_back({id,br,take.x,take.y,w,h}); ans.val+=a[id].v; used[id]++;
            vector<R> nf; nf.reserve(fr.size()+4);
            for(const R& q:fr){
                if(!intersect(q,take)){ nf.push_back(q); continue; }
                // Split every intersected maximal empty rectangle. The resulting rectangles may overlap,
                // but every one remains disjoint from all placed items.
                if(take.x>q.x) nf.push_back({q.x,q.y,take.x-q.x,q.h});
                if(take.x+take.w<q.x+q.w) nf.push_back({take.x+take.w,q.y,q.x+q.w-(take.x+take.w),q.h});
                long long lx=max(q.x,take.x), rx=min(q.x+q.w,take.x+take.w);
                if(take.y>q.y && rx>lx) nf.push_back({lx,q.y,rx-lx,take.y-q.y});
                if(take.y+take.h<q.y+q.h && rx>lx) nf.push_back({lx,take.y+take.h,rx-lx,q.y+q.h-(take.y+take.h)});
            }
            // Remove dominated free rectangles. This is the critical distinction from a skyline:
            // a later item may use any surviving hole, not only a current column top.
            vector<char> dead(nf.size());
            for(int i=0;i<(int)nf.size();i++) if(!dead[i]) for(int j=0;j<(int)nf.size();j++) if(i!=j && !dead[i] && contains(nf[j],nf[i]) && (nf[j].w!=nf[i].w || nf[j].h!=nf[i].h || nf[j].x!=nf[i].x || nf[j].y!=nf[i].y || j<i)) dead[i]=1;
            fr.clear(); fr.reserve(nf.size()); for(int i=0;i<(int)nf.size();i++) if(!dead[i]&&nf[i].w>0&&nf[i].h>0) fr.push_back(nf[i]);
        }
    }
    return ans;
}
int main(){
    JS z; z.eat('{'); long long W=0,H=0; bool allow=false; vector<Item>a;
    while(z.ch()== '"'){
        string key=z.str(); z.eat(':');
        if(key=="bin"){
            z.eat('{'); while(z.ch()!='}') { string k=z.str(); z.eat(':'); if(k=="W")W=z.num(); else if(k=="H")H=z.num(); else allow=z.boolean(); if(z.ch()==',')z.eat(','); } z.eat('}');
        } else {
            z.eat('['); while(z.ch()!=']') { z.eat('{'); Item it; while(z.ch()!='}') { string k=z.str(); z.eat(':'); if(k=="type")it.name=z.str(); else if(k=="w")it.w=z.num(); else if(k=="h")it.h=z.num(); else if(k=="v")it.v=z.num(); else it.lim=z.num(); if(z.ch()==',')z.eat(','); } z.eat('}'); a.push_back(it); if(z.ch()==',')z.eat(','); } z.eat(']');
        }
        if(z.ch()==',') z.eat(','); else break;
    }
    vector<int> den(a.size()); iota(den.begin(),den.end(),0);
    auto density=[&](int i){return (long double)a[i].v/((long double)a[i].w*a[i].h);};
    sort(den.begin(),den.end(),[&](int i,int j){ if(density(i)!=density(j))return density(i)>density(j); return a[i].v>a[j].v; });
    Sol best;
    auto tryit=[&](vector<int> o,int st,int rb){ Sol q=pack(a,W,H,allow,o,st,rb); if(q.val>best.val) best=move(q); };
    // Different free-rectangle fitness functions discriminate placement fragmentation, while value
    // selection remains density-led to match the fractional-area scoring bound.
    for(int st=0;st<3;st++) for(int rb=0;rb<=(allow?1:0);rb++) tryit(den,st,rb);
    vector<int> byval=den; sort(byval.begin(),byval.end(),[&](int i,int j){return a[i].v>a[j].v;});
    tryit(byval,0,0);
    vector<int> small=den; sort(small.begin(),small.end(),[&](int i,int j){ if(density(i)!=density(j))return density(i)>density(j); return a[i].w*a[i].h<a[j].w*a[j].h;});
    tryit(small,1,0);
    // Seed each possible high-value shape once; after it is exhausted, density order closes its holes.
    for(int first=0;first<(int)a.size();first++){ vector<int> o={first}; for(int x:den)if(x!=first)o.push_back(x); tryit(o,0,first&1); }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.p.size();i++){ const P&q=best.p[i]; if(i)cout<<','; cout<<"{\"type\":\""<<a[q.t].name<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.rot<<'}'; }
    cout << "]}";
}
