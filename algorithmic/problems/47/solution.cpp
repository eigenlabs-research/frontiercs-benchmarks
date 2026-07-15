#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v, lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };

// Small JSON reader: input has the fixed schema from the task, but object member order is not assumed.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), istreambuf_iterator<char>()); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p]:0; }
    void eat(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); string r; if(p>=s.size()||s[p++]!='"') return r;
        while(p<s.size() && s[p]!='"') { char c=s[p++]; if(c=='\\' && p<s.size()) { char q=s[p++]; if(q=='n') r+='\n'; else if(q=='t') r+='\t'; else r+=q; } else r+=c; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); int sg=1; if(s[p]=='-') sg=-1,++p; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return x*sg; }
    bool boolean(){ ws(); if(s.compare(p,4,"true")==0){p+=4;return true;} p+=5; return false; }
};

static bool contains(const R&a,const R&b) { return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h; }
static bool intersects(const R&a,const R&b) { return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }

struct Pack {
    int W,H; vector<R> f; vector<P> out;
    Pack(int W,int H):W(W),H(H) { f.push_back({0,0,W,H}); }
    // Return the actual MaxRects fit key used by put.  In particular, area fit
    // is remaining AREA, not the product of the two leftover side lengths.
    // The latter incorrectly gives a perfect score to every flush placement.
    pair<long long,long long> fitScore(int w,int h,int mode) const {
        pair<long long,long long> best={LLONG_MAX,LLONG_MAX};
        for(const R&r:f) if(w<=r.w && h<=r.h) {
            long long a=r.w-w, b=r.h-h, s1,s2;
            if(mode==1) s1=(long long)r.w*r.h-(long long)w*h, s2=min(a,b);
            else if(mode==2) s1=max(a,b), s2=min(a,b);
            else s1=min(a,b), s2=max(a,b);
            best=min(best,{s1,s2});
        }
        return best;
    }
    bool put(int t,int w,int h,int rot,int mode) {
        int bi=-1; auto want=fitScore(w,h,mode);
        if(want.first==LLONG_MAX) return false;
        for(int i=0;i<(int)f.size();++i) if(w<=f[i].w && h<=f[i].h) {
            long long a=f[i].w-w, b=f[i].h-h, s1,s2;
            if(mode==1) s1=(long long)f[i].w*f[i].h-(long long)w*h, s2=min(a,b);
            else if(mode==2) s1=max(a,b), s2=min(a,b);
            else s1=min(a,b), s2=max(a,b);
            if(make_pair(s1,s2)==want) { bi=i; break; }
        }
        if(bi<0) return false;
        R u{f[bi].x,f[bi].y,w,h}; out.push_back({t,u.x,u.y,rot,w,h});
        vector<R> nf; nf.reserve(f.size()+4);
        for(const R&q:f) {
            if(!intersects(q,u)) { nf.push_back(q); continue; }
            if(u.x>q.x) nf.push_back({q.x,q.y,u.x-q.x,q.h});
            if(u.x+u.w<q.x+q.w) nf.push_back({u.x+u.w,q.y,q.x+q.w-u.x-u.w,q.h});
            if(u.y>q.y) nf.push_back({q.x,q.y,q.w,u.y-q.y});
            if(u.y+u.h<q.y+q.h) nf.push_back({q.x,u.y+u.h,q.w,q.y+q.h-u.y-u.h});
        }
        // Removing contained maximal empty rectangles controls fragmentation without changing free space.
        vector<R> clean; clean.reserve(nf.size());
        for(const R&q:nf) {
            if(q.w<=0||q.h<=0) continue;
            bool dead=false;
            for(int i=0;i<(int)clean.size();) {
                if(contains(clean[i],q)) { dead=true; break; }
                if(contains(q,clean[i])) clean[i]=clean.back(),clean.pop_back(); else ++i;
            }
            if(!dead) clean.push_back(q);
        }
        f.swap(clean); return true;
    }
};

int main(){
    Json j; int W=0,H=0; bool allow=false; vector<Item> it;
    j.eat('{');
    while(j.ch()== '"') {
        string key=j.str(); j.eat(':');
        if(key=="bin") {
            j.eat('{'); while(j.ch()!='}') { string k=j.str(); j.eat(':'); if(k=="W")W=(int)j.num(); else if(k=="H")H=(int)j.num(); else allow=j.boolean(); if(j.ch()==',')j.eat(','); } j.eat('}');
        } else if(key=="items") {
            j.eat('['); while(j.ch()!=']') { j.eat('{'); Item z; while(j.ch()!='}') { string k=j.str(); j.eat(':'); if(k=="type")z.id=j.str(); else if(k=="w")z.w=(int)j.num(); else if(k=="h")z.h=(int)j.num(); else if(k=="v")z.v=j.num(); else z.lim=j.num(); if(j.ch()==',')j.eat(','); } j.eat('}'); it.push_back(z); if(j.ch()==',')j.eat(','); } j.eat(']');
        }
        if(j.ch()==',') j.eat(',');
    }
    vector<int> copies;
    for(int i=0;i<(int)it.size();++i) {
        bool fits=(it[i].w<=W&&it[i].h<=H)||(allow&&it[i].h<=W&&it[i].w<=H);
        if(!fits || it[i].v==0) continue;
        long long cap=min(it[i].lim, (long long)W*H/((long long)it[i].w*it[i].h));
        for(long long k=0;k<cap;k++) copies.push_back(i);
    }
    vector<P> best; long long bestv=-1;
    // Discriminating variants: same MaxRects mechanism, but deliberately vary value-density
    // ordering and fit rule. This tests whether loss is from selection order rather than output format.
    const double noise[]={0,.035,.08,.16,.30,.50,.75,1.05,.12,.24,.42,.65};
    clock_t start=clock();
    for(int run=0;run<12;run++) {
        if((double)(clock()-start)/CLOCKS_PER_SEC>.82) break;
        vector<int> a=copies;
        vector<double> pri(it.size());
        for(int q=0;q<(int)it.size();q++) {
            double den=(double)it[q].v/((double)it[q].w*it[q].h);
            unsigned x=(unsigned)(q*1103515245u + run*2654435761u + 12345u);
            double r=((x>>8)&65535)/65535.0-.5;
            // A type-level perturbation keeps copies together, yielding useful rectangular blocks.
            pri[q]=log(den+1e-30)+noise[run]*r;
            if(run==10) pri[q]=log((double)it[q].v+1.0)+.08*r;
            if(run==11) pri[q]=log(den+1e-30)+.20*r + .0005*max(it[q].w,it[q].h);
        }
        sort(a.begin(),a.end(),[&](int x,int y){ if(pri[x]!=pri[y]) return pri[x]>pri[y]; long long ax=(long long)it[x].w*it[x].h, ay=(long long)it[y].w*it[y].h; return ax>ay; });
        Pack pk(W,H);
        for(int q:a) {
            int w=it[q].w,h=it[q].h, rot=0;
            // Select the orientation by precisely the same criterion as the
            // free-rectangle choice.  This rules out the prior mismatch where
            // orientation used short-side fit even on area/long-side runs.
            if(allow && w!=h) {
                auto s0=pk.fitScore(w,h,run%3), s1=pk.fitScore(h,w,run%3);
                if(s1<s0) w=it[q].h,h=it[q].w,rot=1;
            }
            bool placed=pk.put(q,w,h,rot,run%3);
            if(!placed && allow) {
                if(rot) pk.put(q,it[q].w,it[q].h,0,run%3);
                else pk.put(q,it[q].h,it[q].w,1,run%3);
            }
        }
        long long val=0; for(auto&p:pk.out) val+=it[p.t].v;
        if(val>bestv) bestv=val,best.swap(pk.out);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<best.size();++i) { if(i) cout<<','; const P&p=best[i]; cout<<"{\"type\":\""<<it[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}";
}
