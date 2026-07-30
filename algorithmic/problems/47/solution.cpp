#include <bits/stdc++.h>
using namespace std;
using ll = long long;

struct Parser {
    const string &s; size_t p=0;
    Parser(const string &s):s(s){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void take(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); take('"'); string r;
        while(p<s.size() && s[p]!='"') {
            char c=s[p++];
            if(c=='\\' && p<s.size()) {
                char e=s[p++];
                if(e=='n') r+='\n'; else if(e=='r') r+='\r'; else if(e=='t') r+='\t';
                else r+=e;
            } else r+=c;
        }
        if(p<s.size()) ++p;
        return r;
    }
    ll num(){ ws(); bool neg=p<s.size()&&s[p]=='-'; if(neg)++p; ll x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return neg?-x:x; }
    bool boolean(){ ws(); bool x=p<s.size()&&s[p]=='t'; p+=x?4:5; return x; }
};
struct Item { string id; int w=0,h=0,limit=0; ll v=0; };
struct Rect { int x,y,w,h; };
struct Place { int t,x,y,r,w,h; };

static bool intersects(const Rect&a,const Rect&b){
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}
static bool contains(const Rect&a,const Rect&b){
    return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h;
}

struct MaxRects {
    int W,H; vector<Rect> freeR; vector<Place> placed;
    MaxRects(int W,int H):W(W),H(H){ freeR.push_back({0,0,W,H}); }

    bool put(int t, int aw, int ah, bool rotate, int mode) {
        long long best1=LLONG_MAX,best2=LLONG_MAX,best3=LLONG_MAX;
        int bi=-1,bw=0,bh=0,br=0,bx=0,by=0;
        for(int i=0;i<(int)freeR.size();++i) {
            const Rect &f=freeR[i];
            for(int r=0;r<(rotate && aw!=ah ? 2:1);++r) {
                int w=r?ah:aw, h=r?aw:ah;
                if(w>f.w || h>f.h) continue;
                long long a,b,c;
                if(mode==0) { // best short side fit
                    a=min(f.w-w,f.h-h); b=max(f.w-w,f.h-h); c=(ll)f.y*W+f.x;
                } else if(mode==1) { // bottom-left, good for regular grids
                    a=(ll)f.y+h; b=f.x; c=(ll)(f.w-w)*(f.h-h);
                } else if(mode==2) { // best area fit
                    a=(ll)f.w*f.h-(ll)w*h; b=min(f.w-w,f.h-h); c=(ll)f.y*W+f.x;
                } else { // minimize the larger residue first
                    a=max(f.w-w,f.h-h); b=min(f.w-w,f.h-h); c=(ll)f.y*W+f.x;
                }
                if(tuple<ll,ll,ll>(a,b,c)<tuple<ll,ll,ll>(best1,best2,best3)) {
                    best1=a;best2=b;best3=c;bi=i;bw=w;bh=h;br=r;bx=f.x;by=f.y;
                }
            }
        }
        if(bi<0) return false;
        Rect u{bx,by,bw,bh};
        vector<Rect> nf; nf.reserve(freeR.size()+4);
        for(const Rect &f: freeR) {
            if(!intersects(f,u)) { nf.push_back(f); continue; }
            if(u.x>f.x) nf.push_back({f.x,f.y,u.x-f.x,f.h});
            if(u.x+u.w<f.x+f.w) nf.push_back({u.x+u.w,f.y,f.x+f.w-u.x-u.w,f.h});
            if(u.y>f.y) nf.push_back({f.x,f.y,f.w,u.y-f.y});
            if(u.y+u.h<f.y+f.h) nf.push_back({f.x,u.y+u.h,f.w,f.y+f.h-u.y-u.h});
        }
        // Remove free rectangles contained in another one.  The free list is normally small;
        // sorting out duplicates first keeps the quadratic pass cheap.
        sort(nf.begin(),nf.end(),[](const Rect&a,const Rect&b){
            if(a.x!=b.x)return a.x<b.x;if(a.y!=b.y)return a.y<b.y;if(a.w!=b.w)return a.w>b.w;return a.h>b.h;
        });
        nf.erase(unique(nf.begin(),nf.end(),[](const Rect&a,const Rect&b){return a.x==b.x&&a.y==b.y&&a.w==b.w&&a.h==b.h;}),nf.end());
        vector<char> dead(nf.size());
        for(int i=0;i<(int)nf.size();++i) if(!dead[i])
            for(int j=0;j<(int)nf.size();++j) if(i!=j && !dead[i] && contains(nf[j],nf[i])) dead[i]=1;
        freeR.clear(); freeR.reserve(nf.size());
        for(int i=0;i<(int)nf.size();++i) if(!dead[i] && nf[i].w>0 && nf[i].h>0) freeR.push_back(nf[i]);
        placed.push_back({t,bx,by,br,bw,bh});
        return true;
    }
};

static string quoteJSON(const string&s){
    string r="\"";
    for(unsigned char c:s){
        if(c=='"'||c=='\\'){r+='\\';r+=c;}
        else if(c=='\n')r+="\\n"; else if(c=='\r')r+="\\r"; else if(c=='\t')r+="\\t";
        else r+=c;
    }
    return r+'"';
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)),{}); Parser q(input);
    int W=0,H=0; bool allow=false; vector<Item> it;
    q.take('{');
    for(int top=0;top<2;++top){
        if(top)q.take(','); string key=q.str();q.take(':');
        if(key=="bin"){
            q.take('{');for(int j=0;j<3;++j){if(j)q.take(',');string k=q.str();q.take(':');if(k=="W")W=q.num();else if(k=="H")H=q.num();else allow=q.boolean();}q.take('}');
        } else {
            q.take('['); bool first=true;
            while(1){q.ws();if(q.p<input.size()&&input[q.p]==']'){++q.p;break;}if(!first)q.take(',');first=false;q.take('{');Item a;
                for(int j=0;j<5;++j){if(j)q.take(',');string k=q.str();q.take(':');if(k=="type")a.id=q.str();else if(k=="w")a.w=q.num();else if(k=="h")a.h=q.num();else if(k=="v")a.v=q.num();else a.limit=q.num();}
                q.take('}');it.push_back(a);
            }
        }
    }
    const int M=it.size();
    vector<Place> answer; ll answerValue=0;
    mt19937_64 rng(20260732ULL ^ ((uint64_t)W<<32) ^ H);
    auto started=chrono::steady_clock::now();

    // Each run changes both the economic ordering and the geometric placement rule.
    // The first runs are deterministic (important under a tight time limit); later runs
    // mildly perturb density to escape bad fragmentation patterns.
    for(int run=0; run<300; ++run){
        if(run>=8 && chrono::duration<double>(chrono::steady_clock::now()-started).count()>0.78) break;
        vector<pair<long double,int>> order;
        for(int t=0;t<M;++t){
            long double area=(long double)it[t].w*it[t].h;
            long double den=it[t].v/area, key;
            if(run==0||run==1) key=den;
            else if(run==2) key=den*pow(area,0.08L);
            else if(run==3) key=den/pow(area,0.08L);
            else if(run==4) key=(long double)it[t].v/sqrt(area);
            else if(run==5) key=den*(1.0L+0.10L*min(it[t].w,it[t].h)/max(it[t].w,it[t].h));
            else if(run==6) key=den*(1.0L+0.12L/sqrt((long double)it[t].limit));
            else if(run==7) key=den*(1.0L+0.08L*max(it[t].w,it[t].h)/max(W,H));
            else {
                long double noise=((long double)(rng()%1000001)/1000000.0L-0.5L);
                // Alternate conservative and adventurous perturbations rather than
                // letting late runs become essentially random.
                long double spread=0.08L+0.035L*((run-8)%12);
                key=log(max((long double)1e-30,den))+noise*spread;
            }
            order.push_back({-key,t});
        }
        sort(order.begin(),order.end(),[&](auto&a,auto&b){if(a.first!=b.first)return a.first<b.first;return it[a.second].v>it[b.second].v;});
        MaxRects pack(W,H); int mode=run%4;
        if(run>=8 && run%3==2) {
            // Interleave copies of similarly valuable types.  Grouping copies makes
            // clean grids, while interleaving often closes the seams those grids leave.
            vector<pair<long double,int>> seq;
            int copies=0; for(auto&a:it) copies+=a.limit; seq.reserve(copies);
            long double jitter=0.04L+0.015L*((run/3)%7);
            for(int t=0;t<M;++t) for(int n=0;n<it[t].limit;++n) {
                long double den=(long double)it[t].v/((long double)it[t].w*it[t].h);
                long double z=((long double)(rng()%1000001)/1000000.0L-0.5L);
                seq.push_back({-(log(den)+z*jitter),t});
            }
            sort(seq.begin(),seq.end()); vector<char> blocked(M);
            for(auto [dummy,t]:seq) if(!blocked[t] && !pack.put(t,it[t].w,it[t].h,allow,mode)) blocked[t]=1;
        } else {
            for(auto [dummy,t]:order){
                // If one copy no longer fits, no later copy of this type can fit either,
                // since maximal free rectangles only shrink.
                for(int n=0;n<it[t].limit;++n)
                    if(!pack.put(t,it[t].w,it[t].h,allow,mode)) break;
            }
        }
        ll val=0; for(auto&p:pack.placed) val+=it[p.t].v;
        if(val>answerValue){answerValue=val;answer=move(pack.placed);}
    }

    cout<<"{\"placements\":[";
    for(size_t i=0;i<answer.size();++i){
        const Place&p=answer[i];if(i)cout<<',';
        cout<<"{\"type\":"<<quoteJSON(it[p.t].id)<<",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<'}';
    }
    cout<<"]}\n";
}
