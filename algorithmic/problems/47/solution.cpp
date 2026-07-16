#include <bits/stdc++.h>
using namespace std;

struct Item { string name; int w,h; long long v, lim; };
struct R { int x,y,w,h; };
struct P { int id,x,y,rot; };

// The input grammar is fixed JSON; this small parser deliberately accepts arbitrary key order.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){ ws(); string r; ch('"'); while(p<s.size() && s[p]!='"') { if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; } ch('"'); return r; }
    long long num(){ ws(); long long z=0, sg=1; if(s[p]=='-') sg=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p])) z=z*10+s[p++]-'0'; return z*sg; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};

struct Solver {
    int W,H; bool rotate;
    vector<Item> a;
    vector<int> left;
    vector<P> out, best;
    long long val=0,bestval=-1;

    // Put the largest currently useful free box first.  Every placement splits one box
    // into two disjoint boxes, so feasibility follows without expensive intersection tests.
    void run(int mode) {
        left.resize(a.size()); for(int i=0;i<(int)a.size();++i) left[i]=a[i].lim;
        vector<R> boxes; boxes.push_back({0,0,W,H});
        priority_queue<pair<long long,int>> q; q.push({1LL*W*H,0});
        auto addbox = [&](R r) { boxes.push_back(r); q.push({1LL*r.w*r.h,(int)boxes.size()-1}); };
        out.clear(); val=0;
        while(!q.empty()) {
            R r=boxes[q.top().second]; q.pop();
            int bi=-1, br=0; long double bs=-1;
            for(int i=0;i<(int)a.size();++i) if(left[i]) {
                for(int z=0;z<=(rotate && a[i].w!=a[i].h);++z) {
                    int iw=z?a[i].h:a[i].w, ih=z?a[i].w:a[i].h;
                    if(iw>r.w || ih>r.h) continue;
                    long double den=(long double)a[i].v/(iw*1LL*ih);
                    long double used=(long double)(iw*1LL*ih)/(r.w*1LL*r.h);
                    // Density is the governing criterion.  The small fit term changes only
                    // near ties and makes the residual rectangles less fragmented.
                    long double sc=den*(1.0L + (mode%3==0 ? .10L : mode%3==1 ? .035L : .18L)*used);
                    if(mode>=6) sc += (long double)a[i].v/(W*1LL*H)*.025L;
                    if(sc>bs+1e-18L || (fabsl(sc-bs)<1e-18L && iw*ih > (bi<0?0:(br?a[bi].h*a[bi].w:a[bi].w*a[bi].h)))) bs=sc,bi=i,br=z;
                }
            }
            if(bi<0) continue;
            int iw=br?a[bi].h:a[bi].w, ih=br?a[bi].w:a[bi].h;
            out.push_back({bi,r.x,r.y,br}); val+=a[bi].v; --left[bi];
            int rw=r.w-iw, rh=r.h-ih;
            // Both guillotine cuts are valid.  Alternate the cut direction by variant and
            // shape, avoiding a systematic bias toward either long axis.
            bool vertical;
            if(mode%4==0) vertical = rw*ih >= iw*rh;
            else if(mode%4==1) vertical = rw*ih < iw*rh;
            else if(mode%4==2) vertical = r.w>=r.h;
            else vertical = r.w<r.h;
            if(vertical) {
                if(rw) addbox({r.x+iw,r.y,rw,r.h});
                if(rh) addbox({r.x,r.y+ih,iw,rh});
            } else {
                if(rh) addbox({r.x,r.y+ih,r.w,rh});
                if(rw) addbox({r.x+iw,r.y,rw,ih});
            }
        }
        if(val>bestval) bestval=val,best=out;
    }
};

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Json j; j.ws(); j.ch('{');
    int W=0,H=0; bool rot=false; vector<Item> items;
    // Top-level has exactly bin and items, but accept either order.
    for(int top=0;top<2;++top){
        string key=j.str(); j.ch(':');
        if(key=="bin"){
            j.ch('{'); for(int k=0;k<3;++k){ string x=j.str(); j.ch(':'); if(x=="W") W=j.num(); else if(x=="H") H=j.num(); else rot=j.boolean(); if(k<2) j.ch(','); } j.ch('}');
        } else {
            j.ch('['); j.ws(); while(j.p<j.s.size() && j.s[j.p]!=']'){
                j.ch('{'); Item it; for(int k=0;k<5;++k){ string x=j.str(); j.ch(':'); if(x=="type") it.name=j.str(); else if(x=="w") it.w=j.num(); else if(x=="h") it.h=j.num(); else if(x=="v") it.v=j.num(); else it.lim=j.num(); if(k<4) j.ch(','); } j.ch('}'); items.push_back(it); j.ws(); if(j.s[j.p]==',') j.ch(',');
            } j.ch(']');
        }
        if(top==0) j.ch(',');
    }
    Solver s; s.W=W;s.H=H;s.rotate=rot;s.a=items;
    // A compact portfolio changes only residual-cut and tie behavior, retaining the same
    // density-first constructive mechanism in every candidate.
    for(int k=0;k<10;++k) s.run(k);
    cout << "{\"placements\":[";
    for(size_t k=0;k<s.best.size();++k){ if(k) cout<<','; const P&p=s.best[k]; cout<<"{\"type\":\""<<items[p.id].name<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
