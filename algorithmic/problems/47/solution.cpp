#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v; int lim; };
struct Pl { int t,x,y,r; };

// Small JSON reader: the input grammar is deliberately simple, but accepting arbitrary key order
// makes the solver independent of formatting used by the generator.
struct Json {
    string s; size_t p=0;
    Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char ch(){ ws(); return p<s.size()?s[p]:0; }
    void eat(char c){ ws(); if(p<s.size()&&s[p]==c) ++p; }
    string str(){ ws(); string r; if(p>=s.size()||s[p++]!='"') return r; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\'&&p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size())++p; return r; }
    long long num(){ ws(); long long z=0,sg=1; if(s[p]=='-')sg=-1,++p; while(p<s.size()&&isdigit((unsigned char)s[p]))z=z*10+s[p++]-'0'; return z*sg; }
    bool boolean(){ ws(); bool b=s.compare(p,4,"true")==0; p+=b?4:5; return b; }
};

struct Candidate {
    long long val=0; vector<Pl> p;
};
int W,H; bool rotateOK; vector<Item> a;

// A skyline has one height per integer x coordinate.  Every new rectangle is set above the
// maximum of its footprint, so it is non-overlapping without any later geometric repair.
Candidate packSkyline(vector<int> ord, int flavor) {
    vector<int> sky(W,0), used(a.size()); Candidate out;
    for(int id:ord) {
        const Item &q=a[id];
        for(int copy=0; copy<q.lim; ++copy) {
            int bestx=-1,besty=INT_MAX,bestr=0,besttop=INT_MAX;
            long long bestspread=LLONG_MAX;
            vector<long long> pref(W+1); for(int z=0;z<W;z++) pref[z+1]=pref[z]+sky[z];
            int tries=rotateOK && q.w!=q.h ? 2:1;
            for(int r=0;r<tries;r++) {
                int rw=r?q.h:q.w, rh=r?q.w:q.h;
                if(rw>W||rh>H) continue;
                // Sliding maximum gives the floor height at every horizontal location in O(W).
                deque<int> dq;
                for(int x=0;x<W;x++) {
                    while(!dq.empty() && dq.front()<x-rw+1) dq.pop_front();
                    while(!dq.empty() && sky[dq.back()]<=sky[x]) dq.pop_back();
                    dq.push_back(x);
                    if(x+1<rw) continue;
                    int left=x-rw+1, y=sky[dq.front()];
                    if(y+rh>H) continue;
                    int top=y+rh;
                    long long spread=0;
                    // Tie criteria deliberately vary across starts, retaining diverse contours.
                    if(flavor&1) spread=(long long)rw*top-(pref[x+1]-pref[left]);
                    bool take=false;
                    if(y<besty) take=true;
                    else if(y==besty) {
                        if((flavor&2) ? top<besttop : spread<bestspread) take=true;
                        else if(top==besttop && spread==bestspread && left<bestx) take=true;
                    }
                    if(take) bestx=left,besty=y,bestr=r,besttop=top,bestspread=spread;
                }
            }
            if(bestx<0) break; // This type cannot be placed under the current contour.
            int rw=bestr?q.h:q.w, rh=bestr?q.w:q.h;
            for(int x=bestx;x<bestx+rw;x++) sky[x]=besty+rh;
            out.p.push_back({id,bestx,besty,bestr}); out.val+=q.v; ++used[id];
        }
    }
    return out;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Json j; j.eat('{');
    // Top-level input has exactly bin and items, although their order is not assumed here.
    for(int top=0;top<2;top++) {
        if(top) j.eat(','); string key=j.str(); j.eat(':');
        if(key=="bin") {
            j.eat('{'); for(int k=0;k<3;k++){ if(k)j.eat(','); string z=j.str();j.eat(':');
                if(z=="W") W=(int)j.num(); else if(z=="H") H=(int)j.num(); else rotateOK=j.boolean(); }
            j.eat('}');
        } else {
            j.eat('['); bool first=true;
            while(j.ch()!=']') { if(!first)j.eat(','); first=false; j.eat('{'); Item q;
                for(int k=0;k<5;k++){if(k)j.eat(',');string z=j.str();j.eat(':');
                    if(z=="type")q.id=j.str(); else if(z=="w")q.w=(int)j.num(); else if(z=="h")q.h=(int)j.num(); else if(z=="v")q.v=j.num(); else q.lim=(int)j.num();}
                j.eat('}'); a.push_back(q);
            } j.eat(']');
        }
    }
    vector<vector<int>> orders;
    int n=a.size(); vector<int> base(n); iota(base.begin(),base.end(),0);
    // Independent objective orderings are cheap multi-starts over the same geometric mechanism.
    for(int mode=0;mode<8;mode++) { auto o=base;
        sort(o.begin(),o.end(),[&](int x,int y){
            const Item &u=a[x],&v=a[y]; long double du=(long double)u.v/(u.w*u.h), dv=(long double)v.v/(v.w*v.h);
            if(mode==0) {if(du!=dv)return du>dv;}
            else if(mode==1) {if(u.v!=v.v)return u.v>v.v;}
            else if(mode==2) {if(u.h!=v.h)return u.h<v.h; if(du!=dv)return du>dv;}
            else if(mode==3) {if(u.w!=v.w)return u.w<v.w; if(du!=dv)return du>dv;}
            else if(mode==4) {if(max(u.w,u.h)!=max(v.w,v.h))return max(u.w,u.h)<max(v.w,v.h); if(du!=dv)return du>dv;}
            else if(mode==5) {if((long long)u.v*v.w*v.h!=(long long)v.v*u.w*u.h)return (long long)u.v*v.w*v.h>(long long)v.v*u.w*u.h;}
            else if(mode==6) {if(u.v!=v.v)return u.v>v.v; if(u.w*u.h!=v.w*v.h)return u.w*u.h<v.w*v.h;}
            else { unsigned hx=x*1103515245u+12345u,hy=y*1103515245u+12345u; if((hx&1023)!=(hy&1023))return (hx&1023)>(hy&1023);}
            return x<y;
        }); orders.push_back(o);
    }
    Candidate best;
    for(int k=0;k<(int)orders.size();k++) for(int f=0;f<4;f++) { Candidate c=packSkyline(orders[k],f); if(c.val>best.val)best=move(c); }
    cout<<"{\"placements\":[";
    for(size_t k=0;k<best.p.size();k++){ if(k)cout<<','; auto z=best.p[k]; cout<<"{\"type\":\""<<a[z.t].id<<"\",\"x\":"<<z.x<<",\"y\":"<<z.y<<",\"rot\":"<<z.r<<'}'; }
    cout<<"]}\n";
}
