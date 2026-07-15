#include <bits/stdc++.h>
using namespace std;

struct Item { string id; int w,h; long long v, lim; };
struct Placement { int t,x,y,rot; };

// Tiny JSON reader: the input grammar contains only objects, arrays, strings, integers and booleans.
struct Json {
    string s; size_t p=0;
    Json() { s.assign(istreambuf_iterator<char>(cin), {}); }
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    char get(){ ws(); return p<s.size()?s[p++]:0; }
    string str(){
        ws(); string r; if(p>=s.size() || s[p++]!='\"') return r;
        while(p<s.size() && s[p]!='\"') { char c=s[p++]; if(c=='\\' && p<s.size()) c=s[p++]; r+=c; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); bool neg=false; if(s[p]=='-') neg=true,++p; long long x=0; while(p<s.size()&&isdigit((unsigned char)s[p])) x=x*10+s[p++]-'0'; return neg?-x:x; }
    bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; while(p<s.size()&&isalpha((unsigned char)s[p]))++p; return r; }
};

struct Node { int x,w,h; };
struct Result { long long value=0; vector<Placement> p; };

static bool fitAt(const vector<Node>& a, int i, int rw, int rh, int W, int H, int &yy) {
    int left=rw, y=0;
    if(a[i].x+rw>W) return false;
    for(int j=i; left>0; ++j) {
        if(j==(int)a.size()) return false;
        y=max(y,a[j].h); left-=a[j].w;
    }
    yy=y;
    return y+rh<=H;
}
static void raiseSkyline(vector<Node>& a, int x, int rw, int nh) {
    int z=x+rw; vector<Node> b; b.reserve(a.size()+2);
    for(auto q:a) {
        int e=q.x+q.w;
        if(e<=x || q.x>=z) b.push_back(q);
        else {
            if(q.x<x) b.push_back({q.x,x-q.x,q.h});
            if(e>z) b.push_back({z,e-z,q.h});
        }
    }
    b.push_back({x,rw,nh});
    sort(b.begin(),b.end(),[](const Node&A,const Node&B){return A.x<B.x;});
    a.clear();
    for(auto q:b) {
        if(q.w<=0) continue;
        if(!a.empty() && a.back().x+a.back().w==q.x && a.back().h==q.h) a.back().w+=q.w;
        else a.push_back(q);
    }
}

// A skyline is a compact state representation: every candidate is a different
// density/value ordering and is constructed independently, rather than repaired from one packing.
static Result build(const vector<Item>& it, int W, int H, bool allowRot, int trial) {
    int n=it.size(); vector<Node> sky={{0,W,0}}; vector<long long> used(n);
    double md=0,mv=0;
    for(auto &q:it) md=max(md,(double)q.v/q.w/q.h), mv=max(mv,(double)q.v);
    static const double aa[]={.25,.55,.85,1.15,1.6,2.2,3.0};
    static const double bb[]={0,.04,.10,.18,.28};
    static const double cc[]={0,.015,.04,.08,.13};
    double ap=aa[trial%7], bp=bb[(trial/7)%5], cp=cc[(trial/35)%5];
    vector<double> key(n);
    unsigned seed=1234567u+trial*10007u;
    for(int i=0;i<n;i++) {
        seed=seed*1664525u+1013904223u;
        double jitter=.84 + ((seed>>8)&1023)/1023.0*.32;
        double d=((double)it[i].v/it[i].w/it[i].h)/md;
        key[i]=pow(d,ap)*pow((double)it[i].v/mv,bp)*jitter;
    }
    Result out;
    const int CAP=5000;
    while((int)out.p.size()<CAP) {
        double best=-1e100; int bi=-1,bj=-1,br=0,by=0;
        for(int j=0;j<(int)sky.size();j++) for(int i=0;i<n;i++) if(used[i]<it[i].lim) {
            for(int r=0;r<= (allowRot && it[i].w!=it[i].h);r++) {
                int w=r?it[i].h:it[i].w, h=r?it[i].w:it[i].h, y;
                if(!fitAt(sky,j,w,h,W,H,y)) continue;
                // Prefer valuable dense pieces, but break priority ties by keeping the skyline low.
                double score=key[i]-cp*(double)(y+h)/H;
                score-=.0015*(double)(W-(sky[j].x+w))/W;
                if(r && (trial&1)) score+=.0007;
                if(score>best) best=score,bi=i,bj=j,br=r,by=y;
            }
        }
        if(bi<0) break;
        int w=br?it[bi].h:it[bi].w, h=br?it[bi].w:it[bi].h;
        int x=sky[bj].x;
        out.p.push_back({bi,x,by,br}); out.value+=it[bi].v; ++used[bi];
        raiseSkyline(sky,x,w,by+h);
    }
    return out;
}
static string esc(const string& x) { string r; for(char c:x) { if(c=='\"'||c=='\\') r+='\\'; r+=c; } return r; }

int main(){
    Json j; int W=0,H=0; bool rot=false; vector<Item> it;
    if(j.get()!='{') return 0;
    while(true) {
        j.ws(); if(j.p>=j.s.size() || j.s[j.p]=='}') { ++j.p; break; }
        string k=j.str(); j.get();
        if(k=="bin") {
            j.get(); while(true) { j.ws(); if(j.s[j.p]=='}'){++j.p;break;} string z=j.str();j.get();
                if(z=="W") W=(int)j.num(); else if(z=="H") H=(int)j.num(); else if(z=="allow_rotate") rot=j.boolean();
                j.ws(); if(j.s[j.p]==',')++j.p;
            }
        } else if(k=="items") {
            j.get(); while(true) { j.ws(); if(j.s[j.p]==']'){++j.p;break;} j.get(); Item q;
                while(true) { j.ws(); if(j.s[j.p]=='}'){++j.p;break;} string z=j.str();j.get();
                    if(z=="type")q.id=j.str(); else if(z=="w")q.w=(int)j.num(); else if(z=="h")q.h=(int)j.num(); else if(z=="v")q.v=j.num(); else if(z=="limit")q.lim=j.num();
                    j.ws();if(j.s[j.p]==',')++j.p;
                } it.push_back(q); j.ws();if(j.s[j.p]==',')++j.p;
            }
        }
        j.ws(); if(j.p<j.s.size()&&j.s[j.p]==',') ++j.p;
    }
    Result best;
    // 7*5*5 independent orderings, still small because only 8--12 types are inspected per skyline edge.
    for(int t=0;t<175;t++) { Result r=build(it,W,H,rot,t); if(r.value>best.value) best=move(r); }
    cout << "{\"placements\":[";
    for(size_t k=0;k<best.p.size();k++) { auto q=best.p[k]; if(k) cout<<',';
        cout<<"{\"type\":\""<<esc(it[q.t].id)<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.rot<<'}';
    }
    cout<<"]}\n";
}
