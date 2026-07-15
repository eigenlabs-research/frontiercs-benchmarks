#include <bits/stdc++.h>
using namespace std;

struct Item { string name; int w,h,lim; long long v; };
struct R { int x,y,w,h; };
struct Pl { int id,x,y,rot,w,h; };

struct Parser {
    string s; size_t p=0;
    Parser(string z):s(move(z)){}
    void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
    void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
    string str(){
        ws(); string r; if(p<s.size() && s[p]=='"') ++p;
        while(p<s.size() && s[p]!='"') { if(s[p]=='\\' && p+1<s.size()) ++p; r+=s[p++]; }
        if(p<s.size()) ++p; return r;
    }
    long long num(){ ws(); int sg=1; if(s[p]=='-'){sg=-1;++p;} long long z=0; while(p<s.size()&&isdigit((unsigned char)s[p])) z=z*10+s[p++]-'0'; return sg*z; }
    bool boolean(){ ws(); bool z=s.compare(p,4,"true")==0; p+=z?4:5; return z; }
};

static bool contains(const R&a,const R&b) {
    return a.x<=b.x && a.y<=b.y && a.x+a.w>=b.x+b.w && a.y+a.h>=b.y+b.h;
}
static void addFree(vector<R>& a, R q) {
    if(q.w<=0 || q.h<=0) return;
    for(int i=0;i<(int)a.size();) {
        if(contains(a[i],q)) return;
        if(contains(q,a[i])) { a[i]=a.back(); a.pop_back(); }
        else ++i;
    }
    a.push_back(q);
}
static void cut(vector<R>& fr, const R& u) {
    vector<R> nxt; nxt.reserve(fr.size()+8);
    for(const R& f: fr) {
        int ix1=max(f.x,u.x), ix2=min(f.x+f.w,u.x+u.w);
        int iy1=max(f.y,u.y), iy2=min(f.y+f.h,u.y+u.h);
        if(ix1>=ix2 || iy1>=iy2) { addFree(nxt,f); continue; }
        // The four pieces are the usual maximal-rectangles subtraction pieces.
        addFree(nxt,{f.x,f.y,u.x-f.x,f.h});
        addFree(nxt,{u.x+u.w,f.y,f.x+f.w-(u.x+u.w),f.h});
        addFree(nxt,{f.x,f.y,f.w,u.y-f.y});
        addFree(nxt,{f.x,u.y+u.h,f.w,f.y+f.h-(u.y+u.h)});
    }
    fr.swap(nxt);
}

struct Choice { int id=-1, fi=-1, rot=0, w=0,h=0; long long a=0,b=0,c=0; };
static vector<Pl> pack(const vector<Item>& it, int W, int H, bool canrot, int mode) {
    vector<R> fr; fr.push_back({0,0,W,H}); vector<int> used(it.size()); vector<Pl> ans;
    while(true) {
        Choice best;
        for(int k=0;k<(int)it.size();++k) if(used[k]<it[k].lim) {
            Choice own;
            for(int ro=0;ro<=(canrot?1:0);++ro) {
                if(ro && it[k].w==it[k].h) continue;
                int w=ro?it[k].h:it[k].w, h=ro?it[k].w:it[k].h;
                for(int j=0;j<(int)fr.size();++j) if(w<=fr[j].w && h<=fr[j].h) {
                    long long dw=fr[j].w-w, dh=fr[j].h-h;
                    Choice z; z.id=k;z.fi=j;z.rot=ro;z.w=w;z.h=h;
                    z.a=min(dw,dh); z.b=max(dw,dh); z.c=1LL*fr[j].w*fr[j].h-1LL*w*h;
                    auto betterFit=[&](const Choice& x,const Choice& y) {
                        if(y.id<0) return true;
                        if(mode==0) { if(x.a!=y.a)return x.a<y.a; if(x.b!=y.b)return x.b<y.b; }
                        else if(mode==1) { if(x.c!=y.c)return x.c<y.c; if(x.a!=y.a)return x.a<y.a; }
                        else { if(x.b!=y.b)return x.b<y.b; if(x.a!=y.a)return x.a<y.a; }
                        const R& X=fr[x.fi]; const R& Y=fr[y.fi];
                        if(X.y!=Y.y) return X.y<Y.y;
                        return X.x<Y.x;
                    };
                    if(betterFit(z,own)) own=z;
                }
            }
            if(own.id<0) continue;
            bool take=false;
            if(best.id<0) take=true;
            else {
                // Compare value density exactly, avoiding floating point ranking errors.
                __int128 L=(__int128)it[own.id].v*best.w*best.h;
                __int128 Rr=(__int128)it[best.id].v*own.w*own.h;
                if(L!=Rr) take=L>Rr;
                else if(own.c!=best.c) take=own.c<best.c;
                else if(own.a!=best.a) take=own.a<best.a;
            }
            if(take) best=own;
        }
        if(best.id<0) break;
        R f=fr[best.fi], u{f.x,f.y,best.w,best.h};
        ans.push_back({best.id,u.x,u.y,best.rot,best.w,best.h});
        ++used[best.id]; cut(fr,u);
    }
    return ans;
}

static string esc(const string& s) { string r; for(char c:s) { if(c=='"'||c=='\\') r+='\\'; r+=c; } return r; }
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)),{}); Parser q(input);
    vector<Item> items; int W=0,H=0; bool rot=false;
    q.ch('{');
    while(true) {
        q.ws(); if(q.p>=q.s.size()||q.s[q.p]=='}'){q.ch('}');break;}
        string key=q.str(); q.ch(':');
        if(key=="bin") {
            q.ch('{'); while(true){ q.ws(); if(q.s[q.p]=='}'){q.ch('}');break;} string k=q.str();q.ch(':');
                if(k=="W") W=(int)q.num(); else if(k=="H") H=(int)q.num(); else if(k=="allow_rotate") rot=q.boolean();
                q.ws(); if(q.s[q.p]==',')q.ch(','); }
        } else if(key=="items") {
            q.ch('['); while(true){ q.ws(); if(q.s[q.p]==']'){q.ch(']');break;} q.ch('{'); Item z;
                while(true){q.ws();if(q.s[q.p]=='}'){q.ch('}');break;} string k=q.str();q.ch(':');
                    if(k=="type")z.name=q.str(); else if(k=="w")z.w=(int)q.num(); else if(k=="h")z.h=(int)q.num(); else if(k=="v")z.v=q.num(); else if(k=="limit")z.lim=(int)q.num();
                    q.ws();if(q.s[q.p]==',')q.ch(',');}
                items.push_back(z); q.ws();if(q.s[q.p]==',')q.ch(','); }
        }
        q.ws(); if(q.p<q.s.size()&&q.s[q.p]==',')q.ch(',');
    }
    vector<Pl> answer; long long bv=-1;
    // Same density-first mechanism, with independent free-rectangle fit tie-breaks.
    for(int mode=0;mode<3;++mode) {
        vector<Pl> z=pack(items,W,H,rot,mode); long long v=0; for(auto&p:z)v+=items[p.id].v;
        if(v>bv){bv=v;answer.swap(z);}
    }
    cout << "{\"placements\":[";
    for(int i=0;i<(int)answer.size();++i){ if(i)cout<<','; const Pl&p=answer[i];
        cout<<"{\"type\":\""<<esc(items[p.id].name)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
