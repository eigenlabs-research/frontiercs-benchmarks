#include <bits/stdc++.h>
using namespace std;

// Small JSON reader: input is a single, ordinary JSON object as specified by the task.
struct J {
    enum {NUM, STR, BOOL, ARR, OBJ} k;
    long long n=0; string s; bool b=false;
    vector<J> a; map<string,J> o;
};
struct Parser {
    string z; size_t p=0;
    Parser(string x):z(move(x)){}
    void ws(){ while(p<z.size() && isspace((unsigned char)z[p])) ++p; }
    string str(){
        ws(); ++p; string r;
        while(p<z.size() && z[p]!='"') {
            if(z[p]=='\\' && p+1<z.size()) {
                char c=z[++p];
                if(c=='n') r+='\n'; else if(c=='t') r+='\t'; else if(c=='r') r+='\r'; else r+=c;
                ++p;
            } else r+=z[p++];
        }
        if(p<z.size()) ++p; return r;
    }
    J val(){
        ws(); J r;
        if(z[p]=='"'){ r.k=J::STR; r.s=str(); return r; }
        if(z[p]=='{'){
            r.k=J::OBJ; ++p; ws();
            while(p<z.size() && z[p]!='}') { string q=str(); ws(); ++p; r.o[q]=val(); ws(); if(z[p]==',') ++p; ws(); }
            ++p; return r;
        }
        if(z[p]=='['){
            r.k=J::ARR; ++p; ws();
            while(p<z.size() && z[p]!=']') { r.a.push_back(val()); ws(); if(z[p]==',') ++p; ws(); }
            ++p; return r;
        }
        if(z.compare(p,4,"true")==0){p+=4;r.k=J::BOOL;r.b=true;return r;}
        if(z.compare(p,5,"false")==0){p+=5;r.k=J::BOOL;r.b=false;return r;}
        r.k=J::NUM; bool neg=false; if(z[p]=='-'){neg=true;++p;} while(p<z.size()&&isdigit((unsigned char)z[p]))r.n=r.n*10+z[p++]-'0'; if(neg)r.n=-r.n; return r;
    }
};
struct Item { string name; int w,h,lim; long long v; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };

static bool intersects(const R&a,const R&b) {
    return a.x < b.x+b.w && b.x < a.x+a.w && a.y < b.y+b.h && b.y < a.y+a.h;
}
static bool contains(const R&a,const R&b) {
    return b.x>=a.x && b.y>=a.y && b.x+b.w<=a.x+a.w && b.y+b.h<=a.y+a.h;
}

// Standard maximal-rectangles packing.  Free rectangles may overlap, but every insertion
// removes the inserted area from every free rectangle, which preserves feasibility.
struct Pack {
    int W,H,rule; const vector<Item>& it;
    vector<R> freeR; vector<int> used; vector<P> out; long long score=0;
    Pack(int W,int H,int rule,const vector<Item>&it):W(W),H(H),rule(rule),it(it),used(it.size()) { freeR.push_back({0,0,W,H}); }
    bool put(int t, bool rotation) {
        int w=rotation?it[t].h:it[t].w, h=rotation?it[t].w:it[t].h;
        int bi=-1; array<long long,3> best={LLONG_MAX,LLONG_MAX,LLONG_MAX};
        for(int i=0;i<(int)freeR.size();++i) {
            R f=freeR[i]; if(w>f.w||h>f.h) continue;
            long long dx=f.w-w, dy=f.h-h, ar=1LL*f.w*f.h-1LL*w*h;
            array<long long,3> q;
            if(rule==0) q={min(dx,dy),max(dx,dy),ar};       // short-side fit
            else if(rule==1) q={ar,min(dx,dy),max(dx,dy)};  // area fit
            else if(rule==2) q={f.y,f.x,min(dx,dy)};        // bottom-left
            else q={max(dx,dy),min(dx,dy),ar};              // long-side fit
            if(q<best) best=q,bi=i;
        }
        if(bi<0) return false;
        R q={freeR[bi].x,freeR[bi].y,w,h};
        vector<R> nf; nf.reserve(freeR.size()*2+4);
        for(R f:freeR) {
            if(!intersects(f,q)) { nf.push_back(f); continue; }
            if(q.x>f.x) nf.push_back({f.x,f.y,q.x-f.x,f.h});
            if(q.x+q.w<f.x+f.w) nf.push_back({q.x+q.w,f.y,f.x+f.w-(q.x+q.w),f.h});
            if(q.y>f.y) nf.push_back({f.x,f.y,f.w,q.y-f.y});
            if(q.y+q.h<f.y+f.h) nf.push_back({f.x,q.y+q.h,f.w,f.y+f.h-(q.y+q.h)});
        }
        // Remove contained rectangles.  This is the key MaxRects pruning step.
        vector<char> dead(nf.size());
        for(int i=0;i<(int)nf.size();++i) if(!dead[i])
            for(int j=i+1;j<(int)nf.size();++j) if(!dead[j]) {
                if(contains(nf[i],nf[j])) dead[j]=1;
                else if(contains(nf[j],nf[i])) { dead[i]=1; break; }
            }
        freeR.clear(); freeR.reserve(nf.size());
        for(int i=0;i<(int)nf.size();++i) if(!dead[i]&&nf[i].w>0&&nf[i].h>0) freeR.push_back(nf[i]);
        out.push_back({t,q.x,q.y,(int)rotation,w,h}); ++used[t]; score+=it[t].v;
        return true;
    }
};

static string esc(const string&s) {
    string r; for(char c:s) { if(c=='"'||c=='\\')r+='\\'; r+=c; } return r;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string input((istreambuf_iterator<char>(cin)),{}); if(input.empty()) return 0;
    J root=Parser(input).val(); auto &bn=root.o["bin"], &ja=root.o["items"];
    int W=(int)bn.o["W"].n, H=(int)bn.o["H"].n; bool rotOK=bn.o["allow_rotate"].b;
    vector<Item> a;
    for(const J& x:ja.a) a.push_back({x.o.at("type").s,(int)x.o.at("w").n,(int)x.o.at("h").n,(int)x.o.at("limit").n,x.o.at("v").n});
    int n=a.size(); vector<int> base(n); iota(base.begin(),base.end(),0);
    vector<P> answer; long long best=-1;
    // Several order/tie-break and fit-rule choices are cheap diversification of one MaxRects mechanism.
    for(int run=0;run<8;++run) {
        vector<int> ord=base;
        sort(ord.begin(),ord.end(),[&](int i,int j){
            long double di=(long double)a[i].v/(a[i].w*a[i].h), dj=(long double)a[j].v/(a[j].w*a[j].h);
            if(run==4) { if(a[i].v!=a[j].v) return a[i].v>a[j].v; }
            else if(fabsl(di-dj)>1e-18L) return di>dj;
            long long ai=1LL*a[i].w*a[i].h, aj=1LL*a[j].w*a[j].h;
            if(run==1 && ai!=aj) return ai>aj;
            if(run==2 && ai!=aj) return ai<aj;
            if(run==3 && max(a[i].w,a[i].h)!=max(a[j].w,a[j].h)) return max(a[i].w,a[i].h)>max(a[j].w,a[j].h);
            return i<j;
        });
        Pack pk(W,H,run%4,a);
        for(int t:ord) {
            while(pk.used[t]<a[t].lim && pk.out.size()<2500) {
                bool can0=(a[t].w<=W&&a[t].h<=H), can1=rotOK&&a[t].h<=W&&a[t].w<=H;
                if(!can0&&!can1) break;
                // Test both orientations by their actual best placement metric indirectly:
                // try the orientation whose resulting bounding box is less elongated first;
                bool first=false;
                if(can1 && (!can0 || (run&1 ? a[t].h<a[t].w : a[t].h>a[t].w))) first=true;
                if(pk.put(t,first)) continue;
                if(can0&&can1&&pk.put(t,!first)) continue;
                break;
            }
        }
        if(pk.score>best) best=pk.score,answer=move(pk.out);
    }
    cout << "{\"placements\":[";
    for(size_t i=0;i<answer.size();++i) { const P&p=answer[i]; if(i) cout<<',';
        cout<<"{\"type\":\""<<esc(a[p.t].name)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
    cout << "]}\n";
}
