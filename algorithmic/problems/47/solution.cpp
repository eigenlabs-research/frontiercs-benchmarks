#include <bits/stdc++.h>
using namespace std;

struct Tok {
    string s; size_t p=0;
    Tok(string x):s(move(x)){}
    void ws(){while(p<s.size() && isspace((unsigned char)s[p])) ++p;}
    void ch(char c){ws(); if(p<s.size() && s[p]==c) ++p;}
    string str(){ ws(); string r; if(p>=s.size()||s[p++]!='"') return r; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\'&&p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size())++p; return r; }
    long long num(){ws(); int sg=1;if(s[p]=='-')sg=-1,++p; long long x=0;while(p<s.size()&&isdigit((unsigned char)s[p]))x=x*10+s[p++]-'0';return sg*x;}
    bool boolean(){ws(); if(s.compare(p,4,"true")==0){p+=4;return true;}p+=5;return false;}
};
struct Item { string id; int w,h,lim; long long v; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };
int W,H; bool rotation; vector<Item> it;

static bool intersects(const R&a,const R&b){return a.x<b.x+b.w&&b.x<a.x+a.w&&a.y<b.y+b.h&&b.y<a.y+a.h;}
static bool contains(const R&a,const R&b){return b.x>=a.x&&b.y>=a.y&&b.x+b.w<=a.x+a.w&&b.y+b.h<=a.y+a.h;}

vector<P> pack(const vector<int>& order, int mode) {
    vector<R> fr(1,{0,0,W,H}); vector<P> out;
    for(int t:order){
        int bw=0,bh=0,br=0,bi=-1; long long best1=LLONG_MAX,best2=LLONG_MAX,best3=LLONG_MAX;
        for(int i=0;i<(int)fr.size();++i) for(int rot=0;rot<=rotation;++rot){
            int w=rot?it[t].h:it[t].w, h=rot?it[t].w:it[t].h;
            if(w>fr[i].w||h>fr[i].h) continue;
            long long a=fr[i].w-w, b=fr[i].h-h;
            long long q1,q2,q3;
            if(mode==0){q1=min(a,b);q2=max(a,b);q3=a*b;} // tight short-side fit
            else if(mode==1){q1=a*b;q2=min(a,b);q3=max(a,b);} // minimum immediate waste
            else {q1=max(a,b);q2=min(a,b);q3=a*b;} // alternate aspect fit
            if(q1<best1||(q1==best1&&(q2<best2||(q2==best2&&q3<best3)))){
                best1=q1;best2=q2;best3=q3;bi=i;bw=w;bh=h;br=rot;
            }
        }
        if(bi<0) continue;
        R used{fr[bi].x,fr[bi].y,bw,bh};
        out.push_back({t,used.x,used.y,br,bw,bh});
        vector<R> nf; nf.reserve(fr.size()+4);
        for(const R& f:fr){
            if(!intersects(f,used)){nf.push_back(f);continue;}
            // Split every intersecting maximal rectangle.  These pieces may overlap each other,
            // but all remain empty; containment pruning preserves the MaxRects invariant.
            if(used.x>f.x) nf.push_back({f.x,f.y,used.x-f.x,f.h});
            if(used.x+used.w<f.x+f.w) nf.push_back({used.x+used.w,f.y,f.x+f.w-(used.x+used.w),f.h});
            if(used.y>f.y) nf.push_back({f.x,f.y,f.w,used.y-f.y});
            if(used.y+used.h<f.y+f.h) nf.push_back({f.x,used.y+used.h,f.w,f.y+f.h-(used.y+used.h)});
        }
        vector<char> dead(nf.size());
        for(int a=0;a<(int)nf.size();++a) for(int b=0;b<(int)nf.size();++b)
            if(a!=b && contains(nf[b],nf[a]) && (!contains(nf[a],nf[b]) || b<a)) { dead[a]=1; break; }
        fr.clear(); for(int i=0;i<(int)nf.size();++i) if(!dead[i]&&nf[i].w&&nf[i].h) fr.push_back(nf[i]);
        if(out.size()>=3500) break; // protects runtime/output on pathological tiny-item instances
    }
    return out;
}
int main(){
    ios::sync_with_stdio(false);cin.tie(nullptr);
    string s((istreambuf_iterator<char>(cin)),{}); Tok z(s); z.ch('{');
    while(true){ z.ws(); if(z.p>=s.size()||s[z.p]=='}'){z.ch('}');break;} string key=z.str();z.ch(':');
        if(key=="bin") { z.ch('{'); for(int k=0;k<3;k++){string q=z.str();z.ch(':');if(q=="W")W=z.num();else if(q=="H")H=z.num();else rotation=z.boolean();z.ws();if(z.p<s.size()&&s[z.p]==',')z.p++;} z.ch('}'); }
        else { z.ch('['); while(true){z.ws();if(s[z.p]==']'){z.p++;break;}z.ch('{');Item a;for(int k=0;k<5;k++){string q=z.str();z.ch(':');if(q=="type")a.id=z.str();else if(q=="w")a.w=z.num();else if(q=="h")a.h=z.num();else if(q=="v")a.v=z.num();else a.lim=z.num();z.ws();if(s[z.p]==',')z.p++;}z.ch('}');it.push_back(a);z.ws();if(s[z.p]==',')z.p++;} }
        z.ws();if(z.p<s.size()&&s[z.p]==',')z.p++;
    }
    vector<int> copies; long long area=1LL*W*H;
    for(int i=0;i<(int)it.size();++i){int cap=min<long long>(it[i].lim,area/(1LL*it[i].w*it[i].h)+1);for(int k=0;k<cap;k++)copies.push_back(i);}
    vector<P> answer; long long best=-1;
    for(int pass=0;pass<5;pass++){
        vector<int> o=copies;
        sort(o.begin(),o.end(),[&](int a,int b){
            __int128 da=(__int128)it[a].v*it[b].w*it[b].h, db=(__int128)it[b].v*it[a].w*it[a].h;
            if(pass==3){ if(it[a].v!=it[b].v)return it[a].v>it[b].v; }
            else if(pass==4){ if(it[a].w*it[a].h!=it[b].w*it[b].h)return it[a].w*it[a].h<it[b].w*it[b].h; }
            if(da!=db)return da>db;
            if(pass==1 && it[a].h!=it[b].h)return it[a].h>it[b].h;
            if(pass==2 && it[a].w!=it[b].w)return it[a].w>it[b].w;
            return a<b;
        });
        auto q=pack(o,pass%3); long long val=0;for(auto&p:q)val+=it[p.t].v;
        if(val>best)best=val,answer=move(q);
    }
    cout<<"{\"placements\":[";
    for(size_t i=0;i<answer.size();++i){if(i)cout<<',';auto&p=answer[i];cout<<"{\"type\":\""<<it[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}
    cout<<"]}\n";
}
