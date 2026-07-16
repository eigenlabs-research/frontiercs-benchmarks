#include <bits/stdc++.h>
using namespace std;
struct Item { string name; int w,h; long long v,lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };

// Minimal parser for the fixed JSON input schema.
struct Json {
 string s; size_t p=0;
 Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
 void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
 void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
 string str(){ ws(); string r; if(p<s.size()&&s[p]=='"')++p; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\'&&p+1<s.size())++p; r+=s[p++]; } if(p<s.size())++p; return r; }
 long long num(){ ws(); int sign=1; if(p<s.size()&&s[p]=='-')sign=-1,++p; long long r=0; while(p<s.size()&&isdigit((unsigned char)s[p]))r=r*10+s[p++]-'0'; return sign*r; }
 bool boolean(){ ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r; }
};
static inline bool inter(const R&a,const R&b){ return a.x<b.x+b.w && b.x<a.x+a.w && a.y<b.y+b.h && b.y<a.y+a.h; }
static inline bool inside(const R&a,const R&b){ return a.x>=b.x&&a.y>=b.y&&a.x+a.w<=b.x+b.w&&a.y+a.h<=b.y+b.h; }

int main(){
 ios::sync_with_stdio(false); cin.tie(nullptr);
 Json j; vector<Item>a; int W=0,H=0; bool canrot=false;
 j.ch('{');
 for(int z=0;z<2;z++){
  if(z)j.ch(','); string key=j.str(); j.ch(':');
  if(key=="bin"){
   j.ch('{'); for(int q=0;q<3;q++){ if(q)j.ch(','); string k=j.str();j.ch(':'); if(k=="W")W=(int)j.num(); else if(k=="H")H=(int)j.num(); else canrot=j.boolean(); } j.ch('}');
  }else{
   j.ch('['); j.ws(); while(j.p<j.s.size()&&j.s[j.p]!=']'){
    if(!a.empty())j.ch(','); j.ch('{'); Item it;
    for(int q=0;q<5;q++){ if(q)j.ch(','); string k=j.str();j.ch(':'); if(k=="type")it.name=j.str(); else if(k=="w")it.w=(int)j.num(); else if(k=="h")it.h=(int)j.num(); else if(k=="v")it.v=j.num(); else it.lim=j.num(); }
    j.ch('}'); a.push_back(it); j.ws();
   } j.ch(']');
  }
 }
 int n=(int)a.size(); long double md=1,mv=1;
 for(auto &q:a){ md=max(md,(long double)q.v/(q.w*q.h)); mv=max(mv,(long double)q.v); }
 vector<P> best; long long bestval=-1;
 // This is deliberately an ablation of the old guillotine representation: every placement
 // cuts all intersected maximal empty rectangles, so L-shaped residual space remains usable.
 const long double bias[6]={0,.04L,.12L,.30L,.75L,.12L};
 // MaxRects has quadratic dominance cleanup. Bound its representation so all six
 // constructions remain usable under the one-second limit on high-limit inputs.
 constexpr int MAX_RAW_FREE=96;
 for(int run=0;run<6;run++){
  vector<R> free={{0,0,W,H}}; vector<long long> used(n); vector<P> out; long long value=0;
  while(!free.empty()){
   int bf=-1,bt=-1,bro=0,bw=0,bh=0; long double bscore=-1;
   for(int fi=0;fi<(int)free.size();fi++){
    const R& f=free[fi];
    for(int t=0;t<n;t++) if(used[t]<a[t].lim) for(int ro=0;ro<=(canrot&&a[t].w!=a[t].h);ro++){
     int w=ro?a[t].h:a[t].w, h=ro?a[t].w:a[t].h;
     if(w>f.w||h>f.h) continue;
     long double den=(long double)a[t].v/(w*h)/md;
     long double fill=(long double)(w*h)/(f.w*f.h);
     long double shortfit=(long double)min(f.w-w,f.h-h)/max(W,H);
     // Density remains dominant.  Runs vary only the fit pressure, exposing different layouts.
     long double score=den*(1+bias[run]*fill) + (run==5 ? .10L*(long double)a[t].v/mv : 0);
     score+=1e-7L*(1-shortfit);
     if(score>bscore){ bscore=score;bf=fi;bt=t;bro=ro;bw=w;bh=h; }
    }
   }
   if(bf<0) break;
   R placed={free[bf].x,free[bf].y,bw,bh};
   out.push_back({bt,placed.x,placed.y,bro,bw,bh}); used[bt]++; value+=a[bt].v;
   vector<R> next; next.reserve(free.size()+4);
   for(const R& f:free){
    if(!inter(f,placed)){ next.push_back(f); continue; }
    // Each generated rectangle is a subset of f that is disjoint from placed.
    if(placed.x>f.x) next.push_back({f.x,f.y,placed.x-f.x,f.h});
    if(placed.x+placed.w<f.x+f.w) next.push_back({placed.x+placed.w,f.y,f.x+f.w-(placed.x+placed.w),f.h});
    if(placed.y>f.y) next.push_back({f.x,f.y,f.w,placed.y-f.y});
    if(placed.y+placed.h<f.y+f.h) next.push_back({f.x,placed.y+placed.h,f.w,f.y+f.h-(placed.y+placed.h)});
   }
   // Keep a bounded, high-area subset before the quadratic dominance pass.  Every
   // retained rectangle is still empty, so this can affect quality but not validity.
   if((int)next.size()>MAX_RAW_FREE){
    nth_element(next.begin(),next.begin()+MAX_RAW_FREE,next.end(),[](const R& u,const R& v){
     return (long long)u.w*u.h>(long long)v.w*v.h;
    });
    next.resize(MAX_RAW_FREE);
   }
   // Overlap among maximal rectangles is intentional; discard only dominated ones.
   vector<char> dead(next.size());
   for(int i=0;i<(int)next.size();i++) if(!dead[i]) for(int k=0;k<(int)next.size();k++) if(i!=k&&!dead[i]&&inside(next[i],next[k])) dead[i]=1;
   free.clear(); free.reserve(next.size());
   for(int i=0;i<(int)next.size();i++) if(!dead[i]&&next[i].w>0&&next[i].h>0) free.push_back(next[i]);
  }
  if(value>bestval){ bestval=value; best.swap(out); }
 }
 cout<<"{\"placements\":[";
 for(size_t i=0;i<best.size();i++){ if(i)cout<<','; const P&p=best[i]; cout<<"{\"type\":\""<<a[p.t].name<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}'; }
 cout<<"]}\n";
}
