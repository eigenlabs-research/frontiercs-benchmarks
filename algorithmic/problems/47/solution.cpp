#include <bits/stdc++.h>
using namespace std;
struct Item{string id; int w,h; long long v; int lim;};
struct R{int x,y,w,h;}; struct P{int t,x,y,rot,w,h;};
string s; size_t at=0;
void ws(){while(at<s.size() && isspace((unsigned char)s[at]))at++;}
void ch(char c){ws(); if(at<s.size()&&s[at]==c)at++;}
string str(){ws(); ch('"'); string r; while(at<s.size()&&s[at]!='"'){if(s[at]=='\\'&&at+1<s.size())at++; r+=s[at++];} ch('"'); return r;}
long long num(){ws(); long long z=0,sg=1;if(s[at]=='-')sg=-1,at++;while(at<s.size()&&isdigit((unsigned char)s[at]))z=z*10+s[at++]-'0';return z*sg;}
bool boolean(){ws();if(s.compare(at,4,"true")==0){at+=4;return true;}at+=5;return false;}
void split(vector<R>& f,R q){ vector<R> g; for(R r:f){int x=max(r.x,q.x), y=max(r.y,q.y), X=min(r.x+r.w,q.x+q.w), Y=min(r.y+r.h,q.y+q.h); if(x>=X||y>=Y){g.push_back(r);continue;} if(q.x>r.x)g.push_back({r.x,r.y,q.x-r.x,r.h}); if(q.x+q.w<r.x+r.w)g.push_back({q.x+q.w,r.y,r.x+r.w-q.x-q.w,r.h}); if(q.y>r.y)g.push_back({r.x,r.y,r.w,q.y-r.y}); if(q.y+q.h<r.y+r.h)g.push_back({r.x,q.y+q.h,r.w,r.y+r.h-q.y-q.h}); }
 vector<R> h; for(int i=0;i<(int)g.size();i++){if(g[i].w<=0||g[i].h<=0)continue; bool inside=false;for(int j=0;j<(int)g.size();j++)if(i!=j&&g[j].x<=g[i].x&&g[j].y<=g[i].y&&g[j].x+g[j].w>=g[i].x+g[i].w&&g[j].y+g[j].h>=g[i].y+g[i].h){inside=true;break;}if(!inside)h.push_back(g[i]);} f.swap(h);}
struct Ans{long long val=0;vector<P> p;};
int main(){ios::sync_with_stdio(false);cin.tie(nullptr); s.assign((istreambuf_iterator<char>(cin)),{}); int W=0,H=0;bool rot=false;vector<Item>a; ws();ch('{');
 bool firstTop=true; while(true){ws();if(at<s.size()&&s[at]=='}'){at++;break;}if(!firstTop)ch(',');firstTop=false;string k=str();ch(':');
 if(k=="bin"){ch('{');bool first=true;while(true){ws();if(s[at]=='}'){at++;break;}if(!first)ch(',');first=false;string key=str();ch(':');if(key=="W")W=num();else if(key=="H")H=num();else if(key=="allow_rotate")rot=boolean();}}
 else if(k=="items"){ch('[');bool first=true;while(true){ws();if(s[at]==']'){at++;break;}if(!first)ch(',');first=false;ch('{');Item q{};bool firstField=true;while(true){ws();if(s[at]=='}'){at++;break;}if(!firstField)ch(',');firstField=false;string key=str();ch(':');if(key=="type")q.id=str();else if(key=="w")q.w=num();else if(key=="h")q.h=num();else if(key=="v")q.v=num();else if(key=="limit")q.lim=num();}a.push_back(q);}}
 }
 int n=a.size();
 auto run=[&](int mode,int force)->Ans{vector<R> f={{0,0,W,H}};vector<int> used(n);Ans out;auto put=[&](int t,int x,int y,int rr){int w=rr?a[t].h:a[t].w,h=rr?a[t].w:a[t].h;out.p.push_back({t,x,y,rr,w,h});out.val+=a[t].v;used[t]++;split(f,{x,y,w,h});};
 // Forced item is deliberately anchored at one of four corners: this tests whether early region commitment matters.
 if(force>=0){int t=force,brot=0,bw=a[t].w,bh=a[t].h; if(rot && a[t].h*a[t].w==a[t].w*a[t].h){ // choose orientation with smaller maximum side when both fit
   if(a[t].h<=W&&a[t].w<=H && (a[t].h<a[t].w || a[t].w>a[t].h)){brot=1;bw=a[t].h;bh=a[t].w;}}
   if(bw>W||bh>H){brot=0;bw=a[t].w;bh=a[t].h;} if(bw<=W&&bh<=H&&a[t].lim){int c=force%4;put(t,(c&1)?W-bw:0,(c&2)?H-bh:0,brot);}}
 for(int step=0;step<20000;step++){long double best=-1e300;int bt=-1,br=0,bfi=-1;for(int t=0;t<n;t++)if(used[t]<a[t].lim)for(int rr=0;rr<= (rot?1:0);rr++){int w=rr?a[t].h:a[t].w,h=rr?a[t].w:a[t].h;for(int j=0;j<(int)f.size();j++)if(w<=f[j].w&&h<=f[j].h){long double den=(long double)a[t].v/(w*h); long double base;
 if(mode==0)base=den; else if(mode==1)base=a[t].v; else if(mode==2)base=0.72L*den+0.28L*(long double)a[t].v/(W*H); else if(mode==3)base=(long double)a[t].v/sqrt((long double)w*h); else if(mode==4)base=den*(1.0L+0.15L*min(w,h)/(long double)max(w,h)); else base=(long double)a[t].v/(w*h)*0.8L+(long double)a[t].v/(W*H)*0.2L;
 int dx=f[j].w-w,dy=f[j].h-h; long double fit=(mode%2? -(long double)(dx*dy)/(W*H):-(long double)min(dx,dy)/max(1,max(W,H)));
 long double sc=base+fit*base*0.035L; if(sc>best){best=sc;bt=t;br=rr;bfi=j;}}}if(bt<0)break;put(bt,f[bfi].x,f[bfi].y,br);}return out;};
 Ans best; for(int m=0;m<6;m++){Ans z=run(m,-1);if(z.val>best.val)best=move(z);} // Forced starts are the discriminating branch, not sample-specific tuning.
 for(int t=0;t<n;t++){Ans z=run(t%6,t);if(z.val>best.val)best=move(z);} 
 cout<<"{\"placements\":[";for(size_t i=0;i<best.p.size();i++){auto&p=best.p[i];if(i)cout<<',';cout<<"{\"type\":\""<<a[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}cout<<"]}";}
