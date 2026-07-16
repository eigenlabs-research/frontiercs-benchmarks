#include <bits/stdc++.h>
using namespace std;
struct Item { string id; int w,h; long long v,lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot; };

struct Json {
 string s; size_t p=0;
 Json(){ s.assign(istreambuf_iterator<char>(cin),{}); }
 void ws(){ while(p<s.size() && isspace((unsigned char)s[p])) ++p; }
 void ch(char c){ ws(); if(p<s.size() && s[p]==c) ++p; }
 string str(){ ws(); string r; if(p>=s.size()||s[p++]!='"') return r; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\'&&p+1<s.size()) ++p; r+=s[p++]; } if(p<s.size())++p; return r; }
 long long num(){ ws(); int z=1; if(s[p]=='-')z=-1,++p; long long r=0; while(p<s.size()&&isdigit((unsigned char)s[p]))r=r*10+s[p++]-'0'; return z*r; }
 bool boolean(){ ws(); bool r=p<s.size()&&s[p]=='t'; while(p<s.size()&&isalpha((unsigned char)s[p]))++p; return r; }
};
static bool inter(const R&a,const R&b){ return a.x<b.x+b.w&&b.x<a.x+a.w&&a.y<b.y+b.h&&b.y<a.y+a.h; }

vector<P> pack(const vector<Item>& a,int W,int H,bool rotok, vector<int> ord,int mode){
 vector<R> f(1,{0,0,W,H}); vector<P> ans; vector<int> used(a.size());
 const int MAXP=1800, KEEP=80;
 for(int t:ord){
  for(int rep=0;rep<a[t].lim && (int)ans.size()<MAXP;rep++){
   int bi=-1, br=0; long long bs=LLONG_MAX, bt=LLONG_MAX;
   for(int i=0;i<(int)f.size();i++) for(int q=0;q<(rotok?2:1);q++){
    int w=q?a[t].h:a[t].w, h=q?a[t].w:a[t].h;
    if(w>f[i].w||h>f[i].h) continue;
    long long dx=f[i].w-w, dy=f[i].h-h, u=min(dx,dy), v=max(dx,dy), area=dx*dy;
    long long s1,s2;
    if(mode==0){s1=u;s2=v;} else if(mode==1){s1=area;s2=u;}
    else if(mode==2){s1=dy;s2=dx;} else {s1=dx;s2=dy;}
    // tie preference deliberately varies orientation, while every choice stays in one free box
    if(s1<bs || (s1==bs && (s2<bt || (s2==bt && ((q^mode)&1))))){bs=s1;bt=s2;bi=i;br=q;}
   }
   if(bi<0) break;
   int w=br?a[t].h:a[t].w, h=br?a[t].w:a[t].h;
   R z{f[bi].x,f[bi].y,w,h}; ans.push_back({t,z.x,z.y,br}); used[t]++;
   vector<R> nf; nf.reserve(f.size()*2+4);
   for(const R& e:f){
    if(!inter(e,z)){nf.push_back(e);continue;}
    if(z.x>e.x) nf.push_back({e.x,e.y,z.x-e.x,e.h});
    if(z.x+z.w<e.x+e.w) nf.push_back({z.x+z.w,e.y,e.x+e.w-(z.x+z.w),e.h});
    if(z.y>e.y) nf.push_back({e.x,e.y,e.w,z.y-e.y});
    if(z.y+z.h<e.y+e.h) nf.push_back({e.x,z.y+z.h,e.w,e.y+e.h-(z.y+z.h)});
   }
   nf.erase(remove_if(nf.begin(),nf.end(),[](const R&r){return r.w<=0||r.h<=0;}),nf.end());
   vector<char> dead(nf.size());
   for(int i=0;i<(int)nf.size();i++) if(!dead[i]) for(int j=0;j<(int)nf.size();j++) if(i!=j&&!dead[i] && nf[j].x<=nf[i].x&&nf[j].y<=nf[i].y&&nf[j].x+nf[j].w>=nf[i].x+nf[i].w&&nf[j].y+nf[j].h>=nf[i].y+nf[i].h) dead[i]=1;
   f.clear(); for(int i=0;i<(int)nf.size();i++) if(!dead[i]) f.push_back(nf[i]);
   if((int)f.size()>KEEP){ nth_element(f.begin(),f.begin()+KEEP,f.end(),[](const R&A,const R&B){return 1LL*A.w*A.h>1LL*B.w*B.h;}); f.resize(KEEP); }
  }
 }
 return ans;
}
int main(){
 Json j; j.ch('{'); int W=0,H=0; bool allow=false; vector<Item>a;
 for(int top=0;top<2;top++){
  if(top)j.ch(','); string key=j.str();j.ch(':');
  if(key=="bin"){
   j.ch('{');for(int k=0;k<3;k++){if(k)j.ch(',');string q=j.str();j.ch(':');if(q=="W")W=j.num();else if(q=="H")H=j.num();else allow=j.boolean();}j.ch('}');
  }else{
   j.ch('['); bool first=true; while(true){j.ws();if(j.p<j.s.size()&&j.s[j.p]==']'){j.p++;break;}if(!first)j.ch(',');first=false;j.ch('{');Item x;for(int k=0;k<5;k++){if(k)j.ch(',');string q=j.str();j.ch(':');if(q=="type")x.id=j.str();else if(q=="w")x.w=j.num();else if(q=="h")x.h=j.num();else if(q=="v")x.v=j.num();else x.lim=j.num();}j.ch('}');a.push_back(x);}
  }
 } j.ch('}');
 vector<P> best; long long bv=-1; int n=a.size();
 for(int style=0;style<6;style++){
  vector<int> o(n);iota(o.begin(),o.end(),0);
  sort(o.begin(),o.end(),[&](int x,int y){
   long double ax=(long double)a[x].v/(a[x].w*a[x].h), ay=(long double)a[y].v/(a[y].w*a[y].h);
   if(style==0) return ax>ay;
   if(style==1) return a[x].v>a[y].v;
   if(style==2) return a[x].w*a[x].h>a[y].w*a[y].h;
   if(style==3) return a[x].h>a[y].h;
   if(style==4) return a[x].w>a[y].w;
   return ax*(1+min(a[x].lim,50LL)/100.0L)>ay*(1+min(a[y].lim,50LL)/100.0L);
  });
  // Eight bounded starts: density/value get three fit rules; shape starts are ablations.
  int modes = style<2 ? 3 : (style==2 || style==5 ? 1 : 0);
  for(int mode=0;mode<modes;mode++){
   auto z=pack(a,W,H,allow,o,mode);long long v=0;for(auto&p:z)v+=a[p.t].v;
   if(v>bv)bv=v,best.swap(z);
  }
 }
 cout<<"{\"placements\":[";
 for(int i=0;i<(int)best.size();i++){if(i)cout<<',';auto&p=best[i];cout<<"{\"type\":\""<<a[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}
 cout<<"]}";
}
