#include <bits/stdc++.h>
using namespace std;
struct It { string name; long long w,h,v,lim; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot,w,h; };

// Small JSON reader: input is deliberately restricted to the statement's JSON types.
struct JS {
 string s; size_t p=0;
 JS(){ s.assign(istreambuf_iterator<char>(cin),{}); }
 void ws(){while(p<s.size() && isspace((unsigned char)s[p]))p++;}
 void ch(char c){ws(); if(p<s.size()&&s[p]==c)p++;}
 string str(){ws(); string r; ch('"'); while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\'&&p+1<s.size())p++; r+=s[p++]; } if(p<s.size())p++; return r;}
 long long num(){ws(); long long z=0,sg=1;if(s[p]=='-')sg=-1,p++;while(p<s.size()&&isdigit((unsigned char)s[p]))z=z*10+s[p++]-'0';return z*sg;}
 bool boolean(){ws();bool b=s.compare(p,4,"true")==0;p+=b?4:5;return b;}
};
struct Pack {
 int W,H,n,mode; bool rotate; vector<It> a; vector<int> used; vector<R> fr; vector<P> out; long long val=0;
 Pack(int W,int H,bool ro,const vector<It>& a,int m):W(W),H(H),n(a.size()),mode(m),rotate(ro),a(a),used(n),fr{{0,0,W,H}}{}
 void prune(){
  vector<R> q; for(auto r:fr) if(r.w>0&&r.h>0) q.push_back(r); fr.swap(q);
  vector<char> bad(fr.size());
  for(int i=0;i<(int)fr.size();i++) for(int j=0;j<(int)fr.size();j++) if(i!=j && fr[i].x>=fr[j].x&&fr[i].y>=fr[j].y&&fr[i].x+fr[i].w<=fr[j].x+fr[j].w&&fr[i].y+fr[i].h<=fr[j].y+fr[j].h) {bad[i]=1;break;}
  q.clear(); for(int i=0;i<(int)fr.size();i++)if(!bad[i])q.push_back(fr[i]); fr.swap(q);
 }
 void put(P p){
  vector<R> q; int px=p.x,py=p.y,px2=px+p.w,py2=py+p.h;
  for(R f:fr){
   int fx2=f.x+f.w,fy2=f.y+f.h;
   if(px>=fx2||px2<=f.x||py>=fy2||py2<=f.y){q.push_back(f);continue;}
   if(py>f.y) q.push_back({f.x,f.y,f.w,py-f.y});
   if(py2<fy2) q.push_back({f.x,py2,f.w,fy2-py2});
   if(px>f.x) q.push_back({f.x,f.y,px-f.x,f.h});
   if(px2<fx2) q.push_back({px2,f.y,fx2-px2,f.h});
  }
  fr.swap(q); prune(); used[p.t]++; val+=a[p.t].v; out.push_back(p);
 }
 // Return a position for this particular item orientation.  The alternatives only affect
 // fragmentation tie-breaking; item choice remains density-first.
 bool position(int w,int h,int &bx,int &by,long long &u,long long &v){
  bool ok=false;
  for(R f:fr) if(w<=f.w&&h<=f.h){
   long long dx=f.w-w,dy=f.h-h, A=dx*dy+dx*h+dy*w; // free area not occupied in this host
   long long x,y, z1,z2;
   if(mode==0){z1=min(dx,dy);z2=max(dx,dy); x=f.x;y=f.y;}
   else if(mode==1){z1=A;z2=min(dx,dy);x=f.x;y=f.y;}
   else if(mode==2){z1=f.y;z2=f.x; x=f.x;y=f.y;}
   else {z1=max(dx,dy);z2=min(dx,dy);x=f.x;y=f.y;}
   if(!ok || z1<u || (z1==u&&(z2<v || (z2==v&&(y<by||(y==by&&x<bx)))))){ok=true;u=z1;v=z2;bx=x;by=y;}
  }
  return ok;
 }
 void run(){
  while(true){
   bool got=false; P best{}; long long bu=0,bv=0; // density comparison is done by cross multiplication
   for(int i=0;i<n;i++) if(used[i]<a[i].lim) for(int r=0;r<(rotate&&a[i].w!=a[i].h?2:1);r++){
    int w=r?a[i].h:a[i].w, h=r?a[i].w:a[i].h,bx,by; long long u,v;
    if(!position(w,h,bx,by,u,v))continue;
    if(!got){got=true;best={i,bx,by,r,w,h};bu=u;bv=v;continue;}
    __int128 lhs=(__int128)a[i].v*best.w*best.h, rhs=(__int128)a[best.t].v*w*h;
    // Prefer profit density; within equal density use a different maximal-rectangle fit rule.
    if(lhs>rhs || (lhs==rhs && (u<bu || (u==bu&&v<bv)))) {best={i,bx,by,r,w,h};bu=u;bv=v;}
   }
   if(!got)break; put(best);
  }
 }
};
int main(){
 ios::sync_with_stdio(false);cin.tie(nullptr);
 JS j; int W=0,H=0; bool rot=false; vector<It>a;
 j.ch('{'); while(true){j.ws();if(j.p>=j.s.size()||j.s[j.p]=='}'){j.p++;break;} string k=j.str();j.ch(':');
  if(k=="bin"){j.ch('{');while(true){j.ws();if(j.s[j.p]=='}'){j.p++;break;}string q=j.str();j.ch(':');if(q=="W")W=j.num();else if(q=="H")H=j.num();else rot=j.boolean();j.ws();if(j.s[j.p]==',')j.p++;}}
  else {j.ch('[');while(true){j.ws();if(j.s[j.p]==']'){j.p++;break;}j.ch('{');It z;while(true){j.ws();if(j.s[j.p]=='}'){j.p++;break;}string q=j.str();j.ch(':');if(q=="type")z.name=j.str();else if(q=="w")z.w=j.num();else if(q=="h")z.h=j.num();else if(q=="v")z.v=j.num();else z.lim=j.num();j.ws();if(j.s[j.p]==',')j.p++;}a.push_back(z);j.ws();if(j.s[j.p]==',')j.p++;}}
  j.ws();if(j.p<j.s.size()&&j.s[j.p]==',')j.p++;
 }
 Pack *ans=nullptr; for(int m=0;m<4;m++){Pack *q=new Pack(W,H,rot,a,m);q->run();if(!ans||q->val>ans->val){delete ans;ans=q;}else delete q;}
 cout<<"{\"placements\":["; for(int i=0;i<(int)ans->out.size();i++){auto p=ans->out[i];if(i)cout<<',';cout<<"{\"type\":\""<<a[p.t].name<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}cout<<"]}\n";
 delete ans;
}
