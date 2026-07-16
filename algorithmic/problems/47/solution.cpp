#include <bits/stdc++.h>
using namespace std;
struct It { string id; int w,h,lim; long long v; };
struct P { int t,x,y,r; };
struct N { int x,y,w; };
struct JS {
 string s; size_t p=0; JS(){ s.assign(istreambuf_iterator<char>(cin),{}); }
 void ws(){while(p<s.size()&&isspace((unsigned char)s[p]))p++;} void ex(char c){ws(); if(p<s.size()&&s[p]==c)p++;}
 string str(){ws(); ex('"'); string r; while(p<s.size()&&s[p]!='"'){ if(s[p]=='\\'&&p+1<s.size())p++; r+=s[p++]; } ex('"'); return r;}
 long long num(){ws(); long long z=0,sg=1;if(s[p]=='-')sg=-1,p++;while(p<s.size()&&isdigit((unsigned char)s[p]))z=z*10+s[p++]-'0';return z*sg;}
 bool boolean(){ws(); bool r=s[p]=='t'; p+=r?4:5; return r;}
};
struct Pack {
 int W,H; bool rot; vector<It>& a; vector<int> used; vector<N> sky; vector<P> out; long long val=0;
 Pack(int W,int H,bool r,vector<It>& a):W(W),H(H),rot(r),a(a),used(a.size()),sky{{0,0,W}}{}
 bool pos(int w,int h,int &bx,int &by){
  int bestY=INT_MAX,bestX=INT_MAX;
  for(int k=0;k<(int)sky.size();k++) { int x=sky[k].x; if(x+w>W) continue; int y=0, end=x+w;
   for(auto q:sky) if(q.x<end && q.x+q.w>x) y=max(y,q.y);
   if(y+h<=H && (y<bestY || (y==bestY&&x<bestX))) bestY=y,bestX=x;
  }
  if(bestY==INT_MAX)return false; bx=bestX;by=bestY;return true;
 }
 void put(int t,int x,int y,int w,int h,int r){
  vector<N> z; int l=x,rr=x+w;
  for(auto q:sky){ int e=q.x+q.w; if(e<=l||q.x>=rr)z.push_back(q); else {if(q.x<l)z.push_back({q.x,q.y,l-q.x});if(e>rr)z.push_back({rr,q.y,e-rr});} }
  z.push_back({l,y+h,w}); sort(z.begin(),z.end(),[](N A,N B){return A.x<B.x;});
  sky.clear(); for(auto q:z){if(q.w<=0)continue; if(!sky.empty()&&sky.back().y==q.y&&sky.back().x+sky.back().w==q.x)sky.back().w+=q.w;else sky.push_back(q);}
  used[t]++;out.push_back({t,x,y,r});val+=a[t].v;
 }
 bool bestFor(int t,int &x,int &y,int &r,int orientMode){
  int bx,by; bool ok=false; int bw=0,bh=0;
  for(int q=0;q<=(rot&&a[t].w!=a[t].h);q++){int w=q?a[t].h:a[t].w,h=q?a[t].w:a[t].h,xx,yy;if(!pos(w,h,xx,yy))continue;
   // Alternate skyline runs choose either lowest base or lowest resulting top contour.
   if(!ok || (orientMode==0 ? make_tuple(yy,xx,yy+h) : make_tuple(yy+h,yy,xx)) < (orientMode==0 ? make_tuple(by,bx,by+bh) : make_tuple(by+bh,by,bx))){ok=true;bx=xx;by=yy;r=q;bw=w;bh=h;}
  }
  x=bx;y=by;return ok;
 }
};
int main(){
 JS j; int W=0,H=0;bool allow=false;vector<It>a; j.ex('{');
 for(int top=0;top<2;top++){if(top)j.ex(',');string key=j.str();j.ex(':');if(key=="bin"){j.ex('{');for(int k=0;k<3;k++){if(k)j.ex(',');string q=j.str();j.ex(':');if(q=="W")W=j.num();else if(q=="H")H=j.num();else if(q=="allow_rotate")allow=j.boolean();}j.ex('}');}
 else {j.ex('[');bool first=true;while(true){j.ws();if(j.p<j.s.size()&&j.s[j.p]==']'){j.p++;break;}if(!first)j.ex(',');first=false;j.ex('{');It it;for(int k=0;k<5;k++){if(k)j.ex(',');string q=j.str();j.ex(':');if(q=="type")it.id=j.str();else if(q=="w")it.w=j.num();else if(q=="h")it.h=j.num();else if(q=="v")it.v=j.num();else if(q=="limit")it.lim=j.num();}j.ex('}');a.push_back(it);}}
 } j.ex('}');
 vector<int> base(a.size());iota(base.begin(),base.end(),0); vector<P> answer;long long best=-1;
 // Distinct static priority policies plus two contour-aware selection policies.
 for(int mode=0;mode<8;mode++){
  vector<int> ord=base;
  sort(ord.begin(),ord.end(),[&](int i,int k){auto &A=a[i],&B=a[k]; long double ka,kb;
   if(mode==0)ka=(long double)A.v/A.w/A.h,kb=(long double)B.v/B.w/B.h;
   else if(mode==1)ka=A.v,kb=B.v;
   else if(mode==2)ka=(long double)A.v/min(A.w,A.h),kb=(long double)B.v/min(B.w,B.h);
   else if(mode==3)ka=(long double)A.v/max(A.w,A.h),kb=(long double)B.v/max(B.w,B.h);
   else if(mode==4)ka=-(long double)A.w*A.h,kb=-(long double)B.w*B.h;
   else if(mode==5)ka=(long double)A.v/(A.w+A.h),kb=(long double)B.v/(B.w+B.h);
   else if(mode==6)ka=(long double)A.v/A.h,kb=(long double)B.v/B.h;
   else ka=(long double)A.v/A.w,kb=(long double)B.v/B.w;
   if(ka!=kb)return ka>kb;return A.v>B.v;});
  // cyclic shifts avoid every run giving the same early item type absolute priority.
  if(mode>=4&&ord.size())rotate(ord.begin(),ord.begin()+(mode-4)%ord.size(),ord.end());
  Pack p(W,H,allow,a); int cap=6000;
  for(int t:ord) while(p.used[t]<a[t].lim && (int)p.out.size()<cap){int x,y,r;if(!p.bestFor(t,x,y,r,mode&1))break;int ww=r?a[t].h:a[t].w,hh=r?a[t].w:a[t].h;p.put(t,x,y,ww,hh,r);}
  if(p.val>best)best=p.val,answer=p.out;
 }
 // Lowest-contour choice interleaves types; it explicitly tests whether gap filling, not type order, is decisive.
 for(int om=0;om<2;om++){Pack p(W,H,allow,a);int cap=5000;while((int)p.out.size()<cap){int bt=-1,bx=0,by=0,br=0;long double bs=-1;for(int t=0;t<(int)a.size();t++)if(p.used[t]<a[t].lim){int x,y,r;if(!p.bestFor(t,x,y,r,om))continue; long double den=(long double)a[t].v/a[t].w/a[t].h; long double sc=den; if(bt<0 || (om==0?make_tuple(y,x,-sc):make_tuple(y+(r?a[t].w:a[t].h),y,-sc)) < (om==0?make_tuple(by,bx,-bs):make_tuple(by+(br?a[bt].w:a[bt].h),by,-bs))){bt=t;bx=x;by=y;br=r;bs=sc;}}if(bt<0)break;int ww=br?a[bt].h:a[bt].w,hh=br?a[bt].w:a[bt].h;p.put(bt,bx,by,ww,hh,br);}if(p.val>best)best=p.val,answer=p.out;}
 cout<<"{\"placements\":[";for(int i=0;i<(int)answer.size();i++){if(i)cout<<',';auto q=answer[i];cout<<"{\"type\":\""<<a[q.t].id<<"\",\"x\":"<<q.x<<",\"y\":"<<q.y<<",\"rot\":"<<q.r<<'}';}cout<<"]}";
}
