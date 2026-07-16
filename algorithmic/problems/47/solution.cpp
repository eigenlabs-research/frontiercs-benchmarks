#include <bits/stdc++.h>
using namespace std;
struct P { string s; size_t i=0; P(string z):s(move(z)){} void ws(){while(i<s.size()&&isspace((unsigned char)s[i]))i++;} char ch(){ws();return i<s.size()?s[i++]:0;} string str(){ws(); if(i>=s.size()||s[i++]!='"') return ""; string r; while(i<s.size()&&s[i]!='"'){if(s[i]=='\\'&&i+1<s.size())i++; r+=s[i++];} if(i<s.size())i++; return r;} long long num(){ws(); long long x=0,sg=1;if(s[i]=='-')sg=-1,i++;while(i<s.size()&&isdigit((unsigned char)s[i]))x=x*10+s[i++]-'0';return x*sg;} bool boolean(){ws();if(s.compare(i,4,"true")==0){i+=4;return true;}i+=5;return false;} };
struct Item{string id;int w,h,lim;long long v;};
struct R{int x,y,w,h;}; struct Pl{int t,x,y,rot,w,h;};
struct State{vector<R> f; vector<int> used; vector<Pl> p; long long val=0;};
static void esc(const string& s){ for(char c:s){if(c=='"'||c=='\\') cout<<'\\'; cout<<c;} }
int main(){
 ios::sync_with_stdio(false);cin.tie(nullptr);
 string in((istreambuf_iterator<char>(cin)),{}); P q(in); int W=0,H=0;bool allow=false;vector<Item>a;
 q.ch(); while(true){q.ws();if(q.i>=q.s.size()||q.s[q.i]=='}'){q.ch();break;}string k=q.str();q.ch(); if(k=="bin"){q.ch();while(true){q.ws();if(q.s[q.i]=='}'){q.i++;break;}string z=q.str();q.ch();if(z=="W")W=q.num();else if(z=="H")H=q.num();else allow=q.boolean();q.ws();if(q.s[q.i]==',')q.i++;}}else if(k=="items"){q.ch();while(true){q.ws();if(q.s[q.i]==']'){q.i++;break;}q.ch();Item e;while(true){q.ws();if(q.s[q.i]=='}'){q.i++;break;}string z=q.str();q.ch();if(z=="type")e.id=q.str();else if(z=="w")e.w=q.num();else if(z=="h")e.h=q.num();else if(z=="v")e.v=q.num();else e.lim=q.num();q.ws();if(q.s[q.i]==',')q.i++;}a.push_back(e);q.ws();if(q.s[q.i]==',')q.i++;}}q.ws();if(q.i<q.s.size()&&q.s[q.i]==',')q.i++;}
 long long maxv=1;for(auto&e:a)maxv=max(maxv,e.v); long double binarea=(long double)max(1,W)*H;
 State best;
 // Each run is a guillotine tree: leaves are disjoint rectangles, not overlapping MaxRects regions.
 for(int mode=0;mode<7;mode++){
  State s;s.f.push_back({0,0,W,H});s.used.assign(a.size(),0);
  for(int step=0;step<5000;step++){
   long double bs=-1;int bi=-1,bt=-1,bw=0,bh=0,br=0,bc=0,bd=0;
   for(int fi=0;fi<(int)s.f.size();fi++){R z=s.f[fi];
    for(int t=0;t<(int)a.size();t++)if(s.used[t]<a[t].lim&&a[t].v>0){
     for(int ro=0;ro<= (allow&&a[t].w!=a[t].h);ro++){int iw=ro?a[t].h:a[t].w,ih=ro?a[t].w:a[t].h;if(iw>z.w||ih>z.h)continue;
      long double ar=(long double)iw*ih, den=a[t].v/ar;
      long double fit=(long double)ar/(z.w*z.h), global=ar/binarea;
      int rw=z.w-iw,rh=z.h-ih; long double sliver=(long double)(min(rw,ih)+min(rh,iw))/(iw+ih+1.0L);
      long double score;
      if(mode==0) score=den;
      else if(mode==1) score=den*(1+0.55L*global);
      else if(mode==2) score=den*(1+0.8L*fit);
      else if(mode==3) score=den*(1-0.18L*sliver);
      else if(mode==4) score=den+0.13L*a[t].v/maxv;
      else if(mode==5) score=den*(1+0.25L*global-0.10L*sliver);
      else score=den*(1+0.35L*fit+0.20L*global);
      // A tiny deterministic tie perturbation makes runs choose genuinely different trees.
      score+=1e-12L*((t*17+fi*3+mode*11)%19);
      if(score>bs){bs=score;bi=fi;bt=t;bw=iw;bh=ih;br=ro;bc=(mode+fi+t)%4;bd=(mode/2+fi+t)%2;}
     }
    }
   }
   if(bi<0)break; R z=s.f[bi];s.f[bi]=s.f.back();s.f.pop_back();
   bool right=bc&1, top=bc&2; int px=right?z.x+z.w-bw:z.x, py=top?z.y+z.h-bh:z.y;
   s.p.push_back({bt,px,py,br,bw,bh});s.used[bt]++;s.val+=a[bt].v;
   int rw=z.w-bw,rh=z.h-bh;
   auto add=[&](int x,int y,int w,int h){if(w>0&&h>0)s.f.push_back({x,y,w,h});};
   if(bd==0){ // vertical cut first: a full-height side leaf and a leaf above/below the item
    if(!right)add(z.x+bw,z.y,rw,z.h); else add(z.x,z.y,rw,z.h);
    if(!top)add(px,z.y+bh,bw,rh); else add(px,z.y,bw,rh);
   }else{ // horizontal cut first: a full-width leaf and a side leaf beside the item
    if(!top)add(z.x,z.y+bh,z.w,rh); else add(z.x,z.y,z.w,rh);
    if(!right)add(z.x+bw,py,rw,bh); else add(z.x,py,rw,bh);
   }
  }
  if(s.val>best.val)best=move(s);
 }
 cout<<"{\"placements\":[";for(size_t i=0;i<best.p.size();i++){auto&x=best.p[i];if(i)cout<<',';cout<<"{\"type\":\"";esc(a[x.t].id);cout<<"\",\"x\":"<<x.x<<",\"y\":"<<x.y<<",\"rot\":"<<x.rot<<'}';}cout<<"]}\n";
}
