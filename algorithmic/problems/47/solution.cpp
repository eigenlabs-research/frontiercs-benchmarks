#include <bits/stdc++.h>
using namespace std;
struct Item{string id;int w,h,lim;long long v;};
struct R{int x,y,w,h;}; struct P{int t,x,y,rot,w,h;};
struct Json{
 string s; size_t p=0; Json(){s.assign(istreambuf_iterator<char>(cin),{});} void ws(){while(p<s.size()&&isspace((unsigned char)s[p]))p++;} void ch(char c){ws();if(p<s.size()&&s[p]==c)p++;}
 string str(){ws();ch('"');string r;while(p<s.size()&&s[p]!='"'){if(s[p]=='\\'&&p+1<s.size())p++;r+=s[p++];}ch('"');return r;} long long num(){ws();int z=1;if(s[p]=='-')z=-1,p++;long long r=0;while(p<s.size()&&isdigit((unsigned char)s[p]))r=r*10+s[p++]-'0';return z*r;} bool boolean(){ws();bool r=s.compare(p,4,"true")==0;p+=r?4:5;return r;}
};
static bool inside(const R&a,const R&b){return a.x<=b.x&&a.y<=b.y&&a.x+a.w>=b.x+b.w&&a.y+a.h>=b.y+b.h;}
static bool hit(const R&a,const R&b){return a.x<b.x+b.w&&b.x<a.x+a.w&&a.y<b.y+b.h&&b.y<a.y+a.h;}
struct MR{
 int W,H,heur; vector<Item>* a; vector<R> f; vector<int> used; vector<P> out; long long val=0;
 MR(int W,int H,vector<Item>*a,int q):W(W),H(H),heur(q),a(a),f{{0,0,W,H}},used(a->size()){}
 void add(R q){if(q.w<=0||q.h<=0)return;for(int i=0;i<(int)f.size();){if(inside(f[i],q))return;if(inside(q,f[i]))f[i]=f.back(),f.pop_back();else i++;}f.push_back(q);}
 void split(R q){vector<R> old;old.swap(f);for(R z:old){if(!hit(z,q)){add(z);continue;}add({z.x,z.y,q.x-z.x,z.h});add({q.x+q.w,z.y,z.x+z.w-q.x-q.w,z.h});add({z.x,z.y,z.w,q.y-z.y});add({z.x,q.y+q.h,z.w,z.y+z.h-q.y-q.h});}}
 bool put(int t,int rot){Item&z=(*a)[t];int w=rot?z.h:z.w,h=rot?z.w:z.h,bi=-1,u=INT_MAX,v=INT_MAX;long long ba=LLONG_MAX;for(int i=0;i<(int)f.size();i++)if(f[i].w>=w&&f[i].h>=h){int x=f[i].w-w,y=f[i].h-h;int mn=min(x,y),mx=max(x,y);long long ar=1LL*x*y;bool take=heur==0?make_pair(mn,mx)<make_pair(u,v):heur==1?ar<ba:make_pair(mx,mn)<make_pair(v,u);if(take)bi=i,u=mn,v=mx,ba=ar;}if(bi<0)return false;R q{f[bi].x,f[bi].y,w,h};split(q);out.push_back({t,q.x,q.y,rot,w,h});used[t]++;val+=z.v;return true;}
};
// A separate guillotine representation: free rectangles are a disjoint partition, not overlapping maximal rectangles.
struct GP{
 vector<Item>*a; bool allow; int mode; vector<R> leaf; vector<int> used; vector<P> out; long long val=0;
 GP(int W,int H,vector<Item>*q,bool al,int m):a(q),allow(al),mode(m),leaf{{0,0,W,H}},used(q->size()){}
 void run(){
  int n=a->size();
  for(int step=0;step<1200;step++){
   int bl=-1,bt=-1,br=0; long double bs=-1;
   for(int i=0;i<(int)leaf.size();i++) for(int t=0;t<n;t++) if(used[t]<(*a)[t].lim) for(int r=0;r<= (allow&&(*a)[t].w!=(*a)[t].h);r++){
    int w=r?(*a)[t].h:(*a)[t].w,h=r?(*a)[t].w:(*a)[t].h;if(w>leaf[i].w||h>leaf[i].h)continue;
    long double den=(long double)(*a)[t].v/(w*h), fill=(long double)(w*h)/(leaf[i].w*leaf[i].h);
    // Mildly prefer a nearly filled partition; mode changes only deterministic tie/shape pressure.
    long double shape=(mode==0?min(leaf[i].w-w,leaf[i].h-h):mode==1?max(leaf[i].w-w,leaf[i].h-h):abs((leaf[i].w-w)-(leaf[i].h-h)));
    long double jitter=1.0L+(((t*29+i*7+mode*11)%17)-8)*0.006L;
    long double sc=den*jitter*(1.0L+0.18L*sqrtl(fill))/(1.0L+0.00012L*shape);
    if(sc>bs){bs=sc;bl=i;bt=t;br=r;}
   }
   if(bl<0)break;
   R z=leaf[bl];Item&q=(*a)[bt];int w=br?q.h:q.w,h=br?q.w:q.h;
   out.push_back({bt,z.x,z.y,br,w,h});used[bt]++;val+=q.v;
   // Two guillotine cuts; alternate their order to expose horizontally and vertically layered layouts.
   R one,two;if((mode&1)==0){one={z.x+w,z.y,z.w-w,h};two={z.x,z.y+h,z.w,z.h-h};}else{one={z.x+w,z.y,z.w-w,z.h};two={z.x,z.y+h,w,z.h-h};}
   leaf[bl]=leaf.back();leaf.pop_back();if(one.w>0&&one.h>0)leaf.push_back(one);if(two.w>0&&two.h>0)leaf.push_back(two);
  }
 }
};
int main(){
 Json j;j.ch('{');int W=0,H=0;bool allow=false;vector<Item>a; bool firstTop=true;
 // JSON object-member order is not part of the input contract.
 while(true){j.ws();if(j.p<j.s.size()&&j.s[j.p]=='}'){j.p++;break;}if(!firstTop)j.ch(',');firstTop=false;string k=j.str();j.ch(':');
  if(k=="bin"){j.ch('{');bool first=true;while(true){j.ws();if(j.p<j.s.size()&&j.s[j.p]=='}'){j.p++;break;}if(!first)j.ch(',');first=false;string q=j.str();j.ch(':');if(q=="W")W=j.num();else if(q=="H")H=j.num();else if(q=="allow_rotate")allow=j.boolean();}}
  else if(k=="items"){j.ch('[');bool first=1;while(true){j.ws();if(j.p<j.s.size()&&j.s[j.p]==']'){j.p++;break;}if(!first)j.ch(',');first=0;j.ch('{');Item x{};bool fieldFirst=true;while(true){j.ws();if(j.p<j.s.size()&&j.s[j.p]=='}'){j.p++;break;}if(!fieldFirst)j.ch(',');fieldFirst=false;string q=j.str();j.ch(':');if(q=="type")x.id=j.str();else if(q=="w")x.w=j.num();else if(q=="h")x.h=j.num();else if(q=="v")x.v=j.num();else if(q=="limit")x.lim=j.num();}a.push_back(x);}}
 }
 vector<P>best;long long bv=-1;int n=a.size();
 // Incumbent maximal-rectangle family retained as a safe fallback.
 for(int mode=0;mode<10;mode++){vector<int>ord(n);iota(ord.begin(),ord.end(),0);sort(ord.begin(),ord.end(),[&](int i,int k){long double x=(long double)a[i].v/(a[i].w*a[i].h),y=(long double)a[k].v/(a[k].w*a[k].h);if(mode==1)x=a[i].v,y=a[k].v;if(mode>=2){x*=.88L+((i*37+mode*19)%23)/80.0L;y*=.88L+((k*37+mode*19)%23)/80.0L;}return fabsl(x-y)>1e-12?x>y:a[i].v>a[k].v;});MR q(W,H,&a,mode%3);for(int t:ord){int r=0;if(allow&&a[t].w!=a[t].h){bool tall=mode<2?mode==1:((t*17+mode*13)&1)==0;r=tall?(a[t].h<a[t].w):(a[t].h>a[t].w);}while(q.used[t]<a[t].lim&&q.out.size()<3000)if(!q.put(t,r))break;}if(q.val>bv)bv=q.val,best=q.out;}
 // Diversification branch: guillotine leaf packings with per-copy rotation competition.
 for(int m=0;m<2;m++){GP q(W,H,&a,allow,m);q.run();if(q.val>bv)bv=q.val,best=q.out;}
 cout<<"{\"placements\":[";for(size_t i=0;i<best.size();i++){if(i)cout<<',';P&p=best[i];cout<<"{\"type\":\""<<a[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}cout<<"]}\n";
}
