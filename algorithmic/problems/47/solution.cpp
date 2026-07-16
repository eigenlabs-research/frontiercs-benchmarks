#include <bits/stdc++.h>
using namespace std;
struct Item{string id; long long w,h,v,lim;};
struct R{long long x,y,w,h;};
struct P{int t,rot; long long x,y,w,h;};
struct J{
 string s; size_t p=0; J(){s.assign(istreambuf_iterator<char>(cin),{});} void ws(){while(p<s.size()&&isspace((unsigned char)s[p]))p++;} char ch(){ws();return p<s.size()?s[p++]:0;} string str(){ws(); if(p<s.size()&&s[p]=='"')p++; string r; while(p<s.size()&&s[p]!='"'){if(s[p]=='\\'&&p+1<s.size())p++;r+=s[p++];}if(p<s.size())p++;return r;} long long num(){ws();long long z=1,r=0;if(s[p]=='-')z=-1,p++;while(p<s.size()&&isdigit((unsigned char)s[p]))r=r*10+s[p++]-'0';return z*r;} bool boolean(){ws();bool r=s.compare(p,4,"true")==0;p+=r?4:5;return r;}
};
static bool inter(const R&a,const R&b){return a.x<b.x+b.w&&b.x<a.x+a.w&&a.y<b.y+b.h&&b.y<a.y+a.h;}
struct Pack{
 long long W,H; vector<Item>& a; bool allow; vector<R> fr; vector<P> q; vector<long long> used;
 Pack(long long W,long long H,vector<Item>&a,bool al):W(W),H(H),a(a),allow(al),fr{{0,0,W,H}},used(a.size()){}
 void cut(R z){ vector<R> n; for(R f:fr){ if(!inter(f,z)){n.push_back(f);continue;} if(z.x>f.x)n.push_back({f.x,f.y,z.x-f.x,f.h}); if(z.x+z.w<f.x+f.w)n.push_back({z.x+z.w,f.y,f.x+f.w-z.x-z.w,f.h}); if(z.y>f.y)n.push_back({f.x,f.y,f.w,z.y-f.y}); if(z.y+z.h<f.y+f.h)n.push_back({f.x,z.y+z.h,f.w,f.y+f.h-z.y-z.h}); }
  vector<R> k; for(int i=0;i<(int)n.size();i++){if(n[i].w<=0||n[i].h<=0)continue;bool in=false;for(int j=0;j<(int)n.size();j++)if(i!=j&&n[j].x<=n[i].x&&n[j].y<=n[i].y&&n[j].x+n[j].w>=n[i].x+n[i].w&&n[j].y+n[j].h>=n[i].y+n[i].h){in=true;break;}if(!in)k.push_back(n[i]);} fr.swap(k);
 }
 bool one(int t,int mode){ if(used[t]>=a[t].lim)return false; long long bw=0,bh=0;int br=0,bi=-1; pair<long long,long long> best={LLONG_MAX,LLONG_MAX};
  for(int r=0;r<=(allow&&a[t].w!=a[t].h);r++){long long w=r?a[t].h:a[t].w,h=r?a[t].w:a[t].h;for(int i=0;i<(int)fr.size();i++)if(w<=fr[i].w&&h<=fr[i].h){long long dx=fr[i].w-w,dy=fr[i].h-h;pair<long long,long long> z; if(mode==0)z={min(dx,dy),max(dx,dy)};else if(mode==1)z={dx*dy,min(dx,dy)};else z={fr[i].y+ h,fr[i].x+w}; if(z<best){best=z;bi=i;bw=w;bh=h;br=r;}}}
  if(bi<0)return false; R z{fr[bi].x,fr[bi].y,bw,bh};q.push_back({t,br,z.x,z.y,bw,bh});used[t]++;cut(z);return true;
 }
 long long val(){long long z=0;for(auto&p:q)z+=a[p.t].v;return z;}
};
int main(){J j; long long W=0,H=0;bool al=false;vector<Item>a; j.ch();for(int top=0;top<2;top++){if(top)j.ch();string key=j.str();j.ch();if(key=="bin"){j.ch();for(int k=0;k<3;k++){if(k)j.ch();string x=j.str();j.ch();if(x=="W")W=j.num();else if(x=="H")H=j.num();else al=j.boolean();}j.ch();}else{j.ch();bool first=true;while(1){j.ws();if(j.p<j.s.size()&&j.s[j.p]==']'){j.p++;break;}if(!first)j.ch();first=false;j.ch();Item it;for(int k=0;k<5;k++){if(k)j.ch();string x=j.str();j.ch();if(x=="type")it.id=j.str();else if(x=="w")it.w=j.num();else if(x=="h")it.h=j.num();else if(x=="v")it.v=j.num();else it.lim=j.num();}j.ch();a.push_back(it);}}}
 vector<P> ans;long long bv=-1; int n=a.size();
 for(int run=0;run<8;run++){vector<int> ord(n);iota(ord.begin(),ord.end(),0);sort(ord.begin(),ord.end(),[&](int x,int y){__int128 X=(__int128)a[x].v*a[y].w*a[y].h,Y=(__int128)a[y].v*a[x].w*a[x].h;if(X!=Y)return X>Y;if(run%3==0)return a[x].v>a[y].v;if(run%3==1)return a[x].w*a[x].h>a[y].w*a[y].h;return max(a[x].w,a[x].h)>max(a[y].w,a[y].h);}); if(run>=3) rotate(ord.begin(),ord.begin()+(run-2)%n,ord.end());Pack z(W,H,a,al);for(int t:ord)while(z.one(t,run%3));if(z.val()>bv){bv=z.val();ans=z.q;}}
 auto esc=[](const string& s){string r;for(char c:s){if(c=='"'||c=='\\')r+='\\';r+=c;}return r;};
 cout<<"{\"placements\":[";for(size_t i=0;i<ans.size();i++){if(i)cout<<',';auto&p=ans[i];cout<<"{\"type\":\""<<esc(a[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}cout<<"]}\n";
}
