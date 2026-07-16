#include <bits/stdc++.h>
using namespace std;
struct Item { string id; int w,h,lim; long long v; };
struct R { int x,y,w,h; };
struct P { int t,x,y,rot; };
int W,H; bool canrot; vector<Item> a;

static string fieldString(const string& s,const string& k){
    regex q("\\\""+k+"\\\"\\s*:\\s*\\\"([^\\\"]*)\\\""); smatch m;
    return regex_search(s,m,q)?m[1].str():"";
}
static long long fieldNum(const string& s,const string& k){
    regex q("\\\""+k+"\\\"\\s*:\\s*(-?[0-9]+)"); smatch m;
    return regex_search(s,m,q)?stoll(m[1].str()):0;
}
static bool fieldBool(const string& s,const string& k){
    regex q("\\\""+k+"\\\"\\s*:\\s*(true|false)"); smatch m;
    return regex_search(s,m,q) && m[1]=="true";
}
static int closeBrace(const string&s,int p){ int d=0; bool quote=false,esc=false; for(int i=p;i<(int)s.size();i++){char c=s[i]; if(quote){if(esc)esc=false; else if(c=='\\')esc=true; else if(c=='\"')quote=false;} else {if(c=='\"')quote=true; else if(c=='{')d++; else if(c=='}'&&!--d)return i;}} return -1; }
static bool contains(const R& u,const R& v){return v.x>=u.x&&v.y>=u.y&&v.x+v.w<=u.x+u.w&&v.y+v.h<=u.y+u.h;}

struct Pack {
 vector<R> f; vector<P> p; vector<int> used; long long value=0;
 Pack():f(1,{0,0,W,H}),used(a.size()){}
 void prune(){
   vector<char> bad(f.size());
   for(int i=0;i<(int)f.size();i++) for(int j=0;j<(int)f.size();j++) if(i!=j && contains(f[j],f[i])) {bad[i]=1;break;}
   vector<R> z; z.reserve(f.size()); for(int i=0;i<(int)f.size();i++) if(!bad[i]&&f[i].w>0&&f[i].h>0) z.push_back(f[i]); f.swap(z);
 }
 void put(int t,int fi,int rw,int rh,int rot){
   R q=f[fi]; int x=q.x,y=q.y; vector<R> z; z.reserve(f.size()+4);
   for(R r:f){
     int ix=max(r.x,x), iy=max(r.y,y), ax=min(r.x+r.w,x+rw), ay=min(r.y+r.h,y+rh);
     if(ix>=ax||iy>=ay) z.push_back(r);
     else { if(r.x<x)z.push_back({r.x,r.y,x-r.x,r.h}); if(x+rw<r.x+r.w)z.push_back({x+rw,r.y,r.x+r.w-x-rw,r.h}); if(r.y<y)z.push_back({r.x,r.y,r.w,y-r.y}); if(y+rh<r.y+r.h)z.push_back({r.x,y+rh,r.w,r.y+r.h-y-rh}); }
   }
   f.swap(z); prune(); p.push_back({t,x,y,rot}); used[t]++; value+=a[t].v;
 }
};

int main(){
 ios::sync_with_stdio(false);cin.tie(nullptr);
 string s((istreambuf_iterator<char>(cin)),{});
 size_t bp=s.find("\"bin\""); int b0=s.find('{',bp), b1=closeBrace(s,b0); string bin=s.substr(b0,b1-b0+1);
 W=fieldNum(bin,"W"); H=fieldNum(bin,"H"); canrot=fieldBool(bin,"allow_rotate");
 size_t ip=s.find("\"items\""); int at=s.find('[',ip); bool q=false,esc=false; int end=at;
 for(;end<(int)s.size();end++){char c=s[end];if(q){if(esc)esc=false;else if(c=='\\')esc=true;else if(c=='\"')q=false;}else if(c=='\"')q=true;else if(c==']')break;}
 for(int i=at+1;i<end;){ if(s[i]=='{'){int j=closeBrace(s,i);string o=s.substr(i,j-i+1); a.push_back({fieldString(o,"type"),(int)fieldNum(o,"w"),(int)fieldNum(o,"h"),(int)fieldNum(o,"limit"),fieldNum(o,"v")});i=j+1;}else i++; }
 vector<P> best; long long bv=-1; int n=a.size();
 // The portfolios deliberately isolate ordering/selection effects while retaining one geometry engine.
 for(int mode=0;mode<12;mode++){
   vector<int> ord(n); iota(ord.begin(),ord.end(),0);
   auto key=[&](int i){double ar=1.0*a[i].w*a[i].h, d=a[i].v/ar; double asp=(double)max(a[i].w,a[i].h)/min(a[i].w,a[i].h);
     if(mode==0)return d; if(mode==1)return (double)a[i].v; if(mode==2)return d*sqrt(ar); if(mode==3)return d/sqrt(ar);
     if(mode==4)return d*(1+.12*asp); if(mode==5)return d/(1+.12*asp); if(mode==6)return (double)a[i].v/sqrt(ar);
     if(mode==7)return d+1e-8*a[i].v; if(mode==8)return d*(a[i].lim<80?1.15:1); if(mode==9)return (double)a[i].v/ar;
     if(mode==10)return ar; return -ar; };
   sort(ord.begin(),ord.end(),[&](int x,int y){double X=key(x),Y=key(y);if(X!=Y)return X>Y;return x<y;});
   Pack z;
   for(int t:ord) while(z.used[t]<a[t].lim){
     int bi=-1,bw=0,bh=0,br=0; long long bs=LLONG_MAX, ba=LLONG_MAX;
     for(int i=0;i<(int)z.f.size();i++) for(int r=0;r<(canrot?2:1);r++){
       int ww=r?a[t].h:a[t].w, hh=r?a[t].w:a[t].h; R F=z.f[i]; if(ww>F.w||hh>F.h)continue;
       long long ss=min(F.w-ww,F.h-hh), aa=1LL*(F.w-ww)*(F.h-hh);
       // Alternate fit rule and orientation ties, a controlled ablation of placement selection.
       long long primary=(mode&1)?aa:ss, secondary=(mode&1)?ss:aa;
       if(primary<bs||(primary==bs&&secondary<ba)||(primary==bs&&secondary==ba&&r==((mode>>1)&1))){bs=primary;ba=secondary;bi=i;bw=ww;bh=hh;br=r;}
     }
     if(bi<0)break; z.put(t,bi,bw,bh,br);
   }
   if(z.value>bv){bv=z.value;best=z.p;}
 }
 cout<<"{\"placements\":[";
 for(int i=0;i<(int)best.size();i++){if(i)cout<<',';auto&p=best[i];cout<<"{\"type\":\""<<a[p.t].id<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.rot<<'}';}
 cout<<"]}\n";
}
