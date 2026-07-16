#include <bits/stdc++.h>
using namespace std;
struct Item{string id; int w,h,lim; long long v;};
struct P{int t,x,y,r;};
struct JS{
 string s; size_t p=0; JS(){s.assign(istreambuf_iterator<char>(cin),{});} void ws(){while(p<s.size()&&isspace((unsigned char)s[p]))p++;}
 char ch(){ws(); return p<s.size()?s[p]:0;} void eat(char c){ws(); if(p<s.size()&&s[p]==c)p++;}
 string str(){eat('"'); string r; while(p<s.size()&&s[p]!='"'){if(s[p]=='\\'&&p+1<s.size())p++; r+=s[p++];} if(p<s.size())p++; return r;}
 long long num(){ws(); long long z=0; bool n=false;if(s[p]=='-'){n=true;p++;}while(p<s.size()&&isdigit((unsigned char)s[p]))z=z*10+s[p++]-'0';return n?-z:z;}
 bool boolean(){ws(); bool r=s.compare(p,4,"true")==0; p+=r?4:5; return r;}
};
static uint64_t seed=88172645463325252ull;
uint64_t rnd(){seed^=seed<<7;seed^=seed>>9;return seed;}
string esc(const string& a){string r;for(char c:a){if(c=='"'||c=='\\')r+='\\';r+=c;}return r;}
int main(){
 JS q; int W=0,H=0; bool rotok=false; vector<Item> a;
 q.eat('{'); for(int top=0;top<2;top++) { if(top)q.eat(','); string key=q.str();q.eat(':');
  if(key=="bin"){q.eat('{');for(int z=0;z<3;z++){if(z)q.eat(',');string k=q.str();q.eat(':');if(k=="W")W=q.num();else if(k=="H")H=q.num();else rotok=q.boolean();}q.eat('}');}
  else {q.eat('[');bool first=true;while(q.ch()!=']'){if(!first)q.eat(',');first=false;q.eat('{');Item it;for(int z=0;z<5;z++){if(z)q.eat(',');string k=q.str();q.eat(':');if(k=="type")it.id=q.str();else if(k=="w")it.w=q.num();else if(k=="h")it.h=q.num();else if(k=="v")it.v=q.num();else it.lim=q.num();}q.eat('}');a.push_back(it);}q.eat(']');}
 } q.eat('}');
 int n=a.size(); vector<P> best; long long bestv=-1;
 auto start=chrono::steady_clock::now();
 // A trial is a horizontal-band guillotine packing. Each band has one chosen height,
 // but its interior is a value-priority one-dimensional mix of all shorter items.
 auto expired=[&](){return chrono::duration<double,milli>(chrono::steady_clock::now()-start).count()>820;};
 int trials=0;
 while(trials<200 && !expired()){
  trials++; vector<int> used(n); vector<P> out; long long val=0; int y=0;
  vector<double> bias(n); for(int i=0;i<n;i++) bias[i]=0.72 + (rnd()%10001)/10000.0*0.62;
  while(y<H){
   int remH=H-y, pickH=-1; vector<P> pick; vector<int> pickUse; long long pickV=-1; double pickKey=-1;
   // Test all meaningful band heights. This explicit comparison is the distinction from a skyline choice.
   for(int owner=0;owner<n && !expired();owner++) for(int rr=0;rr<(rotok?2:1);rr++){
    if(expired()) break;
    int bh=rr?a[owner].w:a[owner].h; if(bh>remH || used[owner]>=a[owner].lim) continue;
    vector<int> uu=used; vector<P> row; long long rv=0; int x=0;
    while(x<W){
     int bi=-1,br=0; double bs=-1;
     for(int i=0;i<n;i++) if(uu[i]<a[i].lim) for(int r=0;r<(rotok?2:1);r++){
      int iw=r?a[i].h:a[i].w, ih=r?a[i].w:a[i].h;
      if(iw>W-x||ih>bh)continue;
      // Favor density, penalize unused vertical band space, and randomize per type once per trial.
      double sc=(double)a[i].v/(iw*bh)*bias[i];
      sc*=1.0+0.08*(double)(rnd()%1000)/1000.0;
      if(sc>bs){bs=sc;bi=i;br=r;}
     }
     if(bi<0)break;
     int iw=br?a[bi].h:a[bi].w; row.push_back({bi,x,y,br}); x+=iw; uu[bi]++; rv+=a[bi].v;
    }
    if(row.empty())continue;
    double key=(double)rv/bh * (0.94+(rnd()%1001)/1000.0*0.12);
    if(key>pickKey){pickKey=key;pickH=bh;pick.swap(row);pickUse.swap(uu);pickV=rv;}
   }
   if(pickH<0 || expired()) break;
   out.insert(out.end(),pick.begin(),pick.end()); used.swap(pickUse); val+=pickV; y+=pickH;
  }
  if(val>bestv){bestv=val;best.swap(out);}
 }
 cout<<"{\"placements\":[";
 for(size_t i=0;i<best.size();i++){if(i)cout<<',';auto&p=best[i];cout<<"{\"type\":\""<<esc(a[p.t].id)<<"\",\"x\":"<<p.x<<",\"y\":"<<p.y<<",\"rot\":"<<p.r<<'}';}
 cout<<"]}";
}
