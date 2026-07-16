#include <bits/stdc++.h>
using namespace std;
struct Pt { long long x,y; };
int main(){
 ios::sync_with_stdio(false); cin.tie(nullptr);
 int n; if(!(cin>>n)) return 0;
 vector<Pt> p(n); for(auto &z:p) cin>>z.x>>z.y;
 vector<char> prime(n,1); if(n) prime[0]=0; if(n>1) prime[1]=0;
 for(int i=2;i*1LL*i<n;i++) if(prime[i]) for(int j=i*i;j<n;j+=i) prime[j]=0;
 auto dis=[&](int a,int b){ return hypot((double)p[a].x-p[b].x,(double)p[a].y-p[b].y); };
 auto weight=[&](int step,int id){ return (step%10==0 && !prime[id]) ? 1.1 : 1.0; };
 auto cost=[&](const vector<int>& v){
   double s=0; int last=0, step=1;
   for(int x:v){s+=weight(step++,last)*dis(last,x); last=x;}
   return s+weight(step,last)*dis(last,0);
 };
 vector<vector<int>> cand;
 vector<int> f; for(int i=1;i<n;i++) f.push_back(i); cand.push_back(f); reverse(f.begin(),f.end()); cand.push_back(f);
 // Alternating sweeps in equal-population vertical and horizontal strips.
 for(int bands: {4,8,16,32,64}) if(n>2){
   vector<int> v(n-1); iota(v.begin(),v.end(),1);
   int m=v.size(), bs=(m+bands-1)/bands;
   for(int b=0;b<bands;b++){
     int l=b*bs,r=min(m,l+bs); if(l>=r) break;
     sort(v.begin()+l,v.begin()+r,[&](int a,int c){ return p[a].y==p[c].y?a<c:p[a].y<p[c].y; });
     if(b&1) reverse(v.begin()+l,v.begin()+r);
   } cand.push_back(v);
   v.resize(n-1); iota(v.begin(),v.end(),1);
   sort(v.begin(),v.end(),[&](int a,int c){return p[a].y==p[c].y?a<c:p[a].y<p[c].y;});
   for(int b=0;b<bands;b++){
     int l=b*bs,r=min(m,l+bs); if(l>=r) break;
     sort(v.begin()+l,v.begin()+r,[&](int a,int c){ return p[a].x==p[c].x?a<c:p[a].x<p[c].x; });
     if(b&1) reverse(v.begin()+l,v.begin()+r);
   } cand.push_back(v);
 }
 // A space-filling ordering supplies a non-strip alternative.
 long long minx=p[0].x,maxx=p[0].x,miny=p[0].y,maxy=p[0].y;
 for(auto z:p){minx=min(minx,z.x);maxx=max(maxx,z.x);miny=min(miny,z.y);maxy=max(maxy,z.y);}
 auto morton=[&](int id){
   unsigned long long X=(unsigned long long)((__int128)(p[id].x-minx)*((1<<21)-1)/max(1LL,maxx-minx));
   unsigned long long Y=(unsigned long long)((__int128)(p[id].y-miny)*((1<<21)-1)/max(1LL,maxy-miny));
   unsigned long long z=0; for(int b=0;b<21;b++) z|=((X>>b)&1)<<(2*b),z|=((Y>>b)&1)<<(2*b+1); return z;
 };
 f.resize(n-1); iota(f.begin(),f.end(),1); sort(f.begin(),f.end(),[&](int a,int b){auto x=morton(a),y=morton(b);return x==y?a<b:x<y;}); cand.push_back(f); reverse(f.begin(),f.end()); cand.push_back(f);
 // Reversing a geometric cycle preserves its distances but changes which IDs occupy carrot steps.
 int originalCandidates=cand.size();
 for(int i=0;i<originalCandidates;i++){
   vector<int> r=cand[i]; reverse(r.begin(),r.end()); cand.push_back(move(r));
 }
 int best=0; double bc=cost(cand[0]); for(int i=1;i<(int)cand.size();i++){double z=cost(cand[i]);if(z<bc)bc=z,best=i;}
 vector<int> v=cand[best]; int m=v.size();
 // Exact position-aware 2-opt.  Radius 16 is deliberately wider than the prior short repair.
 vector<double> d(n); // edge from route position k to k+1
 auto idat=[&](int pos)->int { return pos==0||pos==n?0:v[pos-1]; };
 for(int k=0;k<n;k++) d[k]=dis(idat(k),idat(k+1));
 for(int pass=0;pass<2;pass++){
   bool any=false;
   for(int i=0;i<m;i++) for(int j=i+1;j<m && j<=i+16;j++){
     int a=i+1,b=j+1; double old=0, neu=0;
     for(int k=a-1;k<=b;k++) old+=d[k]*weight(k+1,idat(k));
     neu+=dis(idat(a-1),idat(b))*weight(a,idat(a-1));
     for(int k=a;k<b;k++) { int oldpos=a+b-k; neu+=d[oldpos-1]*weight(k+1,idat(oldpos)); }
     neu+=dis(idat(a),idat(b+1))*weight(b+1,idat(a));
     if(neu+1e-7<old){
       reverse(v.begin()+i,v.begin()+j+1);
       reverse(d.begin()+a,d.begin()+b);
       d[a-1]=dis(idat(a-1),idat(a)); d[b]=dis(idat(b),idat(b+1)); any=true;
     }
   }
   if(!any) break;
 }
 cout<<n+1<<'\n'<<0<<'\n'; for(int x:v) cout<<x<<'\n'; cout<<0<<'\n';
}
