#include <bits/stdc++.h>
using namespace std;
struct Pt { long long x,y; };
static vector<Pt> p;
static vector<char> prime;
static inline double edge(int a,int b){ return hypot((double)p[a].x-p[b].x,(double)p[a].y-p[b].y); }
static inline double term(const vector<int>& r,int t){
    double z=edge(r[t-1],r[t]);
    return (t%10==0 && !prime[r[t-1]]) ? z*1.1 : z;
}
static double value(const vector<int>& r){ double s=0; for(int t=1;t<(int)r.size();++t)s+=term(r,t); return s; }
static uint64_t hilbert(uint32_t x,uint32_t y, int bits){
    // rot() reflects in the full grid, not the current bit cell.
    const uint32_t side=1u<<bits;
    uint64_t d=0;
    for(uint32_t s=side>>1;s;s>>=1){
        uint32_t rx=(x&s)!=0, ry=(y&s)!=0;
        d+=(uint64_t)s*s*((3*rx)^ry);
        if(!ry){ if(rx){ x=side-1-x; y=side-1-y; } swap(x,y); }
    }
    return d;
}
static void descend(vector<int>& r){
    int n=(int)r.size()-1;
    // Position-aware adjacent exchanges are cheap and also repair carrot slots.
    for(int pass=0;pass<3;pass++){
        bool any=false;
        for(int i=1;i<n;i++){
            int lo=i, hi=min(n,i+2); double old=0; for(int t=lo;t<=hi;t++) old+=term(r,t);
            swap(r[i],r[i+1]); double nw=0; for(int t=lo;t<=hi;t++) nw+=term(r,t);
            if(nw+1e-7<old) any=true; else swap(r[i],r[i+1]);
        }
        if(!any) break;
    }
    // Bounded reversals are an exact small 2-opt neighborhood despite varying multipliers.
    for(int pass=0;pass<2;pass++){
        bool any=false;
        for(int len=3;len<=6;len++) for(int l=1;l+len-1<n;l++){
            int rr=l+len-1; double old=0; for(int t=l;t<=rr+1;t++) old+=term(r,t);
            reverse(r.begin()+l,r.begin()+rr+1); double nw=0; for(int t=l;t<=rr+1;t++) nw+=term(r,t);
            if(nw+1e-7<old) any=true; else reverse(r.begin()+l,r.begin()+rr+1);
        }
        if(!any) break;
    }
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0; p.resize(n);
    for(auto &q:p) cin>>q.x>>q.y;
    prime.assign(n,true); if(n>0)prime[0]=false; if(n>1)prime[1]=false;
    for(int i=2;i*(long long)i<n;i++) if(prime[i]) for(long long j=(long long)i*i;j<n;j+=i) prime[j]=false;
    vector<int> best; double bestv=1e300;
    auto consider=[&](vector<int> a){
        vector<int> r; r.reserve(n+1); r.push_back(0); r.insert(r.end(),a.begin(),a.end()); r.push_back(0);
        double v=value(r); if(v<bestv){bestv=v;best.swap(r);}
    };
    // Preserved family: all dihedral orientations of a correctly scaled Hilbert curve.
    long long xmin=p[0].x,xmax=p[0].x,ymin=p[0].y,ymax=p[0].y;
    for(auto q:p){xmin=min(xmin,q.x);xmax=max(xmax,q.x);ymin=min(ymin,q.y);ymax=max(ymax,q.y);}
    const uint32_t M=(1u<<21)-1;
    vector<uint32_t> xx(n),yy(n);
    for(int i=0;i<n;i++){
        xx[i]=xmax==xmin?M/2:(uint32_t)((__int128)(p[i].x-xmin)*M/(xmax-xmin));
        yy[i]=ymax==ymin?M/2:(uint32_t)((__int128)(p[i].y-ymin)*M/(ymax-ymin));
    }
    for(int o=0;o<8;o++){
        vector<pair<uint64_t,int>> z; z.reserve(max(0,n-1));
        for(int i=1;i<n;i++){
            uint32_t a=xx[i],b=yy[i], u,v;
            switch(o){
              case 0:u=a;v=b;break; case 1:u=M-a;v=b;break; case 2:u=a;v=M-b;break; case 3:u=M-a;v=M-b;break;
              case 4:u=b;v=a;break; case 5:u=M-b;v=a;break; case 6:u=b;v=M-a;break; default:u=M-b;v=M-a;
            }
            z.push_back({hilbert(u,v,21),i});
        }
        sort(z.begin(),z.end()); vector<int>a; for(auto e:z)a.push_back(e.second); consider(a);
    }
    // Diversifying family: x is already ordered, so serpentine strips form a multiscale spatial decomposition.
    int root=max(2,(int)sqrt((double)max(1,n-1)));
    for(int mul: {1,2,4,8}) for(int flip=0;flip<2;flip++){
        int block=max(2,root*mul); vector<int>a; a.reserve(n-1);
        for(int l=1,part=0;l<n;l+=block,part++){
            int r=min(n,l+block); vector<int> b; for(int i=l;i<r;i++)b.push_back(i);
            sort(b.begin(),b.end(),[&](int i,int j){return p[i].y==p[j].y?i<j:p[i].y<p[j].y;});
            if(((part+flip)&1)) reverse(b.begin(),b.end()); a.insert(a.end(),b.begin(),b.end());
        }
        consider(a);
    }
    descend(best);
    cout<<n+1<<'\n'; for(int x:best) cout<<x<<'\n';
}
