#include <bits/stdc++.h>
using namespace std;
struct Pt { long long x,y; };
static vector<Pt> p;
static vector<char> primeid;
static int n;

static inline double d(int a,int b) { return hypot((double)p[a].x-p[b].x,(double)p[a].y-p[b].y); }
static inline double ec(int step0,int a,int b) { // step0 is zero based edge index
    double z=d(a,b);
    if ((step0+1)%10==0 && !primeid[a]) z*=1.1;
    return z;
}
static double cost(const vector<int>& r) {
    double z=0; for(int i=0;i<n;i++) z+=ec(i,r[i],r[i+1]); return z;
}
static vector<int> anchored(vector<int> v) {
    int q=find(v.begin(),v.end(),0)-v.begin();
    vector<int> r; r.reserve(n+1); r.push_back(0);
    for(int k=1;k<n;k++) r.push_back(v[(q+k)%n]);
    r.push_back(0); return r;
}
static unsigned long long hilbert(unsigned x,unsigned y) {
    const unsigned S=1u<<20; unsigned long long ans=0;
    for(unsigned s=S>>1;s;s>>=1) {
        unsigned rx=(x&s)!=0, ry=(y&s)!=0;
        ans += (unsigned long long)s*s*((3*rx)^ry);
        if(!ry) { if(rx) { x=S-1-x; y=S-1-y; } swap(x,y); }
    }
    return ans;
}
static vector<int> horder(int mode) {
    bool sw=mode&1, sx=mode&2, sy=mode&4;
    vector<long long> a(n),b(n); long long amin=LLONG_MAX,amax=LLONG_MIN,bmin=LLONG_MAX,bmax=LLONG_MIN;
    for(int i=0;i<n;i++) {
        long long u=sw?p[i].y:p[i].x, v=sw?p[i].x:p[i].y;
        if(sx) u=-u; if(sy) v=-v;
        a[i]=u;b[i]=v; amin=min(amin,u);amax=max(amax,u);bmin=min(bmin,v);bmax=max(bmax,v);
    }
    const unsigned S=(1u<<20)-1;
    vector<pair<unsigned long long,int>> z; z.reserve(n);
    for(int i=0;i<n;i++) {
        unsigned u=amax==amin?0:(unsigned)((__int128)(a[i]-amin)*S/(amax-amin));
        unsigned v=bmax==bmin?0:(unsigned)((__int128)(b[i]-bmin)*S/(bmax-bmin));
        z.push_back({hilbert(u,v),i});
    }
    sort(z.begin(),z.end()); vector<int> q; for(auto e:z) q.push_back(e.second); return q;
}
// Evaluate an exact weighted 2-opt reversal: reverse r[l+1..rr].
static bool try2opt(vector<int>& r, vector<int>& pos, int l,int rr) {
    if(l<0 || rr>=n || rr<=l+1) return false;
    // Reversal preserves every internal undirected distance.  Only its two
    // boundary edges and the sparse (each tenth) source-dependent weights vary.
    double delta=ec(l,r[l],r[rr])-ec(l,r[l],r[l+1]);
    delta+=ec(rr,r[l+1],r[rr+1])-ec(rr,r[rr],r[rr+1]);
    // Old internal edge i maps to reversed edge l+rr-i, whose source is
    // old r[i+1].  Enumerate only the two relevant residue classes rather
    // than scanning every edge in this bounded interval.
    int lo=l+1;
    int i=lo+(9-lo%10+10)%10; // old step i+1 is a multiple of 10
    for(;i<rr;i+=10) if(!primeid[r[i]]) delta-=0.1*d(r[i],r[i+1]);
    int residue=(l+rr+1)%10; // (l+rr-i)+1 is a multiple of 10
    i=lo+(residue-lo%10+10)%10;
    for(;i<rr;i+=10) if(!primeid[r[i+1]]) delta+=0.1*d(r[i],r[i+1]);
    if(delta >= -1e-7) return false;
    reverse(r.begin()+l+1,r.begin()+rr+1);
    for(int k=l+1;k<=rr;k++) pos[r[k]]=k;
    return true;
}
// Short reversals make the carrot delta exactly affordable, unlike ordinary 2-opt.
static void shortTwoOpt(vector<int>& r) {
    vector<int> pos(n); for(int i=0;i<n;i++) pos[r[i]]=i;
    vector<int> byy(n); iota(byy.begin(),byy.end(),0);
    sort(byy.begin(),byy.end(),[](int a,int b){return p[a].y==p[b].y?a<b:p[a].y<p[b].y;});
    vector<int> yr(n); for(int i=0;i<n;i++) yr[byy[i]]=i;
    const int NEAR=10, SPAN=72;
    // Two sweeps expose reversals made possible by earlier changes.
    for(int pass=0;pass<2;pass++) for(int l=0;l<n-1;l++) {
        int a=r[l];
        for(int dx=-NEAR;dx<=NEAR;dx++) if(dx) {
            int b=a+dx; if(b>=0&&b<n) { int j=pos[b]; if(j>l+1 && j-l<=SPAN) try2opt(r,pos,l,j); }
        }
        int q=yr[a];
        for(int dy=-NEAR;dy<=NEAR;dy++) if(dy) {
            int bq=q+dy; if(bq>=0&&bq<n) { int j=pos[byy[bq]]; if(j>l+1 && j-l<=SPAN) try2opt(r,pos,l,j); }
        }
    }
}
int main(){
    ios::sync_with_stdio(false);cin.tie(nullptr);
    if(!(cin>>n)) return 0; p.resize(n); for(auto &v:p) cin>>v.x>>v.y;
    primeid.assign(n,true); if(n)primeid[0]=false; if(n>1)primeid[1]=false;
    for(int i=2;(long long)i*i<n;i++) if(primeid[i]) for(int j=i*i;j<n;j+=i) primeid[j]=false;
    vector<int> base(n);iota(base.begin(),base.end(),0);
    vector<int> best=anchored(base); double bestc=cost(best);
    reverse(base.begin()+1,base.end()); vector<int> cand=anchored(base); double cc=cost(cand); if(cc<bestc)best=cand,bestc=cc;
    for(int m=0;m<8;m++) { auto q=horder(m); for(int rev=0;rev<2;rev++) { if(rev) reverse(q.begin(),q.end()); cand=anchored(q); cc=cost(cand); if(cc<bestc) best=cand,bestc=cc; if(rev) reverse(q.begin(),q.end()); } }
    // The refinement is monotone under the real checker objective, so the candidate is safe.
    cand=best; shortTwoOpt(cand); cc=cost(cand); if(cc<bestc) best=cand;
    cout<<n+1<<'\n'; for(int v:best) cout<<v<<'\n';
}
