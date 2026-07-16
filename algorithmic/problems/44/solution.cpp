#include <bits/stdc++.h>
using namespace std;
struct Pt { long long x,y; };
static vector<Pt> a;
static vector<char> primeId;
static int n;

// Hilbert index for a 2^21 by 2^21 square.  The transform before this routine
// supplies the eight rotations/reflections of the normalized coordinate plane.
static unsigned long long hilbert(unsigned x, unsigned y) {
    unsigned long long d=0;
    for (unsigned s=1u<<20; s; s>>=1) {
        unsigned rx=(x&s)!=0, ry=(y&s)!=0;
        d |= (unsigned long long)((3*rx)^ry)*s*s;
        if (!ry) {
            if (rx) { x=(2*s-1)-x; y=(2*s-1)-y; }
            swap(x,y);
        }
    }
    return d;
}
static inline double dis(int u,int v) {
    double dx=(double)a[u].x-a[v].x, dy=(double)a[u].y-a[v].y;
    return hypot(dx,dy);
}
static inline double edge(const vector<int>& p,int i) {
    int j=(i+1==n?0:i+1);
    double z=dis(p[i],p[j]);
    if ((i+1)%10==0 && !primeId[p[i]]) z*=1.1;
    return z;
}
static double score(const vector<int>& p) { double z=0; for(int i=0;i<n;i++) z+=edge(p,i); return z; }

// Exact delta for a bounded reversal.  Reversing across a carrot boundary
// changes the IDs supplying internal penalized steps, so score every affected
// directed edge instead of excluding that useful part of the neighborhood.
static bool tryReverse(vector<int>& p,int l,int r) {
    double old=0;
    for(int i=l-1;i<=r;i++) old+=edge(p,i);
    reverse(p.begin()+l,p.begin()+r+1);
    double nw=0;
    for(int i=l-1;i<=r;i++) nw+=edge(p,i);
    if(nw+1e-7<old) return true;
    reverse(p.begin()+l,p.begin()+r+1);
    return false;
}
static bool trySwap(vector<int>& p,int u,int v) {
    if(u==v) return false;
    int lo=max(0,min(u,v)-1), hi=min(n-1,max(u,v)+1);
    double old=0,nw=0;
    for(int i=lo;i<=hi;i++) old+=edge(p,i);
    swap(p[u],p[v]);
    for(int i=lo;i<=hi;i++) nw+=edge(p,i);
    if(nw+1e-7<old) return true;
    swap(p[u],p[v]); return false;
}
static void improve(vector<int>& p) {
    // Include cuts straddling every tenth step; their carrot delta is exact.
    for(int l=1;l<n-1;l++) for(int r=l+1;r<n && r<=l+10;r++) tryReverse(p,l,r);
    // At carrot positions, prefer a nearby prime only when its full affected
    // edge cost pays for the changed spatial adjacency.
    for(int pos=9;pos<n;pos+=10) if(!primeId[p[pos]]) {
        int best=-1; double gain=0;
        for(int q=max(1,pos-36);q<=min(n-1,pos+36);q++) if(primeId[p[q]]) {
            int lo=max(0,min(pos,q)-1), hi=min(n-1,max(pos,q)+1); double before=0;
            for(int i=lo;i<=hi;i++) before+=edge(p,i);
            swap(p[pos],p[q]); double after=0; for(int i=lo;i<=hi;i++) after+=edge(p,i); swap(p[pos],p[q]);
            if(before-after>gain) gain=before-after,best=q;
        }
        if(best>=0) trySwap(p,pos,best);
    }
    for(int l=1;l<n-1;l++) for(int r=l+1;r<n && r<=l+8;r++) tryReverse(p,l,r);
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>n)) return 0; a.resize(n);
    long long xmin=LLONG_MAX,xmax=LLONG_MIN,ymin=LLONG_MAX,ymax=LLONG_MIN;
    for(auto &q:a){cin>>q.x>>q.y; xmin=min(xmin,q.x);xmax=max(xmax,q.x);ymin=min(ymin,q.y);ymax=max(ymax,q.y);}
    primeId.assign(n,true); if(n>0)primeId[0]=false; if(n>1)primeId[1]=false;
    for(int i=2;(long long)i*i<n;i++) if(primeId[i]) for(int j=i*i;j<n;j+=i) primeId[j]=false;
    vector<int> best; double bestCost=1e300;
    vector<pair<double,vector<int>>> seeds;
    const unsigned LIM=(1u<<21)-1;
    for(int ori=0;ori<8;ori++) {
        vector<pair<unsigned long long,int>> key; key.reserve(n);
        for(int i=0;i<n;i++) {
            unsigned X=(xmax==xmin?0:(unsigned)((__int128)(a[i].x-xmin)*LIM/(xmax-xmin)));
            unsigned Y=(ymax==ymin?0:(unsigned)((__int128)(a[i].y-ymin)*LIM/(ymax-ymin)));
            unsigned u=X,v=Y;
            if(ori&1) u=LIM-u; if(ori&2) v=LIM-v; if(ori&4) swap(u,v);
            key.push_back({hilbert(u,v),i});
        }
        sort(key.begin(),key.end()); vector<int> p; p.reserve(n);
        int z=0; while(key[z].second!=0) z++;
        for(int k=0;k<n;k++) p.push_back(key[(z+k)%n].second);
        seeds.push_back({score(p),move(p)});
    }
    sort(seeds.begin(),seeds.end(),[](const auto& A,const auto& B){return A.first<B.first;});
    // Raw curve cost cheaply filters orientations; optimize the two most
    // promising distinct cyclic embeddings under the full carrot objective.
    for(int k=0;k<min(2,(int)seeds.size());k++) {
        vector<int> p=move(seeds[k].second); improve(p); double c=score(p);
        if(c<bestCost) bestCost=c,best.swap(p);
    }
    cout<<n+1<<'\n'; for(int v:best) cout<<v<<'\n'; cout<<0<<'\n';
}
