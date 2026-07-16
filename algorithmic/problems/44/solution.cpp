#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x,y; };
struct Node {
    int id, l=-1, r=-1, par=-1, live=0;
    long long lx,rx,ly,ry;
};

static vector<Pt> p;
static vector<Node> tr;
static vector<int> where;

int build(vector<int>& a, int lo, int hi, int dep, int parent) {
    if (lo>=hi) return -1;
    int mid=(lo+hi)/2, d=dep&1;
    nth_element(a.begin()+lo,a.begin()+mid,a.begin()+hi,[d](int u,int v) {
        return d ? p[u].y<p[v].y : p[u].x<p[v].x;
    });
    int q=(int)tr.size(); tr.push_back({});
    tr[q].id=a[mid]; tr[q].par=parent;
    int left=build(a,lo,mid,dep+1,q);
    int right=build(a,mid+1,hi,dep+1,q);
    tr[q].l=left; tr[q].r=right;
    Node &z=tr[q]; z.live=1; z.lx=z.rx=p[z.id].x; z.ly=z.ry=p[z.id].y;
    for(int c: {z.l,z.r}) if(c!=-1) {
        z.live+=tr[c].live;
        z.lx=min(z.lx,tr[c].lx); z.rx=max(z.rx,tr[c].rx);
        z.ly=min(z.ly,tr[c].ly); z.ry=max(z.ry,tr[c].ry);
    }
    where[z.id]=q;
    return q;
}
static inline long double d2(int a,int b) {
    long double dx=(long double)p[a].x-p[b].x, dy=(long double)p[a].y-p[b].y;
    return dx*dx+dy*dy;
}
static inline long double boxd2(int q,const Node &z) {
    long double dx=0,dy=0;
    if(p[q].x<z.lx) dx=(long double)z.lx-p[q].x; else if(p[q].x>z.rx) dx=(long double)p[q].x-z.rx;
    if(p[q].y<z.ly) dy=(long double)z.ly-p[q].y; else if(p[q].y>z.ry) dy=(long double)p[q].y-z.ry;
    return dx*dx+dy*dy;
}
void nearest(int q,int v,int &ans,long double &best) {
    if(v<0 || !tr[v].live || boxd2(q,tr[v])>=best) return;
    int id=tr[v].id;
    // A deleted node may still have live descendants.
    if(id!=q && where[id]>=0 && d2(q,id)<best) { ans=id; best=d2(q,id); }
    int a=tr[v].l,b=tr[v].r;
    long double da=(a<0?1e100L:boxd2(q,tr[a])), db=(b<0?1e100L:boxd2(q,tr[b]));
    if(da>db) swap(a,b);
    nearest(q,a,ans,best); nearest(q,b,ans,best);
}
void erasePoint(int id) {
    int v=where[id];
    // Marking where negative distinguishes deleted point nodes from live points.
    where[id]=-1;
    while(v!=-1) { --tr[v].live; v=tr[v].par; }
}

static inline uint64_t hilbert(uint32_t x,uint32_t y) {
    uint64_t z=0;
    // Standard xy-to-Hilbert rotation, with a 21-bit square.
    for(uint32_t s=1u<<20; s; s>>=1) {
        uint32_t rx=(x&s)?1:0, ry=(y&s)?1:0;
        z += (uint64_t)s*s*((3*rx)^ry);
        if(!ry) { if(rx) { x=(1u<<21)-1-x; y=(1u<<21)-1-y; } swap(x,y); }
    }
    return z;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    p.resize(n); for(auto &q:p) cin>>q.x>>q.y;
    vector<char> prime(max(2,n),true); prime[0]=false; if(n>1) prime[1]=false;
    for(int i=2;i*1LL*i<n;i++) if(prime[i]) for(int j=i*i;j<n;j+=i) prime[j]=false;
    auto edge=[&](int a,int b){ return hypot((double)p[a].x-p[b].x,(double)p[a].y-p[b].y); };
    auto cost=[&](const vector<int>& r) {
        double s=0; for(int t=1;t<=n;t++) { double w=(t%10==0&&!prime[r[t-1]])?1.1:1.; s+=w*edge(r[t-1],r[t]); } return s;
    };
    vector<int> best; double bestCost=1e300;
    auto consider=[&](vector<int> cyc) {
        int at=find(cyc.begin(),cyc.end(),0)-cyc.begin();
        vector<int> r; r.reserve(n+1); r.push_back(0);
        for(int k=1;k<n;k++) r.push_back(cyc[(at+k)%n]); r.push_back(0);
        double c=cost(r); if(c<bestCost) bestCost=c,best.swap(r);
    };
    // The dynamic kd-tree walk is the primary construction.
    vector<int> ids(n); iota(ids.begin(),ids.end(),0); tr.reserve(n); where.assign(n,-1);
    int root=build(ids,0,n,0,-1); erasePoint(0);
    vector<int> nn; nn.reserve(n); nn.push_back(0); int cur=0;
    for(int k=1;k<n;k++) { int q=-1; long double b=1e100L; nearest(cur,root,q,b); if(q<0) break; nn.push_back(q); erasePoint(q); cur=q; }
    consider(nn); reverse(nn.begin()+1,nn.end()); consider(nn);
    long long minx=p[0].x,maxx=p[0].x,miny=p[0].y,maxy=p[0].y;
    for(auto q:p) minx=min(minx,q.x),maxx=max(maxx,q.x),miny=min(miny,q.y),maxy=max(maxy,q.y);
    const uint64_t S=(1u<<21)-1;
    // Static fallback family protects simple smooth input distributions.
    for(int mode=0;mode<4;mode++) {
        vector<pair<uint64_t,int>> v; v.reserve(n);
        for(int i=0;i<n;i++) {
            uint32_t X=(maxx==minx?0:(uint64_t)(p[i].x-minx)*S/(maxx-minx));
            uint32_t Y=(maxy==miny?0:(uint64_t)(p[i].y-miny)*S/(maxy-miny));
            if(mode&1) X=S-X; if(mode&2) Y=S-Y;
            v.push_back({hilbert(X,Y),i});
        }
        sort(v.begin(),v.end()); vector<int> a; for(auto e:v)a.push_back(e.second);
        consider(a); reverse(a.begin(),a.end()); consider(a);
    }
    vector<int> ay(n); iota(ay.begin(),ay.end(),0);
    sort(ay.begin(),ay.end(),[](int a,int b){return p[a].y<p[b].y;}); consider(ay); reverse(ay.begin(),ay.end()); consider(ay);
    vector<int> ax(n); iota(ax.begin(),ax.end(),0); consider(ax); reverse(ax.begin(),ax.end()); consider(ax);
    // Position-aware short 2-opt: all changed carrot-weighted edges are recomputed exactly.
    int W=(n<=10000?24:6), passes=(n<=10000?2:1);
    for(int pass=0;pass<passes;pass++) for(int l=1;l<n-1;l++) {
        for(int rr=l+1;rr<n && rr<=l+W;rr++) {
            auto val=[&](int pos){ return (pos>=l&&pos<=rr)?best[l+rr-pos]:best[pos]; };
            double old=0, nw=0;
            for(int t=l;t<=rr+1;t++) {
                double w=(t%10==0&&!prime[best[t-1]])?1.1:1.; old+=w*edge(best[t-1],best[t]);
                int a=val(t-1),b=val(t); double w2=(t%10==0&&!prime[a])?1.1:1.; nw+=w2*edge(a,b);
            }
            if(nw+1e-7<old) { reverse(best.begin()+l,best.begin()+rr+1); bestCost+=nw-old; }
        }
    }
    cout<<n+1<<'\n'; for(int x:best) cout<<x<<'\n';
}
