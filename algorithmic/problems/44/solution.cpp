#include <bits/stdc++.h>
using namespace std;

struct Pt { long long x, y; };
struct Node {
    int id, l = -1, r = -1, p = -1, cnt = 0;
    long long lx, rx, ly, ry;
};
int n;
vector<Pt> q;
vector<Node> tr;
vector<int> where;
vector<unsigned char> alive, primeId;

int build(vector<int>& a, int lo, int hi, int dep, int parent) {
    if (lo >= hi) return -1;
    int mid = (lo + hi) / 2;
    int axis = dep & 1;
    nth_element(a.begin()+lo, a.begin()+mid, a.begin()+hi, [&](int u, int v) {
        return axis ? q[u].y < q[v].y : q[u].x < q[v].x;
    });
    int z = (int)tr.size();
    tr.push_back({a[mid], -1, -1, parent, hi-lo, q[a[mid]].x, q[a[mid]].x, q[a[mid]].y, q[a[mid]].y});
    where[a[mid]] = z;
    int l = build(a, lo, mid, dep+1, z), r = build(a, mid+1, hi, dep+1, z);
    tr[z].l = l; tr[z].r = r;
    for (int c : {l,r}) if (c >= 0) {
        tr[z].lx=min(tr[z].lx,tr[c].lx); tr[z].rx=max(tr[z].rx,tr[c].rx);
        tr[z].ly=min(tr[z].ly,tr[c].ly); tr[z].ry=max(tr[z].ry,tr[c].ry);
    }
    return z;
}
inline long double dsq(int a, int b) {
    long double dx=(long double)q[a].x-q[b].x, dy=(long double)q[a].y-q[b].y;
    return dx*dx+dy*dy;
}
inline long double boxd(int id, int b) {
    const Node &v=tr[id]; long double dx=0,dy=0;
    if(q[b].x<v.lx) dx=(long double)v.lx-q[b].x; else if(q[b].x>v.rx) dx=(long double)q[b].x-v.rx;
    if(q[b].y<v.ly) dy=(long double)v.ly-q[b].y; else if(q[b].y>v.ry) dy=(long double)q[b].y-v.ry;
    return dx*dx+dy*dy;
}
int bestId, seen;
long double bestD;
void nearestDfs(int z, int from) {
    if (z < 0 || !tr[z].cnt || seen >= 112 || boxd(z,from) > bestD) return;
    ++seen;
    int u=tr[z].id;
    if (alive[u]) { long double d=dsq(u,from); if(d<bestD || (d==bestD && u<bestId)) bestD=d,bestId=u; }
    int a=tr[z].l,b=tr[z].r;
    if(a>=0 && b>=0 && boxd(b,from)<boxd(a,from)) swap(a,b);
    nearestDfs(a,from); nearestDfs(b,from);
}
int getNearest(int from) {
    bestId=-1; bestD=1e100L; seen=0; nearestDfs(0,from);
    return bestId;
}
void eraseCity(int u) {
    alive[u]=0;
    for(int z=where[u];z>=0;z=tr[z].p) --tr[z].cnt;
}
inline double edgeCost(const vector<int>& a, int t) {
    int u=a[t-1], v=a[t];
    double dx=(double)q[u].x-q[v].x, dy=(double)q[u].y-q[v].y;
    double d=hypot(dx,dy);
    return (t%10==0 && !primeId[u]) ? 1.1*d : d;
}
double cost(const vector<int>& a) { double s=0; for(int t=1;t<=n;t++) s+=edgeCost(a,t); return s; }
void adjacentImprove(vector<int>& a, int rounds) {
    for(int pass=0;pass<rounds;pass++) {
        bool changed=false;
        int start=(pass&1)?n-2:1, stop=(pass&1)?0:n-1;
        int step=(pass&1)?-1:1;
        for(int i=start;i!=stop;i+=step) { // swap route positions i and i+1; both are nonzero cities
            double old=edgeCost(a,i)+edgeCost(a,i+1)+edgeCost(a,i+2);
            swap(a[i],a[i+1]);
            double nw=edgeCost(a,i)+edgeCost(a,i+1)+edgeCost(a,i+2);
            if(nw+1e-7 < old) changed=true;
            else swap(a[i],a[i+1]);
        }
        if(!changed) break;
    }
}
uint64_t morton(long long x, long long y) {
    uint32_t X=(uint32_t)(x+1000000000LL), Y=(uint32_t)(y+1000000000LL);
    uint64_t z=0;
    for(int b=0;b<31;b++) z|=((uint64_t)((X>>b)&1)<<(2*b))|((uint64_t)((Y>>b)&1)<<(2*b+1));
    return z;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>n)) return 0;
    q.resize(n); for(auto &p:q) cin>>p.x>>p.y;
    primeId.assign(n,true); if(n>0)primeId[0]=false; if(n>1)primeId[1]=false;
    for(int i=2;(long long)i*i<n;i++) if(primeId[i]) for(int j=i*i;j<n;j+=i) primeId[j]=false;

    vector<int> ids(n); iota(ids.begin(),ids.end(),0); tr.reserve(n); where.resize(n); alive.assign(n,1);
    build(ids,0,n,0,-1); eraseCity(0);
    vector<int> nn; nn.reserve(n+1); nn.push_back(0); int cur=0;
    for(int k=1;k<n;k++) { int v=getNearest(cur); if(v<0) { for(int j=1;j<n;j++) if(alive[j]) {v=j;break;} } nn.push_back(v); eraseCity(v); cur=v; }
    nn.push_back(0); adjacentImprove(nn,3);

    // A space-filling ordering is a cheap robust fallback when greedy nearest-neighbor gets trapped.
    vector<int> mo(n-1); iota(mo.begin(),mo.end(),1);
    sort(mo.begin(),mo.end(),[](int a,int b){ return morton(q[a].x,q[a].y)<morton(q[b].x,q[b].y); });
    vector<int> spatial; spatial.reserve(n+1); spatial.push_back(0); spatial.insert(spatial.end(),mo.begin(),mo.end()); spatial.push_back(0);
    vector<int> rev=spatial; reverse(rev.begin()+1,rev.end()-1);
    adjacentImprove(spatial,2); adjacentImprove(rev,2);
    vector<int> ans=nn; double bc=cost(ans);
    for(auto *cand : {&spatial,&rev}) { double z=cost(*cand); if(z<bc) bc=z,ans=*cand; }
    cout<<n+1<<'\n'; for(int v:ans) cout<<v<<'\n';
}
