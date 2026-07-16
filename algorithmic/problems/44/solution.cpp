#include <bits/stdc++.h>
using namespace std;

struct Point { long long x,y; };
struct Node {
    int id=-1, left=-1, right=-1, parent=-1, cnt=0, pcnt=0;
    long long xl, xr, yl, yr;
};
static vector<Point> p;
static vector<Node> tr;
static vector<int> where;
static vector<unsigned char> alive, prime;

static int build(vector<int>& a, int lo, int hi, int par) {
    if (lo>=hi) return -1;
    long long xmin=p[a[lo]].x, xmax=xmin, ymin=p[a[lo]].y, ymax=ymin;
    for(int i=lo+1;i<hi;i++) {
        xmin=min(xmin,p[a[i]].x); xmax=max(xmax,p[a[i]].x);
        ymin=min(ymin,p[a[i]].y); ymax=max(ymax,p[a[i]].y);
    }
    int ax=(xmax-xmin < ymax-ymin);
    int mid=(lo+hi)/2;
    nth_element(a.begin()+lo,a.begin()+mid,a.begin()+hi,[&](int u,int v) {
        return ax ? p[u].y<p[v].y : p[u].x<p[v].x;
    });
    int q=(int)tr.size(); tr.push_back({});
    tr[q].id=a[mid]; tr[q].parent=par; tr[q].xl=xmin; tr[q].xr=xmax; tr[q].yl=ymin; tr[q].yr=ymax;
    where[a[mid]]=q;
    tr[q].left=build(a,lo,mid,q); tr[q].right=build(a,mid+1,hi,q);
    tr[q].cnt=hi-lo;
    tr[q].pcnt=0; for(int i=lo;i<hi;i++) tr[q].pcnt+=prime[a[i]];
    return q;
}
static inline double d2box(const Node& z, int id) {
    double dx=0,dy=0;
    if(p[id].x<z.xl) dx=(double)z.xl-p[id].x; else if(p[id].x>z.xr) dx=(double)p[id].x-z.xr;
    if(p[id].y<z.yl) dy=(double)z.yl-p[id].y; else if(p[id].y>z.yr) dy=(double)p[id].y-z.yr;
    return dx*dx+dy*dy;
}
static inline double dd2(int a,int b) { double x=(double)p[a].x-p[b].x, y=(double)p[a].y-p[b].y; return x*x+y*y; }
static void nearestRec(int q,int from,bool onlyPrime,int& best,double& val) {
    if(q<0 || tr[q].cnt==0 || (onlyPrime && tr[q].pcnt==0) || d2box(tr[q],from)>=val) return;
    int v=tr[q].id;
    if(alive[v] && (!onlyPrime || prime[v])) { double z=dd2(from,v); if(z<val) val=z,best=v; }
    int a=tr[q].left,b=tr[q].right;
    if(a>=0 && b>=0 && d2box(tr[b],from)<d2box(tr[a],from)) swap(a,b);
    nearestRec(a,from,onlyPrime,best,val); nearestRec(b,from,onlyPrime,best,val);
}
static int nearest(int from,bool onlyPrime) {
    int ans=-1; double v=numeric_limits<double>::infinity();
    nearestRec(0,from,onlyPrime,ans,v); return ans;
}
static void erasePoint(int v) {
    alive[v]=0;
    for(int q=where[v];q>=0;q=tr[q].parent) { --tr[q].cnt; if(prime[v]) --tr[q].pcnt; }
}
// Coordinate differences are bounded by 2e9, so squared double distances are safe.
static inline double dist(int a,int b) { double x=(double)p[a].x-p[b].x, y=(double)p[a].y-p[b].y; return sqrt(x*x+y*y); }
static inline double edgeCost(const vector<int>& r,int t) {
    double z=dist(r[t-1],r[t]);
    if(t%10==0 && !prime[r[t-1]]) z*=1.1;
    return z;
}
// Reversal changes precisely steps l .. r+1; weights remain tied to their global step.
static bool improveReverse(vector<int>& r,int l,int rr) {
    double before=0, after=0;
    for(int t=l;t<=rr+1;t++) before+=edgeCost(r,t);
    reverse(r.begin()+l,r.begin()+rr+1);
    for(int t=l;t<=rr+1;t++) after+=edgeCost(r,t);
    if(after+1e-7<before) return true;
    reverse(r.begin()+l,r.begin()+rr+1); return false;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    p.resize(n); for(auto &z:p) cin>>z.x>>z.y;
    prime.assign(n,true); if(n>0) prime[0]=false; if(n>1) prime[1]=false;
    for(int i=2;i*(long long)i<n;i++) if(prime[i]) for(long long j=1LL*i*i;j<n;j+=i) prime[j]=false;
    vector<int> a(n); iota(a.begin(),a.end(),0); where.resize(n); tr.reserve(n); build(a,0,n,-1);
    alive.assign(n,1); erasePoint(0);
    vector<int> route; route.reserve(n+1); route.push_back(0); int cur=0;
    for(int pos=1;pos<n;pos++) {
        int ordinary=nearest(cur,false), pick=ordinary;
        // P[pos] is the source of step pos+1, hence this is the only useful prime bias.
        if((pos+1)%10==0) {
            int pp=nearest(cur,true);
            if(pp>=0 && (ordinary<0 || dd2(cur,pp)<=dd2(cur,ordinary)*1.1025)) pick=pp;
        }
        if(pick<0) break;
        route.push_back(pick); erasePoint(pick); cur=pick;
    }
    route.push_back(0);
    // Large variable windows are a broader neighborhood than the incumbent's 2--4 reversals.
    for(int l=1;l<n;l++) for(int len=2;len<=16 && l+len-1<n;len++) improveReverse(route,l,l+len-1);
    // A short cleanup catches moves exposed by a preceding overlapping long reversal.
    for(int l=1;l<n;l++) for(int len=2;len<=6 && l+len-1<n;len++) improveReverse(route,l,l+len-1);
    cout<<n+1<<'\n'; for(int v:route) cout<<v<<'\n';
}
