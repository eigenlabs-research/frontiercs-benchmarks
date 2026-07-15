#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };
struct Node {
    int l=-1, r=-1, parent=-1, lo=0, hi=0, cnt=0;
    double xmin, xmax, ymin, ymax;
};

struct KDTree {
    const vector<Point>& p;
    vector<int> ord, leaf;
    vector<Node> tr;
    static constexpr int B = 20;
    KDTree(const vector<Point>& pp, vector<int> ids) : p(pp), ord(move(ids)), leaf(pp.size(), -1) {
        if (!ord.empty()) build(0, (int)ord.size(), -1);
    }
    int build(int lo, int hi, int parent) {
        int v=tr.size(); tr.push_back(Node());
        Node &n=tr[v]; n.lo=lo; n.hi=hi; n.parent=parent; n.cnt=hi-lo;
        n.xmin=n.ymin=1e100; n.xmax=n.ymax=-1e100;
        for(int i=lo;i<hi;i++) {
            const auto &q=p[ord[i]];
            n.xmin=min(n.xmin,(double)q.x); n.xmax=max(n.xmax,(double)q.x);
            n.ymin=min(n.ymin,(double)q.y); n.ymax=max(n.ymax,(double)q.y);
        }
        if(hi-lo<=B) { for(int i=lo;i<hi;i++) leaf[ord[i]]=v; return v; }
        bool ax=(n.xmax-n.xmin >= n.ymax-n.ymin);
        int mid=(lo+hi)/2;
        nth_element(ord.begin()+lo,ord.begin()+mid,ord.begin()+hi,[&](int a,int b){
            return ax ? p[a].x<p[b].x : p[a].y<p[b].y;
        });
        int L=build(lo,mid,v), R=build(mid,hi,v);
        tr[v].l=L; tr[v].r=R;
        return v;
    }
    double boxd(int v, int a) const {
        const Node &n=tr[v]; double x=p[a].x,y=p[a].y, dx=0,dy=0;
        if(x<n.xmin) dx=n.xmin-x; else if(x>n.xmax) dx=x-n.xmax;
        if(y<n.ymin) dy=n.ymin-y; else if(y>n.ymax) dy=y-n.ymax;
        return dx*dx+dy*dy;
    }
    void queryRec(int v, int a, int &best, double &bd) const {
        const Node &n=tr[v]; if(!n.cnt || boxd(v,a)>=bd) return;
        if(n.l<0) {
            for(int i=n.lo;i<n.hi;i++) { int b=ord[i]; if(leaf[b]>=0) { // leaf[] remains; cnt is authoritative only for this leaf
                double dx=(double)p[a].x-p[b].x, dy=(double)p[a].y-p[b].y, d=dx*dx+dy*dy;
                // A removed point is identified by a negative marker in pointLeaf.
                if(active[b] && d<bd) bd=d,best=b;
            }}
            return;
        }
        int u=n.l,w=n.r; if(boxd(w,a)<boxd(u,a)) swap(u,w);
        queryRec(u,a,best,bd); queryRec(w,a,best,bd);
    }
    vector<char> active;
    void initialize() { active.assign(p.size(),0); for(int x:ord) active[x]=1; }
    int nearest(int a) const { int b=-1; double d=1e300; if(!tr.empty()) queryRec(0,a,b,d); return b; }
    void erase(int a) {
        if(a<0 || !active[a]) return;
        active[a]=0; int v=leaf[a];
        while(v>=0) { --tr[v].cnt; v=tr[v].parent; }
    }
};

static inline bool prime(int x, const vector<char>& isp) { return isp[x]; }
static double edgeCost(const vector<int>& r, int t, const vector<Point>& p, const vector<char>& isp) {
    int a=r[t-1], b=r[t];
    double d=hypot((double)p[a].x-p[b].x,(double)p[a].y-p[b].y);
    return d*((t%10==0 && !isp[a]) ? 1.1 : 1.0);
}
static double cost(const vector<int>& r, const vector<Point>& p, const vector<char>& isp) {
    double z=0; for(int t=1;t<(int)r.size();t++) z+=edgeCost(r,t,p,isp); return z;
}
static void adjacentImprove(vector<int>& r, const vector<Point>& p, const vector<char>& isp) {
    int n=r.size()-1;
    for(int pass=0;pass<5;pass++) {
        bool any=false;
        for(int i=1;i<=n-2;i++) {
            double old=edgeCost(r,i,p,isp)+edgeCost(r,i+1,p,isp)+edgeCost(r,i+2,p,isp);
            swap(r[i],r[i+1]);
            double nw=edgeCost(r,i,p,isp)+edgeCost(r,i+1,p,isp)+edgeCost(r,i+2,p,isp);
            if(nw+1e-8<old) any=true;
            else swap(r[i],r[i+1]);
        }
        if(!any) break;
    }
}

static inline double directEdge(int a, int b, int t, const vector<Point>& p, const vector<char>& isp) {
    double d=hypot((double)p[a].x-p[b].x,(double)p[a].y-p[b].y);
    return d*((t%10==0 && !isp[a]) ? 1.1 : 1.0);
}
// Reverse positions [i+1,j].  Ordinary internal edge lengths cancel in pairs;
// only the two boundary edges and internal tenth edges need to be examined.
static double reverseDelta(const vector<int>& r, int i, int j, const vector<Point>& p, const vector<char>& isp) {
    double old=directEdge(r[i],r[i+1],i+1,p,isp)+directEdge(r[j],r[j+1],j+1,p,isp);
    double nw=directEdge(r[i],r[j],i+1,p,isp)+directEdge(r[i+1],r[j+1],j+1,p,isp);
    int k=((i+2+9)/10)*10;
    for(;k<=j;k+=10) {
        old+=directEdge(r[k-1],r[k],k,p,isp);
        // At new position k, the reversed segment contains these old vertices.
        nw+=directEdge(r[i+j+2-k],r[i+j+1-k],k,p,isp);
    }
    return nw-old;
}
// A bounded, exact 2-opt pass.  A moderately wider x/y candidate graph catches
// crossings which are local geometrically but not in the greedy tour.  Most
// candidates are discarded by an unpenalized boundary screen before paying for
// the exact carrot-aware reversal delta.
static void boundedTwoOpt(vector<int>& r, const vector<Point>& p, const vector<char>& isp) {
    int n=(int)r.size()-1;
    vector<int> byY(n), yrank(n), pos(n);
    iota(byY.begin(),byY.end(),0);
    sort(byY.begin(),byY.end(),[&](int a,int b) { return p[a].y==p[b].y ? p[a].x<p[b].x : p[a].y<p[b].y; });
    for(int k=0;k<n;k++) yrank[byY[k]]=k;
    for(int pass=0;pass<2;pass++) {
        for(int k=0;k<n;k++) pos[r[k]]=k;
        bool any=false;
        for(int i=0;i<=n-2;i++) {
            int a=r[i];
            int bestJ=-1;
            double best=0;
            auto consider=[&](int b) {
                if(b<0 || b>=n) return;
                int j=pos[b];
                if(j<=i+1 || j>=n || j-i>300) return;
                // Ordinary edge crossings are necessary for the geometric
                // improvements sought here.  This cheap screen keeps long
                // reversal penalty accounting off the overwhelmingly bad pairs.
                auto euclid=[&](int u,int v) {
                    return hypot((double)p[u].x-p[v].x,(double)p[u].y-p[v].y);
                };
                if(euclid(r[i],r[j])+euclid(r[i+1],r[j+1]) >=
                   euclid(r[i],r[i+1])+euclid(r[j],r[j+1])-1e-8) return;
                double d=reverseDelta(r,i,j,p,isp);
                if(d<best-1e-7) best=d,bestJ=j;
            };
            for(int z=1;z<=8;z++) { consider(a-z); consider(a+z); }
            int q=yrank[a];
            for(int z=1;z<=8;z++) { if(q>=z) consider(byY[q-z]); if(q+z<n) consider(byY[q+z]); }
            if(bestJ>=0) {
                reverse(r.begin()+i+1,r.begin()+bestJ+1);
                // Only this interval moved, so refresh exactly those inverse entries.
                for(int k=i+1;k<=bestJ;k++) pos[r[k]]=k;
                any=true;
            }
        }
        if(!any) break;
    }
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<Point> p(n); for(auto &q:p) cin>>q.x>>q.y;
    vector<char> isp(max(2,n),true); isp[0]=false; if(n>1) isp[1]=false;
    for(int i=2;i*i<n;i++) if(isp[i]) for(int j=i*i;j<n;j+=i) isp[j]=false;
    vector<int> all, prs; all.reserve(n-1);
    for(int i=1;i<n;i++) { all.push_back(i); if(isp[i]) prs.push_back(i); }
    KDTree general(p,all), primes(p,prs); general.initialize(); primes.initialize();
    vector<int> r; r.reserve(n+1); r.push_back(0); int cur=0;
    while((int)r.size()<n) {
        int a=general.nearest(cur); if(a<0) break;
        // The city put at position 9,19,... becomes the source of the following penalty edge.
        if(((int)r.size()%10)==9 && !primes.tr.empty()) {
            int q=primes.nearest(cur);
            if(q>=0) {
                double da=hypot((double)p[cur].x-p[a].x,(double)p[cur].y-p[a].y);
                double dq=hypot((double)p[cur].x-p[q].x,(double)p[cur].y-p[q].y);
                if(dq <= da*1.06) a=q;
            }
        }
        r.push_back(a); general.erase(a); if(isp[a]) primes.erase(a); cur=a;
    }
    r.push_back(0);
    adjacentImprove(r,p,isp);
    boundedTwoOpt(r,p,isp);
    adjacentImprove(r,p,isp);
    vector<int> base; base.reserve(n+1); for(int i=0;i<n;i++) base.push_back(i); base.push_back(0);
    if(cost(base,p,isp)<cost(r,p,isp)) r.swap(base);
    cout<<n+1<<'\n'; for(int x:r) cout<<x<<'\n';
}
