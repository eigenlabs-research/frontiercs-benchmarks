#include <bits/stdc++.h>
using namespace std;

struct Point { long long x, y; };
struct Node {
    int id, l = -1, r = -1, par = -1, cnt = 0;
    bool live = true;
    long long xl, xr, yl, yr;
};

int n;
vector<Point> p;
vector<Node> tr;
vector<int> where;

void pull(int u) {
    Node &a = tr[u];
    a.cnt = 0;
    bool any = false;
    if (a.live) {
        a.cnt = 1; any = true;
        a.xl = a.xr = p[a.id].x;
        a.yl = a.yr = p[a.id].y;
    }
    for (int v : {a.l, a.r}) if (v != -1 && tr[v].cnt) {
        const Node &b = tr[v];
        if (!any) { a.xl=b.xl; a.xr=b.xr; a.yl=b.yl; a.yr=b.yr; any=true; }
        else {
            a.xl=min(a.xl,b.xl); a.xr=max(a.xr,b.xr);
            a.yl=min(a.yl,b.yl); a.yr=max(a.yr,b.yr);
        }
        a.cnt += b.cnt;
    }
}

int build(vector<int> &a, int lo, int hi, int dep, int parent) {
    if (lo >= hi) return -1;
    int mid = (lo + hi) >> 1, axis = dep & 1;
    nth_element(a.begin()+lo, a.begin()+mid, a.begin()+hi, [axis](int i, int j) {
        return axis ? p[i].y < p[j].y : p[i].x < p[j].x;
    });
    int u = (int)tr.size();
    tr.push_back(Node());
    tr[u].id = a[mid]; tr[u].par = parent;
    where[a[mid]] = u;
    tr[u].l = build(a, lo, mid, dep+1, u);
    tr[u].r = build(a, mid+1, hi, dep+1, u);
    pull(u);
    return u;
}

long double boxDist2(int u, long long x, long long y) {
    const Node &a = tr[u];
    long double dx = 0, dy = 0;
    if (x < a.xl) dx = (long double)a.xl-x;
    else if (x > a.xr) dx = (long double)x-a.xr;
    if (y < a.yl) dy = (long double)a.yl-y;
    else if (y > a.yr) dy = (long double)y-a.yr;
    return dx*dx + dy*dy;
}
void nearest(int u, int from, int &ans, long double &best) {
    if (u == -1 || tr[u].cnt == 0 || boxDist2(u,p[from].x,p[from].y) >= best) return;
    const Node &a = tr[u];
    if (a.live) {
        long double dx=(long double)p[a.id].x-p[from].x, dy=(long double)p[a.id].y-p[from].y;
        long double d=dx*dx+dy*dy;
        if (d < best || (d == best && a.id < ans)) best=d, ans=a.id;
    }
    int q[2]={a.l,a.r};
    if (q[0]!=-1 && q[1]!=-1 && boxDist2(q[1],p[from].x,p[from].y) < boxDist2(q[0],p[from].x,p[from].y)) swap(q[0],q[1]);
    nearest(q[0],from,ans,best); nearest(q[1],from,ans,best);
}
void eraseId(int id) {
    int u=where[id]; tr[u].live=false;
    for (; u!=-1; u=tr[u].par) pull(u);
}

bool isprime(int x, const vector<bool>& prime) { return x >= 2 && prime[x]; }
double routeCost(const vector<int>& r, const vector<bool>& prime) {
    double z=0;
    for (int t=1;t<=n;t++) {
        Point a=p[r[t-1]], b=p[r[t]];
        double d=hypot((double)a.x-b.x,(double)a.y-b.y);
        if (t%10==0 && !isprime(r[t-1],prime)) d*=1.1;
        z+=d;
    }
    return z;
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>n)) return 0;
    p.resize(n); for(auto &q:p) cin>>q.x>>q.y;
    vector<bool> prime(n,true); if(n>0) prime[0]=false; if(n>1) prime[1]=false;
    for(int i=2;i*(long long)i<n;i++) if(prime[i]) for(long long j=(long long)i*i;j<n;j+=i) prime[j]=false;

    vector<int> ids(n); iota(ids.begin(),ids.end(),0);
    tr.reserve(n); where.resize(n);
    int root=build(ids,0,n,0,-1);
    vector<int> greedy; greedy.reserve(n+1); greedy.push_back(0); eraseId(0);
    int cur=0;
    while(tr[root].cnt) {
        int v=-1; long double best=numeric_limits<long double>::infinity();
        // The carrot factor is determined by cur (the edge source), so it is
        // constant during this query; Euclidean nearest is the exact greedy choice.
        nearest(root,cur,v,best);
        greedy.push_back(v); eraseId(v); cur=v;
    }
    greedy.push_back(0);
    vector<int> baseline; baseline.reserve(n+1); for(int i=0;i<n;i++) baseline.push_back(i); baseline.push_back(0);
    const vector<int>& out = routeCost(greedy,prime) < routeCost(baseline,prime) ? greedy : baseline;
    cout << n+1 << '\n'; for(int v:out) cout << v << '\n';
}
