#include <bits/stdc++.h>
using namespace std;

struct Point { long long x,y; };
struct Node {
    long long lx, rx, ly, ry;
    int l=-1, r=-1, lo=0, hi=0, cnt=0;
    bool leaf=false;
};

int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<Point> p(n);
    for (auto &q:p) cin >> q.x >> q.y;
    vector<char> prime(n, true);
    if(n) prime[0]=false;
    if(n>1) prime[1]=false;
    for(int i=2;i*(long long)i<n;i++) if(prime[i])
        for(int j=i*i;j<n;j+=i) prime[j]=false;

    vector<int> a(n), leafOf(n), parent;
    iota(a.begin(),a.end(),0);
    vector<Node> tr; tr.reserve(2*n);
    parent.reserve(2*n);
    function<int(int,int,int)> build = [&](int lo,int hi,int par)->int {
        int v=(int)tr.size(); tr.push_back(Node()); parent.push_back(par);
        Node &nd=tr[v]; nd.lo=lo; nd.hi=hi; nd.cnt=hi-lo;
        nd.lx=nd.ly=LLONG_MAX; nd.rx=nd.ry=LLONG_MIN;
        for(int k=lo;k<hi;k++) {
            auto q=p[a[k]];
            nd.lx=min(nd.lx,q.x); nd.rx=max(nd.rx,q.x);
            nd.ly=min(nd.ly,q.y); nd.ry=max(nd.ry,q.y);
        }
        if(hi-lo<=12) {
            nd.leaf=true;
            for(int k=lo;k<hi;k++) leafOf[a[k]]=v;
        } else {
            bool xsplit=(nd.rx-nd.lx >= nd.ry-nd.ly);
            int mid=(lo+hi)/2;
            nth_element(a.begin()+lo,a.begin()+mid,a.begin()+hi,[&](int u,int w){
                return xsplit ? p[u].x<p[w].x : p[u].y<p[w].y;
            });
            int L=build(lo,mid,v), R=build(mid,hi,v);
            tr[v].l=L; tr[v].r=R;
        }
        return v;
    };
    int root=build(0,n,-1);
    vector<char> used(n,0);
    auto erasePoint = [&](int id) {
        used[id]=1;
        for(int v=leafOf[id];v!=-1;v=parent[v]) --tr[v].cnt;
    };
    auto boxDist = [&](int v, int id) {
        double dx=0,dy=0;
        if(p[id].x<tr[v].lx) dx=(double)tr[v].lx-p[id].x;
        else if(p[id].x>tr[v].rx) dx=(double)p[id].x-tr[v].rx;
        if(p[id].y<tr[v].ly) dy=(double)tr[v].ly-p[id].y;
        else if(p[id].y>tr[v].ry) dy=(double)p[id].y-tr[v].ry;
        return dx*dx+dy*dy;
    };
    auto d2 = [&](int u,int v) { double x=(double)p[u].x-p[v].x, y=(double)p[u].y-p[v].y; return x*x+y*y; };
    erasePoint(0);
    vector<int> route; route.reserve(n+1); route.push_back(0);
    int cur=0;
    for(int taken=1;taken<n;taken++) {
        int best=-1; double bestd=numeric_limits<double>::infinity();
        function<void(int)> query = [&](int v) {
            if(tr[v].cnt==0 || boxDist(v,cur)>=bestd) return;
            if(tr[v].leaf) {
                for(int k=tr[v].lo;k<tr[v].hi;k++) { int z=a[k]; if(!used[z]) {
                    double q=d2(cur,z); if(q<bestd) bestd=q,best=z;
                }}
            } else {
                int L=tr[v].l,R=tr[v].r;
                double dl=boxDist(L,cur), dr=boxDist(R,cur);
                if(dl>dr) { swap(L,R); swap(dl,dr); }
                if(dl<bestd) query(L);
                if(dr<bestd) query(R);
            }
        };
        query(root);
        // best always exists because exactly n-1 active points remain.
        erasePoint(best); route.push_back(best); cur=best;
    }
    route.push_back(0);

    auto edge = [&](int step, int u, int v) {
        double z=hypot((double)p[u].x-p[v].x,(double)p[u].y-p[v].y);
        return z*((step%10==0 && !prime[u]) ? 1.1 : 1.0);
    };
    auto total = [&]() { double s=0; for(int t=1;t<=n;t++) s+=edge(t,route[t-1],route[t]); return s; };
    // Reversing a tour preserves ordinary internal distances; only its periodic carrot
    // sources, plus the two boundary edges, need to be considered.
    auto improve = [&]() {
        bool changed=false;
        for(int i=1;i<=n-2;i++) for(int len=2;len<=24 && i+len-1<=n-1;len++) {
            int j=i+len-1;
            double delta=edge(i,route[i-1],route[j])+edge(j+1,route[i],route[j+1])
                        -edge(i,route[i-1],route[i])-edge(j+1,route[j],route[j+1]);
            // Internal edges have step numbers i+1 through j.  Step i is the
            // left boundary and was already handled above (including its multiplier).
            int first=((i+10)/10)*10;
            for(int t=first;t<=j;t+=10) {
                int k=i+j-t+1;
                double dd=hypot((double)p[route[k-1]].x-p[route[k]].x,
                                (double)p[route[k-1]].y-p[route[k]].y);
                if(!prime[route[t-1]]) delta-=.1*dd;
                if(!prime[route[k]]) delta+=.1*dd;
            }
            if(delta < -1e-7) { reverse(route.begin()+i,route.begin()+j+1); changed=true; }
        }
        return changed;
    };
    // Compare the same geometric cycle in both directions, then make two bounded passes.
    vector<int> rev=route;
    reverse(rev.begin()+1,rev.end()-1);
    double forward=total();
    route=rev;
    double backward=total();
    if(forward<=backward) reverse(route.begin()+1,route.end()-1);
    improve(); improve();
    cout << n+1 << '\n';
    for(int x:route) cout << x << '\n';
}
