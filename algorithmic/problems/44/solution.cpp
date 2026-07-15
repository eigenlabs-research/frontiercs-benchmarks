#include <bits/stdc++.h>
using namespace std;
struct Pt { long long x,y; };

static uint64_t hilbert(unsigned x, unsigned y) {
    uint64_t d=0; const unsigned M=(1u<<21)-1;
    for (unsigned s=1u<<20;s;s>>=1) {
        unsigned rx=(x&s)!=0, ry=(y&s)!=0;
        d+=(uint64_t)s*s*((3*rx)^ry);
        if (!ry) { if(rx) { x=M-x; y=M-y; } swap(x,y); }
    }
    return d;
}

// A deletion-aware kd tree gives an exact nearest unvisited city.  This is a
// deliberately different tour family from the space-filling-curve candidates.
struct KD {
    struct Node { int id, l=-1, r=-1, par=-1, live=0; long long xl,xh,yl,yh; };
    const vector<Pt>& p; vector<int> a, where; vector<Node> tr;
    KD(const vector<Pt>& q):p(q),a(q.size()),where(q.size(),-1) {
        iota(a.begin(),a.end(),0); tr.reserve(q.size()); build(0,(int)a.size(),0,-1);
    }
    int build(int lo,int hi,int dep,int par) {
        if(lo>=hi) return -1;
        int m=(lo+hi)/2, ax=dep&1;
        nth_element(a.begin()+lo,a.begin()+m,a.begin()+hi,[&](int u,int v) {
            return ax ? p[u].y<p[v].y : p[u].x<p[v].x;
        });
        int z=tr.size(); tr.push_back({a[m],-1,-1,par,1,p[a[m]].x,p[a[m]].x,p[a[m]].y,p[a[m]].y});
        where[a[m]]=z;
        tr[z].l=build(lo,m,dep+1,z); tr[z].r=build(m+1,hi,dep+1,z);
        for(int c:{tr[z].l,tr[z].r}) if(c>=0) {
            tr[z].live+=tr[c].live;
            tr[z].xl=min(tr[z].xl,tr[c].xl); tr[z].xh=max(tr[z].xh,tr[c].xh);
            tr[z].yl=min(tr[z].yl,tr[c].yl); tr[z].yh=max(tr[z].yh,tr[c].yh);
        }
        return z;
    }
    long double lower(int z,long long x,long long y) const {
        const Node& q=tr[z]; long double dx=0,dy=0;
        if(x<q.xl) dx=(long double)q.xl-x; else if(x>q.xh) dx=(long double)x-q.xh;
        if(y<q.yl) dy=(long double)q.yl-y; else if(y>q.yh) dy=(long double)y-q.yh;
        return dx*dx+dy*dy;
    }
    void erase(int id) { for(int z=where[id];z>=0;z=tr[z].par) --tr[z].live; }
    void query(int z,long long x,long long y,int& ans,long double& best) const {
        if(z<0 || !tr[z].live || lower(z,x,y)>=best) return;
        const Node& q=tr[z];
        long double dx=(long double)p[q.id].x-x, dy=(long double)p[q.id].y-y, d=dx*dx+dy*dy;
        // A node's own point is live iff it has not been erased; live counts
        // alone do not identify it, so its self point is tested through where.
        // Erased points have where=-2 (set by the constructor below's caller).
        if(where[q.id]>=0 && d<best) best=d,ans=q.id;
        int u=q.l,v=q.r;
        if(u>=0 && v>=0 && lower(v,x,y)<lower(u,x,y)) swap(u,v);
        query(u,x,y,ans,best); query(v,x,y,ans,best);
    }
    int nearest(long long x,long long y) const { int ans=-1; long double b=numeric_limits<long double>::infinity(); query(0,x,y,ans,b); return ans; }
};

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0; vector<Pt> p(n);
    long long minx=LLONG_MAX,maxx=LLONG_MIN,miny=LLONG_MAX,maxy=LLONG_MIN;
    for(auto& q:p){cin>>q.x>>q.y;minx=min(minx,q.x);maxx=max(maxx,q.x);miny=min(miny,q.y);maxy=max(maxy,q.y);}
    vector<char> prime(n,true); prime[0]=false; if(n>1) prime[1]=false;
    for(int i=2;1LL*i*i<n;i++) if(prime[i]) for(int j=i*i;j<n;j+=i) prime[j]=false;
    auto dist=[&](int u,int v){long double dx=(long double)p[u].x-p[v].x,dy=(long double)p[u].y-p[v].y;return sqrtl(dx*dx+dy*dy);};
    auto score=[&](const vector<int>& r){long double z=0;for(int t=1;t<=n;t++)z+=(t%10==0&&!prime[r[t-1]]?1.1L:1)*dist(r[t-1],r[t]);return z;};
    vector<int> best; long double bestv=numeric_limits<long double>::infinity();
    auto consider=[&](vector<int> r){long double z=score(r);if(z<bestv)bestv=z,best=move(r);};
    auto cycle=[&](const vector<int>& c,bool rev){int at=find(c.begin(),c.end(),0)-c.begin();vector<int> r{0};r.reserve(n+1);for(int k=1;k<n;k++)r.push_back(c[(at+(rev?-k:k)+n)%n]);r.push_back(0);consider(move(r));};
    vector<int> mono(n);iota(mono.begin(),mono.end(),0);cycle(mono,0);cycle(mono,1);
    const unsigned LIM=(1u<<21)-1; vector<unsigned> X(n),Y(n);
    for(int i=0;i<n;i++){X[i]=maxx==minx?LIM/2:(unsigned)((unsigned long long)(p[i].x-minx)*LIM/(maxx-minx));Y[i]=maxy==miny?LIM/2:(unsigned)((unsigned long long)(p[i].y-miny)*LIM/(maxy-miny));}
    for(int sw=0;sw<2;sw++)for(int fx=0;fx<2;fx++)for(int fy=0;fy<2;fy++){
        vector<pair<uint64_t,int>> o;o.reserve(n);for(int i=0;i<n;i++){unsigned x=X[i],y=Y[i];if(sw)swap(x,y);if(fx)x=LIM-x;if(fy)y=LIM-y;o.push_back({hilbert(x,y),i});}sort(o.begin(),o.end());vector<int> c(n);for(int i=0;i<n;i++)c[i]=o[i].second;cycle(c,0);cycle(c,1);
    }
    // Exact greedy nearest-neighbor chain from the required depot.
    KD kd(p); vector<int> nn; nn.reserve(n+1); nn.push_back(0); kd.erase(0); kd.where[0]=-2;
    for(int k=1;k<n;k++){int v=kd.nearest(p[nn.back()].x,p[nn.back()].y); if(v<0) break; nn.push_back(v);kd.erase(v);kd.where[v]=-2;} nn.push_back(0); if((int)nn.size()==n+1)consider(move(nn));
    auto edge=[&](int t){return (t%10==0&&!prime[best[t-1]]?1.1L:1)*dist(best[t-1],best[t]);};
    auto extra=[&](int t){return t%10==0&&!prime[best[t-1]]?.1L*dist(best[t-1],best[t]):0.L;};
    // Preserve the incumbent's short exact carrot-aware 2-opt polish.
    for(int pass=0;pass<3;pass++){bool any=false;for(int l=1;l<n;l++)for(int r=l+1;r<=min(n-1,l+7);r++){
        // The changed internal edges are l+1..r; edge l is already a boundary edge.
        long double before=edge(l)+edge(r+1);for(int t=(l/10+1)*10;t<=r;t+=10)before+=extra(t);
        reverse(best.begin()+l,best.begin()+r+1);long double after=edge(l)+edge(r+1);for(int t=(l/10+1)*10;t<=r;t+=10)after+=extra(t);
        if(after+1e-9L<before)bestv+=after-before,any=true;else reverse(best.begin()+l,best.begin()+r+1);
    }if(!any)break;}
    cout<<n+1<<'\n';for(int v:best)cout<<v<<'\n';
}
