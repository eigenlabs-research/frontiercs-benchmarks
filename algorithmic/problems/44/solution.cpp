#include <bits/stdc++.h>
using namespace std;
struct Pt { long long x,y; };
struct Node { int id,l=-1,r=-1; long long lx,rx,ly,ry; };
int n; vector<Pt> a; vector<Node> tr; vector<char> used, primeId;

int build(vector<int>& v, int lo, int hi, int dep) {
    if(lo>=hi) return -1;
    int m=(lo+hi)/2, ax=dep&1;
    nth_element(v.begin()+lo,v.begin()+m,v.begin()+hi,[&](int i,int j){
        return ax ? a[i].y<a[j].y : a[i].x<a[j].x;
    });
    int z=tr.size(); tr.push_back({}); tr[z].id=v[m];
    tr[z].l=build(v,lo,m,dep+1); tr[z].r=build(v,m+1,hi,dep+1);
    Node &q=tr[z]; q.lx=q.rx=a[q.id].x; q.ly=q.ry=a[q.id].y;
    for(int c: {q.l,q.r}) if(c>=0) {
        q.lx=min(q.lx,tr[c].lx); q.rx=max(q.rx,tr[c].rx);
        q.ly=min(q.ly,tr[c].ly); q.ry=max(q.ry,tr[c].ry);
    }
    return z;
}
inline double d2(int i,int j) { return hypot((double)a[i].x-a[j].x,(double)a[i].y-a[j].y); }
inline double boxd(int z,int id) {
    const Node&q=tr[z]; double dx=0,dy=0;
    if(a[id].x<q.lx) dx=(double)q.lx-a[id].x; else if(a[id].x>q.rx) dx=(double)a[id].x-q.rx;
    if(a[id].y<q.ly) dy=(double)q.ly-a[id].y; else if(a[id].y>q.ry) dy=(double)a[id].y-q.ry;
    return dx*dx+dy*dy;
}
void nearRec(int z,int cur,int &best,double &bd) {
    if(z<0 || boxd(z,cur)>=bd*bd) return;
    int id=tr[z].id;
    if(!used[id]) { double q=d2(cur,id); if(q<bd) bd=q,best=id; }
    int L=tr[z].l,R=tr[z].r;
    if(L>=0 && R>=0 && boxd(R,cur)<boxd(L,cur)) swap(L,R);
    nearRec(L,cur,best,bd); nearRec(R,cur,best,bd);
}
vector<int> kdTour() {
    vector<int> ans; ans.reserve(n+1); ans.push_back(0); used.assign(n,0); used[0]=1;
    int cur=0, remain=n-1, root=-1, built=0;
    while(remain) {
        // Rebuild after a quarter of the current tree was consumed: lazy deleted
        // bounding boxes otherwise make a static KD tree degenerate near the end.
        if(root<0 || built-remain > max(256,built/4)) {
            vector<int> v; v.reserve(remain); for(int i=1;i<n;i++) if(!used[i]) v.push_back(i);
            tr.clear(); tr.reserve(v.size()); root=build(v,0,(int)v.size(),0); built=remain;
        }
        int b=-1; double bd=1e300; nearRec(root,cur,b,bd);
        if(b<0) { for(int i=1;i<n;i++) if(!used[i]) { b=i; break; } }
        ans.push_back(b); used[b]=1; cur=b; --remain;
    }
    ans.push_back(0); return ans;
}
inline double mult(int source,int step) { return (step%10==0 && !primeId[source]) ? 1.1 : 1.0; }
double cost(const vector<int>& p) { double z=0; for(int t=1;t<=n;t++) z+=mult(p[t-1],t)*d2(p[t-1],p[t]); return z; }
// Exact delta for a reversal; all internal step multipliers are reconsidered.
double revDelta(const vector<int>& p,int l,int r) {
    double old=0,nw=0;
    for(int t=l;t<=r+1;t++) {
        old+=mult(p[t-1],t)*d2(p[t-1],p[t]);
        int u=(t-1<l ? p[t-1] : (t-1>r ? p[t-1] : p[l+r-(t-1)]));
        int v=(t<l ? p[t] : (t>r ? p[t] : p[l+r-t]));
        nw+=mult(u,t)*d2(u,v);
    }
    return nw-old;
}
void local2opt(vector<int>& p) {
    int W = n<4000 ? 36 : 12;
    for(int pass=0;pass<2;pass++) {
        bool any=false;
        for(int l=1;l<n-1;l++) for(int r=l+1;r<=min(n-1,l+W);r++) {
            if(revDelta(p,l,r)<-1e-7) { reverse(p.begin()+l,p.begin()+r+1); any=true; }
        }
        if(!any) break;
    }
}
double swapDelta(const vector<int>& p,int u,int v) {
    int q[4]={u,u+1,v,v+1}; double old=0,nw=0; sort(q,q+4); int last=-1;
    for(int k=0;k<4;k++) { int t=q[k]; if(t==last||t<1||t>n) continue; last=t;
        int x=p[t-1],y=p[t]; old+=mult(x,t)*d2(x,y);
        if(t-1==u) x=p[v]; else if(t-1==v) x=p[u];
        if(t==u) y=p[v]; else if(t==v) y=p[u];
        nw+=mult(x,t)*d2(x,y);
    } return nw-old;
}
void carrotSwaps(vector<int>& p) {
    // Only a short geometric neighborhood is considered, so prime placement cannot
    // trade an arbitrarily long detour for the 10% discount.
    for(int u=9;u<n;u+=10) if(!primeId[p[u]]) {
        int best=u; double gain=0;
        for(int v=max(1,u-28);v<=min(n-1,u+28);v++) if(primeId[p[v]]) {
            double q=swapDelta(p,u,v); if(q<gain) gain=q,best=v;
        }
        if(best!=u) swap(p[u],p[best]);
    }
}
uint64_t mortonKey(int id,long long minx,long long miny) {
    uint32_t x=(uint32_t)(a[id].x-minx), y=(uint32_t)(a[id].y-miny); uint64_t z=0;
    for(int b=0;b<31;b++) z|=((uint64_t)((x>>b)&1)<<(2*b))|((uint64_t)((y>>b)&1)<<(2*b+1));
    return z;
}
int main(){ ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>n)) return 0; a.resize(n); for(auto &p:a) cin>>p.x>>p.y;
    primeId.assign(n,true); if(n>0) primeId[0]=false; if(n>1) primeId[1]=false;
    for(int i=2;i*1LL*i<n;i++) if(primeId[i]) for(int j=i*i;j<n;j+=i) primeId[j]=false;
    vector<vector<int>> cand;
    auto makeSorted=[&](auto cmp){ vector<int> p(n+1); p[0]=p[n]=0; vector<int> v; for(int i=1;i<n;i++)v.push_back(i); sort(v.begin(),v.end(),cmp); for(int i=1;i<n;i++)p[i]=v[i-1]; cand.push_back(move(p)); };
    makeSorted([&](int i,int j){return a[i].x<a[j].x;});
    makeSorted([&](int i,int j){ if(a[i].y!=a[j].y)return a[i].y<a[j].y; return a[i].x<a[j].x; });
    long long minx=a[0].x,miny=a[0].y; for(auto&p:a) minx=min(minx,p.x),miny=min(miny,p.y);
    makeSorted([&](int i,int j){return mortonKey(i,minx,miny)<mortonKey(j,minx,miny);});
    cand.push_back(kdTour());
    int take=0; double best=cost(cand[0]); for(int i=1;i<(int)cand.size();i++){double z=cost(cand[i]);if(z<best)best=z,take=i;}
    vector<int> ans=move(cand[take]); carrotSwaps(ans); local2opt(ans); carrotSwaps(ans);
    cout<<n+1<<'\n'; for(int x:ans) cout<<x<<'\n';
}
