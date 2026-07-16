#include <bits/stdc++.h>
using namespace std;

struct City { long long x,y; };
static vector<City> a;
static vector<char> primeId;

static inline double edgeCost(int u,int v,int step) {
    double d=hypot((double)a[u].x-a[v].x,(double)a[u].y-a[v].y);
    return (step%10==0 && !primeId[u]) ? 1.1*d : d;
}
static double routeCost(const vector<int>& p) {
    double z=0; int n=(int)p.size()-1;
    for(int i=0;i<n;i++) z+=edgeCost(p[i],p[i+1],i+1);
    return z;
}

// Turn a cyclic geometric ordering into a tour rooted at city zero.
static vector<int> rooted(const vector<int>& q, bool backward) {
    int n=q.size(), at=find(q.begin(),q.end(),0)-q.begin();
    vector<int> p; p.reserve(n+1); p.push_back(0);
    for(int k=1;k<n;k++) p.push_back(q[(at+(backward?-k:k)+n)%n]);
    p.push_back(0); return p;
}
static unsigned long long morton(unsigned x,unsigned y) {
    unsigned long long r=0;
    for(int b=0;b<31;b++) { r|=(unsigned long long)((x>>b)&1)<<(2*b); r|=(unsigned long long)((y>>b)&1)<<(2*b+1); }
    return r;
}
static unsigned long long hilbert(unsigned x,unsigned y) {
    // Standard xy-to-Hilbert conversion; coordinates have at most 31 useful bits.
    unsigned long long d=0;
    for(unsigned s=1u<<30;s;s>>=1) {
        unsigned rx=(x&s)?1:0, ry=(y&s)?1:0;
        d += (unsigned long long)s*(unsigned long long)s*((3*rx)^ry);
        if(!ry) { if(rx) { x=(~x); y=(~y); } swap(x,y); }
    }
    return d;
}

// Exact small-window 2-opt: unlike a usual Euclidean 2-opt test, every changed
// directed edge is rescored at its new global step, including carrot positions.
static void boundedTwoOpt(vector<int>& p) {
    const int n=(int)p.size()-1, W=9;
    for(int pass=0;pass<2;pass++) {
        for(int i=0;i+2<n;i++) {
            bool changed=false;
            for(int len=2;len<=W && i+len<n;len++) {
                int j=i+len;
                double oldc=0,newc=0;
                for(int e=i;e<=j;e++) oldc+=edgeCost(p[e],p[e+1],e+1);
                // After reversal of vertices i+1..j, only this contiguous edge block changes.
                newc+=edgeCost(p[i],p[j],i+1);
                for(int e=i+1;e<j;e++) newc+=edgeCost(p[i+j-e+1],p[i+j-e],e+1);
                newc+=edgeCost(p[i+1],p[j+1],j+1);
                if(newc+1e-8<oldc) {
                    reverse(p.begin()+i+1,p.begin()+j+1);
                    changed=true; break;
                }
            }
            if(changed) continue;
        }
    }
}

// A reversal cannot always extract one outlying city from an otherwise good
// spatial order.  This exact Or-opt pass tests such moves without assuming
// that the ordinary (position-independent) TSP delta is valid under carrots.
static void boundedRelocate(vector<int>& p) {
    const int n=(int)p.size()-1, W=12;
    for(int pass=0; pass<2; ++pass) {
        for(int i=1; i<n; ++i) {
            bool changed=false;
            for(int j=max(0,i-W); j<=min(n-1,i+W); ++j) {
                if(j==i || j==i-1) continue; // inserting after itself/predecessor
                int from=(i<j ? i-1 : j), to=(i<j ? j : i);
                double oldc=0, newc=0;
                for(int e=from; e<=to; ++e) oldc+=edgeCost(p[e],p[e+1],e+1);
                auto after=[&](int pos) {
                    if(i<j) {
                        if(pos>=i && pos<j) return p[pos+1];
                        if(pos==j) return p[i];
                    } else {
                        if(pos==j+1) return p[i];
                        if(pos>=j+2 && pos<=i) return p[pos-1];
                    }
                    return p[pos];
                };
                for(int e=from; e<=to; ++e) newc+=edgeCost(after(e),after(e+1),e+1);
                if(newc+1e-8<oldc) {
                    if(i<j) rotate(p.begin()+i,p.begin()+i+1,p.begin()+j+1);
                    else rotate(p.begin()+j+1,p.begin()+i,p.begin()+i+1);
                    changed=true;
                    break;
                }
            }
            if(changed) continue;
        }
    }
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    a.resize(n); for(auto &c:a) cin>>c.x>>c.y;
    primeId.assign(n,true); if(n>0) primeId[0]=false; if(n>1) primeId[1]=false;
    for(int i=2;i*(long long)i<n;i++) if(primeId[i]) for(int j=i*i;j<n;j+=i) primeId[j]=false;
    vector<vector<int>> orders;
    vector<int> q(n); iota(q.begin(),q.end(),0);
    orders.push_back(q); // input/x sweep, retained as a safe baseline candidate
    auto byY=q; sort(byY.begin(),byY.end(),[&](int i,int j){if(a[i].y!=a[j].y)return a[i].y<a[j].y;return a[i].x<a[j].x;}); orders.push_back(move(byY));
    long long xmin=a[0].x,xmax=a[0].x,ymin=a[0].y,ymax=a[0].y;
    long double cx=0,cy=0; for(auto c:a){xmin=min(xmin,c.x);xmax=max(xmax,c.x);ymin=min(ymin,c.y);ymax=max(ymax,c.y);cx+=c.x;cy+=c.y;} cx/=n;cy/=n;
    auto scaled=[&](long long v,long long lo,long long hi){ return hi==lo?0u:(unsigned)(((long double)(v-lo)*2147483647.0L)/(hi-lo)); };
    vector<unsigned long long> mk(n),hk(n); vector<long double> ang(n);
    for(int i=0;i<n;i++){unsigned X=scaled(a[i].x,xmin,xmax),Y=scaled(a[i].y,ymin,ymax);mk[i]=morton(X,Y);hk[i]=hilbert(X,Y);ang[i]=atan2((long double)a[i].y-cy,(long double)a[i].x-cx);}
    auto mo=q; sort(mo.begin(),mo.end(),[&](int i,int j){return mk[i]!=mk[j]?mk[i]<mk[j]:i<j;}); orders.push_back(move(mo));
    auto hi=q; sort(hi.begin(),hi.end(),[&](int i,int j){return hk[i]!=hk[j]?hk[i]<hk[j]:i<j;}); orders.push_back(move(hi));
    auto ar=q; sort(ar.begin(),ar.end(),[&](int i,int j){return ang[i]!=ang[j]?ang[i]<ang[j]:i<j;}); orders.push_back(move(ar));
    vector<int> best; double bestCost=numeric_limits<double>::infinity();
    for(const auto &o:orders) for(bool rev:{false,true}) { auto p=rooted(o,rev); double c=routeCost(p); if(c<bestCost){bestCost=c;best=move(p);} }
    boundedTwoOpt(best);
    boundedRelocate(best);
    cout<<n+1<<'\n'; for(int v:best) cout<<v<<'\n';
}
