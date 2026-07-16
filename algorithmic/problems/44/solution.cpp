#include <bits/stdc++.h>
using namespace std;

struct Point { long long x,y; };
static vector<Point> p;
static vector<char> primeId;
static int n;

static inline double edgeCost(int a, int b, int step) {
    double dx=(double)p[a].x-p[b].x, dy=(double)p[a].y-p[b].y;
    double z=hypot(dx,dy);
    if (step%10==0 && !primeId[a]) z*=1.1;
    return z;
}
static double cost(const vector<int>& r) {
    double ans=0;
    for(int t=1;t<=n;t++) ans+=edgeCost(r[t-1],r[t],t);
    return ans;
}
static vector<int> anchored(const vector<int>& cyc, bool backwards) {
    int at=find(cyc.begin(),cyc.end(),0)-cyc.begin();
    vector<int> r(n+1); r[0]=r[n]=0;
    for(int k=1;k<n;k++) {
        int q=backwards ? at-k : at+k;
        q%=n; if(q<0) q+=n;
        r[k]=cyc[q];
    }
    return r;
}
static void consider(const vector<int>& cyc, vector<int>& best, double& bestCost) {
    for(bool rev: {false,true}) {
        vector<int> r=anchored(cyc,rev);
        double v=cost(r);
        if(v<bestCost) { bestCost=v; best.swap(r); }
    }
}
// Standard Hilbert index, with coordinates already in [0,2^BITS).
static uint64_t hilbert(uint32_t x,uint32_t y) {
    const uint32_t S=1u<<20; uint64_t d=0;
    for(uint32_t s=S>>1;s;s>>=1) {
        uint32_t rx=(x&s)!=0, ry=(y&s)!=0;
        d+=(uint64_t)s*s*((3*rx)^ry);
        if(!ry) { if(rx) { x=s-1-x; y=s-1-y; } swap(x,y); }
    }
    return d;
}
static void improveShort(vector<int>& r, double& cur) {
    // This deliberately recomputes every changed, globally numbered edge: a usual
    // position-independent 2-opt delta is invalid under the carrot constraint.
    if(n>70000) return;
    static const int lens[]={2,3,4,6,8,12};
    for(int pass=0;pass<2;pass++) {
        bool any=false;
        for(int l=1;l<n;l++) for(int len:lens) {
            int rr=l+len-1; if(rr>=n) continue;
            double oldv=0,newv=0;
            for(int t=l;t<=rr+1;t++) {
                oldv+=edgeCost(r[t-1],r[t],t);
                auto at=[&](int pos) { return (pos>=l && pos<=rr) ? r[l+rr-pos] : r[pos]; };
                newv+=edgeCost(at(t-1),at(t),t);
            }
            if(newv+1e-7<oldv) {
                reverse(r.begin()+l,r.begin()+rr+1);
                cur+=newv-oldv; any=true;
            }
        }
        if(!any) break;
    }
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>n)) return 0;
    p.resize(n);
    for(auto &a:p) cin>>a.x>>a.y;
    primeId.assign(n,true); if(n>0) primeId[0]=false; if(n>1) primeId[1]=false;
    for(int i=2;(long long)i*i<n;i++) if(primeId[i])
        for(int j=i*i;j<n;j+=i) primeId[j]=false;

    vector<int> best; double bestCost=numeric_limits<double>::infinity();
    vector<int> v(n); iota(v.begin(),v.end(),0);
    consider(v,best,bestCost); // strengthened input-order guard

    long long minx=p[0].x,maxx=p[0].x,miny=p[0].y,maxy=p[0].y;
    for(auto a:p) { minx=min(minx,a.x);maxx=max(maxx,a.x);miny=min(miny,a.y);maxy=max(maxy,a.y); }
    vector<pair<uint64_t,int>> hv; hv.reserve(n);
    long double sx=maxx==minx?0:((long double)((1u<<20)-1)/(maxx-minx));
    long double sy=maxy==miny?0:((long double)((1u<<20)-1)/(maxy-miny));
    for(int i=0;i<n;i++) hv.push_back({hilbert((uint32_t)((p[i].x-minx)*sx),(uint32_t)((p[i].y-miny)*sy)),i});
    sort(hv.begin(),hv.end()); for(int i=0;i<n;i++) v[i]=hv[i].second; consider(v,best,bestCost);

    // Independent sweep guards.  The y-band family is the experiment's primary
    // representation; x bands protect strongly vertical or clustered inputs.
    for(int bands: {2,4,8,16,32}) if(bands<n) {
        vector<int> ord(n); iota(ord.begin(),ord.end(),0);
        sort(ord.begin(),ord.end(),[&](int a,int b){return p[a].y==p[b].y?a<b:p[a].y<p[b].y;});
        for(int b=0;b<bands;b++) {
            int L=(long long)b*n/bands, R=(long long)(b+1)*n/bands;
            sort(ord.begin()+L,ord.begin()+R,[&](int a,int c){return p[a].x<p[c].x;});
            if(b&1) reverse(ord.begin()+L,ord.begin()+R);
        }
        consider(ord,best,bestCost);
        iota(ord.begin(),ord.end(),0);
        for(int b=0;b<bands;b++) {
            int L=(long long)b*n/bands, R=(long long)(b+1)*n/bands;
            sort(ord.begin()+L,ord.begin()+R,[&](int a,int c){return p[a].y<p[c].y;});
            if(b&1) reverse(ord.begin()+L,ord.begin()+R);
        }
        consider(ord,best,bestCost);
    }
    improveShort(best,bestCost);
    cout<<n+1<<'\n'; for(int x:best) cout<<x<<'\n';
}
