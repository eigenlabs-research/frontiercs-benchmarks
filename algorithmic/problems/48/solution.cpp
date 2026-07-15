#include <bits/stdc++.h>
using namespace std;

struct P { double x,y,z; };

static vector<P> cubicGrid(int n) {
    int lim = 1;
    while (1LL*lim*lim*lim < n) ++lim;
    int bx=lim, by=lim, bz=lim;
    double best=-1;
    for(int a=1;a<=lim;a++) for(int b=1;b<=lim;b++) for(int c=1;c<=lim;c++) {
        if(1LL*a*b*c < n) continue;
        double r=1.0/(2.0*max(a,max(b,c)));
        if(r>best) best=r,bx=a,by=b,bz=c;
    }
    vector<P> v;
    for(int i=0;i<bx && (int)v.size()<n;i++)
        for(int j=0;j<by && (int)v.size()<n;j++)
            for(int k=0;k<bz && (int)v.size()<n;k++)
                v.push_back({(i+.5)/bx,(j+.5)/by,(k+.5)/bz});
    return v;
}

// FCC lattice represented by integer triples with even coordinate sum.
static vector<P> fccGrid(int n) {
    int cap=1;
    while ((1LL*cap*cap*cap+1)/2 < n) ++cap;
    int ax=cap, ay=cap, az=cap;
    double best=-1;
    for(int a=1;a<=cap;a++) for(int b=1;b<=cap;b++) for(int c=1;c<=cap;c++) {
        long long cnt=0;
        for(int i=0;i<a;i++) for(int j=0;j<b;j++) {
            int first=((i+j)&1);
            if(first<c) cnt+=(c-first+1)/2;
        }
        if(cnt<n) continue;
        int mx=max(a,max(b,c));
        double d=1.0/(1.0+(mx-1)/sqrt(2.0));
        if(d>best+1e-15) best=d,ax=a,ay=b,az=c;
    }
    double h=best/sqrt(2.0);
    vector<array<int,3>> q;
    q.reserve(ax*ay*az/2+2);
    for(int i=0;i<ax;i++) for(int j=0;j<ay;j++) for(int k=0;k<az;k++)
        if(((i+j+k)&1)==0) q.push_back({i,j,k});

    // A deterministic randomized order avoids systematically selecting one corner;
    // all selected FCC points retain the same separation objective.
    uint64_t seed=1469598103934665603ULL ^ (uint64_t)n;
    for(int i=(int)q.size()-1;i>0;i--) {
        seed^=seed<<7; seed^=seed>>9; seed^=seed<<8;
        int j=seed%(i+1);
        swap(q[i],q[j]);
    }
    vector<P> v;
    v.reserve(n);
    for(int t=0;t<n;t++) v.push_back({best*.5+q[t][0]*h, best*.5+q[t][1]*h, best*.5+q[t][2]*h});
    return v;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if(!(cin>>n)) return 0;
    vector<P> ans;
    // Exact small symmetric configurations beat finite clipped lattices.
    if(n==2) {
        // Opposite body-diagonal points: 2r = sqrt(3)(1-2r).
        double r=sqrt(3.0)/(2.0+2.0*sqrt(3.0));
        ans={{r,r,r},{1-r,1-r,1-r}};
    } else if(n==3) {
        double r=sqrt(2.0)/(2.0+2.0*sqrt(2.0));
        ans={{r,r,r},{1-r,1-r,r},{1-r,r,1-r}};
    } else if(n==4) {
        double r=sqrt(2.0)/(2.0+2.0*sqrt(2.0));
        ans={{r,r,r},{r,1-r,1-r},{1-r,r,1-r},{1-r,1-r,r}};
    } else {
        vector<P> a=cubicGrid(n), b=fccGrid(n);
        // FCC has certified radius best/2; compare to the balanced grid value.
        double ra=1;
        for(auto&p:a) ra=min(ra,min({p.x,1-p.x,p.y,1-p.y,p.z,1-p.z}));
        double rb=1;
        for(auto&p:b) rb=min(rb,min({p.x,1-p.x,p.y,1-p.y,p.z,1-p.z}));
        // Pair spacing of FCC is twice its boundary lower bound, as is grid spacing.
        ans = (rb>ra+1e-14 ? b : a);
    }
    cout<<setprecision(17);
    for(const auto&p:ans) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    return 0;
}
