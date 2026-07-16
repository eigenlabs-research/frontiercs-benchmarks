#include <bits/stdc++.h>
using namespace std;
using P = array<int,3>;

struct Candidate {
    vector<P> p;
    double step, gap, radius;
};

static Candidate scaled(vector<P> p, double gap) {
    int lo[3] = {INT_MAX,INT_MAX,INT_MAX}, hi[3] = {INT_MIN,INT_MIN,INT_MIN};
    for (auto a:p) for (int d=0; d<3; ++d) lo[d]=min(lo[d],a[d]), hi[d]=max(hi[d],a[d]);
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double step=1.0/(span+gap);
    return {move(p),step,gap,gap*step*0.5};
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n; if (!(cin>>n)) return 0;

    // Incumbent: a cropped FCC lattice (integer points of even coordinate sum).
    int q=1;
    while ((q*q*q+1LL)/2<n) ++q;
    vector<P> f; f.reserve(n);
    for(int x=0;x<q && (int)f.size()<n;++x)
        for(int y=0;y<q && (int)f.size()<n;++y)
            for(int z=0;z<q && (int)f.size()<n;++z)
                if(!((x+y+z)&1)) f.push_back({x,y,z});
    Candidate best=scaled(move(f),sqrt(2.0));

    // Ablation: BCC is less dense in bulk, but can have a better boundary cost
    // for a partly filled small block.  Coordinates are doubled BCC coordinates.
    q=1;
    auto bccCount=[](int s) { long long c=0; for(int x=0;x<s;++x) for(int y=0;y<s;++y) for(int z=0;z<s;++z)
        if ((x&1)==(y&1) && (y&1)==(z&1)) ++c; return c; };
    while(bccCount(q)<n) ++q;
    vector<P> b; b.reserve(n);
    for(int x=0;x<q && (int)b.size()<n;++x)
        for(int y=0;y<q && (int)b.size()<n;++y)
            for(int z=0;z<q && (int)b.size()<n;++z)
                if ((x&1)==(y&1) && (y&1)==(z&1)) b.push_back({x,y,z});
    // A mixed-parity adjacent pair has BCC's sqrt(3) nearest distance;
    // otherwise the selected points lie in one cubic sublattice, distance >= 2.
    int base=q+2; unordered_set<long long> have; have.reserve(2*n);
    auto id=[&](int x,int y,int z) { return (long long)x*base*base+(long long)y*base+z; };
    for(auto a:b) have.insert(id(a[0],a[1],a[2]));
    bool mixed=false;
    for(auto a:b) for(int dx:{-1,1}) for(int dy:{-1,1}) for(int dz:{-1,1}) {
        int x=a[0]+dx,y=a[1]+dy,z=a[2]+dz;
        if(x>=0&&y>=0&&z>=0&&x<q&&y<q&&z<q&&have.count(id(x,y,z))) mixed=true;
    }
    Candidate bc=scaled(move(b),mixed?sqrt(3.0):2.0);
    if(bc.radius>best.radius) best=move(bc);

    // Keep the stated balanced cubic grid as a guaranteed fallback.
    int m=1; while(1LL*m*m*m<n) ++m;
    double gr=1.0/(2.0*m);
    if(gr>best.radius) {
        best.p.clear(); best.step=1.0/m; best.gap=1.0; best.radius=gr;
        for(int x=0;x<m && (int)best.p.size()<n;++x)
            for(int y=0;y<m && (int)best.p.size()<n;++y)
                for(int z=0;z<m && (int)best.p.size()<n;++z) best.p.push_back({x,y,z});
    }
    cout<<setprecision(17);
    if(best.gap==1.0 && best.step==1.0/m && best.radius==gr) {
        for(auto a:best.p) cout<<(a[0]+.5)/m<<' '<<(a[1]+.5)/m<<' '<<(a[2]+.5)/m<<'\n';
    } else {
        int lo[3]={INT_MAX,INT_MAX,INT_MAX};
        for(auto a:best.p) for(int d=0;d<3;++d) lo[d]=min(lo[d],a[d]);
        double margin=best.radius;
        for(auto a:best.p) cout<<margin+best.step*(a[0]-lo[0])<<' '<<margin+best.step*(a[1]-lo[1])<<' '<<margin+best.step*(a[2]-lo[2])<<'\n';
    }
}
