#include <bits/stdc++.h>
using namespace std;

using P = array<double,3>;
static void emit(const vector<P>& a) {
    cout << setprecision(17);
    for (auto p:a) cout << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
}

int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // A few genuinely better small packings.  They are formulae rather than
    // output tables, and cover the cases where a bulk lattice has ugly edges.
    if (n == 2) {
        double r=sqrt(3.0)/(2.0*(1.0+sqrt(3.0)));
        emit({P{r,r,r},P{1-r,1-r,1-r}}); return 0;
    }
    if (n == 3) {
        // Three mutually touching corner spheres (a robust small-n packing).
        const double r=.25;
        emit({P{r,r,r},P{1-r,r,r},P{r,1-r,r}}); return 0;
    }
    if (n == 4) {
        double r=1.0/(2.0+sqrt(2.0));
        emit({P{r,r,r},P{r,1-r,1-r},P{1-r,r,1-r},P{1-r,1-r,r}}); return 0;
    }
    if (n == 5) {
        const double r=.25;
        emit({P{r,r,r},P{r,1-r,r},P{1-r,r,r},P{1-r,1-r,r},P{.5,.5,1-r}}); return 0;
    }
    if (n == 8) {
        const double r=.25; vector<P> a;
        for(int x=0;x<2;x++) for(int y=0;y<2;y++) for(int z=0;z<2;z++)
            a.push_back({x?1-r:r,y?1-r:r,z?1-r:r});
        emit(a); return 0;
    }
    if (n == 9) { // one BCC cell: eight corners and its body centre
        double r=sqrt(3.0)/(4.0+2.0*sqrt(3.0)); vector<P> a;
        for(int x=0;x<2;x++) for(int y=0;y<2;y++) for(int z=0;z<2;z++)
            a.push_back({x?1-r:r,y?1-r:r,z?1-r:r});
        a.push_back({.5,.5,.5}); emit(a); return 0;
    }

    // Face-centred cubic lattice.  Its nearest-neighbour density is much
    // better than the reference Cartesian grid.  Build a near-cubic finite
    // block, then take a deterministic dispersed subset; finally translate
    // and scale the actual subset, so partially filled final shells do not
    // waste cube width.
    int L=1;
    while (((1LL*L*L*L+1)/2) < n) ++L;
    vector<array<int,3>> v;
    v.reserve((L*L*L+1)/2);
    for(int x=0;x<L;x++) for(int y=0;y<L;y++) for(int z=0;z<L;z++)
        if (((x+y+z)&1)==0) v.push_back({x,y,z});

    // Partial Fisher--Yates gives a reproducible, spatially uncorrelated
    // subset without paying to sort or optimize quadratically.
    uint64_t state=0x9e3779b97f4a7c15ULL ^ (uint64_t)n;
    auto rnd = [&]() { state^=state<<7; state^=state>>9; return state; };
    for(int i=0;i<n;i++) {
        int j=i+(int)(rnd()%(v.size()-i));
        swap(v[i],v[j]);
    }
    int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for(int i=0;i<n;i++) for(int d=0;d<3;d++) lo[d]=min(lo[d],v[i][d]),hi[d]=max(hi[d],v[i][d]);
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double step=1.0/(span+sqrt(2.0)), r=step/sqrt(2.0);
    vector<P> ans; ans.reserve(n);
    for(int i=0;i<n;i++) ans.push_back({r+step*(v[i][0]-lo[0]),r+step*(v[i][1]-lo[1]),r+step*(v[i][2]-lo[2])});
    emit(ans);
}
