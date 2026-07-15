#include <bits/stdc++.h>
using namespace std;

struct Candidate {
    vector<array<int,3>> p;
    double d = 0, value = -1;
};

static Candidate makeLattice(int n, int kind) {
    // kind 0: FCC (x+y+z even), kind 1: BCC (all coordinate parities equal).
    // Integer coordinates are used only as a convenient, exactly representable lattice.
    int L = (int)ceil(cbrt((double)n / 2.0)) + 3;
    vector<array<int,3>> all;
    for (int x=-L; x<=L; ++x) for (int y=-L; y<=L; ++y) for (int z=-L; z<=L; ++z) {
        bool take = (kind == 0) ? (((x+y+z)&1) == 0)
                                : (((x&1) == (y&1)) && ((y&1) == (z&1)));
        if (take) all.push_back({x,y,z});
    }
    // A radial prefix is substantially less elongated than a lexicographic prefix.
    auto h = [](const array<int,3>& a) {
        unsigned long long v = (unsigned)(a[0]*11939 + a[1]*31337 + a[2]*6971);
        return v * 11400714819323198485ull;
    };
    sort(all.begin(), all.end(), [&](const auto& a, const auto& b) {
        int ra=a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
        int rb=b[0]*b[0]+b[1]*b[1]+b[2]*b[2];
        if (ra != rb) return ra < rb;
        return h(a) < h(b);
    });
    Candidate c;
    c.p.assign(all.begin(), all.begin()+n);
    c.d = (kind == 0 ? sqrt(2.0) : sqrt(3.0));
    return c;
}

static void evaluate(Candidate& c) {
    int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for (auto q:c.p) for (int k=0;k<3;k++) lo[k]=min(lo[k],q[k]), hi[k]=max(hi[k],q[k]);
    int e=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    // Scale the lattice into the inner cube [r,1-r]^3.  This is a certified
    // radius: all lattice separations are at least d.
    c.value=c.d/(2.0*(e+c.d));
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin>>n) || n<2) return 0;

    Candidate best;
    // The balanced cubic grid is a useful fallback, particularly for small n.
    int m=(int)ceil(cbrt((double)n));
    int a=m, b=m, c=(n + m*m-1)/(m*m);
    Candidate grid; grid.d=1.0;
    for (int x=0; x<a && (int)grid.p.size()<n; ++x)
        for (int y=0; y<b && (int)grid.p.size()<n; ++y)
            for (int z=0; z<c && (int)grid.p.size()<n; ++z)
                grid.p.push_back({x,y,z});
    evaluate(grid); best=grid;

    for (int kind=0; kind<2; ++kind) {
        Candidate q=makeLattice(n,kind);
        evaluate(q);
        if (q.value > best.value) best=move(q);
    }

    int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for (auto q:best.p) for (int k=0;k<3;k++) lo[k]=min(lo[k],q[k]), hi[k]=max(hi[k],q[k]);
    int e=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double r=best.d/(2.0*(e+best.d));
    double scale=1.0/(e+best.d);
    cout<<setprecision(17);
    for (auto q:best.p) {
        // Dimensions shorter than e simply retain extra empty room at the high face.
        cout << r + scale*(q[0]-lo[0]) << ' '
             << r + scale*(q[1]-lo[1]) << ' '
             << r + scale*(q[2]-lo[2]) << '\n';
    }
}
