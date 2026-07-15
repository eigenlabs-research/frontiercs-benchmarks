#include <bits/stdc++.h>
using namespace std;

struct Candidate {
    vector<array<int,3>> p;
    double nearest;
};

static Candidate cubic(int n) {
    int s = 1;
    while (1LL*s*s*s < n) ++s;
    Candidate a;
    a.nearest = 1.0;
    for (int x=0; x<s && (int)a.p.size()<n; ++x)
        for (int y=0; y<s && (int)a.p.size()<n; ++y)
            for (int z=0; z<s && (int)a.p.size()<n; ++z)
                a.p.push_back({x,y,z});
    return a;
}

// Points with all three coordinates of equal parity form a BCC lattice.
// Its shortest nonzero vector has length sqrt(3).
static Candidate bcc(int n) {
    int q = 0;
    for (;; ++q) {
        long long ev = q/2 + 1, od = (q+1)/2;
        if (ev*ev*ev + od*od*od >= n) break;
    }
    Candidate a;
    a.nearest = sqrt(3.0);
    for (int x=0; x<=q && (int)a.p.size()<n; ++x)
        for (int y=0; y<=q && (int)a.p.size()<n; ++y)
            for (int z=0; z<=q && (int)a.p.size()<n; ++z)
                if ((x&1)==(y&1) && (y&1)==(z&1)) a.p.push_back({x,y,z});
    return a;
}

// Integer points whose coordinate sum is even form an FCC lattice.
static Candidate fcc(int n) {
    int q = 0;
    for (;; ++q) {
        long long side = q + 1;
        long long cnt = 0;
        for (int x=0; x<=q; ++x)
            for (int y=0; y<=q; ++y)
                cnt += (side + (((x+y)&1) ? 0 : 1)) / 2;
        if (cnt >= n) break;
    }
    Candidate a;
    a.nearest = sqrt(2.0);
    for (int x=0; x<=q && (int)a.p.size()<n; ++x)
        for (int y=0; y<=q && (int)a.p.size()<n; ++y)
            for (int z=0; z<=q && (int)a.p.size()<n; ++z)
                if (((x+y+z)&1)==0) a.p.push_back({x,y,z});
    return a;
}

static double radiusBound(const Candidate& a) {
    int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for (auto v:a.p) for (int d=0; d<3; ++d) {
        lo[d]=min(lo[d],v[d]); hi[d]=max(hi[d],v[d]);
    }
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    return a.nearest/(2.0*(span+a.nearest));
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin>>n)) return 0;

    // The two-sphere optimum is the cube diagonal, rather than an FCC edge.
    if (n==2) {
        double t=sqrt(3.0), r=t/(2.0*(1.0+t));
        cout<<setprecision(17);
        cout<<r<<' '<<r<<' '<<r<<'\n';
        cout<<1-r<<' '<<1-r<<' '<<1-r<<'\n';
        return 0;
    }

    Candidate a=fcc(n), b=bcc(n), c=cubic(n);
    const Candidate* use=&a;
    if (radiusBound(b) > radiusBound(*use)) use=&b;
    if (radiusBound(c) > radiusBound(*use)) use=&c;
    int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for (auto v:use->p) for (int d=0; d<3; ++d) {
        lo[d]=min(lo[d],v[d]); hi[d]=max(hi[d],v[d]);
    }
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double r=use->nearest/(2.0*(span+use->nearest));
    double scale=1.0/(span+use->nearest);
    cout<<setprecision(17);
    for (auto v:use->p)
        cout << r+scale*(v[0]-lo[0]) << ' '
             << r+scale*(v[1]-lo[1]) << ' '
             << r+scale*(v[2]-lo[2]) << '\n';
}
