#include <bits/stdc++.h>
using namespace std;

struct Candidate {
    vector<array<double,3>> p;
    double radius = -1;
};

// Points are integer coordinates in a smallest cube containing at least n
// points of the indicated Bravais lattice.  After scaling, the lattice's
// nearest-neighbour distance and the clearance to every face are equal.
static Candidate lattice(int n, int kind) {
    auto ok = [kind](int x, int y, int z) {
        if (kind == 0) return ((x + y + z) & 1) == 0; // FCC
        if (kind == 1) return ((x & 1) == (y & 1) && (y & 1) == (z & 1)); // BCC
        return true; // simple cubic
    };
    int s = 0;
    for (;; ++s) {
        int cnt = 0;
        for (int x=0; x<=s; ++x) for (int y=0; y<=s; ++y) for (int z=0; z<=s; ++z)
            cnt += ok(x,y,z);
        if (cnt >= n) break;
    }
    vector<array<int,3>> q;
    q.reserve(n);
    for (int x=0; x<=s && (int)q.size()<n; ++x)
        for (int y=0; y<=s && (int)q.size()<n; ++y)
            for (int z=0; z<=s && (int)q.size()<n; ++z)
                if (ok(x,y,z)) q.push_back({x,y,z});
    int lo[3]={INT_MAX,INT_MAX,INT_MAX}, hi[3]={INT_MIN,INT_MIN,INT_MIN};
    for (auto a:q) for(int d=0;d<3;++d) lo[d]=min(lo[d],a[d]), hi[d]=max(hi[d],a[d]);
    int span=max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double sep = kind==0 ? sqrt(2.0) : (kind==1 ? sqrt(3.0) : 1.0);
    double h=1.0/(span+sep), r=sep*h*0.5;
    Candidate ans; ans.radius=r; ans.p.reserve(n);
    for(auto a:q) ans.p.push_back({r+(a[0]-lo[0])*h, r+(a[1]-lo[1])*h, r+(a[2]-lo[2])*h});
    return ans;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    Candidate best=lattice(n,0); // retain FCC as the asymptotically dense option
    for (int kind=1; kind<=2; ++kind) {
        Candidate cur=lattice(n,kind);
        if (cur.radius > best.radius) best=move(cur);
    }
    cout << setprecision(17);
    for (auto a:best.p) cout << a[0] << ' ' << a[1] << ' ' << a[2] << '\n';
}
