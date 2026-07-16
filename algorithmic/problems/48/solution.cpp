#include <bits/stdc++.h>
using namespace std;

struct Candidate {
    vector<array<int,3>> p;
    double d = 0, value = -1;
};

// Put a finite piece of a lattice in the cube.  If its shortest lattice
// vector is d and its largest coordinate span is R, this has radius d/2/(R+d).
static Candidate make_lattice(int n, int kind) {
    // kind 0: Z^3, 1: FCC (sum parity fixed), 2: BCC (all parities equal)
    double d = kind == 0 ? 1.0 : (kind == 1 ? sqrt(2.0) : sqrt(3.0));
    for (int R = 0;; ++R) {
        for (int phase = 0; phase < (kind == 1 ? 2 : 1); ++phase) {
            Candidate q; q.d = d;
            for (int x = 0; x <= R && (int)q.p.size() < n; ++x)
                for (int y = 0; y <= R && (int)q.p.size() < n; ++y)
                    for (int z = 0; z <= R && (int)q.p.size() < n; ++z) {
                        bool take;
                        if (kind == 0) take = true;
                        else if (kind == 1) take = ((x + y + z) & 1) == phase;
                        else take = ((x & 1) == (y & 1) && (y & 1) == (z & 1));
                        if (take) q.p.push_back({x,y,z});
                    }
            if ((int)q.p.size() != n) continue;
            // For small finite sets the selected subset can be farther apart
            // than the infinite lattice's nearest-neighbour distance.  Measure
            // that case exactly; for large n the lattice bound avoids O(n^2).
            if (n <= 256) {
                double md2 = 1e100;
                for (int i=0;i<n;i++) for (int j=0;j<i;j++) {
                    double dx=q.p[i][0]-q.p[j][0], dy=q.p[i][1]-q.p[j][1], dz=q.p[i][2]-q.p[j][2];
                    md2=min(md2, dx*dx+dy*dy+dz*dz);
                }
                q.d = sqrt(md2);
            }
            int lo[3] = {INT_MAX,INT_MAX,INT_MAX}, hi[3] = {INT_MIN,INT_MIN,INT_MIN};
            for (auto a : q.p) for (int k=0;k<3;k++) lo[k]=min(lo[k],a[k]), hi[k]=max(hi[k],a[k]);
            int span = max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
            q.value = q.d / (2.0 * (span + q.d));
            return q;
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // This is deliberately a finite-box comparison: asymptotic FCC density
    // alone is not a reliable choice when only a few lattice sites are used.
    Candidate best = make_lattice(n, 1);
    for (int kind : {0, 2}) {
        Candidate q = make_lattice(n, kind);
        if (q.value > best.value) best = move(q);
    }

    int lo[3] = {INT_MAX,INT_MAX,INT_MAX}, hi[3] = {INT_MIN,INT_MIN,INT_MIN};
    for (auto a : best.p) for (int k=0;k<3;k++) lo[k]=min(lo[k],a[k]), hi[k]=max(hi[k],a[k]);
    int largest = max({hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]});
    double scale = 1.0 / (largest + best.d);
    cout << setprecision(17);
    for (auto a : best.p) {
        for (int k=0;k<3;k++) {
            int span = hi[k]-lo[k];
            double v = (a[k]-lo[k])*scale + (1.0-span*scale)*0.5;
            if (k) cout << ' ';
            cout << v;
        }
        cout << '\n';
    }
}
