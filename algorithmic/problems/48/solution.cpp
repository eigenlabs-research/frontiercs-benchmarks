#include <bits/stdc++.h>
using namespace std;

using P = array<double,3>;
struct Candidate { double d = -1; vector<P> p; };

// Hexagonal close packing: triangular layers, with every other layer in the
// triangular holes.  perm chooses which crystal axis meets each cube axis.
static vector<P> hcp(double d, const array<int,3>& perm, int phase, int need,
                     bool stopEarly) {
    const double sy = sqrt(3.0) * .5 * d;
    const double sz = sqrt(2.0 / 3.0) * d;
    const double r = .5 * d;
    double sh = (phase ? .5 * d : 0.0);
    int M = (int)ceil(2.5 / d) + 4;
    vector<P> out;
    if (!stopEarly) out.reserve(need);
    for (int k = -M; k <= M; ++k) {
        double z = k * sz;
        double yoff = (k & 1) ? .5 * sy : 0.0;
        for (int j = -M; j <= M; ++j) {
            double y = j * sy + yoff;
            double xoff = ((j & 1) ? .5 * d : 0.0) + ((k & 1) ? .5 * d : 0.0);
            for (int i = -M; i <= M; ++i) {
                double q[3] = {i * d + xoff, y, z};
                P a;
                for (int t = 0; t < 3; ++t) a[t] = .5 + sh + q[perm[t]];
                if (a[0] >= r && a[0] <= 1-r && a[1] >= r && a[1] <= 1-r &&
                    a[2] >= r && a[2] <= 1-r) {
                    out.push_back(a);
                    if ((int)out.size() >= need && stopEarly) return out;
                }
            }
        }
    }
    return out;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // The two-sphere optimum is along a body diagonal; it is a useful small
    // case which a cube-aligned crystal cannot express.
    if (n == 2) {
        double d = sqrt(3.0) / (1.0 + sqrt(3.0));
        cout << setprecision(17);
        cout << d/2 << ' ' << d/2 << ' ' << d/2 << '\n';
        cout << 1-d/2 << ' ' << 1-d/2 << ' ' << 1-d/2 << '\n';
        return 0;
    }

    Candidate best;
    // Always retain the balanced rectangular-grid lower bound.
    int ba=1, bb=1, bc=n, bm=n;
    for (int a=1; a*a*a<=n; ++a)
        for (int b=a; a*b<=n; ++b) {
            int c=(n+a*b-1)/(a*b), m=max({a,b,c});
            if (m<bm) ba=a,bb=b,bc=c,bm=m;
        }
    double gd=1.0/bm;
    best.d=gd;
    for(int i=0;i<ba && (int)best.p.size()<n;i++)
        for(int j=0;j<bb && (int)best.p.size()<n;j++)
            for(int k=0;k<bc && (int)best.p.size()<n;k++)
                best.p.push_back({(i+.5)/ba,(j+.5)/bb,(k+.5)/bc});

    array<int,3> perm={0,1,2};
    do {
        for (int phase=0; phase<2; ++phase) {
            double lo=0, hi=.8;
            // Count-only binary search distinguishes boundary alignment of
            // each axis permutation and of two lattice translations.
            for(int it=0;it<20;++it) {
                double mid=(lo+hi)*.5;
                if ((int)hcp(mid,perm,phase,n,true).size()>=n) lo=mid;
                else hi=mid;
            }
            double safe=lo*(1.0-1e-9);
            auto pts=hcp(safe,perm,phase,n,true);
            if ((int)pts.size()>=n && safe>best.d) best={safe,move(pts)};
        }
    } while(next_permutation(perm.begin(),perm.end()));

    cout << setprecision(17);
    for (const P &p: best.p) cout << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
}
