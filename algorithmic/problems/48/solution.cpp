#include <bits/stdc++.h>
using namespace std;

using P = array<double, 3>;

static double radiusOf(const vector<P>& p) {
    double r = 1.0;
    for (const auto& a : p)
        for (int d = 0; d < 3; ++d)
            r = min(r, min(a[d], 1.0 - a[d]));
    double md2 = 1e100;
    for (int i = 0; i < (int)p.size(); ++i)
        for (int j = 0; j < i; ++j) {
            double x = p[i][0] - p[j][0];
            double y = p[i][1] - p[j][1];
            double z = p[i][2] - p[j][2];
            md2 = min(md2, x*x + y*y + z*z);
        }
    return min(r, 0.5 * sqrt(md2));
}

// A centered portion of an integer lattice.  parity=0 is FCC (i+j+k even),
// parity=1 is BCC (all three coordinate parities equal).
static vector<P> lattice(int n, bool bcc) {
    int L = 0;
    for (;; ++L) {
        int cnt = 0;
        for (int x = 0; x <= L; ++x)
            for (int y = 0; y <= L; ++y)
                for (int z = 0; z <= L; ++z) {
                    bool ok = bcc ? ((x&1) == (y&1) && (y&1) == (z&1))
                                  : (((x+y+z)&1) == 0);
                    cnt += ok;
                }
        if (cnt >= n) break;
    }
    vector<array<int,3>> q;
    q.reserve(n);
    // Taking successive diagonal layers keeps partial outer boxes reasonably spread.
    for (int layer = 0; layer <= 3*L && (int)q.size() < n; ++layer)
        for (int x = 0; x <= L && (int)q.size() < n; ++x)
            for (int y = 0; y <= L && (int)q.size() < n; ++y)
                for (int z = 0; z <= L && (int)q.size() < n; ++z) {
                    if (x+y+z != layer) continue;
                    bool ok = bcc ? ((x&1) == (y&1) && (y&1) == (z&1))
                                  : (((x+y+z)&1) == 0);
                    if (ok) q.push_back({x,y,z});
                }
    int lo[3] = {q[0][0], q[0][1], q[0][2]}, hi[3] = {lo[0],lo[1],lo[2]};
    for (auto a : q) for (int d = 0; d < 3; ++d) lo[d] = min(lo[d],a[d]), hi[d] = max(hi[d],a[d]);
    int R = max({hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2]});
    double near = bcc ? sqrt(3.0) : sqrt(2.0);
    double scale = 1.0 / (R + near);
    double margin = near * scale * 0.5;
    vector<P> ans;
    ans.reserve(n);
    for (auto a : q) ans.push_back({margin + scale*(a[0]-lo[0]), margin + scale*(a[1]-lo[1]), margin + scale*(a[2]-lo[2])});
    return ans;
}

static vector<P> cartesian(int n) {
    int d = 1;
    while (1LL*d*d*d < n) ++d;
    vector<P> a;
    a.reserve(n);
    for (int x = 0; x < d && (int)a.size() < n; ++x)
        for (int y = 0; y < d && (int)a.size() < n; ++y)
            for (int z = 0; z < d && (int)a.size() < n; ++z)
                a.push_back({(x+.5)/d, (y+.5)/d, (z+.5)/d});
    return a;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    vector<vector<P>> cand;
    cand.push_back(cartesian(n));
    cand.push_back(lattice(n, false));
    cand.push_back(lattice(n, true));
    // Exact low-cardinality simplex arrangements are useful where a lattice cell is coarse.
    if (n == 2) {
        double r = sqrt(3.0)/(2.0*(sqrt(3.0)+1.0));
        cand.push_back({P{r,r,r}, P{1-r,1-r,1-r}});
    } else if (n == 3 || n == 4) {
        double r = sqrt(2.0)/(2.0*(sqrt(2.0)+1.0));
        vector<P> a = {{r,r,r},{r,1-r,1-r},{1-r,r,1-r},{1-r,1-r,r}};
        a.resize(n);
        cand.push_back(a);
    }
    int best = 0;
    double br = radiusOf(cand[0]);
    for (int i = 1; i < (int)cand.size(); ++i) {
        double rr = radiusOf(cand[i]);
        if (rr > br) br = rr, best = i;
    }
    cout << setprecision(17);
    for (const auto& p : cand[best]) cout << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
}
