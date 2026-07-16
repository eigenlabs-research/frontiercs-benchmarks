#include <bits/stdc++.h>
using namespace std;

struct P { double x, y, z; };

// Points with integer-coordinate sum even on a cubic grid are an FCC lattice.
// Its nearest-neighbour distance is sqrt(2)*h, so h=sqrt(2)*r
// makes every pair at least 2r apart.
static vector<P> fcc(double r, int need, double ax, double ay, double az) {
    const double h = sqrt(2.0) * r;
    const double hi = 1.0 - r;
    int lim = (int)ceil((1.0 - 2.0*r) / h) + 3;
    vector<P> v;
    v.reserve(need);
    for (int i = -2; i <= lim; ++i) {
        double x = r + h * (i + ax);
        if (x < r - 1e-13 || x > hi + 1e-13) continue;
        for (int j = -2; j <= lim; ++j) {
            double y = r + h * (j + ay);
            if (y < r - 1e-13 || y > hi + 1e-13) continue;
            for (int k = -2; k <= lim; ++k) {
                if ((i + j + k) & 1) continue;
                double z = r + h * (k + az);
                if (z < r - 1e-13 || z > hi + 1e-13) continue;
                v.push_back({x,y,z});
                if ((int)v.size() == need) return v;
            }
        }
    }
    return v;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // The antipodal body-diagonal pair is better than a cropped lattice for n=2.
    if (n == 2) {
        double r = sqrt(3.0) / (2.0 * (1.0 + sqrt(3.0))) * (1.0 - 1e-12);
        cout << setprecision(17);
        cout << r << ' ' << r << ' ' << r << '\n';
        cout << 1-r << ' ' << 1-r << ' ' << 1-r << '\n';
        return 0;
    }

    const double phases[3] = {0.0, 0.25, 0.5};
    double lo = 0.0, hi = 0.5;
    double ba=0, bb=0, bc=0;
    // A fixed finite set of translations substantially reduces boundary waste.
    for (int it = 0; it < 38; ++it) {
        double mid = (lo + hi) * 0.5;
        bool ok = false;
        double ca=0, cb=0, cc=0;
        for (double a : phases) for (double b : phases) for (double c : phases) {
            if ((int)fcc(mid, n, a, b, c).size() == n) {
                ok = true; ca=a; cb=b; cc=c;
                goto found;
            }
        }
        found:
        if (ok) { lo = mid; ba=ca; bb=cb; bc=cc; }
        else hi = mid;
    }
    // Stay inside all equality constraints despite decimal output rounding.
    vector<P> ans = fcc(lo * (1.0 - 1e-11), n, ba, bb, bc);
    // Small crops can waste an FCC layer; never lose to the general cubic grid.
    int m = 1;
    while (1LL*m*m*m < n) ++m;
    if ((int)ans.size() != n || 0.5 / m >= lo * (1.0 - 1e-11)) {
        ans.clear(); ans.reserve(n);
        for (int i=0; i<m && (int)ans.size()<n; ++i)
            for (int j=0; j<m && (int)ans.size()<n; ++j)
                for (int k=0; k<m && (int)ans.size()<n; ++k)
                    ans.push_back({(i+.5)/m,(j+.5)/m,(k+.5)/m});
    }
    cout << setprecision(17);
    for (const P &p : ans) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
}
