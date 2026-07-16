#include <bits/stdc++.h>
using namespace std;

using P = array<double,3>;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    const double rt2 = sqrt(2.0), rt3 = sqrt(3.0);
    vector<P> best;
    double bestR = -1;
    auto take = [&](vector<P> v, double guaranteed) {
        if ((int)v.size() >= n && guaranteed > bestR) {
            v.resize(n);
            best = move(v);
            bestR = guaranteed;
        }
    };

    // Cubic lattice: useful for the small sizes where a complete FCC shell
    // would contain too many unused sites.
    int m = 1;
    while (1LL*m*m*m < n) ++m;
    {
        double r = 0.5 / m;
        vector<P> v; v.reserve(1LL*m*m*m);
        for (int i=0;i<m;i++) for (int j=0;j<m;j++) for (int k=0;k<m;k++)
            v.push_back({(i+.5)/m, (j+.5)/m, (k+.5)/m});
        take(move(v), r);
    }

    // Face-centred cubic lattice.  Its nearest-neighbour distance is sqrt(2)a.
    int f = 1;
    auto fcount = [](int q) { return (q*q*q + 1) / 2; };
    while (fcount(f) < n) ++f;
    {
        double a = 1.0 / (f - 1 + rt2);
        double r = a / rt2;
        vector<P> v; v.reserve(fcount(f));
        for (int i=0;i<f;i++) for (int j=0;j<f;j++) for (int k=0;k<f;k++)
            if (((i+j+k)&1)==0) v.push_back({r+i*a, r+j*a, r+k*a});
        take(move(v), r);
    }

    // Body-centred cubic is occasionally a better finite-size shell.
    int b = 1;
    auto bcount = [](int q) { return q*q*q + (q-1)*(q-1)*(q-1); };
    while (bcount(b) < n) ++b;
    if (b >= 2) {
        double a = 1.0 / (b - 1 + rt3/2.0);
        double r = rt3*a/4.0;
        vector<P> v; v.reserve(bcount(b));
        for (int i=0;i<b;i++) for (int j=0;j<b;j++) for (int k=0;k<b;k++)
            v.push_back({r+i*a, r+j*a, r+k*a});
        for (int i=0;i<b-1;i++) for (int j=0;j<b-1;j++) for (int k=0;k<b-1;k++)
            v.push_back({r+(i+.5)*a, r+(j+.5)*a, r+(k+.5)*a});
        take(move(v), r);
    }

    // Best binary codes in the cube for the exceptional very small cases.
    if (n == 2 || n == 3 || n == 4) {
        vector<array<int,3>> bits;
        double d;
        if (n == 2) {
            bits = {{{0,0,0},{1,1,1}}}; d=rt3;
        } else if (n == 3) {
            bits = {{{0,0,0},{1,1,0},{1,0,1}}}; d=rt2;
        } else {
            bits = {{{0,0,0},{0,1,1},{1,0,1},{1,1,0}}}; d=rt2;
        }
        double r = d/(2.0*(1.0+d));
        vector<P> v;
        for (auto q: bits) v.push_back({r+q[0]*(1-2*r), r+q[1]*(1-2*r), r+q[2]*(1-2*r)});
        take(move(v), r);
    }

    cout << setprecision(17);
    for (const auto &p: best) cout << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
}
