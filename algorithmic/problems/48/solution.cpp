#include <bits/stdc++.h>
using namespace std;

struct Candidate {
    double r = -1;
    vector<array<double,3>> p;
};

static void consider(Candidate &best, Candidate &&x) {
    if (x.r > best.r) best = move(x);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;
    const double rt2 = sqrt(2.0), rt3 = sqrt(3.0);
    Candidate best;

    // FCC: integer triples of one parity.  Its shortest vector is
    // h*(1,1,0), while the face margin is h/sqrt(2).
    for (int L = 0;; ++L) {
        long long even = 0, odd = 0;
        for (int i = 0; i <= L; ++i)
            for (int j = 0; j <= L; ++j)
                for (int k = 0; k <= L; ++k)
                    ((i + j + k) & 1 ? odd : even)++;
        if (max(even, odd) < n) continue;
        // Select the larger parity class (this also avoids needless shells).
        int parity = even >= odd ? 0 : 1;
        double h = 1.0 / (L + rt2), r = h / rt2;
        Candidate x; x.r = r; x.p.reserve(n);
        for (int i = 0; i <= L && (int)x.p.size() < n; ++i)
            for (int j = 0; j <= L && (int)x.p.size() < n; ++j)
                for (int k = 0; k <= L && (int)x.p.size() < n; ++k)
                    if (((i+j+k)&1) == parity)
                        x.p.push_back({r+h*i, r+h*j, r+h*k});
        consider(best, move(x));
        break;
    }

    // BCC ablation: grid indices with all three parities equal.  In contrast
    // to FCC, its short vector is s*(1,1,1), which can give better shell
    // capacities in a finite cube.
    for (int M = 0;; ++M) {
        long long e = M/2 + 1, o = (M+1)/2;
        long long cap = e*e*e + o*o*o;
        if (cap < n) continue;
        double s = 1.0 / (M + rt3), r = rt3*s*0.5;
        Candidate x; x.r = r; x.p.reserve(n);
        for (int i = 0; i <= M && (int)x.p.size() < n; ++i)
            for (int j = 0; j <= M && (int)x.p.size() < n; ++j)
                for (int k = 0; k <= M && (int)x.p.size() < n; ++k)
                    if ((i&1) == (j&1) && (j&1) == (k&1))
                        x.p.push_back({r+s*i, r+s*j, r+s*k});
        consider(best, move(x));
        break;
    }

    // The balanced cubic grid is useful when a lattice shell is sparse.
    int q = 1;
    while (1LL*q*q*q < n) ++q;
    {
        Candidate x; x.r = 0.5/q; x.p.reserve(n);
        for (int i = 0; i < q && (int)x.p.size() < n; ++i)
            for (int j = 0; j < q && (int)x.p.size() < n; ++j)
                for (int k = 0; k < q && (int)x.p.size() < n; ++k)
                    x.p.push_back({(i+.5)/q, (j+.5)/q, (k+.5)/q});
        consider(best, move(x));
    }

    // Exact vertex codes remove the conspicuous small-n lattice losses.
    if (n >= 3 && n <= 4) {
        double r = rt2 / (2.0 + 2.0*rt2);
        int code[4] = {0,3,5,6};
        Candidate x; x.r = r;
        for (int a = 0; a < n; ++a) {
            int b = code[a];
            x.p.push_back({(b&1) ? 1-r : r, (b&2) ? 1-r : r,
                           (b&4) ? 1-r : r});
        }
        consider(best, move(x));
    }
    if (n >= 5 && n <= 8) {
        Candidate x; x.r = .25;
        for (int b = 0; b < n; ++b)
            x.p.push_back({(b&1) ? .75 : .25, (b&2) ? .75 : .25,
                           (b&4) ? .75 : .25});
        consider(best, move(x));
    }

    cout << setprecision(17);
    for (const auto &a : best.p) cout << a[0] << ' ' << a[1] << ' ' << a[2] << '\n';
}
