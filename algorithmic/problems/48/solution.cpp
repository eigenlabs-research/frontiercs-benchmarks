#include <bits/stdc++.h>
using namespace std;

struct Cand {
    double r = -1.0;
    vector<array<double,3>> p;
};

static void try_raw(Cand &best, int n, const vector<array<double,3>> &raw, double dmin,
                    double cx, double cy, double cz) {
    if ((int)raw.size() < n || n <= 0) return;
    vector<pair<double,int>> ord;
    ord.reserve(raw.size());
    for (int i = 0; i < (int)raw.size(); ++i) {
        double dx = raw[i][0] - cx, dy = raw[i][1] - cy, dz = raw[i][2] - cz;
        ord.push_back({dx*dx + dy*dy + dz*dz, i});
    }
    sort(ord.begin(), ord.end(), [](const auto &a, const auto &b){
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    });

    vector<array<double,3>> sel;
    sel.reserve(n);
    double mn[3] = {1e100,1e100,1e100}, mx[3] = {-1e100,-1e100,-1e100};
    for (int t = 0; t < n; ++t) {
        auto q = raw[ord[t].second];
        sel.push_back(q);
        for (int k = 0; k < 3; ++k) mn[k] = min(mn[k], q[k]), mx[k] = max(mx[k], q[k]);
    }
    double B = max({mx[0]-mn[0], mx[1]-mn[1], mx[2]-mn[2]});
    if (B <= 0) return;
    double r = dmin / (2.0 * (B + dmin));
    if (r <= best.r) return;

    double s = 1.0 / (B + dmin);
    vector<array<double,3>> out;
    out.reserve(n);
    for (auto q: sel) {
        array<double,3> v;
        for (int k = 0; k < 3; ++k) {
            double mid = 0.5 * (mn[k] + mx[k]);
            v[k] = 0.5 + (q[k] - mid) * s;
            if (v[k] < 0 && v[k] > -1e-14) v[k] = 0;
            if (v[k] > 1 && v[k] < 1+1e-14) v[k] = 1;
        }
        out.push_back(v);
    }
    best.r = r;
    best.p.swap(out);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    Cand best;

    // Exact/simple strong small cases.
    if (n == 2) {
        double r = sqrt(3.0) / (2.0 * (1.0 + sqrt(3.0)));
        cout << setprecision(17);
        cout << r << ' ' << r << ' ' << r << "\n";
        cout << 1-r << ' ' << 1-r << ' ' << 1-r << "\n";
        return 0;
    }
    if (n == 3 || n == 4) {
        double r = sqrt(2.0) / (2.0 * (1.0 + sqrt(2.0)));
        vector<array<double,3>> v = {
            {r,r,r}, {r,1-r,1-r}, {1-r,r,1-r}, {1-r,1-r,r}
        };
        cout << setprecision(17);
        for (int i = 0; i < n; ++i) cout << v[i][0] << ' ' << v[i][1] << ' ' << v[i][2] << "\n";
        return 0;
    }
    if (n == 5) {
        // Five-point maximin seed in the unit cube: two opposite diagonal
        // endpoints plus the three cyclic permutations of (3/4,3/4,0).
        // Its raw minimum separation is 3/(2*sqrt(2)), which beats any
        // grid with side spacing 1 after the standard margin scaling.
        const double D = 3.0 / (2.0 * sqrt(2.0));
        const double r = D / (2.0 * (1.0 + D));
        const double s = 1.0 - 2.0 * r;
        vector<array<double,3>> raw = {
            {0.0,0.0,0.0}, {1.0,1.0,1.0},
            {0.75,0.75,0.0}, {0.75,0.0,0.75}, {0.0,0.75,0.75}
        };
        cout << setprecision(17);
        for (auto q: raw) cout << r + s*q[0] << ' ' << r + s*q[1] << ' ' << r + s*q[2] << "\n";
        return 0;
    }

    // Enumerate compact finite clusters from three related lattices.  Points closest
    // to the box centre are retained, then uniformly scaled so lattice contacts and
    // cube-face contacts have the same guaranteed radius.
    int L = 0;
    while ((long long)(L+1)*(L+1)*(L+1) < 8LL*n) ++L;
    L += 6;
    L = min(L, 45); // official n<=4096 needs far less; keeps runtime small.

    for (int a = 1; a <= L; ++a) for (int b = 1; b <= L; ++b) for (int c = 1; c <= L; ++c) {
        int mnD = min({a,b,c}), mxD = max({a,b,c});
        if (mxD - mnD > 3) continue; // cube-like clusters are the useful ones here

        long long cubicCnt = 1LL*(a+1)*(b+1)*(c+1);
        if (cubicCnt >= n && cubicCnt <= 2LL*n + 100) {
            vector<array<double,3>> raw;
            raw.reserve((size_t)cubicCnt);
            for (int i=0;i<=a;i++) for (int j=0;j<=b;j++) for (int k=0;k<=c;k++)
                raw.push_back({(double)i,(double)j,(double)k});
            try_raw(best, n, raw, 1.0, a/2.0, b/2.0, c/2.0);
        }

        long long vol = 1LL*(a+1)*(b+1)*(c+1);
        if (vol/2 + 2 >= n && vol <= 3LL*n + 200) {
            vector<array<double,3>> raw;
            raw.reserve((size_t)(vol/2+2));
            for (int i=0;i<=a;i++) for (int j=0;j<=b;j++) for (int k=0;k<=c;k++)
                if (((i+j+k)&1) == 0) raw.push_back({(double)i,(double)j,(double)k});
            try_raw(best, n, raw, sqrt(2.0), a/2.0, b/2.0, c/2.0);
        }

        if (vol/4 + 4 >= n && vol <= 6LL*n + 500) {
            vector<array<double,3>> raw;
            raw.reserve((size_t)(vol/4+4));
            for (int i=0;i<=a;i++) for (int j=0;j<=b;j++) for (int k=0;k<=c;k++)
                if ( (i&1)==(j&1) && (j&1)==(k&1) ) raw.push_back({(double)i,(double)j,(double)k});
            try_raw(best, n, raw, sqrt(3.0), a/2.0, b/2.0, c/2.0);
        }
    }

    // Safety fallback: balanced cubic grid, also covers any unexpectedly huge n.
    if (best.p.empty()) {
        int q = 1; while (1LL*q*q*q < n) ++q;
        double r = 1.0/(2.0*q);
        cout << setprecision(17);
        int cnt = 0;
        for (int i=0;i<q && cnt<n;i++) for (int j=0;j<q && cnt<n;j++) for (int k=0;k<q && cnt<n;k++,cnt++)
            cout << (i+0.5)/q << ' ' << (j+0.5)/q << ' ' << (k+0.5)/q << "\n";
        return 0;
    }

    cout << setprecision(17);
    for (auto &v: best.p) cout << v[0] << ' ' << v[1] << ' ' << v[2] << "\n";
    return 0;
}
