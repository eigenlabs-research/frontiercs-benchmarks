#include <bits/stdc++.h>
using namespace std;

struct Packing {
    long double r;
    int type; // 0: cubic, 1: FCC, 2: small simplex
};

static __int128 fccCount(long long m) {
    // Integer triples in [0,m]^3 whose coordinate sum is even.
    __int128 a = m + 1;
    __int128 e = (a + 1) / 2, o = a / 2;
    return e * e * e + 3 * e * o * o;
}

static __int128 bccCount(long long m) {
    __int128 a = m + 1;
    __int128 e = (a + 1) / 2, o = a / 2;
    return e * e * e + o * o * o;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    long long n;
    if (!(cin >> n)) return 0;

    Packing best{-1, 0};
    if (n <= 4) {
        long double r = (n == 2) ? sqrtl(3.0L) / (2.0L * (1.0L + sqrtl(3.0L)))
                                 : 1.0L / (2.0L + sqrtl(2.0L));
        best = {r, 2};
    }

    long long t = 1;
    while ((__int128)t * t * t < n) ++t;
    Packing cubic{1.0L / (2.0L * t), 0};
    if (cubic.r > best.r) best = cubic;

    long long m = 0;
    while (fccCount(m) < n) ++m;
    Packing fcc{1.0L / (2.0L + sqrtl(2.0L) * m), 1};
    if (fcc.r > best.r) best = fcc;

    long long bm = 0;
    while (bccCount(bm) < n) ++bm;
    Packing bcc{1.0L / (2.0L + 2.0L * bm / sqrtl(3.0L)), 3};
    if (bcc.r > best.r) best = bcc;

    cout << fixed << setprecision(17);
    if (best.type == 2) {
        long double r = best.r, s = 1.0L - r;
        if (n == 2) {
            cout << (double)r << ' ' << (double)r << ' ' << (double)r << '\n';
            cout << (double)s << ' ' << (double)s << ' ' << (double)s << '\n';
        } else {
            // Vertices of a regular tetrahedron selected from cube corners.
            const int v[4][3] = {{0,0,0},{0,1,1},{1,0,1},{1,1,0}};
            for (int q = 0; q < n; ++q)
                cout << (double)(v[q][0] ? s : r) << ' '
                     << (double)(v[q][1] ? s : r) << ' '
                     << (double)(v[q][2] ? s : r) << '\n';
        }
    } else if (best.type == 0) {
        long long made = 0;
        for (long long i = 0; i < t && made < n; ++i)
            for (long long j = 0; j < t && made < n; ++j)
                for (long long k = 0; k < t && made < n; ++k, ++made)
                    cout << (double)((i + 0.5L) / t) << ' '
                         << (double)((j + 0.5L) / t) << ' '
                         << (double)((k + 0.5L) / t) << '\n';
    } else if (best.type == 1) {
        long double r = best.r, h = sqrtl(2.0L) * r;
        long long made = 0;
        for (long long i = 0; i <= m && made < n; ++i)
            for (long long j = 0; j <= m && made < n; ++j)
                for (long long k = 0; k <= m && made < n; ++k) {
                    if ((i + j + k) & 1) continue;
                    cout << (double)(r + h * i) << ' ' << (double)(r + h * j) << ' '
                         << (double)(r + h * k) << '\n';
                    ++made;
                }
    } else {
        long double r = best.r, h = 2.0L * r / sqrtl(3.0L);
        long long made = 0;
        for (long long i = 0; i <= bm && made < n; ++i)
            for (long long j = 0; j <= bm && made < n; ++j)
                for (long long k = 0; k <= bm && made < n; ++k) {
                    if (!((i & 1) == (j & 1) && (j & 1) == (k & 1))) continue;
                    cout << (double)(r + h * i) << ' ' << (double)(r + h * j) << ' '
                         << (double)(r + h * k) << '\n';
                    ++made;
                }
    }
    return 0;
}
