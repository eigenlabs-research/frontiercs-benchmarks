#include <bits/stdc++.h>
using namespace std;

struct Choice {
    double r = -1.0;
    int kind = 0;                 // 0 = cubic, 1 = FCC, 2 = BCC
    int px = 0, py = 0, pz = 0;
};

static long long oneDim(double r, double a, double shift, double base) {
    double lo = (r - shift) / a - base;
    double hi = (1.0 - r - shift) / a - base;
    long long l = (long long)ceil(lo - 1e-12);
    long long h = (long long)floor(hi + 1e-12);
    return max(0LL, h - l + 1);
}

static long long capacity(double r, int kind, int px, int py, int pz, int P) {
    double a = (kind == 1 ? 2.0 * sqrt(2.0) * r : 4.0 * r / sqrt(3.0));
    double sx = a * px / P, sy = a * py / P, sz = a * pz / P;
    static const double fcc[4][3] = {{0,0,0},{0,.5,.5},{.5,0,.5},{.5,.5,0}};
    static const double bcc[2][3] = {{0,0,0},{.5,.5,.5}};
    const double (*bases)[3] = kind == 1 ? fcc : bcc;
    int nb = kind == 1 ? 4 : 2;
    long long total = 0;
    for (int q = 0; q < nb; ++q) {
        long long x = oneDim(r,a,sx,bases[q][0]);
        long long y = oneDim(r,a,sy,bases[q][1]);
        long long z = oneDim(r,a,sz,bases[q][2]);
        total += x*y*z;
        if (total > 1000000000LL) return total;
    }
    return total;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // A translated Bravais lattice is searched rather than fixing its origin.
    // This makes the finite-cube boundary treatment part of the optimization.
    const int P = 16;
    Choice best;
    int d = 1;
    while (1LL*d*d*d < n) ++d;
    best.r = 0.5 / d;

    for (int kind = 1; kind <= 2; ++kind) {
        for (int px = 0; px < P; ++px) for (int py = 0; py < P; ++py) for (int pz = 0; pz < P; ++pz) {
            double lo = 0.0, hi = 0.5;
            for (int it = 0; it < 43; ++it) {
                double mid = (lo + hi) * .5;
                if (capacity(mid,kind,px,py,pz,P) >= n) lo = mid;
                else hi = mid;
            }
            if (lo > best.r + 1e-12) best = {lo,kind,px,py,pz};
        }
    }

    if (best.kind == 0) {
        int made = 0;
        for (int x=0; x<d && made<n; ++x)
            for (int y=0; y<d && made<n; ++y)
                for (int z=0; z<d && made<n; ++z,++made)
                    cout << setprecision(17) << (x+.5)/d << ' ' << (y+.5)/d << ' ' << (z+.5)/d << '\n';
        return 0;
    }

    // Retreat by a tiny relative amount so roundoff cannot put a lattice point
    // just beyond a face after decimal conversion.
    double r = best.r * (1.0 - 1e-10);
    double a = (best.kind == 1 ? 2.0 * sqrt(2.0) * r : 4.0 * r / sqrt(3.0));
    double sx = a * best.px / P, sy = a * best.py / P, sz = a * best.pz / P;
    static const double fcc[4][3] = {{0,0,0},{0,.5,.5},{.5,0,.5},{.5,.5,0}};
    static const double bcc[2][3] = {{0,0,0},{.5,.5,.5}};
    const double (*bases)[3] = best.kind == 1 ? fcc : bcc;
    int nb = best.kind == 1 ? 4 : 2;
    int made = 0;
    cout << setprecision(17);
    for (int q=0; q<nb && made<n; ++q) {
        long long lx = (long long)ceil((r-sx)/a-bases[q][0]-1e-12);
        long long hx = (long long)floor((1-r-sx)/a-bases[q][0]+1e-12);
        long long ly = (long long)ceil((r-sy)/a-bases[q][1]-1e-12);
        long long hy = (long long)floor((1-r-sy)/a-bases[q][1]+1e-12);
        long long lz = (long long)ceil((r-sz)/a-bases[q][2]-1e-12);
        long long hz = (long long)floor((1-r-sz)/a-bases[q][2]+1e-12);
        for (long long i=lx; i<=hx && made<n; ++i)
            for (long long j=ly; j<=hy && made<n; ++j)
                for (long long k=lz; k<=hz && made<n; ++k) {
                    double x=sx+a*(i+bases[q][0]);
                    double y=sy+a*(j+bases[q][1]);
                    double z=sz+a*(k+bases[q][2]);
                    // The inequalities above leave these safely inside; clamp
                    // only insignificant arithmetic noise at a cube boundary.
                    x=min(1.0,max(0.0,x)); y=min(1.0,max(0.0,y)); z=min(1.0,max(0.0,z));
                    cout << x << ' ' << y << ' ' << z << '\n';
                    ++made;
                }
    }
    return 0;
}
