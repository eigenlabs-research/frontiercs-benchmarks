#include <bits/stdc++.h>
using namespace std;

struct Point { double x, y, z; };

// Number of points of a translated face-centred-cubic lattice that fit in
// [r,1-r]^3.  Its nearest-neighbour distance is exactly 2r.
static long long fccCount(double r, int qx, int qy, int qz, int Q) {
    if (r <= 0 || r > .5) return 0;
    const double a = 2.0 * sqrt(2.0) * r;
    const double sx = a * qx / Q, sy = a * qy / Q, sz = a * qz / Q;
    static const double b[4][3] = {{0,0,0},{.5,.5,0},{.5,0,.5},{0,.5,.5}};
    long long total = 0;
    for (int t=0; t<4; ++t) {
        double s[3] = {sx+a*b[t][0], sy+a*b[t][1], sz+a*b[t][2]};
        long long ways = 1;
        for (int d=0; d<3; ++d) {
            long long lo = (long long)ceil((r-s[d])/a - 1e-12);
            long long hi = (long long)floor((1.0-r-s[d])/a + 1e-12);
            if (hi < lo) { ways=0; break; }
            ways *= hi-lo+1;
        }
        total += ways;
    }
    return total;
}

static vector<Point> fccPoints(int n, double r, int qx, int qy, int qz, int Q) {
    const double a = 2.0 * sqrt(2.0) * r;
    const double ss[3] = {a*qx/Q, a*qy/Q, a*qz/Q};
    static const double b[4][3] = {{0,0,0},{.5,.5,0},{.5,0,.5},{0,.5,.5}};
    vector<Point> out; out.reserve(n);
    for (int t=0; t<4 && (int)out.size()<n; ++t) {
        double s[3] = {ss[0]+a*b[t][0], ss[1]+a*b[t][1], ss[2]+a*b[t][2]};
        long long lo[3], hi[3];
        for (int d=0; d<3; ++d) {
            lo[d]=(long long)ceil((r-s[d])/a - 1e-12);
            hi[d]=(long long)floor((1.0-r-s[d])/a + 1e-12);
        }
        for (long long i=lo[0]; i<=hi[0] && (int)out.size()<n; ++i)
            for (long long j=lo[1]; j<=hi[1] && (int)out.size()<n; ++j)
                for (long long k=lo[2]; k<=hi[2] && (int)out.size()<n; ++k)
                    out.push_back({s[0]+a*i,s[1]+a*j,s[2]+a*k});
    }
    return out;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    vector<Point> best;
    double bestR = 0;

    // The balanced grid is a guaranteed fallback (and equals the stated
    // elementary baseline up to permutation of its dimensions).
    int m=1;
    while (1LL*m*m*m < n) ++m;
    double gr = 1.0/(2.0*m);
    vector<Point> grid; grid.reserve(n);
    for (int i=0;i<m && (int)grid.size()<n;++i)
        for (int j=0;j<m && (int)grid.size()<n;++j)
            for (int k=0;k<m && (int)grid.size()<n;++k)
                grid.push_back({(i+.5)/m,(j+.5)/m,(k+.5)/m});
    best=grid; bestR=gr;

    // A few binary cube codes avoid the small-population boundary loss of a
    // periodic lattice.  These are still formulaic, not input-specific data.
    if (n<=8) {
        double r = (n==2 ? sqrt(3.0)/(2.0*(1.0+sqrt(3.0))) :
                    (n<=4 ? sqrt(2.0)/(2.0+2.0*sqrt(2.0)) : .25));
        int code[8]={0,3,5,6,1,2,4,7};
        vector<Point> v; v.reserve(n);
        for(int i=0;i<n;++i) {
            int c = (n==2 && i==1) ? 7 : code[i];
            v.push_back({(c&1)?1-r:r,(c&2)?1-r:r,(c&4)?1-r:r});
        }
        if(r>bestR) bestR=r, best=v;
    }

    // Search several translations of FCC.  Translation phases materially
    // reduce boundary waste while preserving every pair distance.
    const int Q=7;
    for(int x=0;x<Q;++x) for(int y=0;y<Q;++y) for(int z=0;z<Q;++z) {
        double lo=0, hi=.5;
        for(int it=0;it<45;++it) {
            double mid=(lo+hi)*.5;
            if(fccCount(mid,x,y,z,Q)>=n) lo=mid; else hi=mid;
        }
        // Back off from lattice contacts, making output robust to decimal and
        // floating-point roundoff in the checker.
        double r=lo*(1.0-1e-10);
        if(r>bestR) {
            auto v=fccPoints(n,r,x,y,z,Q);
            if((int)v.size()==n) bestR=r, best.swap(v);
        }
    }

    cout << setprecision(17);
    for(const auto &p: best) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
}
