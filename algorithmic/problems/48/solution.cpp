#include <bits/stdc++.h>
using namespace std;

struct P { int x, y, z; };
struct Candidate {
    vector<P> p;
    double r = -1;
};

// Put an integer lattice subset in the cube.  Its shortest lattice distance
// is d; centering its bounding box makes the limiting face and pair distances
// equal whenever the widest coordinate direction is tight.
static Candidate scale_lattice(vector<P> p, double d) {
    int lx = p[0].x, hx = p[0].x, ly = p[0].y, hy = p[0].y, lz = p[0].z, hz = p[0].z;
    for (const P &q : p) {
        lx=min(lx,q.x); hx=max(hx,q.x); ly=min(ly,q.y); hy=max(hy,q.y); lz=min(lz,q.z); hz=max(hz,q.z);
    }
    double span = max({double(hx-lx), double(hy-ly), double(hz-lz)});
    Candidate ans;
    ans.p = move(p);
    ans.r = 1.0 / (2.0 + 2.0 * span / d);
    return ans;
}

static vector<P> choose_centered(vector<P> a, int n, int K) {
    // A centered prefix avoids the lexicographic thin-slab failure for partial
    // lattice boxes.  The comparison is exact after multiplying coordinates by 2.
    sort(a.begin(), a.end(), [K](const P& u, const P& v) {
        long long ux=2LL*u.x-K, uy=2LL*u.y-K, uz=2LL*u.z-K;
        long long vx=2LL*v.x-K, vy=2LL*v.y-K, vz=2LL*v.z-K;
        long long du=ux*ux+uy*uy+uz*uz, dv=vx*vx+vy*vy+vz*vz;
        if (du != dv) return du < dv;
        if (u.x != v.x) return u.x < v.x;
        if (u.y != v.y) return u.y < v.y;
        return u.z < v.z;
    });
    a.resize(n);
    return a;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    Candidate best;

    // FCC: integer triples with even coordinate sum, nearest-neighbor distance sqrt(2).
    int K = 0;
    while (true) {
        long long cnt = 0;
        for (int x=0; x<=K; ++x) for (int y=0; y<=K; ++y) for (int z=0; z<=K; ++z)
            if (((x+y+z)&1)==0) ++cnt;
        if (cnt >= n) break;
        ++K;
    }
    vector<P> a;
    a.reserve(n + 64);
    for (int x=0; x<=K; ++x) for (int y=0; y<=K; ++y) for (int z=0; z<=K; ++z)
        if (((x+y+z)&1)==0) a.push_back({x,y,z});
    Candidate fcc = scale_lattice(choose_centered(move(a), n, K), sqrt(2.0));
    best = move(fcc);

    // BCC is superior for a few very small finite boxes (notably n=2).
    K = 0;
    while (true) {
        long long cnt = 0;
        for (int x=0; x<=K; ++x) for (int y=0; y<=K; ++y) for (int z=0; z<=K; ++z)
            if ((x&1)==(y&1) && (y&1)==(z&1)) ++cnt;
        if (cnt >= n) break;
        ++K;
    }
    a.clear();
    for (int x=0; x<=K; ++x) for (int y=0; y<=K; ++y) for (int z=0; z<=K; ++z)
        if ((x&1)==(y&1) && (y&1)==(z&1)) a.push_back({x,y,z});
    Candidate bcc = scale_lattice(choose_centered(move(a), n, K), sqrt(3.0));
    bool use_bcc = false;
    if (bcc.r > best.r) { best = move(bcc); use_bcc = true; }

    // The balanced cubic grid is a guaranteed finite-size fallback.
    int m = 1;
    while (1LL*m*m*m < n) ++m;
    Candidate grid;
    grid.r = 1.0/(2.0*m);
    grid.p.reserve(n);
    for (int x=0; x<m && (int)grid.p.size()<n; ++x)
        for (int y=0; y<m && (int)grid.p.size()<n; ++y)
            for (int z=0; z<m && (int)grid.p.size()<n; ++z)
                grid.p.push_back({x,y,z});
    if (grid.r > best.r) {
        cout << setprecision(17);
        for (const P &q : grid.p)
            cout << (q.x+.5)/m << ' ' << (q.y+.5)/m << ' ' << (q.z+.5)/m << '\n';
        return 0;
    }

    int lx=best.p[0].x, hx=lx, ly=best.p[0].y, hy=ly, lz=best.p[0].z, hz=lz;
    for (const P &q : best.p) { lx=min(lx,q.x); hx=max(hx,q.x); ly=min(ly,q.y); hy=max(hy,q.y); lz=min(lz,q.z); hz=max(hz,q.z); }
    // Equal-score ties keep FCC, which is the denser asymptotic construction.
    double d = use_bcc ? sqrt(3.0) : sqrt(2.0);
    double s = 2.0*best.r/d;
    double cx=(lx+hx)*.5, cy=(ly+hy)*.5, cz=(lz+hz)*.5;
    cout << setprecision(17);
    for (const P &q : best.p)
        cout << .5+s*(q.x-cx) << ' ' << .5+s*(q.y-cy) << ' ' << .5+s*(q.z-cz) << '\n';
}
