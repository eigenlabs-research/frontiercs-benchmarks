#include <bits/stdc++.h>
using namespace std;

struct Pt { double x,y,z; };

static double radius_of(const vector<Pt>& p) {
    int n = (int)p.size();
    double r = 1e100;
    for (auto &a: p) {
        r = min(r, min(min(a.x, 1.0-a.x), min(a.y, min(1.0-a.y, min(a.z, 1.0-a.z)))));
    }
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
        double dx=p[i].x-p[j].x, dy=p[i].y-p[j].y, dz=p[i].z-p[j].z;
        r = min(r, 0.5*sqrt(dx*dx+dy*dy+dz*dz));
    }
    return r;
}

static vector<Pt> normalize_lattice(vector<array<int,3>> v) {
    int n = (int)v.size();
    int mn[3] = {INT_MAX,INT_MAX,INT_MAX}, mx[3] = {INT_MIN,INT_MIN,INT_MIN};
    for (auto &a: v) for (int d=0; d<3; ++d) mn[d]=min(mn[d],a[d]), mx[d]=max(mx[d],a[d]);
    double range[3] = {double(mx[0]-mn[0]), double(mx[1]-mn[1]), double(mx[2]-mn[2])};
    double R = max(range[0], max(range[1], range[2]));
    double dmin2 = 1e100;
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
        double dx=v[i][0]-v[j][0], dy=v[i][1]-v[j][1], dz=v[i][2]-v[j][2];
        dmin2 = min(dmin2, dx*dx+dy*dy+dz*dz);
    }
    double dmin = sqrt(dmin2);
    double scale = 1.0 / (R + dmin);
    double margin = 0.5 * dmin * scale;
    vector<Pt> p; p.reserve(n);
    for (auto &a: v) {
        double coord[3];
        for (int d=0; d<3; ++d) {
            double extra = 0.5 * (R - range[d]);
            coord[d] = margin + (double(a[d]-mn[d]) + extra) * scale;
        }
        p.push_back({coord[0],coord[1],coord[2]});
    }
    return p;
}

static vector<Pt> polish_small(vector<Pt> p) {
    int n = (int)p.size();
    if (n < 5 || n > 48) return p;
    vector<Pt> cur = p, best = p;
    double cr = radius_of(cur), br = cr;
    mt19937_64 rng(0x9e3779b97f4a7c15ULL ^ (uint64_t)n * 1000003ULL);
    uniform_real_distribution<double> U(0.0, 1.0);
    normal_distribution<double> G(0.0, 1.0);
    int iters = (n <= 12 ? 180000 : (n <= 24 ? 45000 : 12000));
    for (int it=0; it<iters; ++it) {
        double t = 0.060 * pow(0.001, double(it) / max(1, iters-1));
        int id = (int)(rng() % n);
        Pt old = cur[id];
        double step = 0.20 * pow(0.002, double(it) / max(1, iters-1));
        cur[id].x = min(1.0, max(0.0, cur[id].x + G(rng)*step));
        cur[id].y = min(1.0, max(0.0, cur[id].y + G(rng)*step));
        cur[id].z = min(1.0, max(0.0, cur[id].z + G(rng)*step));
        double nr = radius_of(cur);
        bool accept = nr >= cr || U(rng) < exp((nr - cr) / max(1e-9, t));
        if (accept) {
            cr = nr;
            if (nr > br) br = nr, best = cur;
        } else {
            cur[id] = old;
        }
        if ((it & 4095) == 4095 && cr + 1e-5 < br) cur = best, cr = br;
    }
    return (br > radius_of(p) + 1e-12) ? best : p;
}

static vector<Pt> cubic_candidate(int n) {
    int M = 0;
    while (1LL*(M+1)*(M+1)*(M+1) < n) ++M;
    vector<array<int,3>> all;
    all.reserve((M+1)*(M+1)*(M+1));
    double c = M/2.0;
    for (int i=0;i<=M;i++) for (int j=0;j<=M;j++) for (int k=0;k<=M;k++) all.push_back({i,j,k});
    sort(all.begin(), all.end(), [&](const auto& a, const auto& b){
        double da=(a[0]-c)*(a[0]-c)+(a[1]-c)*(a[1]-c)+(a[2]-c)*(a[2]-c);
        double db=(b[0]-c)*(b[0]-c)+(b[1]-c)*(b[1]-c)+(b[2]-c)*(b[2]-c);
        if (da != db) return da < db;
        return a < b;
    });
    all.resize(n);
    return normalize_lattice(all);
}

static vector<Pt> fcc_candidate(int n) {
    int M = 0;
    auto cnt = [](int M)->long long{
        long long e = M/2 + 1, o = (M+1) - e;
        return e*e*e + 3*e*o*o; // triples with even coordinate sum
    };
    while (cnt(M) < n) ++M;
    vector<array<int,3>> all;
    all.reserve((size_t)cnt(M));
    double c = M/2.0;
    for (int i=0;i<=M;i++) for (int j=0;j<=M;j++) for (int k=0;k<=M;k++) if (((i+j+k)&1)==0) all.push_back({i,j,k});
    sort(all.begin(), all.end(), [&](const auto& a, const auto& b){
        double da=(a[0]-c)*(a[0]-c)+(a[1]-c)*(a[1]-c)+(a[2]-c)*(a[2]-c);
        double db=(b[0]-c)*(b[0]-c)+(b[1]-c)*(b[1]-c)+(b[2]-c)*(b[2]-c);
        if (da != db) return da < db;
        return a < b;
    });
    all.resize(n);
    return normalize_lattice(all);
}

static vector<Pt> bcc_candidate(int n) {
    int M = 0;
    auto cnt = [](int M)->long long{
        long long e = M/2 + 1, o = (M+1) - e;
        return e*e*e + o*o*o; // all coordinates same parity
    };
    while (cnt(M) < n) ++M;
    vector<array<int,3>> all;
    all.reserve((size_t)cnt(M));
    double c = M/2.0;
    for (int i=0;i<=M;i++) for (int j=0;j<=M;j++) for (int k=0;k<=M;k++) if ((i&1)==(j&1) && (i&1)==(k&1)) all.push_back({i,j,k});
    sort(all.begin(), all.end(), [&](const auto& a, const auto& b){
        double da=(a[0]-c)*(a[0]-c)+(a[1]-c)*(a[1]-c)+(a[2]-c)*(a[2]-c);
        double db=(b[0]-c)*(b[0]-c)+(b[1]-c)*(b[1]-c)+(b[2]-c)*(b[2]-c);
        if (da != db) return da < db;
        return a < b;
    });
    all.resize(n);
    return normalize_lattice(all);
}

static vector<Pt> small_candidate(int n) {
    vector<Pt> p;
    if (n == 2) {
        double r = sqrt(3.0)/(2.0+2.0*sqrt(3.0));
        p = {{r,r,r},{1-r,1-r,1-r}};
    } else if (n == 3 || n == 4) {
        double r = sqrt(2.0)/(2.0+2.0*sqrt(2.0));
        vector<Pt> t = {{r,r,r},{r,1-r,1-r},{1-r,r,1-r},{1-r,1-r,r}};
        p.assign(t.begin(), t.begin()+n);
    }
    return p;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    vector<vector<Pt>> cand;
    cand.push_back(cubic_candidate(n));
    cand.push_back(fcc_candidate(n));
    cand.push_back(bcc_candidate(n));
    if (n <= 4) cand.push_back(small_candidate(n));

    int best = 0;
    double br = -1;
    for (int i=0;i<(int)cand.size();++i) {
        if ((int)cand[i].size() != n) continue;
        double r = radius_of(cand[i]);
        if (r > br) br = r, best = i;
    }

    cand[best] = polish_small(cand[best]);

    cout.setf(ios::fixed); cout << setprecision(12);
    for (auto &p: cand[best]) {
        double x=min(1.0,max(0.0,p.x)), y=min(1.0,max(0.0,p.y)), z=min(1.0,max(0.0,p.z));
        cout << x << ' ' << y << ' ' << z << '\n';
    }
    return 0;
}
