// Sphere Packing in a Cube (optimization)
//
// For each n on stdin: output n centers in [0,1]^3 maximizing the common radius
//   r(C) = min( 0.5 * min_{i!=j} ||ci-cj|| ,  min_i dist(ci, cube faces) ).
//
// Method (per n):
//   * Seed with an FCC lattice (densest 3D lattice -> approaches Kepler bound),
//     then run several candidate optimizations with different starts:
//        - FCC-centered subset
//        - cubic grid
//        - maximin (farthest-point) seeding from a random/RNG base
//        - a few randomized perturbations of the best-so-far
//   * Each candidate is relaxed by projection at an increasing target radius R:
//        - clamp into [R, 1-R]  (wall constraint)
//        - separate pairs closer than 2R  (non-overlap)
//     bumping R (warm-started) drives r(C) up to feasibility.
//   * Keep the configuration with the largest achieved r(C).
//   * geom_radius uses a spatial hash (O(n)) for the feasibility gate; an exact
//     O(n^2) value is used only once for the final reported radius.

#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <random>

using namespace std;
using A = array<double, 3>;
using V = vector<A>;

static const double INF = 1e300;

static inline double face_gap(const A& p) {
    double g = p[0];
    if (1.0 - p[0] < g) g = 1.0 - p[0];
    if (p[1] < g) g = p[1];
    if (1.0 - p[1] < g) g = 1.0 - p[1];
    if (p[2] < g) g = p[2];
    if (1.0 - p[2] < g) g = 1.0 - p[2];
    return g;
}

static double geom_radius_exact(const V& c) {
    int n = (int)c.size();
    if (n == 0) return 0.0;
    double fa = INF;
    for (int i = 0; i < n; i++) { double g = face_gap(c[i]); if (g < fa) fa = g; }
    if (n == 1) return fa;
    double m2 = INF;
    for (int i = 0; i < n; i++) {
        const A& a = c[i];
        for (int j = i + 1; j < n; j++) {
            double dx = a[0]-c[j][0], dy = a[1]-c[j][1], dz = a[2]-c[j][2];
            double d2 = dx*dx+dy*dy+dz*dz; if (d2 < m2) m2 = d2;
        }
    }
    double r = min(fa, 0.5*sqrt(m2));
    if (!isfinite(r) || r < 0) r = 0.0;
    return r;
}

struct Hash {
    double cell = 0; int G = 1;
    vector<int> hd, nxt;
    static inline int ic(int v, int hi) { return v < 0 ? 0 : (v > hi ? hi : v); }
    void build(const V& p, double R) {
        int n = (int)p.size();
        cell = max(2.0 * R, 1e-7);
        G = max(1, (int)floor(1.0 / cell) + 1);
        if (G > 4000) G = 4000;
        hd.assign((size_t)G * G * G, -1);
        nxt.assign(n, -1);
        for (int i = 0; i < n; i++) {
            int cx = ic((int)floor(p[i][0]/cell), G-1);
            int cy = ic((int)floor(p[i][1]/cell), G-1);
            int cz = ic((int)floor(p[i][2]/cell), G-1);
            long id = ((long)cx*G + cy)*G + cz;
            nxt[i] = hd[id]; hd[id] = i;
        }
    }
};

static double geom_radius_hash(const V& c, double R) {
    int n = (int)c.size();
    if (n == 0) return 0.0;
    double fa = INF;
    for (int i = 0; i < n; i++) { double g = face_gap(c[i]); if (g < fa) fa = g; }
    if (n == 1) return fa;
    Hash h; h.build(c, R);
    int G = h.G;
    double m2 = INF;
    for (int i = 0; i < n; i++) {
        int cx = Hash::ic((int)floor(c[i][0]/h.cell), G-1);
        int cy = Hash::ic((int)floor(c[i][1]/h.cell), G-1);
        int cz = Hash::ic((int)floor(c[i][2]/h.cell), G-1);
        for (int dx = -1; dx <= 1; dx++) {
            int x = cx+dx; if (x<0||x>=G) continue;
            for (int dy = -1; dy <= 1; dy++) {
                int y = cy+dy; if (y<0||y>=G) continue;
                for (int dz = -1; dz <= 1; dz++) {
                    int z = cz+dz; if (z<0||z>=G) continue;
                    long id = ((long)x*G + y)*G + z;
                    for (int j = h.hd[id]; j != -1; j = h.nxt[j]) {
                        if (j <= i) continue;
                        double dx2 = c[i][0]-c[j][0], dy2 = c[i][1]-c[j][1], dz2 = c[i][2]-c[j][2];
                        double d2 = dx2*dx2+dy2*dy2+dz2*dz2;
                        if (d2 < m2) m2 = d2;
                    }
                }
            }
        }
    }
    double r = min(fa, 0.5*sqrt(m2));
    if (!isfinite(r) || r < 0) r = 0.0;
    return r;
}

// FCC lattice, nearest-n to center (interior-biased). Returns empty if insufficient.
static V fcc_init(int n) {
    V pts;
    double d = pow(sqrt(2.0) / n, 1.0/3.0);
    double a = d * sqrt(2.0);
    int K = (int)floor(1.0 / a) + 1;
    double off[4][3] = {{0,0,0},{0.5,0.5,0},{0.5,0,0.5},{0,0.5,0.5}};
    for (int i = 0; i <= K; i++)
        for (int j = 0; j <= K; j++)
            for (int k = 0; k <= K; k++)
                for (int t = 0; t < 4; t++) {
                    double x = (i + off[t][0]) * a;
                    double y = (j + off[t][1]) * a;
                    double z = (k + off[t][2]) * a;
                    if (x>=-1e-9 && x<=1+1e-9 && y>=-1e-9 && y<=1+1e-9 && z>=-1e-9 && z<=1+1e-9)
                        pts.push_back({x,y,z});
                }
    if ((int)pts.size() < n) return V();
    sort(pts.begin(), pts.end(), [](const A& p, const A& q) {
        double dp = (p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)+(p[2]-0.5)*(p[2]-0.5);
        double dq = (q[0]-0.5)*(q[0]-0.5)+(q[1]-0.5)*(q[1]-0.5)+(q[2]-0.5)*(q[2]-0.5);
        return dp < dq;
    });
    pts.resize(n);
    return pts;
}

static V cubic_init(int n) {
    int m = (int)ceil(cbrt((double)n));
    if (m < 1) m = 1;
    V all;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            for (int k = 0; k < m; k++)
                all.push_back({(i+0.5)/m, (j+0.5)/m, (k+0.5)/m});
    if ((int)all.size() < n) all.resize(n, {0.5,0.5,0.5});
    sort(all.begin(), all.end(), [](const A& p, const A& q) {
        double dp = (p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)+(p[2]-0.5)*(p[2]-0.5);
        double dq = (q[0]-0.5)*(q[0]-0.5)+(q[1]-0.5)*(q[1]-0.5)+(q[2]-0.5)*(q[2]-0.5);
        return dp < dq;
    });
    all.resize(n);
    return all;
}

// Maximin (farthest-point) seeding: pick 1st near center, then each next point
// farthest from the already-chosen set (greedy), within the cube. Good spread.
static V maximin_init(int n, mt19937_64& rng) {
    V pts; pts.reserve(n);
    uniform_real_distribution<double> U(0.05, 0.95);
    pts.push_back({0.5, 0.5, 0.5});
    // initial random pool for first few
    for (int i = 1; i < n; i++) {
        double best_d2 = -1; A best = {0.5,0.5,0.5};
        for (int s = 0; s < 40; s++) {
            A cand = {U(rng), U(rng), U(rng)};
            double md2 = INF;
            for (auto& p : pts) {
                double dx = cand[0]-p[0], dy = cand[1]-p[1], dz = cand[2]-p[2];
                double d2 = dx*dx+dy*dy+dz*dz; if (d2 < md2) md2 = d2;
            }
            if (md2 > best_d2) { best_d2 = md2; best = cand; }
        }
        pts.push_back(best);
    }
    return pts;
}

static void project(V& p, double R, int iters) {
    int n = (int)p.size();
    double twoR = 2.0*R, twoR2 = twoR*twoR;
    for (int it = 0; it < iters; it++) {
        for (int i = 0; i < n; i++)
            for (int d = 0; d < 3; d++) {
                if (p[i][d] < R) p[i][d] = R;
                else if (p[i][d] > 1.0-R) p[i][d] = 1.0-R;
            }
        Hash h; h.build(p, R); int G = h.G;
        for (int i = 0; i < n; i++) {
            int cx = Hash::ic((int)floor(p[i][0]/h.cell), G-1);
            int cy = Hash::ic((int)floor(p[i][1]/h.cell), G-1);
            int cz = Hash::ic((int)floor(p[i][2]/h.cell), G-1);
            for (int dx = -1; dx <= 1; dx++) {
                int x = cx+dx; if (x<0||x>=G) continue;
                for (int dy = -1; dy <= 1; dy++) {
                    int y = cy+dy; if (y<0||y>=G) continue;
                    for (int dz = -1; dz <= 1; dz++) {
                        int z = cz+dz; if (z<0||z>=G) continue;
                        long id = ((long)x*G + y)*G + z;
                        for (int j = h.hd[id]; j != -1; j = h.nxt[j]) {
                            if (j <= i) continue;
                            double ex = p[i][0]-p[j][0], ey = p[i][1]-p[j][1], ez = p[i][2]-p[j][2];
                            double d2 = ex*ex+ey*ey+ez*ez;
                            if (d2 < twoR2) {
                                if (d2 < 1e-18) {
                                    p[i][0] += R*0.5; p[j][0] -= R*0.5;
                                    if (p[i][0] > 1.0-R) p[i][0] = 1.0-R;
                                    if (p[j][0] < R) p[j][0] = R;
                                    continue;
                                }
                                double d = sqrt(d2); double push = (twoR - d)*0.5;
                                double ux = ex/d, uy = ey/d, uz = ez/d;
                                p[i][0]+=ux*push; p[i][1]+=uy*push; p[i][2]+=uz*push;
                                p[j][0]-=ux*push; p[j][1]-=uy*push; p[j][2]-=uz*push;
                            }
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < n; i++)
        for (int d = 0; d < 3; d++) {
            if (p[i][d] < R) p[i][d] = R; else if (p[i][d] > 1.0-R) p[i][d] = 1.0-R;
        }
}

// Relax a config upward to maximum feasible radius, using at most `budget` seconds.
static pair<V,double> optimize(const V& start, double budget, int max_iter) {
    V pts = start;
    double R = geom_radius_exact(pts);
    auto t0 = chrono::steady_clock::now();
    auto elapsed = [&]() { return chrono::duration<double>(chrono::steady_clock::now()-t0).count(); };
    for (int t = 0; t < max_iter; t++) {
        if (elapsed() > budget) break;
        double Rtry = R + (0.5 - R)*0.05 + 1e-6;
        V cand = pts; project(cand, Rtry, 20);
        double r = geom_radius_hash(cand, Rtry);
        if (r >= R + 1e-12) { pts = cand; R = r; continue; }
        double Rtry2 = R + (0.5 - R)*0.015 + 1e-6;
        V cand2 = pts; project(cand2, Rtry2, 20);
        double r2 = geom_radius_hash(cand2, Rtry2);
        if (r2 >= R + 1e-12) { pts = cand2; R = r2; continue; }
        // tiny step
        double Rtry3 = R + 1e-4;
        V cand3 = pts; project(cand3, Rtry3, 20);
        double r3 = geom_radius_hash(cand3, Rtry3);
        if (r3 >= R + 1e-12) { pts = cand3; R = r3; continue; }
        break;
    }
    project(pts, R, 30);
    return {pts, geom_radius_exact(pts)};
}

// Corner-biased start for very small n (FCC subsets underperform for tiny n).
static V corner_init(int n) {
    // place points near cube corners / face centers in a spread-out pattern
    V cand;
    double c[8][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1}};
    if (n <= 8) {
        for (int i = 0; i < n; i++) cand.push_back({c[i][0], c[i][1], c[i][2]});
        return cand;
    }
    // fallback: cubic grid
    return cubic_init(n);
}

static void solve_one(int n) {
    if (n <= 0) return;
    if (n == 1) { printf("0.5000000000 0.5000000000 0.5000000000\n"); return; }
    double total_budget = 0.80;
    mt19937_64 rng(12345ULL + (unsigned long)n * 2654435761ULL);

    V best = cubic_init(n);
    double bestR = geom_radius_exact(best);

    // candidate starts
    vector<V> starts;
    V f = fcc_init(n); if (!f.empty()) starts.push_back(f);
    starts.push_back(cubic_init(n));
    if (n <= 8) starts.push_back(corner_init(n));
    starts.push_back(maximin_init(n, rng));

    double per = total_budget * 0.5 / (double)starts.size();
    for (auto& s : starts) {
        auto [p, r] = optimize(s, per, 4000);
        if (r > bestR + 1e-12) { best = p; bestR = r; }
    }
    // random-restart perturbations of the best, time-boxed (scaled down for large n)
    double restart_budget = total_budget * 0.45;
    if (n > 1000) restart_budget = total_budget * 0.25;
    if (n > 2000) restart_budget = total_budget * 0.15;
    auto t0b = chrono::steady_clock::now();
    while (chrono::duration<double>(chrono::steady_clock::now()-t0b).count() < restart_budget) {
        V cand = best;
        normal_distribution<double> N(0.0, 0.02);
        for (auto& pt : cand) {
            pt[0] = min(1.0, max(0.0, pt[0] + N(rng)));
            pt[1] = min(1.0, max(0.0, pt[1] + N(rng)));
            pt[2] = min(1.0, max(0.0, pt[2] + N(rng)));
        }
        double R0 = bestR;
        for (auto& pt : cand) for (int d=0;d<3;d++){ if(pt[d]<R0)pt[d]=R0; else if(pt[d]>1-R0)pt[d]=1-R0; }
        auto [p, r] = optimize(cand, total_budget * 0.06, 2500);
        if (r > bestR + 1e-12) { best = p; bestR = r; }
    }

    for (int i = 0; i < n; i++) printf("%.10f %.10f %.10f\n", best[i][0], best[i][1], best[i][2]);
}

int main() {
    int n;
    while (scanf("%d", &n) == 1) solve_one(n);
    return 0;
}
