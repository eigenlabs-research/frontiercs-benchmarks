#include <bits/stdc++.h>
using namespace std;

struct P { double x, y, z; };

// This is r^2 for precisely the quantity evaluated by the checker.
static double quality2(const vector<P>& a) {
    double v = 1.0;
    int n = (int)a.size();
    for (const P& p : a) {
        double f = min(min(p.x, 1.0-p.x), min(p.y, 1.0-p.y));
        f = min(f, min(p.z, 1.0-p.z));
        v = min(v, f*f);
    }
    for (int i=0; i<n; ++i) for (int j=0; j<i; ++j) {
        double dx=a[i].x-a[j].x, dy=a[i].y-a[j].y, dz=a[i].z-a[j].z;
        v=min(v, 0.25*(dx*dx+dy*dy+dz*dz));
    }
    return v;
}

static vector<P> cube_seed(int n) {
    int q=1;
    while (q*q*q<n) ++q;
    vector<P> a; a.reserve(n);
    for (int z=0; z<q && (int)a.size()<n; ++z)
        for (int y=0; y<q && (int)a.size()<n; ++y)
            for (int x=0; x<q && (int)a.size()<n; ++x)
                a.push_back({(x+.5)/q,(y+.5)/q,(z+.5)/q});
    return a;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    // The two-sphere optimum is useful both as a correctness anchor and to
    // avoid asking a stochastic search to discover a one-dimensional fact.
    if (n == 2) {
        double r=sqrt(3.0)/(2.0*(1.0+sqrt(3.0)));
        cout << setprecision(17);
        cout << r << ' ' << r << ' ' << r << '\n';
        cout << 1-r << ' ' << 1-r << ' ' << 1-r << '\n';
        return 0;
    }

    vector<P> seed=cube_seed(n);
    // A regular tetrahedron is an exact feasible seed for 3 and 4 points.
    // It is supplied to the same max-min search, rather than selected among
    // lattice families.
    if (n==3 || n==4) {
        double r=1.0/(2.0+sqrt(2.0));
        double h=1.0-r;
        seed={{r,r,r},{r,h,h},{h,r,h}};
        if (n==4) seed.push_back({h,h,r});
    }
    vector<P> answer=seed;
    double best=quality2(answer);

    // Exact bottleneck annealing is deliberately restricted to the range in
    // which evaluating every pair is cheap.  Larger instances keep the safe
    // balanced seed rather than risking a time-limit-dependent construction.
    if (n <= 32) {
        int starts = n<=8 ? 6 : (n<=16 ? 4 : 2);
        int rounds = n<=8 ? 35000 : 22000;
        mt19937_64 rng(0x9e3779b97f4a7c15ULL ^ (unsigned)n*0xbf58476d1ce4e5b9ULL);
        uniform_real_distribution<double> uni(0.0,1.0);
        for (int s=0; s<starts; ++s) {
            vector<P> cur = (s==0 ? seed : vector<P>(n));
            if (s) for (P& p:cur) {
                // A modest initial margin avoids spending most trials repairing
                // points nearly on a face.
                p={.06+.88*uni(rng), .06+.88*uni(rng), .06+.88*uni(rng)};
            }
            double now=quality2(cur);
            if (now>best) best=now, answer=cur;
            for (int it=0; it<rounds; ++it) {
                int k=(int)(rng()%n);
                P old=cur[k], trial=old;
                double t=(double)it/(rounds-1);
                double step=.22*(1.0-t)+.0025*t;
                if ((rng()%101)==0) {
                    trial={uni(rng),uni(rng),uni(rng)};
                } else {
                    trial.x=min(1.0,max(0.0,trial.x+(2*uni(rng)-1)*step));
                    trial.y=min(1.0,max(0.0,trial.y+(2*uni(rng)-1)*step));
                    trial.z=min(1.0,max(0.0,trial.z+(2*uni(rng)-1)*step));
                }
                cur[k]=trial;
                double nxt=quality2(cur);
                // Temperature is in radius-squared units, matching quality2.
                double temp=.008*(1.0-t)+.000003;
                bool take=nxt>=now || uni(rng)<exp((nxt-now)/temp);
                if (take) now=nxt; else cur[k]=old;
                if (now>best) best=now, answer=cur;
            }
        }
    }

    cout << setprecision(17);
    for (const P& p:answer) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
}
