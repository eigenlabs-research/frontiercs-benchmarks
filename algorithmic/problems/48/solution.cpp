#include <bits/stdc++.h>
using namespace std;

struct Point { double x, y, z; };

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n;
    if (!(cin >> n)) return 0;

    vector<Point> ans;
    // Small codes in the vertices of a cube give substantially better boundary
    // behavior than a truncated bulk lattice.
    if (n == 2) {
        double r = sqrt(3.0) / (2.0 * (1.0 + sqrt(3.0)));
        ans = {{r,r,r}, {1-r,1-r,1-r}};
    } else if (n == 3 || n == 4) {
        double r = sqrt(2.0) / (2.0 * (1.0 + sqrt(2.0)));
        int v[4][3] = {{0,0,0},{0,1,1},{1,0,1},{1,1,0}};
        for (int q=0; q<n; ++q)
            ans.push_back({v[q][0] ? 1-r : r, v[q][1] ? 1-r : r, v[q][2] ? 1-r : r});
    } else if (n <= 8) {
        // Any subset of cube vertices; an edge is tight at radius 1/4.
        double r = .25;
        for (int q=0; q<n; ++q)
            ans.push_back({(q&1)?1-r:r, (q&2)?1-r:r, (q&4)?1-r:r});
    } else {
        // Candidate A: an FCC lattice.  Integer triples of even sum, divided
        // by sqrt(2), have nearest-neighbour distance one.
        int D = 0;
        while ((((long long)(D+1)*(D+1)*(D+1) + 1) / 2) < n) ++D;
        double fccR = 1.0 / (2.0 * (1.0 + D / sqrt(2.0)));

        // Candidate B is precisely the baseline-style cubic grid.  Selecting
        // a prefix cannot reduce its certified face margin or pair distance.
        int M = 1;
        while ((long long)M*M*M < n) ++M;
        double gridR = 1.0 / (2.0 * M);

        if (fccR > gridR) {
            double step = 1.0 / (1.0 + D / sqrt(2.0));
            double r = step * .5;
            for (int i=0; i<=D && (int)ans.size()<n; ++i)
                for (int j=0; j<=D && (int)ans.size()<n; ++j)
                    for (int k=0; k<=D && (int)ans.size()<n; ++k)
                        if (((i+j+k)&1)==0)
                            ans.push_back({r + step*i/sqrt(2.0),
                                           r + step*j/sqrt(2.0),
                                           r + step*k/sqrt(2.0)});
        } else {
            double step = 1.0 / M;
            for (int i=0; i<M && (int)ans.size()<n; ++i)
                for (int j=0; j<M && (int)ans.size()<n; ++j)
                    for (int k=0; k<M && (int)ans.size()<n; ++k)
                        ans.push_back({(i+.5)*step, (j+.5)*step, (k+.5)*step});
        }
    }
    cout << setprecision(17);
    for (const auto &p : ans) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
}
